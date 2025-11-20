from qgis.PyQt.QtCore import QVariant
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFile,
    QgsProcessingParameterFeatureSink,
    QgsFeatureSink,
    QgsFields,
    QgsField,
    QgsFeature,
    QgsPointXY,
    QgsGeometry,
    QgsWkbTypes,
    QgsProcessingException,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsProject,
)
import math


class ApplyZenTemplateType2_WithVectors(QgsProcessingAlgorithm):

    INPUT = "INPUT"
    TEMPLATE = "TEMPLATE"
    OUTPUT_POINTS = "OUTPUT_POINTS"
    OUTPUT_VECTORS = "OUTPUT_VECTORS"

    def createInstance(self):
        return ApplyZenTemplateType2_WithVectors()

    def name(self):
        return "apply_template_to_waypoints"

    def displayName(self):
        return "2 - Apply Template to Waypoints"

    def group(self):
        return "Survey Tools"

    def groupId(self):
        return "survey_tools"

    def shortHelpString(self):
        return (
            "Template TYPE 2 avec offsets locaux X/Y/Z :\n"
            "- X = le long de la ligne (m)\n"
            "- Y = perpendiculaire à la ligne (m, + à gauche)\n"
            "- Z = vertical\n\n"
            "La couche de setups doit contenir STN, line_id, et line_az ou az (deg).\n"
            "Le script tourne la géométrie selon l'azimut et calcule :\n"
            "- dist = X, offset_y = Y, alongLineLabel = STN + round(X)."
        )

    # -------- TEMPLATE PARSER (TYPE 2) -------- #

    def parse_template(self, filepath):
        channels = []
        current = None

        with open(filepath, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("//"):
                    continue

                if line.startswith("<CH>"):
                    current = {}
                    continue
                if line.startswith("</CH>"):
                    if current:
                        channels.append(current)
                    current = None
                    continue

                if current is not None and "=" in line:
                    k, v = [p.strip() for p in line.split("=", 1)]
                    current[k] = v

        if not channels:
            raise QgsProcessingException("Aucun bloc <CH> trouvé dans le template.")

        if not any(("CH.NAME1" in ch) or ("CH.NAME2" in ch) for ch in channels):
            raise QgsProcessingException(
                "Template incompatible : TYPE 2 requis (au moins un CH.NAME1/CH.NAME2)."
            )

        return channels

    def parse_offset(self, s):
        try:
            X, Y, Z = map(float, s.split(":"))
            return X, Y, Z
        except Exception:
            return 0.0, 0.0, 0.0

    # -------- PARAMETERS -------- #

    def initAlgorithm(self, config=None):

        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT, "Stations (waypoints ZEN)", [QgsProcessing.TypeVectorPoint]
            )
        )

        self.addParameter(
            QgsProcessingParameterFile(
                self.TEMPLATE, "ZEN Template Type 2 (.stt)", extension="stt"
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_POINTS,
                "Electrodes / antennes (points UTM)",
                QgsProcessing.TypeVectorPoint,
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT_VECTORS,
                "Dipôles (vecteurs UTM)",
                QgsProcessing.TypeVectorLine,
            )
        )

    # -------- MAIN PROCESS -------- #

    def processAlgorithm(self, parameters, context, feedback):

        layer = self.parameterAsSource(parameters, self.INPUT, context)
        template_path = self.parameterAsFile(parameters, self.TEMPLATE, context)

        if layer is None:
            raise QgsProcessingException("Invalid input layer.")

        channels = self.parse_template(template_path)

        # Champs d'entrée
        field_names = [f.name() for f in layer.fields()]
        has_stn = "STN" in field_names
        has_line_id = "line_id" in field_names
        has_line_az = "line_az" in field_names
        has_az = "az" in field_names

        if not (has_line_az or has_az):
            raise QgsProcessingException(
                "La couche d'entrée doit contenir un champ 'line_az' ou 'az' (azimut de ligne en degrés)."
            )

        # CRS → UTM
        src_crs = layer.sourceCrs()
        feats = list(layer.getFeatures())
        if not feats:
            raise QgsProcessingException("Layer is empty.")

        first = feats[0].geometry().asPoint()

        if src_crs.isGeographic():
            lon = first.x()
            zone = int((lon + 180.0) / 6.0) + 1
            utm_epsg = 32600 + zone
        else:
            try:
                utm_epsg = int(src_crs.authid().split(":")[1])
            except Exception:
                raise QgsProcessingException(
                    "Impossible de déterminer l'EPSG UTM depuis le CRS source."
                )

        utm_crs = QgsCoordinateReferenceSystem(f"EPSG:{utm_epsg}")
        ct_to_utm = QgsCoordinateTransform(src_crs, utm_crs, QgsProject.instance())

        # Groupage par line_id pour chainage
        rows_by_line = {}

        for ft in feats:
            geom = ft.geometry()
            if geom.isEmpty():
                continue

            pt = geom.asPoint()
            utm = ct_to_utm.transform(pt)

            x, y = utm.x(), utm.y()
            stn = int(ft["STN"]) if has_stn and ft["STN"] is not None else int(ft.id())
            line_id = int(ft["line_id"]) if has_line_id and ft["line_id"] is not None else 0

            rows_by_line.setdefault(line_id, []).append((stn, x, y))

        if not rows_by_line:
            raise QgsProcessingException("No station with geometry.")

        chainage = {}
        for lid, rows in rows_by_line.items():
            rows.sort(key=lambda r: r[0])
            dist_acc = 0.0
            prev = None
            for stn, x, y in rows:
                if prev is None:
                    chainage[(lid, stn)] = 0.0
                else:
                    dx = x - prev[1]
                    dy = y - prev[2]
                    dist_acc += math.sqrt(dx * dx + dy * dy)
                    chainage[(lid, stn)] = dist_acc
                prev = (stn, x, y)

        # Champs de sortie - points
        p_fields = QgsFields()
        p_fields.append(QgsField("line_id", QVariant.Int))
        p_fields.append(QgsField("station", QVariant.Int))
        p_fields.append(QgsField("cmp", QVariant.String))
        p_fields.append(QgsField("ch_idx", QVariant.Int))
        p_fields.append(QgsField("name", QVariant.String))
        p_fields.append(QgsField("label", QVariant.String))
        p_fields.append(QgsField("end", QVariant.String))
        p_fields.append(QgsField("offset_z", QVariant.Double))
        p_fields.append(QgsField("ch_az", QVariant.Double))
        p_fields.append(QgsField("ch_incl", QVariant.Double))
        p_fields.append(QgsField("ant_num", QVariant.Int))
        p_fields.append(QgsField("stn_s", QVariant.Double))
        p_fields.append(QgsField("stn_E", QVariant.Double))
        p_fields.append(QgsField("stn_N", QVariant.Double))
        p_fields.append(QgsField("E", QVariant.Double))
        p_fields.append(QgsField("N", QVariant.Double))
        p_fields.append(QgsField("dist", QVariant.Double))        # X le long de la ligne
        p_fields.append(QgsField("offset_y", QVariant.Double))    # Y perpendiculaire
        p_fields.append(QgsField("alongLineLabel", QVariant.Int)) # STN + round(X)

        # Champs de sortie - vecteurs
        v_fields = QgsFields()
        v_fields.append(QgsField("line_id", QVariant.Int))
        v_fields.append(QgsField("station", QVariant.Int))
        v_fields.append(QgsField("cmp", QVariant.String))
        v_fields.append(QgsField("ch_idx", QVariant.Int))
        v_fields.append(QgsField("label1", QVariant.String))
        v_fields.append(QgsField("label2", QVariant.String))
        v_fields.append(QgsField("dE", QVariant.Double))
        v_fields.append(QgsField("dN", QVariant.Double))
        v_fields.append(QgsField("length", QVariant.Double))
        v_fields.append(QgsField("stn_s", QVariant.Double))

        points_sink, points_id = self.parameterAsSink(
            parameters, self.OUTPUT_POINTS, context,
            p_fields, QgsWkbTypes.Point, utm_crs
        )

        vectors_sink, vectors_id = self.parameterAsSink(
            parameters, self.OUTPUT_VECTORS, context,
            v_fields, QgsWkbTypes.LineString, utm_crs
        )

        # MAIN LOOP
        total = max(1, len(feats))

        for i, feat in enumerate(feats):

            if feedback.isCanceled():
                break

            geom = feat.geometry()
            if geom.isEmpty():
                continue

            pt = geom.asPoint()
            utm_pt = ct_to_utm.transform(pt)
            x0, y0 = utm_pt.x(), utm_pt.y()

            stn = int(feat["STN"]) if has_stn and feat["STN"] is not None else int(feat.id())
            line_id = int(feat["line_id"]) if has_line_id and feat["line_id"] is not None else 0
            stn_s = chainage.get((line_id, stn), 0.0)

            # Azimut de ligne (deg) -> dirE, dirN
            try:
                if has_line_az and feat["line_az"] is not None:
                    line_az = float(feat["line_az"])
                else:
                    line_az = float(feat["az"])
            except Exception:
                raise QgsProcessingException(
                    f"Valeur d'azimut invalide pour la station {stn}."
                )

            az_rad = math.radians(line_az)
            dirE = math.sin(az_rad)
            dirN = math.cos(az_rad)
            # perpendiculaire gauche (+90°)
            perpE = -dirN
            perpN = dirE

            for ch in channels:

                cmp_ = ch.get("CH.CMP", "")
                idx = int(ch.get("CH.INDEX", "0"))
                name1 = ch.get("CH.NAME1", "")
                name2 = ch.get("CH.NAME2", "")
                ch_az = float(ch.get("CH.AZIMUTH", 0) or 0)
                ch_incl = float(ch.get("CH.INCL", 0) or 0)

                has1 = False
                has2 = False
                E1 = N1 = E2 = N2 = None

                # -------- END 1 (local X/Y) -------- #
                if "CH.OFFSET.XYZ1" in ch:

                    X1, Y1, Z1 = self.parse_offset(ch["CH.OFFSET.XYZ1"])

                    # Géométrie en UTM
                    E1 = x0 + X1 * dirE + Y1 * perpE
                    N1 = y0 + X1 * dirN + Y1 * perpN

                    dist = X1
                    offy = Y1
                    along_val = int(stn + int(round(dist)))

                    f1 = QgsFeature(p_fields)
                    f1.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(E1, N1)))

                    f1["line_id"] = line_id
                    f1["station"] = stn
                    f1["cmp"] = cmp_
                    f1["ch_idx"] = idx
                    f1["name"] = f"{cmp_}{idx}_1"
                    f1["label"] = name1
                    f1["end"] = "1"
                    f1["offset_z"] = Z1
                    f1["ch_az"] = ch_az
                    f1["ch_incl"] = ch_incl
                    f1["ant_num"] = None
                    f1["stn_s"] = stn_s
                    f1["stn_E"] = x0
                    f1["stn_N"] = y0
                    f1["E"] = E1
                    f1["N"] = N1
                    f1["dist"] = dist
                    f1["offset_y"] = offy
                    f1["alongLineLabel"] = along_val

                    points_sink.addFeature(f1, QgsFeatureSink.FastInsert)
                    has1 = True

                # -------- END 2 (local X/Y) -------- #
                if "CH.OFFSET.XYZ2" in ch:

                    X2, Y2, Z2 = self.parse_offset(ch["CH.OFFSET.XYZ2"])

                    E2 = x0 + X2 * dirE + Y2 * perpE
                    N2 = y0 + X2 * dirN + Y2 * perpN

                    dist = X2
                    offy = Y2
                    along_val = int(stn + int(round(dist)))

                    f2 = QgsFeature(p_fields)
                    f2.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(E2, N2)))

                    f2["line_id"] = line_id
                    f2["station"] = stn
                    f2["cmp"] = cmp_
                    f2["ch_idx"] = idx
                    f2["name"] = f"{cmp_}{idx}_2"
                    f2["label"] = name2
                    f2["end"] = "2"
                    f2["offset_z"] = Z2
                    f2["ch_az"] = ch_az
                    f2["ch_incl"] = ch_incl
                    f2["ant_num"] = None
                    f2["stn_s"] = stn_s
                    f2["stn_E"] = x0
                    f2["stn_N"] = y0
                    f2["E"] = E2
                    f2["N"] = N2
                    f2["dist"] = dist
                    f2["offset_y"] = offy
                    f2["alongLineLabel"] = along_val

                    points_sink.addFeature(f2, QgsFeatureSink.FastInsert)
                    has2 = True

                # -------- VECTOR -------- #
                if has1 and has2 and E1 is not None and E2 is not None:
                    dE = E2 - E1
                    dN = N2 - N1
                    length = math.sqrt(dE * dE + dN * dN)

                    fv = QgsFeature(v_fields)
                    fv.setGeometry(
                        QgsGeometry.fromPolylineXY(
                            [QgsPointXY(E1, N1), QgsPointXY(E2, N2)]
                        )
                    )

                    fv["line_id"] = line_id
                    fv["station"] = stn
                    fv["cmp"] = cmp_
                    fv["ch_idx"] = idx
                    fv["label1"] = name1
                    fv["label2"] = name2
                    fv["dE"] = dE
                    fv["dN"] = dN
                    fv["length"] = length
                    fv["stn_s"] = stn_s

                    vectors_sink.addFeature(fv, QgsFeatureSink.FastInsert)

            feedback.setProgress(int(100 * (i + 1) / total))

        return {
            self.OUTPUT_POINTS: points_id,
            self.OUTPUT_VECTORS: vectors_id,
        }