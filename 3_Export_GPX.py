from qgis.PyQt.QtCore import QVariant
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterEnum,
    QgsProcessingException,
    QgsProject,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsPointXY,
)
import math


class ExportGPXFromZENStations(QgsProcessingAlgorithm):

    INPUT = "INPUT"
    OUTPUT = "OUTPUT"
    MODE = "MODE"
    NAME_STYLE = "NAME_STYLE"

    def createInstance(self):
        return ExportGPXFromZENStations()

    def name(self):
        return "export_gpx"

    def displayName(self):
        return "3 - Export Waypoints to GPX"

    def group(self):
        return "Survey Tools"

    def groupId(self):
        return "survey_tools"

    def shortHelpString(self):
        return (
            "Export ZEN data to GPX waypoints.\n\n"
            "Input layer should be the ZEN points (electrodes/antennas) layer, "
            "with at least:\n"
            "  - station (INT)  [setup ID]\n"
            "  - line_id (INT)  [optional but recommended]\n"
            "  - cmp (CMP: Ex, Ey, Hx, Hy, Hz, Tx, Ref, ...)\n"
            "  - ch_idx (INT)   [channel index]\n"
            "  - end (\"1\" or \"2\") [dipole endpoints]\n"
            "  - label (e.g. Ex1-, Ex2+, Hx4, ...)\n"
            "  - stn_E, stn_N (UTM coordinates of the setup)\n"
            "  - E, N (UTM coordinates of the actual point)\n"
            "  - alongLineLabel (DOUBLE) = chainage along the line\n"
            "  - offset_y (DOUBLE)       = local Y offset (m, +Y / -Y)\n\n"
            "Modes:\n"
            "  1) Setups (one waypoint per setup station)\n"
            "  2) All points (one waypoint per physical location, using alongLineLabel)\n\n"
            "Side naming (ALL_POINTS mode):\n"
            "  1) Dir + distance with line prefix   (L{line}_S{stn}_N50, SW120, ...)\n"
            "  2) Dir + distance without line prefix (S{stn}_N50, SW120, ...)\n"
            "  3) S{stn}_{label} (ignore line_id in name), direction only in description\n"
            "     Direction tag is added only if |offset_y| > ~0 (off-line points)."
        )

    # ---------------------------------------------------------
    # Parameters
    # ---------------------------------------------------------

    def initAlgorithm(self, config=None):

        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                "Electrodes / antennas (points UTM)",
                [QgsProcessing.TypeVectorPoint],
            )
        )

        self.addParameter(
            QgsProcessingParameterFileDestination(
                self.OUTPUT,
                "Output GPX file",
                "GPX files (*.gpx);;All files (*.*)",
            )
        )

        # Default = 1 => "All points (each location)"
        self.addParameter(
            QgsProcessingParameterEnum(
                self.MODE,
                "Mode",
                options=["Setups (one per station)", "All points (each location)"],
                defaultValue=1,
            )
        )

        self.addParameter(
            QgsProcessingParameterEnum(
                self.NAME_STYLE,
                "Side naming (ALL_POINTS mode)",
                options=[
                    "Dir + distance with line prefix (L{line}_S{stn}_N50)",
                    "Dir + distance without line prefix (S{stn}_N50)",
                    "S{stn}_{label} (ignore line_id in name)",
                ],
                defaultValue=0,
            )
        )

    # ---------------------------------------------------------
    # Utils
    # ---------------------------------------------------------

    def _bearing_deg(self, dx, dy):
        """
        Bearing (0–360°) from North, clockwise, given dx = dE, dy = dN.
        """
        ang = math.degrees(math.atan2(dx, dy))
        if ang < 0:
            ang += 360.0
        return ang

    def _dir8(self, bearing_deg):
        """
        Convert bearing in degrees to one of 8 cardinal directions.
        """
        dirs = ["N", "NE", "E", "SE", "S", "SW", "W", "NW"]
        idx = int((bearing_deg + 22.5) // 45) % 8
        return dirs[idx]

    # ---------------------------------------------------------
    # Main
    # ---------------------------------------------------------

    def processAlgorithm(self, parameters, context, feedback):

        source = self.parameterAsSource(parameters, self.INPUT, context)
        out_path = self.parameterAsFileOutput(parameters, self.OUTPUT, context)
        mode_idx = self.parameterAsEnum(parameters, self.MODE, context)
        name_style_idx = self.parameterAsEnum(parameters, self.NAME_STYLE, context)

        mode = ["SETUPS", "ALL_POINTS"][mode_idx]
        # STYLE1: L{line}_S{stn}_DirDist
        # STYLE2: S{stn}_DirDist
        # S_LABEL: S{stn}_{label}
        name_style = ["STYLE1", "STYLE2", "S_LABEL"][name_style_idx]

        if source is None:
            raise QgsProcessingException("Invalid input layer.")

        field_names = [f.name() for f in source.fields()]

        has_line_id = "line_id" in field_names
        has_station = "station" in field_names
        has_cmp = "cmp" in field_names
        has_ch_idx = "ch_idx" in field_names
        has_end = "end" in field_names
        has_label = "label" in field_names
        has_stnE = "stn_E" in field_names
        has_stnN = "stn_N" in field_names
        has_EN = ("E" in field_names) and ("N" in field_names)
        has_al = "alongLineLabel" in field_names
        has_offy = "offset_y" in field_names

        if not has_station:
            raise QgsProcessingException("Input layer must have a 'station' field.")

        if mode == "ALL_POINTS":
            missing = []
            if not has_al:
                missing.append("alongLineLabel")
            if not has_offy:
                missing.append("offset_y")
            if not has_stnE or not has_stnN:
                missing.append("stn_E/stn_N")
            if not has_EN:
                missing.append("E/N")
            if missing:
                raise QgsProcessingException(
                    "ALL_POINTS mode requires fields: " + ", ".join(missing)
                )

        # CRS transform: input -> WGS84 for GPX
        src_crs = source.sourceCrs()
        wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")
        ct_to_wgs = QgsCoordinateTransform(src_crs, wgs84, QgsProject.instance())

        feats = list(source.getFeatures())
        total = max(1, len(feats))

        # -----------------------------------------------------
        # MODE 1: Setups (one waypoint per setup station)
        # -----------------------------------------------------
        if mode == "SETUPS":

            stations = {}

            for i, f in enumerate(feats):

                if feedback.isCanceled():
                    break

                line_id = int(f["line_id"]) if has_line_id and f["line_id"] is not None else None
                stn = int(f["station"]) if f["station"] is not None else int(f.id())
                cmp_val = str(f["cmp"]) if has_cmp and f["cmp"] is not None else ""

                # Station coordinates: prefer stn_E/stn_N, else geometry
                if has_stnE and has_stnN and f["stn_E"] is not None and f["stn_N"] is not None:
                    x = float(f["stn_E"])
                    y = float(f["stn_N"])
                else:
                    geom = f.geometry()
                    if geom.isEmpty():
                        continue
                    pt = geom.asPoint()
                    x = pt.x()
                    y = pt.y()

                key = (line_id, stn)

                if key not in stations:
                    stations[key] = {
                        "line_id": line_id,
                        "station": stn,
                        "x": x,
                        "y": y,
                        "cmps": set(),
                    }

                if cmp_val:
                    stations[key]["cmps"].add(cmp_val)

                feedback.setProgress(int(100 * i / total))

            # sort by line_id, then station
            def sort_key(item):
                (line_id, stn) = item[0]
                li = line_id if line_id is not None else -999999
                return (li, stn)

            items = sorted(stations.items(), key=sort_key)

            with open(out_path, "w", encoding="utf-8") as f:
                f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
                f.write(
                    '<gpx version="1.1" creator="QGIS ZEN Tools" '
                    'xmlns="http://www.topografix.com/GPX/1/1" '
                    'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                    'xsi:schemaLocation="http://www.topografix.com/GPX/1/1 '
                    'http://www.topografix.com/GPX/1/1/gpx.xsd">\n'
                )

                for (line_id, stn), data in items:

                    if feedback.isCanceled():
                        break

                    x = data["x"]
                    y = data["y"]
                    cmps = sorted(list(data["cmps"]))

                    # Transform UTM -> WGS84
                    wgs_pt = ct_to_wgs.transform(QgsPointXY(x, y))
                    lon = wgs_pt.x()
                    lat = wgs_pt.y()

                    if line_id is not None:
                        name = f"L{line_id}_S{stn}"
                        line_str = str(line_id)
                    else:
                        name = f"S{stn}"
                        line_str = "N/A"

                    cmp_str = ",".join(cmps) if cmps else "N/A"
                    desc = f"Line={line_str}, Stn={stn}, CMP={cmp_str}"

                    f.write(f'  <wpt lat="{lat:.7f}" lon="{lon:.7f}">\n')
                    f.write(f"    <name>{name}</name>\n")
                    f.write(f"    <desc>{desc}</desc>\n")
                    f.write("  </wpt>\n")

                f.write("</gpx>\n")

            feedback.pushInfo(f"GPX (setups only) written to: {out_path}")
            return {self.OUTPUT: out_path}

        # -----------------------------------------------------
        # MODE 2: All points (one waypoint per location)
        # -----------------------------------------------------

        # Group features by physical location (E,N) – rounded to mm
        groups = {}
        for i, ftr in enumerate(feats):

            if feedback.isCanceled():
                break

            line_id = int(ftr["line_id"]) if has_line_id and ftr["line_id"] is not None else None

            # chainage along line (this is our "station" here)
            if has_al and ftr["alongLineLabel"] is not None:
                stn_al = float(ftr["alongLineLabel"])
            else:
                stn_al = float(ftr["station"]) if has_station and ftr["station"] is not None else float(ftr.id())

            # CMP, ch_idx, end, label
            cmp_val = str(ftr["cmp"]) if has_cmp and ftr["cmp"] is not None else ""
            ch_idx = int(ftr["ch_idx"]) if has_ch_idx and ftr["ch_idx"] is not None else -1
            end = str(ftr["end"]) if has_end and ftr["end"] is not None else ""
            label_str = str(ftr["label"]) if has_label and ftr["label"] is not None else ""

            # lateral offset in local Y (used to decide if we add direction tag)
            off_y = float(ftr["offset_y"]) if has_offy and ftr["offset_y"] is not None else 0.0

            # setup/station coordinates
            stn_x = float(ftr["stn_E"])
            stn_y = float(ftr["stn_N"])

            # point coordinates (actual electrode/antenna)
            x = float(ftr["E"])
            y = float(ftr["N"])

            # Key: physical location (rounded to 1 mm)
            key = (round(x, 3), round(y, 3))

            if key not in groups:
                groups[key] = {
                    "line_id": line_id,
                    "station_al": stn_al,
                    "x": x,
                    "y": y,
                    "stn_x": stn_x,
                    "stn_y": stn_y,
                    "offset_y": off_y,
                    "cmps": [],
                    "chs": [],
                    "ends": [],
                    "labels": [],
                }

            g = groups[key]
            g["cmps"].append(cmp_val)
            g["chs"].append(ch_idx)
            g["ends"].append(end)
            g["labels"].append(label_str)

            feedback.setProgress(int(100 * i / total))

        # Now write one GPX waypoint per group
        with open(out_path, "w", encoding="utf-8") as f:

            f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
            f.write(
                '<gpx version="1.1" creator="QGIS ZEN Tools" '
                'xmlns="http://www.topografix.com/GPX/1/1" '
                'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" '
                'xsi:schemaLocation="http://www.topografix.com/GPX/1/1 '
                'http://www.topografix.com/GPX/1/1/gpx.xsd">\n'
            )

            # Sort groups by line_id, then along-line coordinate
            def gkey(item):
                (_, data) = item
                li = data["line_id"] if data["line_id"] is not None else -999999
                return (li, data["station_al"])

            for (loc_key, data) in sorted(groups.items(), key=gkey):

                if feedback.isCanceled():
                    break

                line_id = data["line_id"]
                stn_al = data["station_al"]   # chainage
                stn_display = int(round(stn_al))

                x = data["x"]
                y = data["y"]
                stn_x = data["stn_x"]
                stn_y = data["stn_y"]
                off_y = data["offset_y"]

                cmps = [c for c in data["cmps"] if c]
                labels = [lb for lb in data["labels"] if lb]
                chs = [c for c in data["chs"] if c >= 0]
                ends = [e for e in data["ends"] if e]

                # Transform UTM -> WGS84
                wgs_pt = ct_to_wgs.transform(QgsPointXY(x, y))
                lon = wgs_pt.x()
                lat = wgs_pt.y()

                # Base name logic:
                # STYLE1: L{line}_S{stn}  (if line_id present) else S{stn}
                # STYLE2: S{stn}          (no line prefix even if line_id exists)
                # S_LABEL: S{stn}
                if name_style == "S_LABEL":
                    base_name = f"S{stn_display}"
                elif name_style == "STYLE2":
                    base_name = f"S{stn_display}"
                else:  # STYLE1
                    if line_id is not None:
                        base_name = f"L{line_id}_S{stn_display}"
                    else:
                        base_name = f"S{stn_display}"

                line_str = str(line_id) if line_id is not None else "N/A"

                name = base_name
                dir_label = ""
                dist2d = 0.0

                # Vector from setup to point (always computed)
                dx = x - stn_x  # East
                dy = y - stn_y  # North
                dist2d = math.hypot(dx, dy)

                eps_dist = 0.1  # m
                eps_y = 0.1     # m threshold for "off-line"

                # Only add a compass + distance tag if the point is OFF the line
                if abs(off_y) > eps_y and dist2d > eps_dist:
                    bearing = self._bearing_deg(dx, dy)
                    card = self._dir8(bearing)
                    dist_m = int(round(dist2d))
                    dir_label = f"{card}{dist_m}"

                    if name_style in ("STYLE1", "STYLE2"):
                        name = f"{base_name}_{dir_label}"

                # For S_LABEL style, append primary label to name
                if name_style == "S_LABEL":
                    primary_label = labels[0] if labels else ""
                    if primary_label:
                        name = f"{base_name}_{primary_label}"

                # Build description with all labels/CMPs + direction + offset_y
                desc_parts = [
                    f"Line={line_str}",
                    f"Stn={stn_display}",
                    f"Chainage={stn_al:.3f} m",
                ]
                if cmps:
                    desc_parts.append("CMPs=" + ",".join(sorted(set(cmps))))
                if labels:
                    desc_parts.append("Labels=" + ",".join(labels))
                if chs:
                    desc_parts.append("Ch=" + ",".join(str(c) for c in sorted(set(chs))))
                if ends:
                    desc_parts.append("Ends=" + ",".join(sorted(set(ends))))
                if dir_label:
                    desc_parts.append(f"Dir={dir_label}")
                    desc_parts.append(f"d2D={dist2d:.3f} m")
                desc_parts.append(f"offset_y={off_y:.3f} m")

                desc = ", ".join(desc_parts)

                f.write(f'  <wpt lat="{lat:.7f}" lon="{lon:.7f}">\n')
                f.write(f"    <name>{name}</name>\n")
                f.write(f"    <desc>{desc}</desc>\n")
                f.write("  </wpt>\n")

            f.write("</gpx>\n")

        feedback.pushInfo(f"GPX (all points, one per location) written to: {out_path}")
        return {self.OUTPUT: out_path}