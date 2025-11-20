# -*- coding: utf-8 -*-
from qgis.core import (
    QgsProcessing, QgsProcessingAlgorithm,
    QgsProcessingParameterFeatureSource, QgsProcessingParameterField,
    QgsProcessingParameterRasterLayer, QgsProcessingParameterFileDestination,
    QgsProcessingParameterBoolean,
    QgsProcessingException, QgsFeatureRequest, QgsWkbTypes,
    QgsCoordinateTransform, QgsCoordinateReferenceSystem, QgsProject,
    QgsGeometry, QgsPointXY
)
import json, time
from urllib.request import urlopen
from urllib.parse import urlencode

class ExportStationFileAlgorithm(QgsProcessingAlgorithm):
    PARAM_POINTS='POINTS'; PARAM_NAME_FIELD='NAME_FIELD'; PARAM_ELEV_FIELD='ELEV_FIELD'
    PARAM_DEM='DEM'; PARAM_OUTPUT='OUTPUT'; PARAM_USE_API='USE_API'

    def createInstance(self): return ExportStationFileAlgorithm()
    def name(self): return 'export_station_file'
    def displayName(self): return '5 - Export waypoints to Station file (.stn/.wpt)'
    def group(self): return 'Survey Tools'
    def groupId(self): return 'survey_tools'
    def shortHelpString(self):
        return ("Exports points to a station file (UTM East/North/Elev).\n"
                "Elevation priority: attribute field > DEM sample > (optional) Open-Meteo API > 0.0")

    def initAlgorithm(self, config=None):
        self.addParameter(QgsProcessingParameterFeatureSource(self.PARAM_POINTS,'Waypoint layer (points)',[QgsProcessing.TypeVectorPoint]))
        self.addParameter(QgsProcessingParameterField(self.PARAM_NAME_FIELD,'Name field (optional — uses feature id if empty)',
                                                     parentLayerParameterName=self.PARAM_POINTS, optional=True))
        self.addParameter(QgsProcessingParameterField(self.PARAM_ELEV_FIELD,'Elevation field (optional if DEM/API used)',
                                                     parentLayerParameterName=self.PARAM_POINTS, type=QgsProcessingParameterField.Numeric, optional=True))
        self.addParameter(QgsProcessingParameterRasterLayer(self.PARAM_DEM,'DEM for elevation sampling (optional)', optional=True))
        self.addParameter(QgsProcessingParameterBoolean(self.PARAM_USE_API,'Use Open-Meteo API if no Elevation field/DEM', defaultValue=True))
        self.addParameter(QgsProcessingParameterFileDestination(self.PARAM_OUTPUT,'Station file',fileFilter='Station/WPT (*.stn *.wpt *.txt *.csv)'))

    def _is_utm_epsg(self, epsg): return 32601<=epsg<=32660 or 32701<=epsg<=32760
    def _pick_utm(self, src_crs, src):
        wgs84 = QgsCoordinateReferenceSystem('EPSG:4326')
        tf = QgsCoordinateTransform(src_crs,wgs84,QgsProject.instance()) if src_crs.isValid() and src_crs!=wgs84 else None
        c = src.sourceExtent().center();  c = tf.transform(c) if tf else c
        lon, lat = c.x(), c.y()
        zone = int((lon+180)/6)+1; north = lat>=0
        epsg = (32600 if north else 32700)+zone
        return QgsCoordinateReferenceSystem.fromEpsgId(epsg), epsg, zone,('N' if north else 'S')

    def _fetch_open_meteo(self, lats, lons, feedback, retries=3):
        qs = {"latitude": ",".join(f"{x:.8f}" for x in lats),
              "longitude": ",".join(f"{x:.8f}" for x in lons)}
        url = "https://api.open-meteo.com/v1/elevation?" + urlencode(qs)
        last = None
        for _ in range(retries):
            try:
                with urlopen(url, timeout=30) as r:
                    data = json.loads(r.read().decode("utf-8"))
                    return data.get("elevation", [])
            except Exception as e:
                last = e; time.sleep(1)
        raise QgsProcessingException(f"Open-Meteo fetch failed: {last}")

    def processAlgorithm(self, p, context, feedback):
        src = self.parameterAsSource(p, self.PARAM_POINTS, context)
        if src is None: raise QgsProcessingException('Invalid point layer.')
        if QgsWkbTypes.geometryType(src.wkbType()) != QgsWkbTypes.PointGeometry:
            raise QgsProcessingException('Input must be a point layer.')

        name_field = self.parameterAsString(p, self.PARAM_NAME_FIELD, context)
        elev_field = self.parameterAsString(p, self.PARAM_ELEV_FIELD, context)
        dem = self.parameterAsRasterLayer(p, self.PARAM_DEM, context)
        use_api = self.parameterAsBool(p, self.PARAM_USE_API, context)
        out_path = self.parameterAsFileOutput(p, self.PARAM_OUTPUT, context)

        src_crs = src.sourceCrs(); epsg_src = src_crs.postgisSrid()
        if src_crs.isValid() and self._is_utm_epsg(epsg_src):
            tgt_crs, epsg_tgt = src_crs, epsg_src
            zone = (epsg_src-32600) if epsg_src<32700 else (epsg_src-32700)
            hemi = 'N' if epsg_src<32700 else 'S'
        else:
            tgt_crs, epsg_tgt, zone, hemi = self._pick_utm(src_crs, src)

        tf_to_tgt = QgsCoordinateTransform(src_crs, tgt_crs, QgsProject.instance())

        # Check DEM is numeric
        has_dem = False
        if dem:
            try:
                test_tf = QgsCoordinateTransform(src_crs, dem.crs(), QgsProject.instance())
                test_pt = test_tf.transform(src.sourceExtent().center())
                idres = dem.dataProvider().identify(test_pt, 0)
                has_dem = idres.isValid()
                if not has_dem:
                    feedback.reportError("Selected DEM does not return numeric values (XYZ/WMS tiles won’t work).")
            except Exception:
                has_dem = False

        # Header
        with open(out_path,'w',encoding='utf-8') as f:
            f.write("$WPT.Datum=WGS84\n")
            f.write("$WPT.Proj=UTM\n")
            f.write(f"$WPT.UTMZone={zone}{hemi}\n")
            f.write(f"$WPT.EPSG={epsg_tgt}\n")
            f.write("STATION, EASTING, NORTHING, ELEVATION\n")

            feats = list(src.getFeatures(QgsFeatureRequest()))
            total=len(feats); step=max(1,total//20) if total else 1

            # Precompute for API if needed
            need_api = (not elev_field or elev_field not in src.fields().names()) and not has_dem and use_api
            if need_api:
                # Build WGS84 arrays
                wgs84 = QgsCoordinateReferenceSystem('EPSG:4326')
                tf_to_wgs = QgsCoordinateTransform(src_crs, wgs84, QgsProject.instance())
                lats=[]; lons=[]
                for feat in feats:
                    pt = feat.geometry().asPoint() if not feat.geometry().isMultipart() else feat.geometry().asMultiPoint()[0]
                    wpt = tf_to_wgs.transform(QgsPointXY(pt))
                    lons.append(float(wpt.x())); lats.append(float(wpt.y()))
                # Fetch in chunks of 100
                elev_api=[]
                for i in range(0, total, 100):
                    if feedback.isCanceled(): break
                    sub = self._fetch_open_meteo(lats[i:i+100], lons[i:i+100], feedback)
                    if len(sub)!=(min(100,total-i)):
                        raise QgsProcessingException("Open-Meteo returned a different count than requested.")
                    elev_api.extend(sub)

            for i, feat in enumerate(feats):
                if feedback.isCanceled(): break
                if i % step == 0 and total: feedback.setProgress(int(100*i/total))

                geom = feat.geometry()
                if not geom or geom.isEmpty(): continue
                pt = geom.asPoint() if not geom.isMultipart() else geom.asMultiPoint()[0]

                name = (str(feat[name_field]) if name_field and name_field in src.fields().names() and feat[name_field] is not None else str(feat.id()))
                pt_utm = tf_to_tgt.transform(QgsPointXY(pt))
                east, north = float(pt_utm.x()), float(pt_utm.y())

                elev = None
                # 1) field
                if elev_field and elev_field in src.fields().names():
                    try:
                        v = feat[elev_field]
                        elev = float(v) if v is not None else None
                    except Exception:
                        elev = None
                # 2) DEM
                if elev is None and has_dem:
                    tf_to_dem = QgsCoordinateTransform(src_crs, dem.crs(), QgsProject.instance())
                    ptd = tf_to_dem.transform(QgsPointXY(pt))
                    idres = dem.dataProvider().identify(ptd, 0)
                    if idres.isValid() and idres.results():
                        elev = float(list(idres.results().values())[0])
                # 3) API
                if elev is None and use_api:
                    elev = float(elev_api[i]) if need_api else None
                # 4) default
                if elev is None: elev = 0.0

                f.write(f"{name}, {east:.3f}, {north:.3f}, {elev:.2f}\n")

        feedback.pushInfo(f"Written: {out_path}")
        return {self.PARAM_OUTPUT: out_path}

def classFactory(): return ExportStationFileAlgorithm()
