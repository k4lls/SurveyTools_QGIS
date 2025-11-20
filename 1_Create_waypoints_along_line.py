from qgis.PyQt.QtCore import QVariant
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFeatureSink,
    QgsProcessingParameterNumber,
    QgsFeatureSink,
    QgsFields,
    QgsField,
    QgsFeature,
    QgsGeometry,
    QgsWkbTypes,
    QgsPointXY,
    QgsProcessingException,
    QgsPalLayerSettings,
    QgsTextFormat,
    QgsTextBufferSettings,
    QgsVectorLayerSimpleLabeling,
)
import math


class BuildZenSetupsSimple(QgsProcessingAlgorithm):

    INPUT = "INPUT"
    NUM_SETUPS = "NUM_SETUPS"
    SPACING = "SPACING"
    OFFSET0 = "OFFSET0"
    LINE_AZ = "LINE_AZ"
    STN_START = "STN_START"
    STN_STEP = "STN_STEP"
    OUTPUT = "OUTPUT"

    def createInstance(self):
        return BuildZenSetupsSimple()

    def name(self):
        return "build_waypoints_along_line"

    def displayName(self):
        return "1 - Create Waypoints along a line"

    def group(self):
        return "Survey Tools"

    def groupId(self):
        return "survey_tools"

    def shortHelpString(self):
        return (
            "Generate N ZEN setups starting from each waypoint.\n\n"
            "Each input point defines one line of setups.\n\n"
            "Parameters:\n"
            "- Number of setups (N)\n"
            "- Spacing between setups (m)\n"
            "- Offset of first setup from waypoint (m)\n"
            "- Azimuth (deg from North clockwise)\n"
            "- First STN and STN step\n\n"
            "Output fields:\n"
            "  line_id : id of waypoint source (one line per waypoint)\n"
            "  STN     : station id (used for labels)\n"
            "  stn_s   : chainage from waypoint (m)\n"
            "  az      : azimuth used (deg)\n"
        )

    def initAlgorithm(self, config=None):

        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT,
                "Start waypoints (points, projected CRS)",
                [QgsProcessing.TypeVectorPoint],
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.NUM_SETUPS,
                "Number of setups",
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=10,
                minValue=1,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.SPACING,
                "Spacing between setups (m)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=100.0,
                minValue=0.0,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.OFFSET0,
                "First setup offset from waypoint (m)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0,
                minValue=0.0,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.LINE_AZ,
                "Line azimuth (deg from North clockwise)",
                type=QgsProcessingParameterNumber.Double,
                defaultValue=0.0,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.STN_START,
                "First STN",
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=0,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.STN_STEP,
                "STN step",
                type=QgsProcessingParameterNumber.Integer,
                defaultValue=1,
            )
        )

        self.addParameter(
            QgsProcessingParameterFeatureSink(
                self.OUTPUT,
                "Generated setups",
                QgsProcessing.TypeVectorPoint,
            )
        )

    def processAlgorithm(self, parameters, context, feedback):

        source = self.parameterAsSource(parameters, self.INPUT, context)
        num_setups = self.parameterAsInt(parameters, self.NUM_SETUPS, context)
        spacing = self.parameterAsDouble(parameters, self.SPACING, context)
        offset0 = self.parameterAsDouble(parameters, self.OFFSET0, context)
        az = self.parameterAsDouble(parameters, self.LINE_AZ, context)
        stn_start_param = self.parameterAsInt(parameters, self.STN_START, context)
        stn_step = self.parameterAsInt(parameters, self.STN_STEP, context)

        if source is None:
            raise QgsProcessingException("Invalid input waypoint layer.")
        if num_setups < 1:
            raise QgsProcessingException("Number of setups must be >= 1.")
        if spacing < 0:
            raise QgsProcessingException("Spacing must be >= 0.")
        if offset0 < 0:
            raise QgsProcessingException("First setup offset must be >= 0.")

        crs = source.sourceCrs()
        if crs.isGeographic():
            raise QgsProcessingException("Reproject input layer to UTM (projected CRS) first.")

        feats = list(source.getFeatures())
        if not feats:
            raise QgsProcessingException("Waypoint layer is empty.")

        # Direction vector from azimuth (North=0°, East=90°)
        az_rad = math.radians(az)
        dirE = math.sin(az_rad)
        dirN = math.cos(az_rad)

        # Output fields
        fields = QgsFields()
        fields.append(QgsField("line_id", QVariant.Int))
        fields.append(QgsField("STN", QVariant.Int))
        fields.append(QgsField("stn_s", QVariant.Double))
        fields.append(QgsField("az", QVariant.Double))

        sink, out_id = self.parameterAsSink(
            parameters,
            self.OUTPUT,
            context,
            fields,
            QgsWkbTypes.Point,
            crs,
        )

        # Generate setups
        for feat in feats:
            geom = feat.geometry()
            if geom.isEmpty():
                continue

            line_id = int(feat.id())
            pt = geom.asPoint()
            x0, y0 = pt.x(), pt.y()

            current_stn = stn_start_param

            for k in range(num_setups):

                s = offset0 + k * spacing
                x = x0 + s * dirE
                y = y0 + s * dirN

                f_out = QgsFeature(fields)
                f_out.setGeometry(QgsGeometry.fromPointXY(QgsPointXY(x, y)))

                f_out["line_id"] = line_id
                f_out["STN"] = current_stn
                f_out["stn_s"] = float(s)
                f_out["az"] = float(az)

                sink.addFeature(f_out, QgsFeatureSink.FastInsert)

                current_stn += stn_step

        # --- Enable STN labeling on the output layer ---
        out_layer = context.getMapLayer(out_id)
        if out_layer is not None:
            settings = QgsPalLayerSettings()
            settings.fieldName = "STN"
            settings.isExpression = False
            settings.enabled = True

            text_format = QgsTextFormat()
            text_format.setSize(9)          # size in points
            buffer = QgsTextBufferSettings()
            buffer.setEnabled(True)
            buffer.setSize(1)
            text_format.setBuffer(buffer)
            settings.setFormat(text_format)

            labeling = QgsVectorLayerSimpleLabeling(settings)
            out_layer.setLabeling(labeling)
            out_layer.setLabelsEnabled(True)
            out_layer.triggerRepaint()

        return {self.OUTPUT: out_id}