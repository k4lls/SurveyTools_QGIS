from qgis.PyQt.QtCore import QVariant
from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterFileDestination,
    QgsProcessingParameterEnum,
    QgsProject,
    QgsCoordinateReferenceSystem,
    QgsCoordinateTransform,
    QgsPointXY,
    QgsProcessingException,
)
import math


class ExportRXCFromPoints(QgsProcessingAlgorithm):

    INPUT = "INPUT"
    OUTPUT = "OUTPUT"
    MODE = "MODE"

    def createInstance(self):
        return ExportRXCFromPoints()

    def name(self):
        return "export_rxc"

    def displayName(self):
        return "4 - Export RXC (MT/IP)"

    def group(self):
        return "Survey Tools"

    def groupId(self):
        return "survey_tools"

    def shortHelpString(self):
        return (
            "Exports an RXC file from ZEN points.\n\n"
            "Required fields:\n"
            "  line_id, station, cmp, ch_idx, end,\n"
            "  stn_E, stn_N, E, N,\n"
            "  stn_s, alongLineLabel,\n"
            "  ch_az, ch_incl, ant_num.\n\n"
            "MT mode:\n"
            "  - Includes Ex, Ey, Hx, Hy, Hz (excludes Tx/Ref)\n"
            "  - Ch.Stn for Ex = rounded midpoint of electrodes\n"
            "  - Ch.Stn for Ey = rounded alongLineLabel of end1\n"
            "  - Ch.Stn for H-fields = station\n"
            "  - Outputs Ant#, Azm, Incl\n"
            "  - Does NOT output SX/SY/SZ\n\n"
            "IP mode (locked, do not modify):\n"
            "  - Only Ex and Tx/Ref\n"
            "  - Rx.Stn = station for Ex, 999 for Tx/Ref\n"
            "  - Ch.Stn = Name1 = negative electrode along-line value\n"
            "  - Name2 = positive electrode along-line value\n"
            "  - Zen.Chn = ch_idx for Ex, 1 for Tx/Ref\n"
            "  - Outputs SX/SY/SZ but NO Ant# and NO Incl."
        )

    # ---------------------------------------------------------
    # Helper functions
    # ---------------------------------------------------------

    def compute_azimuth(self, E1, N1, E2, N2):
        """Azimuth (degrees) from point 1 to point 2, north clockwise."""
        dE = E2 - E1
        dN = N2 - N1
        ang = math.degrees(math.atan2(dE, dN))
        return ang + 360.0 if ang < 0 else ang

    def fmt(self, val, nd=3):
        """Format float."""
        return f"{val:.{nd}f}"

    # ---------------------------------------------------------
    # Processing parameters
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
                "Output RXC file",
                "RXC files (*.rxc);;All files (*.*)",
            )
        )

        self.addParameter(
            QgsProcessingParameterEnum(
                self.MODE,
                "Mode",
                options=["MT", "IP"],
                defaultValue=0,
            )
        )

    # ---------------------------------------------------------
    # Main algorithm
    # ---------------------------------------------------------

    def processAlgorithm(self, parameters, context, feedback):

        source = self.parameterAsSource(parameters, self.INPUT, context)
        out_path = self.parameterAsFileOutput(parameters, self.OUTPUT, context)
        mode = ["MT", "IP"][self.parameterAsEnum(parameters, self.MODE, context)]

        if source is None:
            raise QgsProcessingException("Invalid input layer.")

        # Required fields
        required = [
            "line_id",
            "station", "cmp", "ch_idx", "end",
            "stn_E", "stn_N", "E", "N",
            "stn_s", "alongLineLabel",
            "ch_az", "ch_incl", "ant_num",
        ]
        fields = [f.name() for f in source.fields()]
        missing = [f for f in required if f not in fields]
        if missing:
            raise QgsProcessingException("Missing fields: " + ", ".join(missing))

        # Coordinate transform to GPS
        src_crs = source.sourceCrs()
        wgs84 = QgsCoordinateReferenceSystem("EPSG:4326")
        ct_to_wgs = QgsCoordinateTransform(src_crs, wgs84, QgsProject.instance())

        # ---------------------------------------------------------
        # Group features by channel (line_id, station, cmp, ch_idx)
        # ---------------------------------------------------------

        feats = list(source.getFeatures())
        groups = {}

        for f in feats:

            cmp_ = str(f["cmp"])
            if cmp_ == "Off":
                continue

            line_id = int(f["line_id"])
            stn     = int(f["station"])
            idx     = int(f["ch_idx"])
            end     = str(f["end"])

            key = (line_id, stn, cmp_, idx)

            if key not in groups:
                groups[key] = {
                    "line_id": line_id,
                    "station": stn,
                    "cmp": cmp_,
                    "idx": idx,
                    "stn_E": float(f["stn_E"]),
                    "stn_N": float(f["stn_N"]),
                    "stn_s": float(f["stn_s"]),
                    "along1": None,
                    "along2": None,
                    "p1": None,
                    "p2": None,
                    "ch_az": float(f["ch_az"]) if f["ch_az"] is not None else 0.0,
                    "ch_incl": float(f["ch_incl"]) if f["ch_incl"] is not None else 0.0,
                    "ant_num": f["ant_num"],
                }

            E = float(f["E"])
            N = float(f["N"])
            along = f["alongLineLabel"]

            if end == "1":
                groups[key]["p1"] = (E, N)
                groups[key]["along1"] = along

            elif end == "2":
                groups[key]["p2"] = (E, N)
                groups[key]["along2"] = along

        # Convert dict to list
        items = sorted(groups.items())

        # ---------------------------------------------------------
        # MODE IP (locked, do not modify)
        # ---------------------------------------------------------

        if mode == "IP":

            # Keep only Ex and Tx/Ref
            items = [
                (k, d)
                for (k, d) in items
                if d["cmp"] in ("Ex", "Tx", "Ref")
            ]

            # Sort by Line.Name → Ex first → Tx/Ref → stn_s → ch_idx
            items = sorted(
                items,
                key=lambda x: (
                    x[1]["line_id"],
                    0 if x[1]["cmp"] == "Ex" else 1,
                    x[1]["stn_s"],
                    x[1]["idx"],
                )
            )

            header = (
                "Rx.Stn,Ch.Stn,Zen.Chn,Ch.Cmp,Line.Name,Azm,"
                "Name1,East1,North1,"
                "Name2,East2,North2,"
                "Name0,East0,North0,"
                "GPS.Lat,GPS.Lon,"
                "SX1,SY1,SZ1,SX2,SY2,SZ2\n"
            )

            with open(out_path, "w") as f:
                f.write(header)

                for (line_id, stn, cmp_, idx), d in items:

                    stn_E = d["stn_E"]
                    stn_N = d["stn_N"]
                    along1 = d["along1"]
                    along2 = d["along2"]
                    p1 = d["p1"]
                    p2 = d["p2"]

                    # GPS position
                    gps_xy = ct_to_wgs.transform(QgsPointXY(stn_E, stn_N))
                    gps_lon = gps_xy.x()
                    gps_lat = gps_xy.y()

                    # Rx.Stn / Ch.Stn
                    rx_stn = 999 if cmp_ in ("Tx", "Ref") else stn
                    ch_stn = str(along1) if along1 is not None else ""

                    # Zen.Chn
                    zen_chn = 1 if cmp_ in ("Tx", "Ref") else idx

                    # Base receiver
                    name0 = str(stn)
                    east0 = self.fmt(stn_E)
                    north0 = self.fmt(stn_N)

                    name1 = name2 = ""
                    east1 = north1 = east2 = north2 = ""
                    azm = 0.0

                    # Dipoles Ex/Tx/Ref
                    if p1 is not None and p2 is not None:
                        E1, N1 = p1
                        E2, N2 = p2

                        name1 = str(along1)
                        name2 = str(along2)

                        east1 = self.fmt(E1)
                        north1 = self.fmt(N1)
                        east2 = self.fmt(E2)
                        north2 = self.fmt(N2)

                        azm = self.compute_azimuth(E1, N1, E2, N2)

                    # SX/SY/SZ
                    sx1 = name1
                    sx2 = name2

                    row = [
                        str(rx_stn),
                        str(ch_stn),
                        str(zen_chn),
                        cmp_,
                        str(line_id),
                        self.fmt(azm, 2),
                        name1, east1, north1,
                        name2, east2, north2,
                        name0, east0, north0,
                        self.fmt(gps_lat, 7),
                        self.fmt(gps_lon, 7),
                        sx1, "0", "0", sx2, "0", "0",
                    ]
                    f.write(",".join(row) + "\n")

            return {self.OUTPUT: out_path}

        # ---------------------------------------------------------
        # MODE MT (rounding logic)
        # ---------------------------------------------------------

        # Keep all but Tx/Ref
        items = [
            (k, d)
            for (k, d) in items
            if d["cmp"] not in ("Tx", "Ref")
        ]

        # Sort by line → station → cmp → idx
        items = sorted(
            items,
            key=lambda x: (
                x[1]["line_id"],
                x[1]["station"],
                x[1]["cmp"],
                x[1]["idx"],
            )
        )

        header = (
            "Rx.Stn,Ch.Stn,Zen.Chn,Ch.Cmp,Line.Name,Ant#,Azm,Incl,"
            "Name1,East1,North1,"
            "Name2,East2,North2,"
            "Name0,East0,North0,"
            "GPS.Lat,GPS.Lon\n"
        )

        with open(out_path, "w") as f:
            f.write(header)

            for (line_id, stn, cmp_, idx), d in items:

                stn_E = d["stn_E"]
                stn_N = d["stn_N"]
                along1 = d["along1"]
                along2 = d["along2"]
                p1 = d["p1"]
                p2 = d["p2"]
                ch_az = d["ch_az"]
                ch_incl = d["ch_incl"]
                ant_raw = d["ant_num"]

                # GPS
                gps_xy = ct_to_wgs.transform(QgsPointXY(stn_E, stn_N))
                gps_lon = gps_xy.x()
                gps_lat = gps_xy.y()

                # Receiver station
                rx_stn = stn

                # -------------------------
                # Channel station (Ch.Stn)
                # -------------------------

                if cmp_ == "Ex" and along1 is not None and along2 is not None:
                    # Midpoint → rounded integer
                    try:
                        a1 = float(along1)
                        a2 = float(along2)
                        ch_stn = str(int(round((a1 + a2) / 2)))
                    except Exception:
                        ch_stn = str(stn)

                elif cmp_.startswith("E") and along1 is not None:
                    # Ey etc → round along1
                    try:
                        ch_stn = str(int(round(float(along1))))
                    except Exception:
                        ch_stn = str(stn)

                else:
                    # H-fields use station
                    ch_stn = str(stn)

                # Zen.Chn
                zen_chn = idx

                # Ant#
                if ant_raw is None or (isinstance(ant_raw, QVariant) and ant_raw.isNull()):
                    ant_num = "####"
                else:
                    try:
                        ant_num = str(int(float(ant_raw)))
                    except Exception:
                        ant_num = "####"

                # Receiver base
                name0 = str(stn)
                east0 = self.fmt(stn_E)
                north0 = self.fmt(stn_N)

                name1 = name2 = ""
                east1 = north1 = east2 = north2 = ""
                azm = 0.0
                incl = 0.0

                # E-field (Ex/Ey)
                if p1 is not None and p2 is not None and cmp_.startswith("E"):
                    E1, N1 = p1
                    E2, N2 = p2

                    # Rounded along-line names
                    if along1 is not None:
                        name1 = str(int(round(float(along1))))
                    if along2 is not None:
                        name2 = str(int(round(float(along2))))

                    east1 = self.fmt(E1)
                    north1 = self.fmt(N1)
                    east2 = self.fmt(E2)
                    north2 = self.fmt(N2)

                    azm = self.compute_azimuth(E1, N1, E2, N2)
                    incl = 0.0

                # H-field
                elif cmp_.startswith("H"):
                    name1 = name2 = str(stn)
                    east1 = east2 = self.fmt(stn_E)
                    north1 = north2 = self.fmt(stn_N)
                    azm = ch_az
                    incl = ch_incl

                row = [
                    str(rx_stn),
                    str(ch_stn),
                    str(zen_chn),
                    cmp_,
                    str(line_id),
                    ant_num,
                    self.fmt(azm, 2),
                    self.fmt(incl, 2),
                    name1, east1, north1,
                    name2, east2, north2,
                    name0, east0, north0,
                    self.fmt(gps_lat, 7),
                    self.fmt(gps_lon, 7),
                ]
                f.write(",".join(row) + "\n")

        return {self.OUTPUT: out_path}