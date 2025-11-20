# Survey Setup Workflow — Documentation

This document describes the complete workflow for generating **ZEN setups**, **GPX navigation files**, **RXC acquisition files**, and **STN station files**, using a combination of:

- **ZEN Template Builder** (web)
- **QGIS + custom tools** (Python Processing scripts)

The objective is to ensure a **fully reproducible**, field-friendly workflow for survey planning and data acquisition.

# Workflow Summary

1. Draw line start points in QGIS  
2. Create Setup Waypoints along lines
3. Build template in web builder → save `.stt`  
4. Run **Apply Template to Setup Waypoints** → electrodes/antennas  
5. Export:
   - **GPX** for field navigation  
   - **RXC (MT/IP)** for ZEN  
6. Export **STN** if needed for inversion

---

# 1. ZEN Template Builder (Web Interface)
**URL:** https://zonge-international.github.io/SetupTemplateBuilder/

This is where the geometry of a ZEN receiver setup is defined.

## What the builder does
- Choose number of channels (1–16)
- Assign CMP (Ex, Ey, Hx, Hy, Hz, Tx)
- Electrode offsets (+/- X and Y)
- Automatic handling of Z+ (Up/Down)
- Live graphical preview
- Export a Type-2 **`.stt`** template containing:
  ```
  <TEMPLATE>
      TEMPLATE.NAME=...
      RX.XAZIMUTH=90
      RX.ZPOSITIVE=UP or DOWN
      <CH>...</CH>
  </TEMPLATE>
  ```

This `.stt` file is used by QGIS to place electrodes and antennas in UTM.

---

# 2. Script 1 — Create Setup Waypoints Along a Line
**File:** `1_Create_waypoints_along_line.py`

Generates evenly spaced setup waypoints along a line.

### Input
- One point per line
- Number of setups
- Line azimuth
- Setup spacing
- Offset for first point
- First STN and increment

### Output layer
- `line_id`
- `STN`
- `stn_s` (chainage)
- `az` (used azimuth)

---

# 3. Script 2 — Apply Template to Waypoints
**File:** `2_Apply_Template_to_Waypoints.py`

Generate a new template of real electrode and antenna coordinates in **UTM** at each setup waypoints.

### Output layers
1. **Electrode/Antenna points**
   - `cmp` (Ex/Ey/Hx/Hy/Hz/Tx)
   - `label` (Ex1-, Ex1+, etc.)
   - `E`, `N` (actual electrode UTM)
   - `stn_E`, `stn_N` (setup location)
   - `offset_y`, `dist`, `alongLineLabel`
   - `ch_az`, `ch_incl`, `ant_num`

2. **Dipole vectors**

This is the geometric backbone of the workflow.

---

# 4. Script 3 — Export GPX
**File:** `3_Export_GPX.py`

Generates a GPX file for navigation in the field.

### Modes
- **Setups only**
- **All points (default)** – one waypoint per physical location

### Naming
- **On-line:**  
  `L0_S300`
- **Off-line:**  
  `L0_S300_N50`, `L0_S300_SW120`, etc.

Styles:
1. With line prefix  
2. Without line prefix  
3. `S{stn}_{label}`

---

# 5. Script 4 — Export RXC (MT/IP)
**File:** `4_Export_RXC.py`

Generates ZEN-compatible RXC files.

## IP mode
- Ex and Tx/Ref only  
- `Rx.Stn = station`, Tx = 999  
- `Ch.Stn = alongLineLabel`  
- Possibility to inlcude Ant#  
- SX/SY/SZ geometry included

## MT mode
- Ex, Ey, Hx, Hy, Hz  
- No Tx/Ref  
- Ex → midpoint  
- Ey → station of negative  
- H-field includes `Incl`, `Azm`, `Ant#`

---

# 6. Script 5 — Export STN File
**File:** `5_Export_StnFile.py`

Generates a `.stn` file for modeling/inversion.

### Elevation sources:
1. Elevation field  
2. DEM sampling  
3. Open-Meteo API  
4. Else → `0.0`

---

# 7. Benefits

- Consistent geometry from planning → acquisition → inversion  
- Fully reproducible  
- Fast and minimize humain errors
- Field-friendly naming (`N50`, `SE120`, etc.)  

