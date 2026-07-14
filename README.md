# ESN Simulation Code – Abaqus and MATLAB Scripts for Electrospun Fiber Network Generation, Post-Processing, and Kinematic Analysis

## Overview
The article "A continuum mechanics based framework to quantify local kinematics in electrospun networks for mechanobiological investigations"  introduces a numerical and continuum-mechanics-based workflow to quantify local surface kinematics in electrospun fiber networks (ESN). The repository provides the simulation and post-processing scripts used to generate the numerical ESN models, extract nodal coordinate data, render virtual ESN image sequences, and compute the local kinematic quantities reported in the manuscript.

This repository can be cited useing DOI XXXXXXX.

### Authors

- **Joël Zimmerli**
- **Jonas Hofmann**
- **Barbara Röhrnbauer**

Affiliation:

**ZHAW School of Engineering, IMES Institute of Mechanical Systems, Technikumstrasse 71, 8401 Winterthur, Switzerland**

This repository contains the simulation, post-processing, and analysis code used in the study, which introduces a continuum-mechanics-based framework for quantifying local kinematics in electrospun fiber networks (ESN):

1. **`A_ESN_builder.py`** (Abaqus/Python) – generates parametric representative volume elements (RVEs) of electrospun fiber networks, meshes them, applies crosslinks and periodic boundary conditions (PBC), and sets up an Abaqus job for finite element (FE) analysis.
2. **`B_Script_Coord_20Frames.py`** (Abaqus/Python) – post-processes the resulting Abaqus output database (ODB) files by extracting nodal coordinates over all analysis frames, providing the raw kinematic data required for strain-field analysis of the simulated ESN.
3. **`C_Import_COORD.m`** (MATLAB) – imports the nodal coordinate text file exported by `B_Script_Coord_20Frames.py` into a MATLAB table `COORD` and saves it as a `.mat` file, ready for use by `D_Kinematic_Fingerprint_Analysis.mlx` and `E_RenderESN.mlx`.
4. **`D_Kinematic_Fingerprint_Analysis.mlx`** (MATLAB) – the main evaluation script: computes local deformation gradients, Green-Lagrange strains, local rotation, and the deformation type exponent from the facet centers of the ESN, and assembles the resulting "kinematic fingerprint" of the network.
5. **`E_RenderESN.mlx`** (MATLAB) – renders the ESN fiber network as grayscale line-drawing images for each analyzed frame, based on the exported nodal coordinates.

Scripts 1 and 2 must be run within an **Abaqus/CAE Python (kernel) environment**; they rely exclusively on Abaqus-internal modules and are not intended to run as standalone Python scripts. Scripts 3–5 are run in **MATLAB**. `C_Import_COORD.m` must be run first, converting the `_COORD.txt` output of script 2 into the `.mat` file used as input by `D_Kinematic_Fingerprint_Analysis.mlx` and `E_RenderESN.mlx` (see below).

### Sample data / quick start

To allow the analysis and rendering scripts to be run without first executing the full Abaqus simulation pipeline, this repository additionally includes one example dataset from the manuscript:

- `R_200x6xlsm5_plast_viso0_9_UniAx_0_COORD.mat` – example `.mat` file (table `COORD`), already imported via `C_Import_COORD.m`
- `2024-05-26_R_200x6xlsm5_plast_viso0_9_UniAx_0_fiberstackorder.txt` – the matching fiber stacking order file

Using these two files, `D_Kinematic_Fingerprint_Analysis.mlx` and `E_RenderESN.mlx` (Sections 4 and 5) can be run directly, without needing Abaqus or having to run `A_ESN_builder.py` / `B_Script_Coord_20Frames.py` / `C_Import_COORD.m` beforehand. Simply point `coordinateFile` / `load(...)` and `fiberStackOrderFile` / `fibre_stack_filename` in these two scripts to the corresponding filenames above.

Also included is **`Matlab_009_createPlots_FEM_Matlab_Q4_reference.m`**, a plotting helper function required by `D_Kinematic_Fingerprint_Analysis.mlx` to generate the averaged kinematic field plots (see Section 4). It must be located on the MATLAB path (e.g. in the same folder as `D_Kinematic_Fingerprint_Analysis.mlx`) whenever `plotKinematicFields = true`.

---

## 1. `A_ESN_builder.py` – ESN Network Generation

### Purpose

Generates fully parametric ESN RVE models, ready for FE simulation, including:

- Randomized fiber geometry 
- Beam-element meshing (B22 elements, circular cross-section)
- Elastic-plastic fiber material definition
- Crosslinks between neighboring fibers (tie constraints)
- Periodic boundary conditions (PBC) emulating an infinitely extended network
- Load-case-specific boundary conditions and a submission-ready Abaqus job

### Requirements

- **Abaqus/CAE** with the Abaqus Scripting Interface (Python 2 / `mbcs` encoding), Version 2022,
- Only Abaqus-internal modules are used (`part`, `material`, `section`, `assembly`, `step`, `interaction`, `load`, `mesh`, `job`, `sketch`, `visualization`, `abaqus`, `abaqusConstants`, `customKernel`) plus the standard libraries `math`, `random`, `datetime`, `operator`
- No external Python environment or additional packages required — the script runs **only inside Abaqus**

### How to run

```
abaqus cae noGUI=A_ESN_builder.py
```

or within Abaqus: File->Run Script


A new model (`Mdb()`) is created for every generated configuration, to avoid conflicts on repeated execution. The corresponding `.cae` file and Abaqus input file (`.inp`) are written to the working directory.

### Key parameters (set at the top of the script)

| Parameter | Meaning |
|---|---|
| `n_odb` | Number of ESN realizations to generate per configuration |
| `b_ns`, `t_ns` | Lists of RVE widths / thicknesses to generate [µm] |
| `l_s_mean` | Mean segment length between two crosslinks [µm] (governs the interaction range `Int_t`) |
| `eps11`, `eps12` | Applied strain / shear for the respective load case |
| `LoadCase` | `'Uniaxial'`, `'PlanarTension'`, or `'SimpleShear'` |
| `youngs`, `poisson`, `plasticity` | Material properties of the individual fiber |
| `d_f`, `por` | Fiber diameter and target porosity of the network (govern the generated total fiber length) |
| `Alg` | Degree of isotropy (RandomWalk only): 1 = isotropic, → 0 = increasingly aligned |
| `ang_av/-stdev/-min/-max` | Statistics of the segment-to-segment angle changes  |
| `l_step`, `res_fp` | Numerical resolution parameters (spline interpolation, segment length, meshing/partitioning resolution) |

### Workflow

1. **Parameter and material definition** – fiber and beam cross-section properties, load step (`StaticStep` with `nlgeom=ON` and stabilization), field/history output requests.
2. **Fiber generation** 
   - *RandomWalk*: fibers are generated step-by-step with a varying angle; when a fiber crosses an RVE boundary, it is trimmed exactly at the edge and continued on the opposite edge (direct PBC coupling).
   It runs until the total fiber length `TotalFibreLength`, derived from the target porosity, is reached.
3. **Meshing** – each fiber is recursively bisected until the element edge length falls below `res_fp`, followed by beam-element meshing (B22) and definition of a node set of potential crosslink nodes (`InternalNodes`).
4. **Crosslink generation** – a randomly shuffled fiber order (`FLS`) is used to create tie constraints between neighboring fibers within the calculated interaction range `Int_t` (derived from target density / mean segment length using a literature-based relation).
5. **Periodic boundary conditions (PBC)** – coupling equations (`m.Equation`) between opposite RVE edges, referenced to three reference points (`Origin`, `X-Dir`, `Y-Dir`).
6. **Load case & job** – depending on `LoadCase`, the reference-point boundary conditions are set within the load step, a description containing all relevant metrics (fiber count, fiber length, isotropy coefficient, etc.) is generated, and an Abaqus job including input file is written.


### Output files (per generated configuration)

- `<date>_[S|R]_<b_n>_<t_n>...<LoadCase>_<i>.cae` – Abaqus model database
- `<...>.inp` – Abaqus input file (via `writeInput`)
- `<...>_fiberstackorder.txt` – fiber order used for crosslink generation, plus the interaction range used


---

## 2. `B_Script_Coord_20Frames.py` – Nodal Coordinate Extraction

### Purpose

To analyze the strain fields of a simulated ESN, the mesh node coordinates must be exported from the FE result file. This script opens one or more Abaqus output database (ODB) files and exports the nodal coordinates from all 21 frames (reference frame 0 plus 20 subsequent frames) of the first analysis step to a plain-text file, providing the raw nodal displacement data used in the subsequent kinematic / strain-field analysis (e.g. deformation gradient, Green-Lagrange strain, "kinematic fingerprint" computation).

### Requirements

- **Abaqus/CAE** with the Abaqus Scripting Interface (Python 2 / `mbcs` encoding)
- Uses only Abaqus-internal modules (`abaqus`, `abaqusConstants`, `odbAccess`)
- Must be run from within Abaqus, with completed `.odb` result files (as produced by `A_ESN_builder.py` followed by an Abaqus/Standard analysis run) available in the working directory

### How to run

```
abaqus python B_Script_Coord_20Frames.py

or within Abaqus: File->Run Script
```

### Configuration

At the top of the script, the `files` list must be set to the ODB filenames to be processed:

```python
files = ['2024-05-26_R_200x6xlsm5_plast_viso0_91_UniAx_1.odb',
         '2024-05-27_R_200x6xlsm5_plast_viso0_92_UniAx_2.odb']
```

Each listed ODB file is opened, processed, and closed in turn.

### Workflow

1. For each ODB file, opens the database and identifies the first analysis step.
2. Accesses the `COORD` (nodal coordinate) field output for frames 0–20 (`frame_numbers = range(21)`).
3. Iterates over every part instance in the root assembly and, for each instance, builds a lookup of nodal coordinates per frame.
4. For every node, writes one line containing:
   - part instance name and node label
   - reference (frame 0) x/y coordinates
   - x/y coordinates for each of the 21 frames
5. Writes the result to a text file named `<odb_basename>_COORD.txt`, with a header row describing all columns.

### Output file format

One `.txt` file per input ODB, comma-separated, with the header:

```
Part_Instance, Node_ID, x_ref, y_ref, Frame_0_X, Frame_0_Y, Frame_1_X, Frame_1_Y, ..., Frame_20_X, Frame_20_Y
```

Each subsequent row corresponds to one mesh node, giving its reference position and its position at every extracted frame — the basis for the strain/kinematic-fingerprint analysis described in the manuscript.

---

## 3. `C_Import_COORD.m` – Importing the Nodal Coordinate File into MATLAB

### Purpose

Converts the plain-text nodal coordinate file exported by `B_Script_Coord_20Frames.py` (`<...>_COORD.txt`) into a MATLAB table `COORD` and saves it as a `.mat` file. This `.mat` file is the required input for both `D_Kinematic_Fingerprint_Analysis.mlx` (Section 4) and `E_RenderESN.mlx` (Section 5), so `C_Import_COORD.m` must be run once for every simulation before either of those scripts is used.

The `.mat` file included in this repository (`R_200x6xlsm5_plast_viso0_9_UniAx_0_COORD.mat`, see "Sample data / quick start" above) was generated with this script and can be used directly if no new simulation needs to be imported.

### Requirements

- MATLAB (tested as a plain `.m` script)
- No toolboxes beyond base MATLAB are required
- Input: the `<...>_COORD.txt` file produced by `B_Script_Coord_20Frames.py`, located in the current working directory

### Configuration

The name of the coordinate file (without the `.txt` extension) is set at the top of the script:

```matlab
name = '2024-09-13_R_40x1xlsm5_plast_viso0_49_UniAx_0_COORD';
```

### Workflow

1. Clears the workspace and command window.
2. Detects the import options for the `<name>.txt` file (`detectImportOptions`) and forces all variables to type `double`.
3. Configures the `Part_Instance` column to strip non-numeric characters and thousands separators (`,`), so the Abaqus part-instance label is imported as a numeric value consistent with the rest of the table.
4. Reads the file into the table `COORD` using `readtable`.
5. Removes the first three rows of the imported table (`COORD(1:3, :) = []`), which correspond to header/metadata rows not needed for the analysis.
6. Saves `COORD` to a `.mat` file, using the input filename with the leading date stamp removed (`name(1,12:end)`), e.g. `R_40x1xlsm5_plast_viso0_49_UniAx_0_COORD.mat`.

### Output file

- `<name_without_date_prefix>.mat` – MATLAB file containing the table `COORD`, ready to be loaded by `D_Kinematic_Fingerprint_Analysis.mlx` and `E_RenderESN.mlx`.


---

## 4. `D_Kinematic_Fingerprint_Analysis.mlx` – Kinematic Fingerprint Evaluation

### Purpose

This is the main MATLAB analysis script of the framework. It processes the nodal coordinates exported from Abaqus (`B_Script_Coord_20Frames.py`, imported via `C_Import_COORD.m`) to compute the local kinematic quantities described in the manuscript — the deformation gradient **F**, the Green-Lagrange strain components **E₁₁, E₂₂, E₁₂**, the local rotation angle **R**, and the deformation type exponent **m** — and assembles them into the statistical **kinematic fingerprint** of the electrospun network.

### Requirements

- MATLAB (tested as a MATLAB Live Script, `.mlx`)
- No toolboxes beyond base MATLAB are required
- Input: a `<...>_COORD.mat` file (table `COORD`, containing at least the columns `Part_Instance`, `Node_ID`, `Frame_0_X`, `Frame_0_Y`, `Frame_20_X`, `Frame_20_Y`, produced by `C_Import_COORD.m`) and the matching `<...>_fiberstackorder.txt` file
- **`Matlab_009_createPlots_FEM_Matlab_Q4_reference.m`** (included in this repository) – plotting helper function, required on the MATLAB path whenever `plotKinematicFields = true`; it generates the averaged kinematic field plots described in step 9 below

### Configuration (User settings section)

```matlab
workingDirectory   = '...';                                        % folder containing the input files
coordinateFile     = 'R_200x6xlsm5_plast_viso0_95_UniAx_0_COORDV1.mat';
fiberStackOrderFile = '2024-05-26_R_200x6xlsm5_plast_viso0_95_UniAx_0_fiberstackorder.txt';

fiberFractionUpper = 1.0;   % upper bound of the fiber-depth fraction to include
fiberFractionLower = 0.5;   % lower bound (e.g. 0.8–1.0 selects the top 20% of fibers)

facetWidths    = [10];      % facet width(s) b_q to analyze [µm]
networkWidths  = [200];     % network width(s) b_n to analyze [µm]

plotFacetSubdivision   = true;
plotKinematicFields    = true;
plotFingerprintHistograms = true;
```

Several `facetWidths` and `networkWidths` values can be supplied as vectors to sweep over multiple facet/network size combinations in a single run.

To run the script directly with the sample data included in this repository, set:

```matlab
coordinateFile      = 'R_200x6xlsm5_plast_viso0_9_UniAx_0_COORD.mat';
fiberStackOrderFile = '2024-05-26_R_200x6xlsm5_plast_viso0_9_UniAx_0_fiberstackorder.txt';
```

### Workflow

1. **Load input data** – loads the `COORD` table and validates that all required columns (`Frame_0_X/Y`, `Frame_20_X/Y`, `Part_Instance`, `Node_ID`) are present.
2. **Region cropping** – for each requested network width `b_n`, crops the nodal point cloud to the corresponding centered square region.
3. **Fiber-depth selection** – reads the fiber stacking order and keeps only the nodes belonging to the fiber fraction defined by `fiberFractionLower`/`fiberFractionUpper` (depth-resolved analysis).
4. **Unique node IDs** – combines `Part_Instance` and `Node_ID` into a single unique node identifier, since Abaqus node labels may repeat across part instances.
5. **Facet grid & facet centers** – subdivides the analysis region into a grid of square facets of width `b_q` (representative area elements, RAE), assigns nodes to facets, and computes the geometric facet center as the arithmetic mean of all node coordinates within each facet, both in the reference (frame 0) and current (frame 20) configuration.
6. **Local Q4-based kinematics** – treats each 2×2 block of neighboring facet centers as a bilinear (Q4) element and, at each of its four corners, computes:
   - the deformation gradient **F = I + ∇u** from the Q4 shape-function derivatives and the corner displacements,
   - the Green-Lagrange strain **E = ½(FᵀF − I)** (components E₁₁, E₂₂, E₁₂),
   - the local rotation angle via polar decomposition of **F** (`polar_rotation_2D`),
   - the deformation type exponent **m**, derived from the eigenvalues of the right Cauchy-Green tensor **C = FᵀF** (`deformation_type_exponent`; distinguishes uniaxial/equibiaxial/planar-type local deformation).

   These four **non-averaged** corner states per element form the basis of the kinematic fingerprint.
7. **Averaged fields for visualization** – **F** is averaged across all element corners contributing to a given facet center; E, R, and m for the plotted fields are then computed from this averaged **F** (not averaged individually), consistent with the manuscript's methodology.
8. **Kinematic fingerprint assembly** – the non-averaged local E₁₁, E₂₂, E₁₂, R, and m values are flattened into fingerprint vectors, invalid (NaN) states are removed, and a fingerprint table plus summary statistics (mean, standard deviation) are compiled.
9. **Plotting** (optional, controlled by the `plot...` flags) – facet subdivision and displacement field, averaged kinematic fields (via `Matlab_009_createPlots_FEM_Matlab_Q4_reference.m`), and histograms of the fingerprint distributions (E₁₁, E₂₂, E₁₂, R, m).
10. **Results storage** – all results per network-width/facet-width combination are stored in the `resultsBySetting` cell array (containing raw, averaged, and fingerprint data), and convenience variables (`E11_vec`, `E22_vec`, `E12_vec`, `R_vec`, `R_rad_vec`, `m_vec`, `fingerprint_table`, `stats`) are exposed for the last analyzed setting.

### Local functions (defined at the end of the script)

- `q4_shape_function_derivatives(xi, eta)` – derivatives of the bilinear Q4 shape functions with respect to natural coordinates.
- `polar_rotation_2D(F)` – rotation tensor from the polar decomposition of a 2D deformation gradient (via SVD), with correction for numerical reflections.
- `deformation_type_exponent(F)` – computes the deformation type exponent **m** from the invariants of **C = FᵀF**, with dedicated handling of the reference/undeformed state and the equibiaxial case.

### Output

- MATLAB workspace variables and the `resultsBySetting` cell array, containing per-setting facet centers, displacement fields, non-averaged local kinematic quantities, averaged fields for visualization, and the fingerprint table/statistics.
- Figures (if enabled): facet subdivision plot, averaged kinematic field plots (via `Matlab_009_createPlots_FEM_Matlab_Q4_reference.m`), and kinematic fingerprint histograms.

---

## 5. `E_RenderESN.mlx` – Rendering ESN Images

### Purpose

After creating the `.mat` file containing the `COORD` table (via `C_Import_COORD.m`, see Section 3), `E_RenderESN.mlx` renders greyscale line-drawing images of the ESN fibre network.

By default, all 21 deformation states, from frame 0 to frame 20, are rendered. The frame range can be changed in the following loop:

```matlab
for frame_index = 0:20
```

### Requirements

* MATLAB, tested as a MATLAB Live Script (`.mlx`)
* No toolboxes beyond base MATLAB are required
* A `<...>_COORD.mat` file containing the table `COORD` (produced by `C_Import_COORD.m`)
* The corresponding `<...>_fiberstackorder.txt` file

Both input files must originate from the same Abaqus simulation generated using `A_ESN_builder.py` and processed using `B_Script_Coord_20Frames.py`.

### Configuration

The two required input files are specified at the beginning of the script:

```matlab
load('R_40x1xlsm5_plast_viso0_49_UniAx_0_COORD.mat');

fibre_stack_filename = ...
    '2024-09-13_R_40x1xlsm5_plast_viso0_49_UniAx_0_fiberstackorder.txt';
```

To run the script directly with the sample data included in this repository, set:

```matlab
load('R_200x6xlsm5_plast_viso0_9_UniAx_0_COORD.mat');

fibre_stack_filename = ...
    '2024-05-26_R_200x6xlsm5_plast_viso0_9_UniAx_0_fiberstackorder.txt';
```

Note that in this case `membrane_width_um` (see below) must be set to `200`, matching the network width of the sample data.

The following rendering parameters can then be adjusted:

| Parameter                | Meaning                                                                                                                  |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------ |
| `fibre_fraction_to_keep` | Fraction of fibres retained from the top of the fibre stack. For example, `0.8` retains the uppermost 80% of the fibres. |
| `membrane_width_um`      | Membrane width in µm. This value must correspond to the actual width of the rendered ESN.                                |
| `darkest_grey_value`     | Grey value of the lowest retained fibre, where `0` is black and `1` is white.                                            |
| `brightest_grey_value`   | Grey value of the uppermost fibre, where `0` is black and `1` is white.                                                  |

Example:

```matlab
fibre_fraction_to_keep = 1.0;
membrane_width_um = 40;

darkest_grey_value = 0.10;
brightest_grey_value = 0.85;
```

To preserve the correct physical fibre thickness in the rendered image, the fibre diameter must correspond to `d_f` used in `A_ESN_builder.py`.

The default value is:

```matlab
fibre_diameter_mm = 0.403;
```

The diameter is converted from millimetres to MATLAB line-width points using:

```matlab
millimetres_per_point = 25.4 / 72;
line_width_points = fibre_diameter_mm / millimetres_per_point;
```

### Workflow

1. The script loads the `COORD` table and reads the fibre stacking order from the corresponding `_fiberstackorder.txt` file.

2. The fibre stacking order is reversed so that the uppermost fibre is listed first.

3. The number of fibres to retain is calculated from `fibre_fraction_to_keep`. Only nodes belonging to the selected upper fraction of the fibre stack are retained.

4. Each fibre path is bisected, and only the first half of its coordinate points is retained for rendering.

5. A linearly varying grey value is assigned to each retained fibre. The gradient extends from `brightest_grey_value` for the uppermost fibre to `darkest_grey_value` for the lowest fibre.

6. When `fibre_fraction_to_keep < 1`, the darkest grey value is adjusted to compensate for the omitted lower fibres and preserve the intended depth-dependent brightness gradient.

7. For each selected frame, the fibres are plotted from the lowest to the uppermost stacking position on a dark-grey background.

8. The line width is scaled according to the physical fibre diameter, while the plot dimensions and axis limits are defined by `membrane_width_um`.

9. Each frame is exported as a PNG image at 500 dpi.

### Output Files

By default, the following files are generated:

```text
fibre_network_frame_00.png
fibre_network_frame_01.png
...
fibre_network_frame_20.png
```

Each image represents the ESN fibre network at one analysed deformation state.

---

## Included example data files

To make the analysis and rendering scripts (Sections 4 and 5) runnable out of the box, this repository includes one example dataset (network width 200 µm, uniaxial loading), consisting of:

| File | Produced by | Used by |
|---|---|---|
| `R_200x6xlsm5_plast_viso0_9_UniAx_0_COORD.mat` | `C_Import_COORD.m` | `D_Kinematic_Fingerprint_Analysis.mlx`, `E_RenderESN.mlx` |
| `2024-05-26_R_200x6xlsm5_plast_viso0_9_UniAx_0_fiberstackorder.txt` | `A_ESN_builder.py` | `D_Kinematic_Fingerprint_Analysis.mlx`, `E_RenderESN.mlx` |
| `Matlab_009_createPlots_FEM_Matlab_Q4_reference.m` | – (standalone plotting helper) | `D_Kinematic_Fingerprint_Analysis.mlx` (when `plotKinematicFields = true`) |

These files allow `D_Kinematic_Fingerprint_Analysis.mlx` and `E_RenderESN.mlx` to be run and inspected without needing Abaqus or access to `A_ESN_builder.py` / `B_Script_Coord_20Frames.py` / `C_Import_COORD.m`.

---

## Relation to the manuscript

`A_ESN_builder.py` generates the ESN RVE realizations that are subsequently subjected to virtual mechanical loading (uniaxial, planar tension, simple shear) in Abaqus. `B_Script_Coord_20Frames.py` extracts the resulting nodal coordinate histories from the simulation output. `C_Import_COORD.m` imports these coordinate histories into MATLAB. `D_Kinematic_Fingerprint_Analysis.mlx` performs the core continuum-mechanics-based evaluation described in the manuscript — RAE-based facet-center kinematics and the resulting "kinematic fingerprint" (E₁₁, E₂₂, E₁₂, R, m) — that forms the main quantitative result of the study. `E_RenderESN.mlx` renders the coordinate histories into grayscale ESN images for visual inspection and figures.
