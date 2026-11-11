# TLS-ERT-validation: TLS-informed ERT inversion of underground cavities, code and derived data

Companion repository of the paper **"Where Does Geometry Help? A Non-Circular Validation of TLS-Informed ERT Inversion of Underground Cavities"** (Grassi S., Morreale G., Imposa S., 2026).

It contains the MATLAB scripts, the constraint geometries (box files), the inverted models and the evaluation outputs used in the paper, for two case studies in eastern Sicily: the **Micio Conti lava tube** and the **Chiusazza karst cave**.

Repository: https://github.com/g-Sabrina87/TLS-ERT-validation

## 1. What the study does

At each site, three ERT inversions are compared on the same apparent-resistivity dataset, mesh and settings:

| Model | Constraint in the starting model |
|---|---|
| **U** | none (unconstrained) |
| **C** | TLS geometry of the whole cavity (positive control) |
| **H** | TLS geometry of the constraint region only; a spatial sector is held out |

The held-out sector is evaluated against the TLS geometry of the full cavity, so that the evaluation is non-circular. Evaluation is carried out on mesh elements, excluding the neighborhood of the constraint points and limiting the set to the vicinity of the cavity. Null constraints (translated copies of the constraint cloud) are used as a negative control. Agreement is measured with rank AUC (primary), average precision and IoU, with a paired spatial block bootstrap. Two further held-out sectors per site are evaluated as a sensitivity analysis. The quick test (`quick_test.m`, Section 8) reproduces the main result at Micio Conti in a few seconds.

## 2. Software

- MATLAB R2020b or later. The code was developed and tested with MATLAB R2025b (25.2, Update 3). Base MATLAB: `alphaShape`/`inShape`, `inpolygon`, `readmatrix`, `fminsearch`, `prctile` (in R2025b; in older releases `prctile` may belong to the Statistics and Machine Learning Toolbox). Statistics and Machine Learning Toolbox: `knnsearch`, `kmeans`, `pdist`/`pdist2`, `silhouette`, `corr`, `boxplot`, `rangesearch`. Image Processing Toolbox: `bwconncomp` (connected components of the cavity cluster). Toolbox membership was checked with `which` on R2025b.
- ERTLab64 - 1.2.0 (Geostudi Astier s.r.l.), commercial, version as reported by the software. The inversions were run in ERTLab64; its code is not distributed. The inputs for ERTLab64 (box files) and the inverted models are provided, so that the evaluation scripts and the quick test do not require ERTLab64.

## 3. Data

| Item | Content | Availability |
|---|---|---|
| ERT data | apparent resistivity datasets of both sites (1,938 quadrupoles used at Micio Conti, 3,535 at Chiusazza) | not deposited; available from the corresponding author on reasonable request |
| Electrode coordinates | electrode positions of both sites | this repository |
| TLS point clouds | cropped to the ERT area (full cloud and constraint-only cloud) | this repository |
| Constraint geometries | box files (`xmin xmax ymin ymax zmin zmax rho`) of C and H at each site (the exact inputs of the inversions of the paper); the boxes of the nulls and of the additional held-out sectors are not stored and are regenerated (see below) | this repository (C and H); regenerated (nulls, additional sectors) |
| Inverted models | element centroid coordinates and resistivity for U, C, H and the nulls | this repository |
| Inversion reports | ERTLab64 run reports | not deposited for all runs because of their size |

Folders of the repository:
- `Model_TLS/`: TLS point clouds of the two caves, full and cropped to the ERT area, electrode coordinates.
- `Constraint geometries/`: box files of C and H at each site, as used in the inversions.
- `5_Models_ERT_U_C_H_Null/`: inverted ERT models (U, C, H and the null constraints) of both sites.
- the numbered MATLAB scripts listed in Section 5, `quick_test.m` (Section 8) and `LICENSE`, in the root folder.

Original survey data and TLS: Grassi, S., Morreale, G., Giudice, G., Messina, D., Quattrocchi, G., Monforte, P., & Imposa, S. (2026). Advanced geophysical data analysis and 3D surveys for the detection and characterization of volcanic and karstic caves. Remote Sensing Applications: Society and Environment, 42, 102101.

## 4. Parameters used in the paper

| | Micio Conti | Chiusazza |
|---|---|---|
| Alpha-shape radius α | 1.35 m | 2.50 m |
| Imposed resistivity (C, H, nulls) | 75,000 Ω·m | 20,000 Ω·m |
| Voxel size for the boxes | 0.5 m | 1.0 m |
| Classification threshold (IoU) | 48,107 Ω·m | 5,513 Ω·m |
| Exclusion distance from constraint points | 1.5 m | 1.0 m |
| Maximum distance from the TLS cloud (evaluation set) | 3 m | 3 m |
| TLS points in the constraint region / held-out sector | 21,966 (33%) / 43,776 | 30,195 (about 60%) / 20,510 |
| Evaluation elements (inside the cavity) | 355 (68) | 11,758 (2,200) |
| Number of nulls | 30 (second series, main); 20 (first series, supplementary) | 20 |
| Position of the nulls | second series: at least 6 m from the full TLS, vertical shift within ±4 m, margin of 8 m from the lateral edges, cloud inside the mesh; first series: unconstrained (12 of 20 touch the cavity) | unconstrained (all 20 touch the cavity) |
| Assumed data error (apparent resistivity) | 2% of the value + 1e-4 Ω·m | 2% of the value + 1e-4 Ω·m |
| Spatial block sizes (bootstrap, x = y; vertical size 2 m) | 3, 5, 7, 10 m | 3, 5, 7, 10, 15, 20 m |
| Bootstrap replicates | 2,000 (seed 1) | 2,000 (seed 1) |

## 5. Workflow and scripts

The scripts are numbered in the order of the workflow. File names follow the convention `<Site>_cropped_ERTarea.txt` (full TLS cloud), `<Site>_constraint_only.txt`, `<Site>_heldout.txt`, `<Site>_holdout_region.mat` (structure `HELDOUT` with the held-out rectangle), `Model_<Site>_U.txt`, `Model_<Site>_C_rho<value>.txt`, `Model_<Site>_H_rho<value>.txt`, with `<Site>` = `Micio Conti` or `Chiusazza`. Most have a site switch (`site = 'Micio Conti'` or `'Chiusazza'`) and a settings block at the top where file names and parameters are set. All scripts are run from the folder that contains the input files. The numbers without scripts correspond to data folders and to the inversions.

**1 — Preparation of the TLS cloud**
- `1_crop_tls_to_ert_area.m`: crops the TLS cloud to the area investigated by the ERT.

**2 — Choice of the alpha-shape radius**
- `2_1_1_MicioConti_alphaRadius_test.m`, `2_1_2_Chiusazza_alphaRadius_test.m`: elbow method and geometric index J(α) at each site.
- `2_2_J_alpha_weight_formalization.m`: normalized weights of the index and sensitivity of the choice to them.

**3 — Hold-out partition and constraint geometry**
- `3_1_split_tls_holdout.m`: splits the TLS cloud (switch `site`) into constraint and held-out parts and saves the held-out polygon, the Z range and the held-out rectangle `HELDOUT` (bounding box of the held-out points). At Micio Conti the held-out region is an XY polygon with Z between 268.50 and 278.93 m; the upper limit lies above the top of the cloud (277.35 m), so there is no upper limit.
- `3_2_4_2_cavity_to_rect_anomalies.m`: builds the box file for ERTLab64 from a cloud (cleaning, alpha shape, voxelization, merge into boxes). The switch `site` selects the parameters of the site and `useHoldout` selects H (boxes inside the held-out region removed) or C (full geometry, from the full cloud). It is used for the geometries of C and H (step 3.2) and for the null constraints (step 4.2).
- `3_3_1_Chiusazza_extra_holdouts.m`, `3_3_2_Micio_extra_holdouts.m`: additional held-out sectors (two strips at Chiusazza, two slabs at Micio Conti).

**Regenerating the box files.** The boxes are a deterministic function of the constraint cloud and of the parameters in the script: `3_2_4_2_cavity_to_rect_anomalies.m` (`useHoldout = true` for H, `false` for C) generates the box files from the TLS clouds; the C and H files used in the inversions of the paper are provided here. The script regenerates H at Micio Conti exactly (166 boxes, 126.0 m³); for the other three files it gives a close but not identical geometry. The files provided here, not the regenerated ones, are the inputs of the inversions of the paper. The boxes of the nulls are obtained by rebuilding the translated clouds with `4_1_make_null_constraint_from_TLS.m` (`mode = 'reproduce'`) and running `3_2_4_2_cavity_to_rect_anomalies.m` with `nullMode = true`; the boxes of the additional held-out sectors by `3_3_1`/`3_3_2` (`modo = 'scrivi'`) followed by `3_2_4_2_cavity_to_rect_anomalies.m` inside each folder.

**4 — Null constraints**
- `4_1_make_null_constraint_from_TLS.m`: translated copies of the constraint cloud. `mode = 'reproduce'` rebuilds the nulls used in the paper from the stored translations (`NULL_summary_<Site>.csv`: 30 nulls at Micio Conti, 20 at Chiusazza); `mode = 'generate'` draws new random translations (Micio Conti: no point closer than 6 m to the full TLS, vertical shift within ±4 m, margin of 8 m from the lateral edges of the model, cloud inside the mesh; Chiusazza: no constraint on the position).
- Boxes of the nulls (step 4.2): `3_2_4_2_cavity_to_rect_anomalies.m` with `nullMode = true` (Micio Conti) processes the 30 translated clouds in a loop, with the clipping area set to the extent of the model mesh instead of the ERT quadrilateral (the nulls lie partly outside the quadrilateral), and writes `cavity_anomalies_NULL_NN_75000.txt`. The regenerated nulls have 163–179 boxes and 124.5–130.0 m³, as in the paper.
- `4_3_check_null_boxes.m`: checks on number, volume and domain of the boxes (Micio Conti).
- `4_4_null_position_check.m`: position of the nulls with respect to the held-out sector, the cavity and the electrodes.

**5 — Inversions in ERTLab64 (not scripted)**
Models U, C, H and the nulls are inverted with the same settings (Section 4, Section 6 and the paper). Output in `5_Models_ERT_U_C_H_Null/`: model files with element centroids and resistivity.

**6 — Classification threshold**
- `6_1_calibrate_threshold_constraint_only.m`: threshold from the constraint-region geometry only (used for the IoU of the hold-out evaluation).
- `6_2_kmeans_threshold_constrained.m`: K-means clustering of the resistivity values of the model with the full geometry (number of clusters by elbow and silhouette).
- `6_3_threshold_FINAL.m`: threshold by volume matching with the K-means cluster and with the TLS geometry.
- `6_4_threshold_sensitivity.m`: sensitivity of the threshold to the imposed resistivity.
The estimates of 6_2–6_4 use the full geometry (circular), are descriptive and are not used for the evaluation.

**7 — Evaluation**
- `7_1_bootstrap_UH_spatial_block_buffer.m`: paired spatial block bootstrap; AUC, IoU, AP; buffer and distance rules.
- `7_2_null_diagnostics_buffer.m`: evaluation of U, H and the nulls on the same elements, for several exclusion radii and distances from the TLS (radius 0 and no distance limit: evaluation on the whole held-out rectangle).
- `7_3_1_null_distance_vs_auc_spearman_Micio_Conti.m`: Spearman rank correlation between the distance of the nulls from the cavity and the absolute AUC difference from U (Micio Conti, second series).
- `7_3_2_null_overlap_vs_auc_Chiusazza.m`: overlap of the nulls with the TLS cloud (share of points within 3 m) and with the held-out sector, and its relation with AUC, average precision and mean IoU (Spearman and Pearson coefficients, descriptive; Chiusazza).
- `7_4_1_iou_threshold_sweep_MicioConti.m`, `7_4_2_iou_threshold_sweep_Chiusazza.m`: IoU and average precision over a range of thresholds, for H, U and the nulls.
- `7_5_alpha_sensitivity_and_distances.m`: sensitivity to the α of the ground truth and to the evaluation distances.
- `7_6_variogram_blocksize.m`: autocorrelation range of log10 resistivity, to compare with the block sizes.

**8 — Behaviour of the constraint**
- `8_1_attenuation_check.m`, `8_2_attenuation_check_fractions.m`: how much of the imposed value is retained in the constrained region (median; fraction of elements at 25, 50, 75 and 100% of the imposed value).
- `8_3_constraint_influence_decay.m`: change produced by the constraint as a function of the distance from it.
- `8_4_1_check_assigned_site_Chiusazza.m`, `8_4_2_check_assigned_vs_evaluated_Micio_Conti.m`, `8_5_constraint_mesh_discretization_check.m`: assignment of the boxes to mesh elements and relation with the evaluation set.
- `8_6_saturation_curve_fit.m`: χ² as a function of the imposed resistivity (saturation curve).

**9–10 — Figures**

| Figure of the paper | Script |
|---|---|
| Figure 1, plan view of the two sites | `10_1_make_figure1_plan_view.m` |
| Figure 2, conceptual scheme of the U/C/H protocol | drawn separately, no script |
| Figure B1 [Supplementary Material], χ² as a function of the imposed resistivity | `10_2_make_figureB1_chi2_vs_rho.m` |
| Figure 3, sensitivity to the block size | `10_3_make_figure3_block_size_sensitivity.m` |
| Figure 4, null-constraint control | `10_4_make_figure4_null_control.m` |
| Figure 5, plan-view maps of U, C−U and H−U | `9_make_figure5_models_UCH.m` |

The figure numbers in the script names are those of the paper.

## 6. Notes

- The inversions stop at the first iteration with a reduced χ² below 1. Reduced χ² is the final weighted data residual divided by the number of quadrupoles and is compared only between inversions with the same assumed error and the same sequence of iterations. At Chiusazza, 17 of the 20 nulls follow the same sequence as U, C and H; the other three are excluded from the comparisons of χ².
- The reduced χ² values reported in the paper are the final weighted residual of each inversion divided by the number of quadrupoles, as printed in the ERTLab64 run report. The run reports of all inversions are not deposited because of their size; those of U, C and H are available on request.
- The nulls are translated copies of one constraint shape and not independent replicates. The counts derived from them are descriptive.
- Scripts are given as used for the paper; the file names in the settings blocks correspond to the model and cloud files of each run.

## 7. License and citation

Code: MIT License (see `LICENSE`). Data: CC BY 4.0. The data licence applies to the files deposited here.

To cite this repository: Grassi S., Morreale G., Imposa S., TTLS-ERT-validation: TLS-informed ERT inversion of underground cavities, code and derived data. GitHub, 2026, https://github.com/g-Sabrina87/TLS-ERT-validation. Please cite also the paper above.

## 8. Quick test

`quick_test.m` checks that the repository reproduces the main result of the paper at Micio Conti without ERTLab64 and in a few seconds. It splits the TLS cloud into constraint region and held-out sector, builds the main evaluation set, computes the rank AUC of the inverted models U, H and C, runs the paired spatial block bootstrap (blocks of 3 x 3 x 2 m, 2,000 replicates) and compares the results with the values of the paper.

**How to run.** Keep the repository folders as downloaded and run in MATLAB (R2020b or later; tested with R2025b; the Statistics and Machine Learning Toolbox is needed for `prctile` in older releases):

```matlab
cd <path to the repository>
quick_test
```

The script searches the files it needs (`Micio Conti_cropped_ERTarea.txt`, `Model_Micio Conti_U.txt`, `Model_Micio Conti_H_rho75000.txt`, `Model_Micio Conti_C_rho75000.txt`) anywhere under its own folder, so no path has to be edited.

**Expected output.**

| Quantity | Value |
|---|---|
| TLS points: constraint region / held-out sector | 21,966 (33.4%) / 43,776 |
| Evaluation elements (inside the cavity) | 355 (68) |
| AUC of U / H / C | 0.614 / 0.614 / 0.807 |
| Delta AUC H - U, 95% interval | 0.000 [-0.006, 0.006] |
| Delta AUC C - U, 95% interval | +0.192 [0.123, 0.269] |

The intervals are those of Table D.1 of the Supplementary Material for blocks of 3 m (the bootstrap uses `rng(1)`). The script prints PASS or FAIL for each comparison and ends with `All checks passed`.
