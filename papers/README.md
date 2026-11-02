# Papers

## In this folder
| File | Reference | Used for |
|---|---|---|
| `Assignment_2_ME_319_Control_Systems_S3.pdf` | ME 319 Assignment 2 spec | Problem statement |
| `deGroot2024_T-MPC_2401.06021v2.pdf` | O. de Groot, L. Ferranti, D. Gavrila, J. Alonso-Mora, *Topology-Driven Parallel Trajectory Optimization in Dynamic Environments*, arXiv:2401.06021v2, 2024 | Parallel distinct hypotheses (homotopy classes), consistent selection (Eq. 12) |
| `Williams2017_MPPI_InformationTheoreticMPC.pdf` | G. Williams et al., *Information Theoretic MPC for Model-Based Reinforcement Learning*, ICRA 2017 | Sampling-based rollouts (MPPI) for P(hit) estimation |
| `Singer1970_Optimal_Tracking_Filter_Manned_Maneuvering.pdf` | R. A. Singer, *Estimating Optimal Tracking Filter Performance for Manned Maneuvering Targets*, IEEE TAES 6(4), 1970 | Correlated-acceleration model, Q from a_max |
| `Li_Jilkov2003_Survey_Part1_Dynamic_Models.pdf` | X. R. Li, V. P. Jilkov, *Survey of Maneuvering Target Tracking. Part I: Dynamic Models*, IEEE TAES 39(4), 2003 | Motion-model choice (CV/CA/Singer/jerk/CT) |
| `Li_Jilkov2005_Survey_Part5_Multiple_Model.pdf` | X. R. Li, V. P. Jilkov, *Survey of Maneuvering Target Tracking. Part V: Multiple-Model Methods*, IEEE TAES 41(4), 2005 | Multiple-model / IMM theory |
| `Mazor_et_al1998_IMM_Survey.pdf` | E. Mazor, A. Averbuch, Y. Bar-Shalom, J. Dayan, *Interacting Multiple Model Methods in Target Tracking: A Survey*, IEEE TAES 34(1), 1998 | IMM algorithm |
| `Pitre_Jilkov_Li2005_Comparative_MM_Algorithms.pdf` | R. R. Pitre, V. P. Jilkov, X. R. Li, *A Comparative Study of Multiple-Model Algorithms for Maneuvering Target Tracking*, SPIE 5809, 2005 | Justification for IMM-class cost/accuracy |
| `Survey_Manoeuvring_Target_Tracking_arXiv1503.07828.pdf` | *A Survey of Manoeuvring Target Tracking Methods*, arXiv:1503.07828 | Singer, CSM, coasting/divergence |
| `Davies2013_Integrated_Target_Tracking_Weapon_Guidance_thesis.pdf` | J. M. Davies, *Integrated Target Tracking and Weapon Guidance*, PhD thesis, Univ. of Liverpool, 2013 | IMM design for weapon guidance |
| `IMM_Warship_Gun_FCS_koreascience.pdf` | IMM tracking filter for a warship gun fire-control system (CV/Singer/CTR) | Prior art: IMM in fire control |
| `Patent_US4794235_NonLinear_Prediction_Gun_FC.pdf` | US 4,794,235, *Non-linear prediction for gun fire control systems* | Prior art: LOS-frame Kalman prediction |
| `Lerro_BarShalom1993_Debiased_Converted_Measurements_vs_EKF.pdf` | D. Lerro, Y. Bar-Shalom, *Tracking with Debiased Consistent Converted Measurements versus EKF*, IEEE TAES 29(3), 1993 | ⚠ **Incomplete**: browser screenshot of page 1 only. Re-download with the IEEE viewer's download button. |

## Still to get (blocked to scripted download; fetch via browser / IITB proxy)
| Reference | Link | Used for |
|---|---|---|
| K. Mehrotra, P. R. Mahapatra, *A Jerk Model for Tracking Highly Maneuvering Targets*, IEEE TAES 33(4):1094–1105, 1997 | https://eprints.iisc.ac.in/2710/ (IISc copy) · https://doi.org/10.1109/7.624345 | Order-4 (jerk) estimator |
| K. D. Glazebrook, A. R. Washburn, *Shoot-Look-Shoot: A Review and Extension*, Operations Research, 2004 | https://apps.dtic.mil/sti/tr/pdf/ADA487449.pdf | Firing-policy prior art |
| Y. Pryluk, T. Shima, O. Golan, *Shoot-Shoot-Look for an Air Defense System* | https://cris.iucc.ac.il/en/publications/shoot-shoot-look-for-an-air-defense-system/ | Firing without kill assessment |
| DRDC, probability of raid annihilation / firing policy report | https://publications.gc.ca/collections/collection_2017/rddc-drdc/D69-35-2014-eng.pdf | Firing-policy prior art |
| NPS, *VLS Missile Mix, Firing Policy, and Deterrence* | https://core.ac.uk/download/528270655.pdf | Firing-policy prior art |
| DTIC ADA063704, anti-aircraft gun fire control (Kalman target-state estimator) | https://apps.dtic.mil/sti/pdfs/ADA063704.pdf | Prior art |

## Web references (not papers)
- Phalanx CIWS rate of fire / bursts: https://navweaps.com/Weapons/WNUS_Phalanx.php
- Laser-weapon magazine/dwell limits: https://www.laserwars.net/p/laser-weapon-infinite-magazine-myth
