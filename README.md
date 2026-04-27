# Spatial Analysis of Oncogene Co-expressing Cells in DLBCL

This repository contains the analysis code accompanying the manuscript:

**"The spatial organization of oncogene co-expressing cells influences survival 
and T cell exhaustion in lymphoma"**

Sridhar S, Ong C, Ringler VVL, Gupta K, et al.

---

## Related Repositories

- **PatternExtract** — automated pipeline for spatial point pattern generation 
  from multiplexed immunofluorescence images:  
  https://github.com/shrutisridhar99/PatternExtract

- **scDEL Calculator** — web-based tool for 2-year relapse risk prediction 
  in DLBCL from routine IHC:  
  https://github.com/Kakarrxt/scDEL_calculator  
  https://scdel-calc.streamlit.app

---

## Dependencies

### R packages
| Package | Version | Use |
|---|---|---|
| spatstat | 3.0 | Spatial point process modelling |
| survival | 3.5 | Survival analysis |
| survminer | 0.4.9 | Kaplan-Meier visualisation |
| contsurvplot | — | Continuous survival heatmaps |
| limma | — | DSP differential expression |
| GSVA | — | Gene set scoring |
| CellChat | 1.6.1 | Cell-cell communication |
| Seurat | 4.3.0 | Single-cell RNA-seq analysis |
| pROC | — | ROC curves and DeLong test |
| irr | — | ICC computation |
| pheatmap | — | Heatmap visualisation |
| dplyr, ggplot2, patchwork | — | Data manipulation and plotting |

### Python packages
| Package | Version | Use |
|---|---|---|
| scanpy | — | Single-cell analysis |
| anndata | — | Data format |
| numpy, pandas | — | Data manipulation |

---

## Data Availability

- **BCA cohort RNA-seq**: European Genome-phenome Archive (EGA), 
  Study ID EGAD00001003783
- **GEP cohorts**: Gene Expression Omnibus (GEO) — accession numbers 
  GSE117556, GSE125966, GSE31312, GSE10846, GSE87371, GSE32918, GSE98588
- **Reddy et al. cohort**: EGA, Study ID EGAS00001002606
- **Schmitz et al. cohort**: dbGaP, accession phs001444.v2.p1
- **POLARIX trial data**: analysed in collaboration with F. Hoffmann-La Roche Ltd.; 
  not publicly available
- **NUH PhenoCycler-Fusion data**: available upon reasonable request

---

## Usage Notes

- All file paths in scripts are set to local directories and must be updated 
  to reflect your data location before running
- The POLARIX biomarker analysis is not included in this repository as the 
  data are restricted; methods are described in the supplementary methods
- SPACEc cell phenotyping pipeline follows the standard workflow available at 
  https://github.com/yuqiyuqitan/SPACEc.git

---

## Citation

If you use this code, please cite:

---

## Contact

For questions regarding the code, please contact:  
Shruti Sridhar — shruti.sridhar@u.nus.edu  
Anand D. Jeyasekharan — csiadj@nus.edu.sg
