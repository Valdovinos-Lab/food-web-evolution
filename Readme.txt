Readme file for Generalist Evolution Project

	- Latest update 9-24-2021
	- Overview of code and data used in preparation of manuscript submitted to Scientific Reports, Sept. 2021 (J. R. Morris, K. A. Allhoff, F. S. Valdovinos)

Contents:
	1) Eco-evolutionary food web (C code)
	2) Data processing and analysis (R code)
	3) Summary data file documentation (R output)
	4) R version and package info


1) Eco-Evolutionary Food Web Model (C code)

The code can be run using the following commands below. See comments in MODEL_CLEAN.c file for details.

  Create executable file:
	gcc MODEL_CLEAN.c -std=c99 -lm -lgsl -lgslcblas -ffast-math -o MODEL_CLEAN.out

  Run simulation:
	./MODEL_CLEAN.out res_CLEAN 12345 0.05 1 0.5 0.2 0.4 0.4 0.3 1


2) Data Processing and Analysis Code (R code)

Two files are included for all data processing and analysis performed in the manuscript.

	ANALYSIS_ZSWEEP_HIST_FINAL.R
	ANALYSIS_ZSWEEP_COM_FINAL.R

Both files include code for processing raw data output from simulations, but also summary data generated in R. Statical analysis and figure generation code is included in both files. Import summary data files included here to run this analysis. See R code files for further details.


3) Data Summary Files (from R output)

All data summary files used in statical and visual output for the manuscript are included here. File names and description are listed below. See data analysis code for more details. All files are located in the folder Z_SWEEP_DATA_SUMMARY_FINAL

ALLBINDsub.csv	Species traits, biomass, etc. (subset of replicates)
ALLTROPsub.csv	Species trophic position, etc. (subset of replicates)
CBsum.csv	Community biomass (SD)
FINALFRBINDSAMPLE.csv	Realized feeding range for example simulation
GSLSEXTALL.csv	Generalist/specialist lifespan (mean)
MIFRALL.csv	Feeding range (mean)
MIFRALLmed.csv	Feeding range (median)
RBMsum.csv	Resource biomass (SD)
REALFRSUM.csv	Realized feeding range (mean)
SLOPEALL.csv	Lifespan slope coefficients
SLOPEEXAMPLE.csv	Lifespan to feeding range for z sweep replicate example
TOSsum.csv	Species turnover (mean)


4) R Version and Package Details for Data Processing and Analysis

Output from R sessionInfo()

R version 4.1.0 (2021-05-18)Platform: x86_64-apple-darwin17.0 (64-bit)Running under: macOS Big Sur 11.5.2Matrix products: defaultLAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylibRandom number generation: RNG:     Mersenne-Twister  Normal:  Inversion  Sample:  Rounding  locale:[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8attached base packages:[1] stats     graphics  grDevices utils     datasets  methods   base     other attached packages: [1] mgcv_1.8-36     nlme_3.1-152    rsq_2.2         pryr_0.1.4      plotrix_3.8-1   forcats_0.5.1   stringr_1.4.0   dplyr_1.0.7     purrr_0.3.4     readr_1.4.0     tidyr_1.1.3    [12] tibble_3.1.2    tidyverse_1.3.1 reshape2_1.4.4  gridExtra_2.3   ggthemes_4.2.4  ggplot2_3.3.5   Rmisc_1.5       plyr_1.8.6      lattice_0.20-44loaded via a namespace (and not attached): [1] Rcpp_1.0.6       lubridate_1.7.10 digest_0.6.27    assertthat_0.2.1 utf8_1.2.1       R6_2.5.0         cellranger_1.1.0 backports_1.2.1  reprex_2.0.0     httr_1.4.2       pillar_1.6.1    [12] rlang_0.4.11     readxl_1.3.1     rstudioapi_0.13  minqa_1.2.4      nloptr_1.2.2.2   Matrix_1.3-4     labeling_0.4.2   splines_4.1.0    lme4_1.1-27.1    munsell_0.5.0    broom_0.7.8     [23] compiler_4.1.0   Deriv_4.1.3      modelr_0.1.8     pkgconfig_2.0.3  tidyselect_1.1.1 codetools_0.2-18 fansi_0.5.0      crayon_1.4.1     dbplyr_2.1.1     withr_2.4.2      MASS_7.3-54     [34] grid_4.1.0       jsonlite_1.7.2   gtable_0.3.0     lifecycle_1.0.0  DBI_1.1.1        magrittr_2.0.1   scales_1.1.1     cli_3.0.0        stringi_1.6.2    farver_2.1.0     fs_1.5.0        [45] xml2_1.3.2       ellipsis_0.3.2   generics_0.1.0   vctrs_0.3.8      boot_1.3-28      tools_4.1.0      glue_1.4.2       hms_1.1.0        colorspace_2.0-2 rvest_1.0.0      haven_2.4.1   
