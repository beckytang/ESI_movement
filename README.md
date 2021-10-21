# Species interactions and movement: modeling environmental effects on community dynamics

## Data: 

* derived_ebird_data.Rda: data were derived from eBird (1) as follows:
  *  Obtained and zero-filled data for the species listed in the manuscript main text using R package auk (2)
  * Filtered data for spatial and temporal windows as described in the manuscript main text
  * For each areal unit BLOB as described in the manuscript, obtained an aggregated count per area for each species using the method as described in the Supplementary Information S6. These are denoted as 'Ypp' in the data frame.

  * Environmental covariates of dominant land cover type, Enhanced Vegetation Index (EVI), and habitat edge length are also included in the data frame

* raw_data.Rda: list containing for both regions (Research Triangle and Arlington) lists of posterior samples of the following parameters. These data were used to generate the tables and figures in the main text:
  *  $\alpha$
  *  $\rho$
  *  $\delta$
  *  $\beta$
  *  $\sigma^{2}_{eta}$
  *  $\sigma^{2}_{gamma}$

## Code: 

* sampler.R : contains code to fit the model as described in the manuscript, along with definitions of modeler-defined inputs
* functions.R: R file containing functions called in the sampler
* cpp_fns.cpp: C++ file containing functions called in the sampler, in addition to functions.R 

## 




(1) Sullivan, B.L., C.L. Wood, M.J. Iliff, R.E. Bonney, D. Fink, and S. Kelling. 2009. eBird: a citizen-based bird observation network in the biological sciences.  Biological Conservation 142: 2282-2292.

(2) Matthew Strimas-Mackey, Eliot Miller, and Wesley Hochachka (2018). auk: eBird Data Extraction and Processing with AWK. R package version 0.3.0. https://cornelllabofornithology.github.io/auk/
