# Data Generating

We use the deconvolution benchmark pipeline (Cobos et al,2020) to generate the simulation data. "helper_functions1.R" is adopted from https://github.com/favilaco/deconv_benchmark with some minor adjustments.


"getdata_pbmc.R" is used to generate the simulation data from the PBMC single cell data, including the bulk samples, cell-type specific gene expression profile and the ground true proportion.

Similar code is also used to generate the simulation data from Brain and Tumor single cell data.

"readydata_pbmc5.Rdata" is the simulation data used for deconvolution and results analysis.
