# Introduction


# Datasets


# How to use the code
Use `main.m` for model estimation and simulation.

Use `main_plot.m` to plot figures in the paper.

Use `main_cf.m` to calculate the estimated transition rates in the following papers
  * Hanewald, K., Li, H., & Shao, A. (2019). Modelling multi-state health transitions in China: A generalised linear model with time trends. Annals of Actuarial Science, 13(1), 145-165. doi:10.1017/S1748499518000167

# Acknowledgement
The following functions are from external sources.
* `artransform.m` is adapted from the Fortran code written by Jouni Helske.
  * https://github.com/helske/KFAS/blob/master/src/artransform.f95
* `logsumexp.m` is written by Tom Minka (c) Microsoft Corporation.
  * https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/39653/versions/1/previews/Shuang_EP_Distributed/utils/logsumexp.m/index.html
* `./plotFiles/ciplot.m` is written by Raymond Reynolds 
  * Raymond Reynolds (2021). Plot confidence intervals (https://www.mathworks.com/matlabcentral/fileexchange/13103-plot-confidence-intervals), MATLAB Central File Exchange. Retrieved February 14, 2021.
* `./plotFiles/rgb.m` is written by Chad A. Greene
  * Chad Greene (2021). Intuitive RGB color values from XKCD (https://www.mathworks.com/matlabcentral/fileexchange/46872-intuitive-rgb-color-values-from-xkcd), MATLAB Central File Exchange. Retrieved February 14, 2021.
