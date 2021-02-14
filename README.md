# Introduction

The code accompanies the paper "Functional Disability with Systematic Trends and Uncertainty: A Comparison between China and the U.S." by Fu, Sherris and Xu. We also provide the datasets (saved in `clhls_transit.csv` and `rndhrs_transit.csv`) used for model estimation.

# Datasets

We use the Chinese Longitudinal Healthy Longevity Survey (CLHLS) and the U.S. Health and Retirement Study (HRS) to estimate the model parameters.


# How to use the code
Use `main.m` for model estimation and simulation. Note that to estimate the trend model, one needs to estimate the static model first since the initial values used in estimating the trend model depend on the estimated parameters of the static model. Similarly, the trend model needs to be estimated before estimating the frailty model. 

Use `main_plot.m` to plot figures in the paper.

Use `main_cf.m` to calculate the estimated transition rates in the following papers
  * Zixi Li, Adam W. Shao & Michael Sherris (2017) The impact of systematic trend and uncertainty on mortality and disability in a multistate latent factor model for transition rates, *North American Actuarial Journal*, 21(4), 594-610, [DOI: 10.1080/10920277.2017.1330157](https://doi.org/10.1080/10920277.2017.1330157)
  * Hanewald, K., Li, H., & Shao, A. (2019). Modelling multi-state health transitions in China: A generalised linear model with time trends. *Annals of Actuarial Science*, 13(1), 145-165. [doi:10.1017/S1748499518000167](https://www.cambridge.org/core/journals/annals-of-actuarial-science/article/modelling-multistate-health-transitions-in-china-a-generalised-linear-model-with-time-trends/93135D3F07A86F260D5F0B9A0B991634)
  * Michael Sherris & Pengyu Wei (2020) A multi-state model of functional disability and health status in the presence of systematic trend and uncertainty, *North American Actuarial Journal*, [DOI: 10.1080/10920277.2019.1708755](https://doi.org/10.1080/10920277.2019.1708755)

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
