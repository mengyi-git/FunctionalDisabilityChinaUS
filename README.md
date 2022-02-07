# Introduction

The code accompanies the paper [Functional Disability with Systematic Trends and Uncertainty: A Comparison between China and the U.S.] (https://doi.org/10.1017/S1748499521000233) by Fu, Sherris and Xu. We also provide the datasets used for model estimation.


# How to use the code
Use `main.m` for model estimation, simulation, and the likelihood ratio test. Note that to estimate the trend model, one needs to estimate the static model first since the initial values used in estimating the trend model depend on the estimated parameters of the static model. Similarly, the trend model needs to be estimated before estimating the frailty model. 

Use `main_plot.m` to plot figures in the paper.

Use `main_cf.m` to calculate the estimated transition rates in the following papers
  * Zixi Li, Adam W. Shao & Michael Sherris (2017) The impact of systematic trend and uncertainty on mortality and disability in a multistate latent factor model for transition rates, *North American Actuarial Journal*, 21(4), 594-610, [DOI: 10.1080/10920277.2017.1330157](https://doi.org/10.1080/10920277.2017.1330157)
  * Hanewald, K., Li, H., & Shao, A. (2019). Modelling multi-state health transitions in China: A generalised linear model with time trends. *Annals of Actuarial Science*, 13(1), 145-165. [doi:10.1017/S1748499518000167](https://www.cambridge.org/core/journals/annals-of-actuarial-science/article/modelling-multistate-health-transitions-in-china-a-generalised-linear-model-with-time-trends/93135D3F07A86F260D5F0B9A0B991634)
  * Michael Sherris & Pengyu Wei (2021) A multi-state model of functional disability and health status in the presence of systematic trend and uncertainty, *North American Actuarial Journal*, 25(1), 17-39. [DOI: 10.1080/10920277.2019.1708755](https://doi.org/10.1080/10920277.2019.1708755)


# Datasets

We use the Chinese Longitudinal Healthy Longevity Survey (CLHLS) and the U.S. Health and Retirement Study (HRS) to estimate the model parameters. The CLHLS data is downloaded from https://doi.org/10.3886/ICPSR36692.v1. It requires some data cleaning before it can be used. Please see https://github.com/mengyi-git/clhls_data_clean for the data cleaning procedure. The HRS data is downloaded from https://hrsdata.isr.umich.edu/data-products/rand. It has been cleaned by the RAND Center for the Study of Aging.

The datasets are saved in `clhls_transit.csv` and `rndhrs_transit.csv`, respectively. The meaning of the variables in each dataset is displayed below.

| Variable Name | Variable Meaning                                                                                   |        CLHLS       |         HRS        |
|---------------|----------------------------------------------------------------------------------------------------|:------------------:|:------------------:|
| `ID`          | Individual identifier                                                                              | :heavy_check_mark: |                    |
| `HHIDPN`      | Individual identifier                                                                              |                    | :heavy_check_mark: |
| `RAFEMALE`    | Gender (=1 if female)                                                                              | :heavy_check_mark: | :heavy_check_mark: |
| `JOINURBAN`   | Residence when joining the survey (=1 if urban, =0 if rural)                                       | :heavy_check_mark: |                    |
| `JOINWAVE`    | Wave when joining the survey                                                                       | :heavy_check_mark: |                    |
| `HACOHORT`    | Sample cohort                                                                                      |                    | :heavy_check_mark: |
| `RxHSTATE`    | Health   state at time t<sub>i</sub>                                                               | :heavy_check_mark: | :heavy_check_mark: |
| `RxHSTATE2`   | Health   state at time t<sub>i+1</sub>                                                             | :heavy_check_mark: | :heavy_check_mark: |
| `TIME`        | Time t<sub>i</sub>                                                                                 | :heavy_check_mark: | :heavy_check_mark: |
| `TAU`         | Duration (in year) between time t<sub>i</sub> and   t<sub>i+1</sub>                                | :heavy_check_mark: | :heavy_check_mark: |
| `RxAGE`       | Age last birthday at time t<sub>i</sub>                                                            | :heavy_check_mark: | :heavy_check_mark: |
| `RxAGE2`      | Age last birthday at time t<sub>i+1</sub>                                                          | :heavy_check_mark: | :heavy_check_mark: |
| `RxAGETRS`    | Age at which the transition occurs (=`RxAGE2` if transition occurs, =-1   if no transition occurs) | :heavy_check_mark: | :heavy_check_mark: |
| `TIMETRS`     | Time at which the transition occurs (=-1 if no transition occurs)                                  | :heavy_check_mark: | :heavy_check_mark: |
| `Y_S1`        | Transition indicator (=1 if transition type 1 is observed)                                         | :heavy_check_mark: | :heavy_check_mark: |
| `Y_S2`        | Transition indicator (=1 if transition type 2 is observed)                                         | :heavy_check_mark: | :heavy_check_mark: |
| `Y_S3`        | Transition indicator (=1 if transition type 3 is observed)                                         | :heavy_check_mark: | :heavy_check_mark: |
| `Y_S4`        | Transition indicator (=1 if transition type 4 is observed)                                         | :heavy_check_mark: | :heavy_check_mark: |
| `R_H`         | Exposure indicator (=1 if in the healthy state)                                                    | :heavy_check_mark: | :heavy_check_mark: |
| `R_D`         | Exposure indicator (=1 if in the disabled state)                                                   | :heavy_check_mark: | :heavy_check_mark: |

Transition type definition
  * 1: healthy to disabled
  * 2: disabled to healthy
  * 3: healthy to dead
  * 4: disabled to dead

Exposure indicator
  * If `R_H=1`, then the individual is exposed to the risk of transition types 1 and 3.
  * If `R_D=1`, then the individual is exposed to the risk of transition types 2 and 4.

Time t<sub>i</sub> is the latest of
  * date of reaching age x
  * date of joining the survey
  * date of observing health transitions.

Time t<sub>i+1</sub> is the earliest of
  * date of reaching age x+1
  * date of exiting the survey
  * date of observing health transitions.

Note that the definition of t<sub>i</sub> is different from that of t<sub>j</sub> in Fu et al. (2021).

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


# Questions and comments
If you have trouble running the code or have better ideas to improve the code, please [log an issue](https://github.com/mengyi-git/FunctionalDisabilityChinaUS/issues).
