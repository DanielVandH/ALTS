# ALTS

This repository contains the R code required to reproduce the results in our paper <...>, which has abstract:

> Standard methods for forecasting electricity loads are not robust to cyberattacks on electricity demand
data, potentially leading to severe consequences such as major economic loss or a system blackout.
Methods are required that can handle forecasting under these conditions and detect outliers that would
otherwise go unnoticed. The key challenge is to remove as many outliers as possible while maintaining
enough clean data to use in the regression. In this paper we investigate robust approaches with data-driven tuning parameters, and in particular, present an adaptive trimmed regression method that can
better detect outliers and provide improved forecasts. In general, data-driven approaches perform much
better than their fixed tuning parameter counterparts. Recommendations for future work are provided.

The main contribution of this code is the extension of Bacher's ALTS [Bacher et al. (2016)](10.1109/ICASSP.2016.7472513) for performing least trimmed square regression with a data-driven tuning parameter. We apply this to the problem of forecasting electricity demand under cyberattacks. We also discuss other data-driven tuning parameters and compare them to fixed tuning parameter methods, and argue that data-driven methods outperform their fixed tuning parameter counterparts. 

The steps to reproduce the results in the paper are given below, and all the contents of this repository are described following this description. 

# Steps to Reproduce the Paper

The main two scripts that perform the work for the results in our paper are `Paper_Code_Data.R` and `Paper_Code_Figure.R`. (If you would only like to reproduce the figures, or just inspect the data sets without rerunning the time consuming algorithms in `Paper_Code_Data.R`, you can skip straight to the description of `Paper_Code_Figure.R` in the next paragraph which reloads the files below.) If you would like to reproduce all the experiments used in the paper, you can run all the code in `Paper_Code_Data.R`. The following files will be saved:
- `data/Simulation/Fixedq.RData`: In this file, we store the results from our simulation study fixed at $q = 1.35$ and for $p = 0.5, 0.55, \ldots, 1$. The name of the variable stored is `results_fixed_q`.
- `data/Simulation/Varyingq.RData`: In this file, we store the results from our simulation study for $q = 1, 1.7, \ldots, 1.35, 1.42, 1.49$, and $p = 0.55,0.6,\ldots,0.90, 0.95$. The name of the variable stored is `results_varying_q`.
- `data/CaseStudy/EconomicLoss.RData`: Results from our economic loss case study. The name of the variable stored is `results_economic_loss`.
- `data/CaseStudy/EconomicLossSmallMu.RData`: Results from our economic loss case study with small $\mu$ and large $\sigma$. The name of the variable stored is `results_economic_loss_small_mu`.
- `data/CaseStudy/SimulationBlackout.RData`: Results from our system blackout case study. The name of the variable stored is `results_system_blackout`.
- `data/CaseStudy/RampLambda.RData`: Results from our ramp attack case study using the scale parameter formulation. The name of the variable stored is `results_ramp_lambda`.
- `data/CaseStudy/RampGamma.RData`: Results from our ramp attack case study using the $\gamma$ formulation. The name of the variable stored is `results_ramp_gamma`.

Once these files are obtained, you can proceed onto `Paper_Code_Figure.R`, which we now describe. 

The code in `Paper_Code_Figure.R` takes all the data from `Paper_Code_Data.R` and makes the figures and tables given in the paper:
- Lines 64--99: This part of the code creates the tables that are used in the paper, in particular Table 2 (we combine the two tables in this code into the single table). This table compares our new ALTS to Bacher's ALTS.
- Lines 100--150: This part of the code reproduces the figure in Appendix C of the paper, which shows how our method's results depend on the choice of the hyperparameter $q$. 
- Lines 261--297: This code produces Figure 1 of the paper, showing how the different cyberattack templates may affect a given time series.
- Lines 298--338: This code produces Figure 2 of the paper, showing the results for an economic loss scenario under a range of parameter values.
- Lines 338--389: Similar to the above, except for Figure 3 and for the economic loss results with small $\mu$.
- Lines 390--429: Similar to the above, except for Figure 4 and for the system blackout results.
- Lines 430--469: Similar to the above, except for Figure 5 and for the ramp attack results under the scale parameter formulation.
- Lines 470--509: Similar to the above, except for Figure 6 and for the ramp attack results under the $\gamma$ formulation.
- Lines 510--516: This prints Table D.3 of the paper, which is the tabulated version of Figure 2.
- Lines 517--523: This prints Table D.4 of the paper, which is the tabulated version of Figure 3.
- Lines 524--529: This prints Table D.5 of the paper, which is the tabulated version of Figure 4.
- Lines 530--535: This prints Table D.6 of the paper, which is the tabulated version of Figure 5.
- Lines 536--543: This prints Table D.7 of the paper, which is the tabulated version of Figure 6.

# Contents

The repository contains the following files. Note that all of the functions mentioned are well documented within each script, should you require more information.

## Scripts of Functions

- _ALTS\_Functions.R_: This script contains functions for our new ALTS algorithm. In particular, we define the following functions:
    1. `modified_mad`: Computes our modified MAD estimator.
    2. `p_estimator`: Computes $p$, the estimated for the proportion of clean data, using our new method.
    3. `weights_calculation`: Computes the weights for the LTS algorithm.
    4. `ALTS`: Our new ALTS method.
    5. `ALTS_Bacher`: Bacher's ALTS, from [Bacher et al. (2016)](10.1109/ICASSP.2016.7472513).
- _Case\_Study\_Functions.R_: This script contains functions used for performing the case study in our paper. We define the following functions:
    1. `random_attack`: Randomly attacks given data.
    2. `ramp_attack`: Ramp attacks given data.
    3. `random_attack_ALTS`: Randomly attacks given data and then fits many models.
    4. `random_attack_all_simulations`: Randomly attacks given data and then fits many models, repeated many times.
    5. `random_attack_all_simulations_grid_vals`: Randomly attacks given data and then fits many models, repeated many times, and for many combinations of provided parameter values.
    6. `process_all_results`: Processes the results from the list of results from the simulation functions.
    7. `ramp_attack_ALTS`: Ramp attacks given data and then fits many models.
    8. `ramp_attack_all_simulations`: Ramp attacks given data and then fits many models, repeated many times.
    9. `ramp_attack_all_simulations_grid_vals`: Ramp attacks given data and then fits many models, repeated many times, and for many combinations of provided parameter values.
    10. `ramp_attack_all_simulations_grid_vals_2`: Ramp attacks given data and then fits many models, repeated many times, and for many combinations of provided parameter values, using the $\gamma = 1 + \ell\lambda/2$ formulation rather than the $\lambda$ formulation as in `ramp_attack_all_simulations_grid_vals`.
    11. `widen_results_and_latexify`: This widens the results for the case study data and then obtains the $\LaTeX$ code for the corresponding table.
- _Simulation\_Functions.R_:
    1. `%notin%`: An infix operator defining the negation of `%in%`.
    2. `sim_fnc_1`: Performs a single iteration for the simulation study.
    3. `process_results`: Processes the results from a number of simulations from the simulation study.
    4. `main_sim_fnc`: For a single $p$, perform the simulation study.
    5. `complete_sim_fnc`: For a vector $\boldsymbol p$, perform the simulation study for each $p_i$.
    6. `complete_sim_fnc_vary_q`: For a vector $\boldsymbol p$ and vector $\boldsymbol q$, perform the simulation study for each $p_i$ and $q_i$.

## Data

We provide many data files in this repository, either as Excel spreadsheets or R data files. These are given in the data folder. We list the contents of this data folder below, separated by type. The simulation study data and case study data files discussed below are produced by `Paper_Code_Data.R`.

### data/GEFCom2012

The data for the case study comes from the Global Energy Forecasting Competititon (GEFCom2012) from [Hong et al. (2014)](https://doi.org/10.1016/j.ijforecast.2013.07.001). This data is given as two Excel files:

- `Load_history.csv`: This data gives the hourly load history at 20 zones.
- `Temperature_history.csv`: This data gives the hourly temperature history at 11 weather stations.
More information on these datasets is given by [Hong et al. (2014)](https://doi.org/10.1016/j.ijforecast.2013.07.001), and in this blogpost by [Hong (2016)](http://blog.drhongtao.com/2016/07/gefcom2012-load-forecasting-data.html). We process this code into a convenient format for our application in the `Paper_Code_Data.R` under the "Read in the case study data" heading.

### data/Simulation

The results from our simulation study are saved as R data files. As described in our paper, this simulation study concerns the model 

$$
y_i = \beta_0+\beta_1x_{1i}+\beta_2x_{2i}+\beta_3x_{3i}+\epsilon_i,\quad i=1,2,\ldots,n, 
$$

where $\beta_0=-1.3$, $\beta_1=2.0$, $\beta_2=1.7$, and $\beta_3=-3.0$. We use $n=2000$ and take $x_1\sim \mathcal U(-1, 1)$, $x_2 \sim \mathcal N(0, 1)$, and $x_3 \sim \mathcal U(0, 1)$. For the residuals, we suppose that their distribution is given as the mixture

$$
\epsilon_i \sim \underbrace{p\mathcal N(0,\sigma_1)}\_{\text{clean data}} + \underbrace{(1-p)\mathcal N(0,\sigma_2)}\_{\text{contaminated data}},
$$

so that $100(1-p)\%$ of the data are outliers, and we take $\sigma_1=0.1$ and $\sigma_2=1.3$.

The simulation study is performed in `Paper_Code_Data.R` under the "Simulation study" heading, producing `data/Simulation/Fixedq.RData` and `data/Simulation/Varyingq.RData`. Please refer to the _Steps to Reproduce the Paper_ section above for a description of these files.

### data/CaseStudy

In addition to the GEFCom2012 data from [Hong et al. (2014)](https://doi.org/10.1016/j.ijforecast.2013.07.001) for our case study, we save the results that we show in the paper. The case study is performed in the sections after the simulation study section in `Paper_Code_Data.R`, and we produce those files listed in the _Steps to Reproduce the Paper_ section above. See our paper for more details.

## Figures 

The figures folder contains all the figures used in the paper.
