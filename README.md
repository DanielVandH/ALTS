# ALTS

This repository contains the R code required to reproduce the results in our paper <...>. The steps to reproduce the results in the paper are given below, and all the contents of this repository are described following this description.

# Steps to Reproduce the Paper

The main two scripts that perform the work for the results in our paper are `Paper_Code_Data.R` and `Paper_Code_Figure.R`. (If you would only like to reproduce the figures, or just inspect the data sets without rerunning the time consuming algorithms in `Paper_Code_Data.R`, you can skip straight to the description of `Paper_Code_Figure.R` in the next paragraph which reloads the files below.) If you would like to reproduce all the experiments used in the paper, you can run all the code in `Paper_Code_Data.R`. The following files will be saved:
- `SimulationStudy_Fixedq.RData`: In this file, we store the results from our simulation study fixed at <!-- $q = 1.35$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=q%20%3D%201.35"> and for <!-- $p = 0.5,0.51,\ldots,1$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=p%20%3D%200.5%2C0.51%2C%5Cldots%2C1">. The name of the variable stored is `results_fixed_q`.
- `SimulationStudy_Varyingq.RData`: In this file, we store the results from our simulation study for <!-- $q = 1,1.7,\ldots,1.35,1.42,1.49$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=q%20%3D%201%2C1.7%2C%5Cldots%2C1.35%2C1.42%2C1.49"> and <!-- $p = 0.55,0.6,\ldots,0.90,0.95$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=p%20%3D%200.55%2C0.6%2C%5Cldots%2C0.90%2C0.95">. The name of the variable stored is `results_varying_q`.
- `SimulationStudy_CVq.RData`: In this file, we store the results from our simulation study with <!-- $q$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=q"> chosen using cross-validation and for <!-- $p = 0.55,0.6,\ldots,0.90,0.95$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=p%20%3D%200.55%2C0.6%2C%5Cldots%2C0.90%2C0.95">. The name of the variable stored is `results_cv`.
- `CaseStudy_EconomicLoss.RData`: Results from our economic loss case study. The name of the variable stored is `results_economic_loss`.
- `CaseStudy_SimulationBlackout.RData`: Results from our system blackout case study. The name of the variable stored is `results_system_blackout`.
- `CaseStudy_RampLambda.RData`: Results from our ramp attack case study using the scale parameter formulation. The name of the variable stored is `results_ramp_lambda`.
- `CaseStudy_RampGamma.RData`: Results from our ramp attack case study using the <!-- $\gamma$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cgamma"> formulation. The name of the variable stored is `results_ramp_gamma`.

Once these files are obtained, you can proceed onto `Paper_Code_Figure.R`, which we now describe. 

The code in `Paper_Code_Figure.R` takes all the data from `Paper_Code_Data.R` and makes the figures given in the below. The following lines produce the corresponding figures:
- _Figure 1_: Lines

# Contents

The repository contains the following files. Note that all of the functions mentioned are well documented within each strip, should you require more information.

## Scripts of Functions

- _ALTS\_Functions.R_: This script contains functions for our new ALTS algorithm. In particular, we define the following functions:
    1. `modified_mad`: Computes our modified MAD estimator.
    2. `p_estimator`: Computes <!-- $p$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=p">, the estimated for the proportion of clean data, using our new method.
    3. `weights_calculation`: Computes the weights for the LTS algorithm.
    4. `ALTS`: Our new ALTS method.
    5. `ALTS_CV`: Our new ALTS method, with <!-- $q$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=q"> selected using cross-validation.
    6. `ALTS_Bacher`: Bacher's ALTS, from [Bacher et al. (2016)](10.1109/ICASSP.2016.7472513).
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
    10. `ramp_attack_all_simulations_grid_vals_2`: Ramp attacks given data and then fits many models, repeated many times, and for many combinations of provided parameter values, using the <!-- $\gamma=1+\ell\lambda/2$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cgamma%3D1%2B%5Cell%5Clambda%2F2"> formulation rather than the <!-- $\lambda$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Clambda"> formulation as in `ramp_attack_all_simulations_grid_vals`.
- _Simulation\_Functions.R_:
    1. `%notin%`: An infix operator defining the negation of `%in%`.
    2. `sim_fnc_1`: Performs a single iteration for the simulation study.
    3. `process_results`: Processes the results from a number of simulations from the simulation study.
    4. `main_sim_fnc`: For a single <!-- $p$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=p">, perform the simulation study.
    5. `complete_sim_fnc`: For a vector <!-- $\mathbf{p}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7Bp%7D">, perform the simulation study for each <!-- $p_i$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=p_i">.
    6. `complete_sim_fnc_vary_q`: For a vector <!-- $\mathbf{p}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7Bp%7D"> and vector <!-- $\mathbf{q}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7Bq%7D">, perform the simulation study for each <!-- $p_i$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=p_i"> and <!-- $q_j$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=q_j">.
    7. `sim_fnc_1_cv`: Performs a single iteration for the cross-validation simulation study.
    8. `main_sim_fnc_cv`: For a single <!-- $p$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=p">, perform the cross-validation simulation study.
    9. `complete_sim_fnc_cv`: For a vector <!-- $\mathbf{p}$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7Bp%7D">, perform the cross-validation simulation study for each <!-- $p_i$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=p_i">.

## Data

We provide many data files in this repository, either as Excel spreadsheets or R data files. These are given in the data folder. We list the contents of this data folder below, separated by type. The simulation study data and case study data files discussed below are produced by `Paper_Code_Data.R`.

### GEFCom2012

The data for the case study comes from the Global Energy Forecasting Competititon (GEFCom2012) from [Hong et al. (2014)](https://doi.org/10.1016/j.ijforecast.2013.07.001). This data is given as two Excel files:

- `Load_history.csv`: This data gives the hourly load history at 20 zones.
- `Temperature_history.csv`: This data gives the hourly temperature history at 11 weather stations.
More information on these datasets is given by [Hong et al. (2014)](https://doi.org/10.1016/j.ijforecast.2013.07.001), and in this blogpost by [Hong (2016)](http://blog.drhongtao.com/2016/07/gefcom2012-load-forecasting-data.html). We process this code into a convenient format for our application in the `Paper_Code_Data.R` under the "Read in the case study data" heading.

### Simulation Study Data

The results from our simulation study are saved as R data files. As described in our paper, this simulation study concerns the model 
<!-- $$ 
y_i = \beta_0 + \beta_1x_{1i} + \beta_2x_{2i} + \beta_3x_{3i} + \epsilon_i, \quad i=1,\ldots,n, 
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=y_i%20%3D%20%5Cbeta_0%20%2B%20%5Cbeta_1x_%7B1i%7D%20%2B%20%5Cbeta_2x_%7B2i%7D%20%2B%20%5Cbeta_3x_%7B3i%7D%20%2B%20%5Cepsilon_i%2C%20%5Cquad%20i%3D1%2C%5Cldots%2Cn%2C%20%0D"></div>

where <!-- $\beta_0 = -1.3$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cbeta_0%20%3D%20-1.3">, <!-- $\beta_1 = 2.0$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cbeta_1%20%3D%202.0">, <!-- $\beta_2 = 1.7$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cbeta_2%20%3D%201.7">, and <!-- $\beta_3 = -3.0$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cbeta_3%20%3D%20-3.0">. We use <!-- $n = 2000$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=n%20%3D%202000"> and take <!-- $x_1 \sim \mathcal U(-1, 1)$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=x_1%20%5Csim%20%5Cmathcal%20U(-1%2C%201)">, <!-- $x_2 \sim \mathcal N(0, 1)$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=x_2%20%5Csim%20%5Cmathcal%20N(0%2C%201)">, <!-- $x_3 \sim \mathcal U(0, 1)$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=x_3%20%5Csim%20%5Cmathcal%20U(0%2C%201)">. For the residuals, we suppose that their distribution is given as the mixture 
<!-- $$ 
\epsilon_i \sim \underbrace{p\mathcal N(0, \sigma_1)}_{\text{clean data}} + \underbrace{(1-p)\mathcal N(0, \sigma_2)}_{\text{contaminated data}}, 
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cepsilon_i%20%5Csim%20%5Cunderbrace%7Bp%5Cmathcal%20N(0%2C%20%5Csigma_1)%7D_%7B%5Ctext%7Bclean%20data%7D%7D%20%2B%20%5Cunderbrace%7B(1-p)%5Cmathcal%20N(0%2C%20%5Csigma_2)%7D_%7B%5Ctext%7Bcontaminated%20data%7D%7D%2C%20%0D"></div>

so that <!-- $100(1-p)\%$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=100(1-p)%5C%25"> of the data are outliers, and we take <!-- $\sigma_1 = 0.1$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Csigma_1%20%3D%200.1"> and <!-- $\sigma_2=1.3$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Csigma_2%3D1.3">. The simulation study is performed in `Paper_Code_Data.R` under the "Simulation study" heading, producing the following files:

- `SimulationStudy_Fixedq.RData`: In this file, we store the results from our simulation study fixed at <!-- $q = 1.35$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=q%20%3D%201.35"> and for <!-- $p = 0.5,0.51,\ldots,1$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=p%20%3D%200.5%2C0.51%2C%5Cldots%2C1">. The name of the variable stored is `results_fixed_q`.
- `SimulationStudy_Varyingq.RData`: In this file, we store the results from our simulation study for <!-- $q = 1,1.7,\ldots,1.35,1.42,1.49$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=q%20%3D%201%2C1.7%2C%5Cldots%2C1.35%2C1.42%2C1.49"> and <!-- $p = 0.55,0.6,\ldots,0.90,0.95$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=p%20%3D%200.55%2C0.6%2C%5Cldots%2C0.90%2C0.95">. The name of the variable stored is `results_varying_q`.
- `SimulationStudy_CVq.RData`: In this file, we store the results from our simulation study with <!-- $q$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=q"> chosen using cross-validation and for <!-- $p = 0.55,0.6,\ldots,0.90,0.95$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=p%20%3D%200.55%2C0.6%2C%5Cldots%2C0.90%2C0.95">. The name of the variable stored is `results_cv`.

### Case Study Data

In addition to the GEFCom2012 data from [Hong et al. (2014)](https://doi.org/10.1016/j.ijforecast.2013.07.001) for our case study, we save the results that we show in the paper. See our paper for more details. The case study is performed in the sections after the simulation study section in `Paper_Code_Data.R`, and we produce the following files:
- `CaseStudy_EconomicLoss.RData`: Results from our economic loss case study. The name of the variable stored is `results_economic_loss`.
- `CaseStudy_SimulationBlackout.RData`: Results from our system blackout case study. The name of the variable stored is `results_system_blackout`.
- `CaseStudy_RampLambda.RData`: Results from our ramp attack case study using the scale parameter formulation. The name of the variable stored is `results_ramp_lambda`.
- `CaseStudy_RampGamma.RData`: Results from our ramp attack case study using the <!-- $\gamma$ --> <img style="transform: translateY(0.1em); background: white;" src="https://render.githubusercontent.com/render/math?math=%5Cgamma"> formulation. The name of the variable stored is `results_ramp_gamma`.

## Figures 

The final file included in this repository is the `Paper_Code_Figure.R`. This script reads in all the previously discussed data files and produces the figures in the paper.
