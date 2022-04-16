# ALTS

This repository contains the R code required to reproduce the results in our paper <...>. 

# Contents

The repository contains the following files. Note that all of the functions mentioned are well documented within each strip, should you require more information.

## Scripts of Functions

- _ALTS\_Functions.R_: This script contains functions for our new ALTS algorithm. In particular, we define the following functions:
    1. `modified_mad`: Computes our modified MAD estimator.
    2. `p_estimator`: Computes $p$, the estimated for the proportion of clean data, using our new method.
    3. `weights_calculation`: Computes the weights for the LTS algorithm.
    4. `ALTS`: Our new ALTS method.
    5. `ALTS_CV`: Our new ALTS method, with $q$ selected using cross-validation.
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
    10. `ramp_attack_all_simulations_grid_vals_2`: Ramp attacks given data and then fits many models, repeated many times, and for many combinations of provided parameter values, using the $\gamma=1+\ell\lambda/2$ formulation rather than the $\lambda$ formulation as in `ramp_attack_all_simulations_grid_vals`.
- _Simulation\_Functions.R_:
    1. `%notin%`: An infix operator defining the negation of `%in%`.
    2. `sim_fnc_1`: Performs a single iteration for the simulation study.
    3. `process_results`: Processes the results from a number of simulations from the simulation study.
    4. `main_sim_fnc`: For a single $p$, perform the simulation study.
    5. `complete_sim_fnc`: For a vector $\mathbf{p}$, perform the simulation study for each $p_i$.
    6. `complete_sim_fnc_vary_q`: For a vector $\mathbf{p}$ and vector $\mathbf{q}$, perform the simulation study for each $p_i$ and $q_j$.
    7. `sim_fnc_1_cv`: Performs a single iteration for the cross-validation simulation study.
    8. `main_sim_fnc_cv`: For a single $p$, perform the cross-validation simulation study.
    9. `complete_sim_fnc_cv`: For a vector $\mathbf{p}$, perform the cross-validation simulation study for each $p_i$.

## Data

We provide many data files in this repository, either as Excel spreadsheets or R data files. These are given in the data folder. We list the contents of this data folder below, separated by type. The simulation study data and case study data files discussed below are produced by `Paper_Code_Data.R`.

### GEFCom2012

The data for the case study comes from the Global Energy Forecasting Competititon (GEFCom2012) from [Hong et al. (2014)](https://doi.org/10.1016/j.ijforecast.2013.07.001). This data is given as two Excel files:
- `Load_history.csv`: This data gives the hourly load history at 20 zones.
- `Temperature_history.csv`: This data gives the hourly temperature history at 11 weather stations.
More information on these datasets is given by [Hong et al. (2014)](https://doi.org/10.1016/j.ijforecast.2013.07.001), and in this blogpost by [Hong (2016)](http://blog.drhongtao.com/2016/07/gefcom2012-load-forecasting-data.html). We process this code into a convenient format for our application in the `Paper_Code_Data.R` under the "Read in the case study data" heading.

### Simulation Study Data

The results from our simulation study are saved as R data files. As described in our paper, this simulation study concerns the model 
<<<<<<< HEAD
<!-- $$ 
y_i = \beta_0 + \beta_1x_{1i} + \beta_2x_{2i} + \beta_3x_{3i} + \epsilon_i, \quad i=1,\ldots,n, 
$$ --> 

<div align="center"><img style="background: white;" src="https://render.githubusercontent.com/render/math?math=y_i%20%3D%20%5Cbeta_0%20%2B%20%5Cbeta_1x_%7B1i%7D%20%2B%20%5Cbeta_2x_%7B2i%7D%20%2B%20%5Cbeta_3x_%7B3i%7D%20%2B%20%5Cepsilon_i%2C%20%5Cquad%20i%3D1%2C%5Cldots%2Cn%2C%20%0D"></div>
=======
<img src="https://bit.ly/3xxzrth" align="center" border="0" alt="y_i = \beta_0 + \beta_1x_{1i} + \beta_2x_{2i} + \beta_3x_{3i} + \epsilon_i, \quad i=1,\ldots,n," width="406" height="19" />

>>>>>>> 03acc1cbfe70e4d032b26594ae857a086d37f09b
where $\beta_0 = -1.3$, $\beta_1 = 2.0$, $\beta_2 = 1.7$, and $\beta_3 = -3.0$. We use $n = 2000$ and take $x_1 \sim \mathcal U(-1, 1)$, $x_2 \sim \mathcal N(0, 1)$, $x_3 \sim \mathcal U(0, 1)$. For the residuals, we suppose that their distribution is given as the mixture $$ \epsilon_i \sim \underbrace{p\mathcal N(0, \sigma_1)}_{\textnormal{clean data}} + \underbrace{(1-p)\mathcal N(0, \sigma_2)}_{\textnormal{contaminated data}}, $$
so that $100(1-p)\%$ of the data are outliers, and we take $\sigma_1 = 0.1$ and $\sigma_2=1.3$. The simulation study is performed in `Paper_Code_Data.R` under the "Simulation study" heading", producing the following files:
- `SimulationStudy_Fixedq.R`: In this file,
