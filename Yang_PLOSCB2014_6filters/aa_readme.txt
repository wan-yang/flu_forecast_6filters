## This package contains the code used for the manuscript: 
## Comparison of filtering methods for the modeling and retrospective forecasting of influenza epidemics (PLoS Compute Biol)
## by Wan Yang, Alicia Karspeck, and Jeffrey Shaman, 2014

## updated on Nov 3, 2014 (fixed error in calculating the standard deviations in PF, MIF, and pMCMC
## updated on July 27, 2015 (fixed error in numerical integration)

## We compared the performance of six filters for forecasting influenza outbreak in U.S. cities
## NOTE:  this comparison was done specifically for an influenza forecast system,
##        NOT a general comparison of the filter algorithms.
## The example data are Google Flu Trend ILI+ for New York City (NYC) 
##     from Week 40, 2003 to Week 47, 2012
## Climatology specific humidity data of NYC are used in the humidity-forcing SIRS model

## To run the forecast systems:
## Store the data files (.csv files) in a directory (dir_home_data)
## Store the functions (.R files) in a directory (dir_home_code)
## Specify the filter (fn.id=1,2…, 6), season (e.g.,season='11-12'), number of weeks for training (ntrn=4,…,26), 
##    number of particles/ensemble members (num_ens, 300 for the ensemble filters and 10,000 for the particle filters in our study) etc. in 'run_forecast.R'
## Run the forecast using the script 'run_forecast.R'

## Have fun!


