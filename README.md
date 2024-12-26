# CoxPH_Analysis_Simple
This code will take a dataset and run multiple univariate analyses. Then, it will take significant univariate analysis covariates and run a multivariate analysis.

## A Couple Notes:
- Significance can be modified but the default is p-values < 0.1.
- The code only accepts datasets that have columns in structure: Survival Time, Survival Events, Covariates. Otherwise, you will need to create a vector of your covariates.


## Sources:
Another Automated Cox Proportional Hazards Function that inspired this code: [Link](https://www.sthda.com/english/wiki/cox-proportional-hazards-model)
