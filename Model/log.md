
**Abundance**<br>
**Poisson**

poisson.stan
Max count per site
One covariate: age


**ZIP**

Simple abundance as a ZIP, no covariates: runs
Simple abundance as a ZIP, with age: runs
Simple abundance as a ZIP, with age & size: does not run
- Size is not normally distributed

Simple abundance as a ZIP, with age & size_log: runs!
Simple abundance as a ZIP, with age, size_log, percent_spruce: does not converge

**Occupancy, no detection**

with intercept, lat, long: runs

**Occupancy, with detection**
