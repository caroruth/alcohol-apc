
import delimited "alcohol deaths.csv", encoding(UTF-8) clear 

drop if singleyearage > 74

**** ie ****

* Using Newton-Raphson optimization (the default in Stata) and scale(x2); (scale parameter = Pearson chi-squared / residual degrees of freedom)
apc_ie deaths, age(singleyearage) period(year) family(poisson) link(log) exposure(population) scale(x2)

* Using IRLS optimization (the default in S-Plus) and scale(dev); (scale parameter = deviance / residual degrees of freedom)
apc_ie deaths, age(singleyearage) period(year) family(poisson) link(log) exposure(population) scale(dev) irls




**** cglim ****

* Using Newton-Raphson optimization (the default in Stata) and scale(x2); (scale parameter = Pearson chi-squared / residual degrees of freedom)
apc_cglim deaths, age(singleyearage) period(year) agepfx("_a") periodpfx("_p") cohortpfx("_c") family(poisson) link(log) exposure(population) scale(x2) constraint("a17=a16")

* Using IRLS optimization (the default in S-Plus) and scale(dev); (scale parameter = deviance / residual degrees of freedom)
apc_cglim deaths, age(singleyearage) period(year) agepfx("_a") periodpfx("_p") cohortpfx("_c") family(poisson) link(log) exposure(population) scale(dev) irls constraint("a17=a16")

