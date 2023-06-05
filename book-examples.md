# Replication


In this article, the functionality of the commands **xtsf1g** and **xtsf2g** is
showcased. At the same time, the code provided below can be used to replicate the results in the chapter:

Stochastic frontier analysis in Stata: using existing and coding new
commands in “Handbook of Research Methods and Applications in Empirical
Microeconomics” (edited by Nigar Hashimzade and Michael A. Thornton),
2021, Chapter 17, ***Edward Elgar Publishing***, [DOI
<img src="man/figures/doi.png"  width="12" height="12">](https://doi.org/10.4337/9781788976480.00027)

Throughout the article, the subset of the banking data (which can be found [here](https://github.com/OlegBadunenko/sf-panel-1st-2nd-gen/tree/main/data)) will be used:

``` stata
use ../../banks00_07, clear
```

# Define specifications

Here are the specification and formulas for the first derivatives of the cost function with respect to input prices and outputs to check monotonicity assumptions and compute returns to scale are defined:

``` stata
global spetech "lny1 lny2 lnw1 trend c.half#(c.lny1#c.lny1 c.lny2#c.lny2 c.lnw1#c.lnw1) c.lny1#(c.lny2 c.lnw1) c.lny2#c.lnw1 c.trend#(c.lny1 c.lny2 c.lnw1 c.trend#c.half)"
global monotonW1 = "_b[lnw1] + _b[c.half#c.lnw1#c.lnw1] *lnw1 + _b[c.lny1#c.lnw1] *lny1 + _b[c.lny2#c.lnw1] *lny2 + _b[c.trend#c.lnw1] *trend"
global monotonY1 = "_b[lny1] + _b[c.half#c.lny1#c.lny1] *lny1 + _b[c.lny1#c.lny2] *lny2 + _b[c.lny1#c.lnw1] *lnw1 + _b[c.trend#c.lny1] *trend"
global monotonY2 = "_b[lny2] + _b[c.half#c.lny2#c.lny2] *lny2 + _b[c.lny1#c.lny2] *lny1 + _b[c.lny2#c.lnw1] *lnw1 + _b[c.trend#c.lny2] *trend"
global year_c = 2001
global itermax = 1000
global decimalopt "cformat(%9.3f) pformat(%5.2f) sformat(%8.2f)"
```

# No determinants of inefficiency

In the first 8 models the specification models don't include determinants of inefficiency.

The inefficiency follow the $u_{it} = \beta(t) u_{i}$. In the 4 models below, $\beta(t)$ is defined as follows:

1. $\beta(t) = 1$ (time-invariant inefficiency)
2. $\beta(t) = ( 1 + \exp(\gamma t + \delta t^2) )^-1$ (time-variant inefficiency)
3. $\beta(t) = 1 + \gamma (t-T_i) + \delta (t-T_i)^2$ (time-variant inefficiency)
4. $\beta(t) = \exp( -\gamma (t-T_i) )$ (time-variant inefficiency)

# Models 1-4: Normal-half normal

Here are the models with normal-half normal assumptions.

## Model 1 (1st generation, time-invariant inefficiency)

We start by the simplest one where $u_{it} = u_{i}$.

Code for Model 1 is detailed. For the rest of the models, fewer details
are shown

``` stata
timer clear 1
timer on 1
xtsf1g lnc $spetech if year > $year_c, cost iterate($itermax) $decimalopt
timer off 1
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M1
quietly predictnl M1_mW1 = $monotonW1 if e(sample)
quietly predictnl M1_mW2 = 1-M1_mW1 if e(sample)
quietly predictnl M1_eY1 = $monotonY1 if e(sample)
quietly predictnl M1_eY2 = $monotonY2 if e(sample)
quietly generate double M1_rts = 1/(M1_eY1 + M1_eY2)
tabstat M1_*, stat(min p25 mean p75 max sd) format(%9.4f)

. timer clear 1

. timer on 1

. xtsf1g lnc $spetech if year > $year_c, cost iterate($itermax) $decimalopt


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =  685.89191  
iter 1:  Log-likelihood =  688.21588  
iter 2:  Log-likelihood =  1222.8627  
iter 3:  Log-likelihood =  1428.5402  
iter 4:  Log-likelihood =  1430.2777  
iter 5:  Log-likelihood =  1430.2788  
iter 6:  Log-likelihood =  1430.2788  

Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.9987
 Adj R-squared    = 0.9984
 AIC              = -2.3400
 BIC              = -2.3033
 Root MSE         = 0.1871
-----------------------------

 The first-generation estimator of
 the cost stochastic frontier model for panel data,
 where effciency is time-invariant

--------------------------------------------------------------------------------------
                 lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
---------------------+----------------------------------------------------------------
                lny1 |      0.397      0.142     2.79    0.01        0.118       0.675
                lny2 |     -0.619      0.399    -1.55    0.12       -1.401       0.163
                lnw1 |      0.385      0.166     2.32    0.02        0.059       0.711
               trend |     -0.561      0.040   -14.00    0.00       -0.640      -0.483
                     |
c.half#c.lny1#c.lny1 |      0.049      0.005    10.29    0.00        0.040       0.059
                     |
c.half#c.lny2#c.lny2 |      0.178      0.032     5.54    0.00        0.115       0.241
                     |
c.half#c.lnw1#c.lnw1 |     -0.069      0.014    -4.96    0.00       -0.096      -0.042
                     |
       c.lny1#c.lny2 |     -0.067      0.011    -5.86    0.00       -0.089      -0.045
                     |
       c.lny1#c.lnw1 |     -0.003      0.006    -0.45    0.65       -0.013       0.008
                     |
       c.lny2#c.lnw1 |     -0.009      0.013    -0.68    0.50       -0.035       0.017
                     |
      c.trend#c.lny1 |      0.005      0.002     3.11    0.00        0.002       0.009
                     |
      c.trend#c.lny2 |      0.012      0.003     3.95    0.00        0.006       0.018
                     |
      c.trend#c.lnw1 |     -0.007      0.003    -2.60    0.01       -0.012      -0.002
                     |
     c.trend#c.trend#|
              c.half |      0.076      0.002    42.70    0.00        0.073       0.080
                     |
               _cons |      5.383      2.642     2.04    0.04        0.205      10.561
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
               _cons |     -4.422      0.032  -138.37    0.00       -4.485      -4.359
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
               _cons |     -2.522      0.084   -30.11    0.00       -2.686      -2.358
--------------------------------------------------------------------------------------

. timer off 1

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -2.3400425

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -2.3033274

. eststo M1

. quietly predictnl M1_mW1 = $monotonW1 if e(sample)

. quietly predictnl M1_mW2 = 1-M1_mW1 if e(sample)

. quietly predictnl M1_eY1 = $monotonY1 if e(sample)

. quietly predictnl M1_eY2 = $monotonY2 if e(sample)

. quietly generate double M1_rts = 1/(M1_eY1 + M1_eY2)

. tabstat M1_*, stat(min p25 mean p75 max sd) format(%9.4f)

   Stats |    M1_mW1    M1_mW2    M1_eY1    M1_eY2    M1_rts
---------+--------------------------------------------------
     Min |   -0.1766    0.8766   -0.1627    0.3993    0.8636
     p25 |   -0.0165    0.9660    0.1332    0.6781    1.0478
    Mean |    0.0067    0.9933    0.1615    0.7475    1.1062
     p75 |    0.0340    1.0165    0.1942    0.8131    1.1583
     Max |    0.1234    1.1766    0.3200    1.2801    1.4827
      SD |    0.0408    0.0408    0.0537    0.1074    0.0824
------------------------------------------------------------
```

## Model 2 (2nd generation, time-varying inefficiency using the Kumbhakar (1990) specification)

This specification is by [Kumbhakar (1990)](https://doi.org/10.1016/0304-4076(90)90055-X)

``` stata
timer clear 2
timer on 2
xtsf2g lnc $spetech if year > $year_c, cost iterate($itermax) $decimalopt
timer off 2
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M2
quietly predictnl M2_mW1 = $monotonW1 if e(sample)
quietly predictnl M2_mW2 = 1-M2_mW1 if e(sample)
quietly predictnl M2_eY1 = $monotonY1 if e(sample)
quietly predictnl M2_eY2 = $monotonY2 if e(sample)
quietly generate double M2_rts = 1/(M2_eY1 + M2_eY2)
tabstat M2_*, stat(min p25 mean p75 max sd) format(%9.4f)

. timer clear 2

. timer on 2

. xtsf2g lnc $spetech if year > $year_c, cost iterate($itermax) $decimalopt


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =  716.15322  (not concave)
iter 1:  Log-likelihood =   1171.048  (not concave)
iter 2:  Log-likelihood =   1300.934  (not concave)
iter 3:  Log-likelihood =   1372.054  
iter 4:  Log-likelihood =  1444.2813  
iter 5:  Log-likelihood =   1454.265  
iter 6:  Log-likelihood =  1480.6447  
iter 7:  Log-likelihood =  1489.7689  
iter 8:  Log-likelihood =  1491.8923  
iter 9:  Log-likelihood =  1492.2628  
iter 10: Log-likelihood =  1492.3095  
iter 11: Log-likelihood =  1492.3115  
iter 12: Log-likelihood =  1492.3115  

Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.9986
 Adj R-squared    = 0.9983
 AIC              = -2.3091
 BIC              = -2.2724
 Root MSE         = 0.1900
-----------------------------

 The second-generation estimator of
 the cost stochastic frontier model for panel data,
 where effciency is time-varying

--------------------------------------------------------------------------------------
                 lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
---------------------+----------------------------------------------------------------
                lny1 |      0.429      0.140     3.06    0.00        0.154       0.703
                lny2 |     -0.483      0.393    -1.23    0.22       -1.254       0.287
                lnw1 |      0.497      0.162     3.06    0.00        0.179       0.815
               trend |     -0.656      0.040   -16.25    0.00       -0.735      -0.577
                     |
c.half#c.lny1#c.lny1 |      0.047      0.005    10.17    0.00        0.038       0.056
                     |
c.half#c.lny2#c.lny2 |      0.171      0.031     5.45    0.00        0.110       0.233
                     |
c.half#c.lnw1#c.lnw1 |     -0.062      0.013    -4.60    0.00       -0.088      -0.035
                     |
       c.lny1#c.lny2 |     -0.069      0.011    -6.14    0.00       -0.090      -0.047
                     |
       c.lny1#c.lnw1 |      0.000      0.005     0.02    0.98       -0.010       0.011
                     |
       c.lny2#c.lnw1 |     -0.024      0.013    -1.80    0.07       -0.049       0.002
                     |
      c.trend#c.lny1 |      0.005      0.002     3.23    0.00        0.002       0.009
                     |
      c.trend#c.lny2 |      0.013      0.003     4.36    0.00        0.007       0.019
                     |
      c.trend#c.lnw1 |     -0.006      0.003    -2.29    0.02       -0.011      -0.001
                     |
     c.trend#c.trend#|
              c.half |      0.094      0.003    32.03    0.00        0.088       0.099
                     |
               _cons |      4.497      2.617     1.72    0.09       -0.633       9.627
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
               _cons |     -4.498      0.032  -139.00    0.00       -4.561      -4.434
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
               _cons |     -2.315      0.089   -26.04    0.00       -2.489      -2.141
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |     -5.000      1.898    -2.63    0.01       -8.721      -1.280
               delta |      0.819      0.315     2.60    0.01        0.203       1.436
--------------------------------------------------------------------------------------
Note, the function of inefficiency change over time is
beta[t] = ( 1 + exp(gamma*t + delta*t^2) )^-1: 
 the larger is the beta[t] the larger is the inefficiency

. timer off 2

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -2.3090757

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -2.2723607

. eststo M2

. quietly predictnl M2_mW1 = $monotonW1 if e(sample)

. quietly predictnl M2_mW2 = 1-M2_mW1 if e(sample)

. quietly predictnl M2_eY1 = $monotonY1 if e(sample)

. quietly predictnl M2_eY2 = $monotonY2 if e(sample)

. quietly generate double M2_rts = 1/(M2_eY1 + M2_eY2)

. tabstat M2_*, stat(min p25 mean p75 max sd) format(%9.4f)

   Stats |    M2_mW1    M2_mW2    M2_eY1    M2_eY2    M2_rts
---------+--------------------------------------------------
     Min |   -0.1651    0.8692   -0.1501    0.4225    0.8582
     p25 |   -0.0112    0.9664    0.1276    0.6782    1.0534
    Mean |    0.0091    0.9909    0.1561    0.7474    1.1129
     p75 |    0.0336    1.0112    0.1882    0.8137    1.1644
     Max |    0.1308    1.1651    0.3162    1.2823    1.4384
      SD |    0.0371    0.0371    0.0525    0.1077    0.0834
------------------------------------------------------------
```

## Model 3 (2nd generation, time-varying inefficiency using the modified Kumbhakar (1990) specification)


``` stata
timer clear 3
timer on 3
xtsf2g lnc $spetech if year > $year_c, model(K1990modified) cost iterate($itermax) $decimalopt
timer off 3
predict M3_te, te_jlms_mean
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M3
quietly predictnl M3_mW1 = $monotonW1 if e(sample)
quietly predictnl M3_mW2 = 1-M3_mW1 if e(sample)
quietly predictnl M3_eY1 = $monotonY1 if e(sample)
quietly predictnl M3_eY2 = $monotonY2 if e(sample)
generate double M3_rts = 1/(M3_eY1 + M3_eY2)
tabstat M3_*, stat(min p25 mean p75 max sd) format(%9.4f)

. timer clear 3

. timer on 3

. xtsf2g lnc $spetech if year > $year_c, model(K1990modified) cost iterate($itermax) $decimalopt


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =  685.89191  
iter 1:  Log-likelihood =  1134.4188  (not concave)
iter 2:  Log-likelihood =  1295.1247  (not concave)
iter 3:  Log-likelihood =  1327.2692  (not concave)
iter 4:  Log-likelihood =  1337.2619  (not concave)
iter 5:  Log-likelihood =  1346.7248  (not concave)
iter 6:  Log-likelihood =  1353.5166  (not concave)
iter 7:  Log-likelihood =   1360.746  (not concave)
iter 8:  Log-likelihood =  1377.9014  
iter 9:  Log-likelihood =  1419.4682  
iter 10: Log-likelihood =  1443.3139  
iter 11: Log-likelihood =  1443.6859  
iter 12: Log-likelihood =  1443.6866  
iter 13: Log-likelihood =  1443.6866  

Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.9987
 Adj R-squared    = 0.9984
 AIC              = -2.3335
 BIC              = -2.2968
 Root MSE         = 0.1877
-----------------------------

 The second-generation estimator of
 the cost stochastic frontier model for panel data,
 where effciency is time-varying

--------------------------------------------------------------------------------------
                 lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
---------------------+----------------------------------------------------------------
                lny1 |      0.475      0.142     3.35    0.00        0.197       0.753
                lny2 |     -0.563      0.394    -1.43    0.15       -1.336       0.210
                lnw1 |      0.393      0.165     2.39    0.02        0.071       0.716
               trend |     -0.592      0.041   -14.62    0.00       -0.672      -0.513
                     |
c.half#c.lny1#c.lny1 |      0.047      0.005    10.12    0.00        0.038       0.057
                     |
c.half#c.lny2#c.lny2 |      0.180      0.032     5.68    0.00        0.118       0.242
                     |
c.half#c.lnw1#c.lnw1 |     -0.066      0.014    -4.82    0.00       -0.093      -0.039
                     |
       c.lny1#c.lny2 |     -0.073      0.011    -6.42    0.00       -0.095      -0.050
                     |
       c.lny1#c.lnw1 |      0.000      0.005     0.01    0.99       -0.011       0.011
                     |
       c.lny2#c.lnw1 |     -0.013      0.013    -0.98    0.33       -0.039       0.013
                     |
      c.trend#c.lny1 |      0.005      0.002     2.75    0.01        0.001       0.008
                     |
      c.trend#c.lny2 |      0.013      0.003     4.26    0.00        0.007       0.019
                     |
      c.trend#c.lnw1 |     -0.006      0.003    -2.35    0.02       -0.011      -0.001
                     |
     c.trend#c.trend#|
              c.half |      0.082      0.002    35.30    0.00        0.078       0.087
                     |
               _cons |      4.705      2.631     1.79    0.07       -0.450       9.861
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
               _cons |     -4.435      0.032  -138.89    0.00       -4.498      -4.373
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
               _cons |     -2.757      0.100   -27.63    0.00       -2.953      -2.561
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |     -0.125      0.028    -4.54    0.00       -0.179      -0.071
               delta |     -0.020      0.005    -3.77    0.00       -0.030      -0.009
--------------------------------------------------------------------------------------
Note, the function of inefficiency change over time is
beta[t] = 1 + gamma*(t-T_i) + delta*(t-T_i)^2: 
 the larger is the beta[t] the larger is the inefficiency

. timer off 3

. predict M3_te, te_jlms_mean
 (872 missing values generated)

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -2.3335014

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -2.2967864

. eststo M3

. quietly predictnl M3_mW1 = $monotonW1 if e(sample)

. quietly predictnl M3_mW2 = 1-M3_mW1 if e(sample)

. quietly predictnl M3_eY1 = $monotonY1 if e(sample)

. quietly predictnl M3_eY2 = $monotonY2 if e(sample)

. generate double M3_rts = 1/(M3_eY1 + M3_eY2)
(872 missing values generated)

. tabstat M3_*, stat(min p25 mean p75 max sd) format(%9.4f)

   Stats |     M3_te    M3_mW1    M3_mW2    M3_eY1    M3_eY2    M3_rts
---------+------------------------------------------------------------
     Min |    1.0071   -0.1710    0.8828   -0.1547    0.4048    0.8383
     p25 |    1.1348   -0.0143    0.9658    0.1307    0.6817    1.0432
    Mean |    1.2764    0.0076    0.9924    0.1597    0.7534    1.1011
     p75 |    1.3506    0.0342    1.0143    0.1928    0.8209    1.1513
     Max |    2.8145    0.1172    1.1710    0.3257    1.3157    1.4546
      SD |    0.2080    0.0391    0.0391    0.0540    0.1114    0.0818
----------------------------------------------------------------------
```

## Model 4 (2nd generation, time-varying inefficiency using the BC specification)

``` stata
timer clear 4
timer on 4
xtsf2g lnc $spetech if year > $year_c, model(BC1992) cost iterate($itermax) $decimalopt
timer off 4
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M4
quietly predictnl M4_mW1 = $monotonW1 if e(sample)
quietly predictnl M4_mW2 = 1-M4_mW1 if e(sample)
quietly predictnl M4_eY1 = $monotonY1 if e(sample)
quietly predictnl M4_eY2 = $monotonY2 if e(sample)
generate double M4_rts = 1/(M4_eY1 + M4_eY2)
tabstat M4_*, stat(min p25 mean p75 max sd) format(%9.4f)

. timer clear 4

. timer on 4

. xtsf2g lnc $spetech if year > $year_c, model(BC1992) cost iterate($itermax) $decimalopt


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =  685.89191  
iter 1:  Log-likelihood =  1198.9156  
iter 2:  Log-likelihood =  1400.7342  
iter 3:  Log-likelihood =  1434.2974  
iter 4:  Log-likelihood =  1435.2554  
iter 5:  Log-likelihood =  1435.2562  
iter 6:  Log-likelihood =  1435.2562  

Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.9987
 Adj R-squared    = 0.9984
 AIC              = -2.3349
 BIC              = -2.2982
 Root MSE         = 0.1876
-----------------------------

 The second-generation estimator of
 the cost stochastic frontier model for panel data,
 where effciency is time-varying

--------------------------------------------------------------------------------------
                 lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
---------------------+----------------------------------------------------------------
                lny1 |      0.459      0.143     3.22    0.00        0.180       0.739
                lny2 |     -0.477      0.399    -1.20    0.23       -1.258       0.304
                lnw1 |      0.417      0.165     2.52    0.01        0.092       0.741
               trend |     -0.556      0.040   -13.99    0.00       -0.634      -0.478
                     |
c.half#c.lny1#c.lny1 |      0.048      0.005    10.18    0.00        0.039       0.057
                     |
c.half#c.lny2#c.lny2 |      0.171      0.032     5.37    0.00        0.109       0.234
                     |
c.half#c.lnw1#c.lnw1 |     -0.066      0.014    -4.79    0.00       -0.093      -0.039
                     |
       c.lny1#c.lny2 |     -0.072      0.011    -6.30    0.00       -0.094      -0.049
                     |
       c.lny1#c.lnw1 |     -0.001      0.005    -0.23    0.82       -0.012       0.009
                     |
       c.lny2#c.lnw1 |     -0.014      0.013    -1.04    0.30       -0.040       0.012
                     |
      c.trend#c.lny1 |      0.005      0.002     3.04    0.00        0.002       0.009
                     |
      c.trend#c.lny2 |      0.012      0.003     4.06    0.00        0.006       0.018
                     |
      c.trend#c.lnw1 |     -0.006      0.003    -2.45    0.01       -0.011      -0.001
                     |
     c.trend#c.trend#|
              c.half |      0.076      0.002    42.52    0.00        0.073       0.080
                     |
               _cons |      4.190      2.657     1.58    0.11       -1.018       9.397
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
               _cons |     -4.429      0.032  -138.49    0.00       -4.492      -4.366
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
               _cons |     -2.638      0.092   -28.53    0.00       -2.819      -2.457
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |      0.029      0.009     3.16    0.00        0.011       0.048
--------------------------------------------------------------------------------------
Note, the function of inefficiency change over time is
beta[t] = exp( -gamma*(t-T_i) ): 
 the larger is the beta[t] the larger is the inefficiency

. timer off 4

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -2.3348946

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -2.2981796

. eststo M4

. quietly predictnl M4_mW1 = $monotonW1 if e(sample)

. quietly predictnl M4_mW2 = 1-M4_mW1 if e(sample)

. quietly predictnl M4_eY1 = $monotonY1 if e(sample)

. quietly predictnl M4_eY2 = $monotonY2 if e(sample)

. generate double M4_rts = 1/(M4_eY1 + M4_eY2)
(872 missing values generated)

. tabstat M4_*, stat(min p25 mean p75 max sd) format(%9.4f)

   Stats |    M4_mW1    M4_mW2    M4_eY1    M4_eY2    M4_rts
---------+--------------------------------------------------
     Min |   -0.1709    0.8793   -0.1587    0.4184    0.8521
     p25 |   -0.0141    0.9655    0.1305    0.6800    1.0495
    Mean |    0.0081    0.9919    0.1601    0.7494    1.1050
     p75 |    0.0345    1.0141    0.1932    0.8137    1.1541
     Max |    0.1207    1.1709    0.3255    1.2968    1.4286
      SD |    0.0391    0.0391    0.0542    0.1077    0.0784
------------------------------------------------------------
```

## Results of Models 1-4

Use the *estout* command for this. Note that parameters $\delta$ and $\gamma$ have different interpretation in models 2-4.

``` stata
estout M1 M2 M3 M4, ///
 cells(b(star fmt(%9.4f)) se(par)) ///
 stats(AIC BIC shat RSS ll N sumTi, ///
 labels("AIC" "BIC" "\$\hat\sigma\$" "RSS" "log-likelihood" "\$N\$" "\$\sum T_{i}\$") ///
 fmt(%9.4f %9.4f %9.4f %9.2f %9.2f %9.0f %9.0f)) ///
 starlevels(* 0.10 ** 0.05 *** 0.01) ///
 varlabels(_cons Constant ) ///
 substitute("_ " "Frontier") ///
 legend label collabels(none) mlabels(none) replace

------------------------------------------------------------------------------------
Frontier                                                                                  
lny1                       0.3966***       0.4286***       0.4747***       0.4592***
                         (0.1422)        (0.1403)        (0.1418)        (0.1427)   
lny2                      -0.6190         -0.4833         -0.5632         -0.4770   
                         (0.3989)        (0.3933)        (0.3945)        (0.3985)   
lnw1                       0.3852**        0.4970***       0.3933**        0.4166** 
                         (0.1663)        (0.1625)        (0.1646)        (0.1654)   
trend                     -0.5613***      -0.6561***      -0.5925***      -0.5562***
                         (0.0401)        (0.0404)        (0.0405)        (0.0398)   
half # lny1 # lny1         0.0495***       0.0468***       0.0475***       0.0482***
                         (0.0048)        (0.0046)        (0.0047)        (0.0047)   
half # lny2 # lny2         0.1783***       0.1713***       0.1796***       0.1713***
                         (0.0322)        (0.0314)        (0.0316)        (0.0319)   
half # lnw1 # lnw1        -0.0687***      -0.0616***      -0.0664***      -0.0659***
                         (0.0139)        (0.0134)        (0.0138)        (0.0138)   
lny1 # lny2               -0.0669***      -0.0686***      -0.0727***      -0.0718***
                         (0.0114)        (0.0112)        (0.0113)        (0.0114)   
lny1 # lnw1               -0.0025          0.0001          0.0000         -0.0013   
                         (0.0055)        (0.0053)        (0.0055)        (0.0055)   
lny2 # lnw1               -0.0091         -0.0237*        -0.0130         -0.0139   
                         (0.0134)        (0.0131)        (0.0133)        (0.0133)   
trend # lny1               0.0054***       0.0054***       0.0048***       0.0053***
                         (0.0017)        (0.0017)        (0.0017)        (0.0017)   
trend # lny2               0.0122***       0.0131***       0.0130***       0.0125***
                         (0.0031)        (0.0030)        (0.0031)        (0.0031)   
trend # lnw1              -0.0068***      -0.0058**       -0.0061**       -0.0064** 
                         (0.0026)        (0.0025)        (0.0026)        (0.0026)   
trend # trend # half       0.0765***       0.0937***       0.0823***       0.0762***
                         (0.0018)        (0.0029)        (0.0023)        (0.0018)   
Constant                   5.3832**        4.4970*         4.7054*         4.1896   
                         (2.6418)        (2.6172)        (2.6305)        (2.6571)   
------------------------------------------------------------------------------------
ln[var(vit)]                                                                        
Constant                  -4.4220***      -4.4975***      -4.4352***      -4.4290***
                         (0.0320)        (0.0324)        (0.0319)        (0.0320)   
------------------------------------------------------------------------------------
ln[var(ui)]                                                                         
Constant                  -2.5216***      -2.3149***      -2.7570***      -2.6381***
                         (0.0837)        (0.0889)        (0.0998)        (0.0925)   
------------------------------------------------------------------------------------
beta[t]                                                                             
gamma                                     -5.0001***      -0.1252***       0.0295***
                                         (1.8983)        (0.0276)        (0.0093)   
delta                                      0.8193***      -0.0197***                
                                         (0.3145)        (0.0052)                   
------------------------------------------------------------------------------------
AIC                       -2.3400         -2.3091         -2.3335         -2.3349   
BIC                       -2.3033         -2.2724         -2.2968         -2.2982   
$\hat\sigma$               0.1871          0.1900          0.1877          0.1876   
RSS                        228.96          250.63          227.45          229.70   
log-likelihood            1430.28         1492.31         1443.69         1435.26   
$N$                           500             500             500             500   
$\sum T_{i}$                 2546            2546            2546            2546   
------------------------------------------------------------------------------------
* p<0.10, ** p<0.05, *** p<0.01
```

## Execution Speed (in seconds), Models 1-4

The estimation is extremely fast. Time in seconds:

``` stata
. timer list 1
   1:      0.79 /        1 =       0.7850

. timer list 2
   2:      0.77 /        1 =       0.7740

. timer list 3
   3:      0.74 /        1 =       0.7450

. timer list 4
   4:      0.49 /        1 =       0.4930
```

# Models 5-8: Normal-truncated normal

Here are the models with normal-truncated normal assumptions.

## Model 5

``` stata
timer clear 5
timer on 5
xtsf1g lnc $spetech if year > $year_c, cost iterate($itermax) distribution(t) $decimalopt
timer off 5
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M5
quietly predictnl M5_mW5 = $monotonW1 if e(sample)
quietly predictnl M5_mW2 = 1-M5_mW5 if e(sample)
quietly predictnl M5_eY5 = $monotonY1 if e(sample)
quietly predictnl M5_eY2 = $monotonY2 if e(sample)
quietly generate double M5_rts = 1/(M5_eY5 + M5_eY2)
tabstat M5_*, stat(min p25 mean p75 max sd) format(%9.4f)

. timer clear 5

. timer on 5

. xtsf1g lnc $spetech if year > $year_c, cost iterate($itermax) distribution(t) $decimalopt


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =  499.06419  (not concave)
iter 1:  Log-likelihood =  1248.0912  
iter 2:  Log-likelihood =  1377.9019  
iter 3:  Log-likelihood =  1425.5041  
iter 4:  Log-likelihood =  1473.1026  
iter 5:  Log-likelihood =  1489.2684  
iter 6:  Log-likelihood =  1496.1943  
iter 7:  Log-likelihood =  1496.4395  
iter 8:  Log-likelihood =  1496.4627  
iter 9:  Log-likelihood =  1496.4649  
iter 10: Log-likelihood =  1496.4652  
iter 11: Log-likelihood =  1496.4652  

Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.9905
 Adj R-squared    = 0.9881
 AIC              = -2.3421
 BIC              = -2.3054
 Root MSE         = 0.1869
-----------------------------

 The first-generation estimator of
 the cost stochastic frontier model for panel data,
 where effciency is time-invariant

--------------------------------------------------------------------------------------
                 lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
---------------------+----------------------------------------------------------------
                lny1 |      0.208      0.141     1.47    0.14       -0.070       0.485
                lny2 |     -1.735      0.398    -4.36    0.00       -2.516      -0.955
                lnw1 |      0.009      0.178     0.05    0.96       -0.341       0.358
               trend |     -0.515      0.039   -13.32    0.00       -0.591      -0.439
                     |
c.half#c.lny1#c.lny1 |      0.048      0.005    10.06    0.00        0.038       0.057
                     |
c.half#c.lny2#c.lny2 |      0.257      0.032     7.95    0.00        0.194       0.321
                     |
c.half#c.lnw1#c.lnw1 |     -0.030      0.015    -2.02    0.04       -0.059      -0.001
                     |
       c.lny1#c.lny2 |     -0.050      0.011    -4.46    0.00       -0.072      -0.028
                     |
       c.lny1#c.lnw1 |      0.000      0.005     0.04    0.97       -0.010       0.011
                     |
       c.lny2#c.lnw1 |      0.014      0.014     0.97    0.33       -0.014       0.041
                     |
      c.trend#c.lny1 |      0.004      0.002     2.64    0.01        0.001       0.008
                     |
      c.trend#c.lny2 |      0.009      0.003     3.15    0.00        0.004       0.015
                     |
      c.trend#c.lnw1 |     -0.009      0.003    -3.37    0.00       -0.014      -0.004
                     |
     c.trend#c.trend#|
              c.half |      0.077      0.002    44.66    0.00        0.073       0.080
                     |
               _cons |     12.540      3.188     3.93    0.00        6.292      18.787
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
               _cons |     -4.512      0.031  -143.66    0.00       -4.574      -4.451
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
               _cons |     -3.646      0.070   -51.84    0.00       -3.784      -3.508
---------------------+----------------------------------------------------------------
E[ui|z]              |
               _cons |      0.797      1.831     0.44    0.66       -2.792       4.386
--------------------------------------------------------------------------------------

. timer off 5

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -2.3420887

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -2.3053737

. eststo M5

. quietly predictnl M5_mW5 = $monotonW1 if e(sample)

. quietly predictnl M5_mW2 = 1-M5_mW5 if e(sample)

. quietly predictnl M5_eY5 = $monotonY1 if e(sample)

. quietly predictnl M5_eY2 = $monotonY2 if e(sample)

. quietly generate double M5_rts = 1/(M5_eY5 + M5_eY2)

. tabstat M5_*, stat(min p25 mean p75 max sd) format(%9.4f)

   Stats |    M5_mW5    M5_mW2    M5_eY5    M5_eY2    M5_rts
---------+--------------------------------------------------
     Min |   -0.0684    0.9204   -0.1494    0.2310    0.8575
     p25 |    0.0091    0.9590    0.1253    0.6525    1.0330
    Mean |    0.0239    0.9761    0.1482    0.7411    1.1409
     p75 |    0.0410    0.9909    0.1778    0.8319    1.2252
     Max |    0.0796    1.0684    0.2884    1.2766    2.1699
      SD |    0.0234    0.0234    0.0481    0.1322    0.1430
------------------------------------------------------------
```

## Model 6 (2nd generation, time-varying inefficiency using the Kumbhakar (1990) specification)

``` stata
timer clear 6
timer on 6
xtsf2g lnc $spetech if year > $year_c, cost iterate($itermax) distribution(t) $decimalopt
timer off 6
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M6
quietly predictnl M6_mW1 = $monotonW1 if e(sample)
quietly predictnl M6_mW2 = 1-M6_mW1 if e(sample)
quietly predictnl M6_eY1 = $monotonY1 if e(sample)
quietly predictnl M6_eY2 = $monotonY2 if e(sample)
quietly generate double M6_rts = 1/(M6_eY1 + M6_eY2)
tabstat M6_*, stat(min p25 mean p75 max sd) format(%9.4f)

. timer clear 6

. timer on 6

. xtsf2g lnc $spetech if year > $year_c, cost iterate($itermax) distribution(t) $decimalopt


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =   649.6318  (not concave)
iter 1:  Log-likelihood =  1241.5285  (not concave)
iter 2:  Log-likelihood =  1322.5145  (not concave)
iter 3:  Log-likelihood =   1355.508  (not concave)
iter 4:  Log-likelihood =   1392.307  (not concave)
iter 5:  Log-likelihood =  1421.9537  (not concave)
iter 6:  Log-likelihood =  1433.7357  
iter 7:  Log-likelihood =  1450.4961  (not concave)
iter 8:  Log-likelihood =  1483.9436  (not concave)
iter 9:  Log-likelihood =  1488.1369  
iter 10: Log-likelihood =  1521.1974  
iter 11: Log-likelihood =  1554.7545  
...
iter 77: Log-likelihood =  1623.0228  
iter 78: Log-likelihood =  1623.0229  
iter 79: Log-likelihood =   1623.023  
iter 80: Log-likelihood =   1623.023  

Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = -0.7572
 Adj R-squared    = -1.2019
 AIC              = 1.1384
 BIC              = 1.1752
 Root MSE         = 1.0652
-----------------------------

 The second-generation estimator of
 the cost stochastic frontier model for panel data,
 where effciency is time-varying

--------------------------------------------------------------------------------------
                 lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
---------------------+----------------------------------------------------------------
                lny1 |      0.178      0.136     1.31    0.19       -0.089       0.444
                lny2 |     -1.712      0.382    -4.48    0.00       -2.461      -0.963
                lnw1 |      0.073      0.169     0.43    0.66       -0.258       0.405
               trend |     -4.351      1.593    -2.73    0.01       -7.475      -1.228
                     |
c.half#c.lny1#c.lny1 |      0.048      0.004    10.82    0.00        0.039       0.057
                     |
c.half#c.lny2#c.lny2 |      0.255      0.031     8.28    0.00        0.195       0.316
                     |
c.half#c.lnw1#c.lnw1 |     -0.023      0.014    -1.68    0.09       -0.051       0.004
                     |
       c.lny1#c.lny2 |     -0.048      0.011    -4.46    0.00       -0.069      -0.027
                     |
       c.lny1#c.lnw1 |      0.001      0.005     0.27    0.79       -0.008       0.011
                     |
       c.lny2#c.lnw1 |      0.005      0.013     0.39    0.70       -0.021       0.032
                     |
      c.trend#c.lny1 |      0.005      0.002     2.94    0.00        0.002       0.008
                     |
      c.trend#c.lny2 |      0.010      0.003     3.55    0.00        0.004       0.015
                     |
      c.trend#c.lnw1 |     -0.009      0.002    -3.69    0.00       -0.013      -0.004
                     |
     c.trend#c.trend#|
              c.half |      0.831      0.315     2.64    0.01        0.214       1.448
                     |
               _cons |     10.760      3.428     3.14    0.00        4.042      17.478
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
               _cons |     -4.638      0.031  -147.71    0.00       -4.699      -4.576
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
               _cons |     -2.781      0.122   -22.89    0.00       -3.019      -2.543
---------------------+----------------------------------------------------------------
E[ui|z]              |
               _cons |     17.054     10.102     1.69    0.09       -2.745      36.852
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |     -0.594      0.129    -4.61    0.00       -0.847      -0.341
               delta |      0.096      0.021     4.69    0.00        0.056       0.136
--------------------------------------------------------------------------------------
Note, the function of inefficiency change over time is
beta[t] = ( 1 + exp(gamma*t + delta*t^2) )^-1: 
 the larger is the beta[t] the larger is the inefficiency

. timer off 6

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  1.1384354

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  1.1751504

. eststo M6

. quietly predictnl M6_mW1 = $monotonW1 if e(sample)

. quietly predictnl M6_mW2 = 1-M6_mW1 if e(sample)

. quietly predictnl M6_eY1 = $monotonY1 if e(sample)

. quietly predictnl M6_eY2 = $monotonY2 if e(sample)

. quietly generate double M6_rts = 1/(M6_eY1 + M6_eY2)

. tabstat M6_*, stat(min p25 mean p75 max sd) format(%9.4f)

   Stats |    M6_mW1    M6_mW2    M6_eY1    M6_eY2    M6_rts
---------+--------------------------------------------------
     Min |   -0.0483    0.9271   -0.1466    0.2418    0.8629
     p25 |    0.0119    0.9603    0.1268    0.6539    1.0308
    Mean |    0.0251    0.9749    0.1495    0.7411    1.1395
     p75 |    0.0397    0.9881    0.1792    0.8326    1.2245
     Max |    0.0729    1.0483    0.2913    1.2670    2.1383
      SD |    0.0199    0.0199    0.0480    0.1315    0.1436
------------------------------------------------------------
```

## Model 7 (2nd generation, time-varying inefficiency using the modified Kumbhakar (1990) specification)


``` stata
timer clear 7
timer on 7
xtsf2g lnc $spetech if year > $year_c, model(K1990modified) cost iterate($itermax) distribution(t) $decimalopt
timer off 7
predict M7_te, te_jlms_mean
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M7
quietly predictnl M7_mW1 = $monotonW1 if e(sample)
quietly predictnl M7_mW2 = 1-M7_mW1 if e(sample)
quietly predictnl M7_eY1 = $monotonY1 if e(sample)
quietly predictnl M7_eY2 = $monotonY2 if e(sample)
generate double M7_rts = 1/(M7_eY1 + M7_eY2)
tabstat M7_*, stat(min p25 mean p75 max sd) format(%9.4f)

. timer on 7

. xtsf2g lnc $spetech if year > $year_c, model(K1990modified) cost iterate($itermax) distribution(t)  $decimalopt


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =  499.06419  (not concave)
iter 1:  Log-likelihood =  1239.3495  
iter 2:  Log-likelihood =  1314.4895  (not concave)
iter 3:  Log-likelihood =  1436.7572  (not concave)
iter 4:  Log-likelihood =  1447.0897  
iter 5:  Log-likelihood =  1471.5866  
iter 6:  Log-likelihood =  1495.7501  
iter 7:  Log-likelihood =  1500.0676  
iter 8:  Log-likelihood =  1501.2273  
iter 9:  Log-likelihood =  1501.3099  
iter 10: Log-likelihood =  1501.3116  
iter 11: Log-likelihood =  1501.3116  

Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.9936
 Adj R-squared    = 0.9920
 AIC              = -2.3338
 BIC              = -2.2971
 Root MSE         = 0.1877
-----------------------------

 The second-generation estimator of
 the cost stochastic frontier model for panel data,
 where effciency is time-varying

--------------------------------------------------------------------------------------
                 lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
---------------------+----------------------------------------------------------------
                lny1 |      0.231      0.141     1.63    0.10       -0.046       0.508
                lny2 |     -1.787      0.398    -4.48    0.00       -2.568      -1.006
                lnw1 |      0.042      0.179     0.23    0.82       -0.309       0.392
               trend |     -0.535      0.041   -13.16    0.00       -0.614      -0.455
                     |
c.half#c.lny1#c.lny1 |      0.047      0.005     9.99    0.00        0.038       0.056
                     |
c.half#c.lny2#c.lny2 |      0.265      0.032     8.19    0.00        0.201       0.328
                     |
c.half#c.lnw1#c.lnw1 |     -0.030      0.015    -2.03    0.04       -0.059      -0.001
                     |
       c.lny1#c.lny2 |     -0.052      0.011    -4.61    0.00       -0.073      -0.030
                     |
       c.lny1#c.lnw1 |      0.001      0.005     0.23    0.82       -0.009       0.012
                     |
       c.lny2#c.lnw1 |      0.010      0.014     0.69    0.49       -0.018       0.037
                     |
      c.trend#c.lny1 |      0.004      0.002     2.53    0.01        0.001       0.008
                     |
      c.trend#c.lny2 |      0.010      0.003     3.31    0.00        0.004       0.016
                     |
      c.trend#c.lnw1 |     -0.008      0.003    -3.33    0.00       -0.013      -0.003
                     |
     c.trend#c.trend#|
              c.half |      0.081      0.003    31.08    0.00        0.076       0.086
                     |
               _cons |     12.816      2.624     4.88    0.00        7.673      17.958
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
               _cons |     -4.519      0.031  -143.72    0.00       -4.580      -4.457
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
               _cons |     -3.728      0.077   -48.57    0.00       -3.878      -3.577
---------------------+----------------------------------------------------------------
E[ui|z]              |
               _cons |      0.619      0.084     7.38    0.00        0.455       0.784
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |     -0.038      0.014    -2.80    0.01       -0.065      -0.011
               delta |     -0.005      0.002    -2.10    0.04       -0.009      -0.000
--------------------------------------------------------------------------------------
Note, the function of inefficiency change over time is
beta[t] = 1 + gamma*(t-T_i) + delta*(t-T_i)^2: 
 the larger is the beta[t] the larger is the inefficiency

. timer off 7

. predict M7_te, te_jlms_mean
 (872 missing values generated)

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -2.3337688

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -2.2970538

. eststo M7

. quietly predictnl M7_mW1 = $monotonW1 if e(sample)

. quietly predictnl M7_mW2 = 1-M7_mW1 if e(sample)

. quietly predictnl M7_eY1 = $monotonY1 if e(sample)

. quietly predictnl M7_eY2 = $monotonY2 if e(sample)

. generate double M7_rts = 1/(M7_eY1 + M7_eY2)
(872 missing values generated)

. tabstat M7_*, stat(min p25 mean p75 max sd) format(%9.4f)

   Stats |     M7_te    M7_mW1    M7_mW2    M7_eY1    M7_eY2    M7_rts
---------+------------------------------------------------------------
     Min |    1.0519   -0.0628    0.9206   -0.1460    0.2248    0.8420
     p25 |    1.7413    0.0097    0.9597    0.1235    0.6545    1.0278
    Mean |    1.9336    0.0241    0.9759    0.1469    0.7455    1.1378
     p75 |    2.0620    0.0403    0.9903    0.1761    0.8388    1.2234
     Max |    3.8968    0.0794    1.0628    0.2888    1.2980    2.1940
      SD |    0.3052    0.0228    0.0228    0.0478    0.1363    0.1466
----------------------------------------------------------------------
```

## Model 8 (2nd generation, time-varying inefficiency using the BC specification)

``` stata
timer clear 8
timer on 8
xtsf2g lnc $spetech if year > $year_c, model(BC1992) cost iterate($itermax) distribution(t) $decimalopt
timer off 8
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M8
quietly predictnl M8_mW1 = $monotonW1 if e(sample)
quietly predictnl M8_mW2 = 1-M8_mW1 if e(sample)
quietly predictnl M8_eY1 = $monotonY1 if e(sample)
quietly predictnl M8_eY2 = $monotonY2 if e(sample)
generate double M8_rts = 1/(M8_eY1 + M8_eY2)
tabstat M8_*, stat(min p25 mean p75 max sd) format(%9.4f)

. timer clear 8

. timer on 8

. xtsf2g lnc $spetech if year > $year_c, model(BC1992) cost iterate($itermax) distribution(t)  $decimalopt


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =  499.06419  (not concave)
iter 1:  Log-likelihood =  1245.8926  
iter 2:  Log-likelihood =  1374.1779  
iter 3:  Log-likelihood =  1433.7164  
iter 4:  Log-likelihood =  1479.5883  
iter 5:  Log-likelihood =  1493.3381  
iter 6:  Log-likelihood =  1497.8862  
iter 7:  Log-likelihood =   1498.517  
iter 8:  Log-likelihood =  1498.5502  
iter 9:  Log-likelihood =  1498.5507  
iter 10: Log-likelihood =  1498.5507  

Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.9927
 Adj R-squared    = 0.9909
 AIC              = -2.3343
 BIC              = -2.2976
 Root MSE         = 0.1876
-----------------------------

 The second-generation estimator of
 the cost stochastic frontier model for panel data,
 where effciency is time-varying

--------------------------------------------------------------------------------------
                 lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
---------------------+----------------------------------------------------------------
                lny1 |      0.231      0.142     1.63    0.10       -0.046       0.509
                lny2 |     -1.714      0.397    -4.31    0.00       -2.493      -0.936
                lnw1 |      0.025      0.178     0.14    0.89       -0.324       0.375
               trend |     -0.505      0.039   -13.02    0.00       -0.581      -0.429
                     |
c.half#c.lny1#c.lny1 |      0.047      0.005     9.99    0.00        0.038       0.056
                     |
c.half#c.lny2#c.lny2 |      0.258      0.032     8.01    0.00        0.195       0.321
                     |
c.half#c.lnw1#c.lnw1 |     -0.028      0.015    -1.91    0.06       -0.057       0.001
                     |
       c.lny1#c.lny2 |     -0.052      0.011    -4.62    0.00       -0.074      -0.030
                     |
       c.lny1#c.lnw1 |      0.001      0.005     0.15    0.88       -0.010       0.011
                     |
       c.lny2#c.lnw1 |      0.011      0.014     0.80    0.42       -0.016       0.039
                     |
      c.trend#c.lny1 |      0.005      0.002     2.73    0.01        0.001       0.008
                     |
      c.trend#c.lny2 |      0.009      0.003     3.15    0.00        0.004       0.015
                     |
      c.trend#c.lnw1 |     -0.009      0.003    -3.44    0.00       -0.014      -0.004
                     |
     c.trend#c.trend#|
              c.half |      0.076      0.002    44.49    0.00        0.073       0.080
                     |
               _cons |     12.319      2.616     4.71    0.00        7.192      17.445
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
               _cons |     -4.516      0.031  -143.61    0.00       -4.578      -4.455
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
               _cons |     -3.697      0.075   -49.51    0.00       -3.843      -3.550
---------------------+----------------------------------------------------------------
E[ui|z]              |
               _cons |      0.674      0.111     6.09    0.00        0.457       0.891
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |      0.014      0.007     2.02    0.04        0.000       0.027
--------------------------------------------------------------------------------------
Note, the function of inefficiency change over time is
beta[t] = exp( -gamma*(t-T_i) ): 
 the larger is the beta[t] the larger is the inefficiency

. timer off 8

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -2.3343331

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -2.297618

. eststo M8

. quietly predictnl M8_mW1 = $monotonW1 if e(sample)

. quietly predictnl M8_mW2 = 1-M8_mW1 if e(sample)

. quietly predictnl M8_eY1 = $monotonY1 if e(sample)

. quietly predictnl M8_eY2 = $monotonY2 if e(sample)

. generate double M8_rts = 1/(M8_eY1 + M8_eY2)
(872 missing values generated)

. tabstat M8_*, stat(min p25 mean p75 max sd) format(%9.4f)

   Stats |    M8_mW1    M8_mW2    M8_eY1    M8_eY2    M8_rts
---------+--------------------------------------------------
     Min |   -0.0618    0.9213   -0.1467    0.2353    0.8485
     p25 |    0.0105    0.9591    0.1239    0.6552    1.0313
    Mean |    0.0247    0.9753    0.1472    0.7442    1.1380
     p75 |    0.0409    0.9895    0.1768    0.8348    1.2215
     Max |    0.0787    1.0618    0.2894    1.2883    2.1433
      SD |    0.0225    0.0225    0.0480    0.1331    0.1417
------------------------------------------------------------
```

## Results of Models 5-8

Use the *estout* command for this. Note that parameters $\delta$ and $\gamma$ have different interpretation in models 6-8.

``` stata
estout M5 M6 M7 M8, ///
 cells(b(star fmt(%9.4f)) se(par)) ///
 stats(AIC BIC shat RSS ll N sumTi, ///
 labels("AIC" "BIC" "\$\hat\sigma\$" "RSS" "log-likelihood" "\$N\$" "\$\sum T_{i}\$") ///
 fmt(%9.4f %9.4f %9.4f %9.2f %9.2f %9.0f %9.0f)) ///
 starlevels(* 0.10 ** 0.05 *** 0.01) ///
 varlabels(_cons Constant ) ///
 substitute("_ " "Frontier") ///
 legend label collabels(none) mlabels(none) replace

------------------------------------------------------------------------------------
Frontier                                                                                  
lny1                       0.2077          0.1778          0.2311          0.2311   
                         (0.1414)        (0.1360)        (0.1415)        (0.1416)   
lny2                      -1.7352***      -1.7119***      -1.7872***      -1.7144***
                         (0.3981)        (0.3822)        (0.3985)        (0.3973)   
lnw1                       0.0086          0.0735          0.0417          0.0254   
                         (0.1781)        (0.1693)        (0.1789)        (0.1782)   
trend                     -0.5151***      -4.3514***      -0.5345***      -0.5053***
                         (0.0387)        (1.5935)        (0.0406)        (0.0388)   
half # lny1 # lny1         0.0476***       0.0480***       0.0468***       0.0469***
                         (0.0047)        (0.0044)        (0.0047)        (0.0047)   
half # lny2 # lny2         0.2571***       0.2553***       0.2645***       0.2578***
                         (0.0324)        (0.0308)        (0.0323)        (0.0322)   
half # lnw1 # lnw1        -0.0298**       -0.0234*        -0.0300**       -0.0282*  
                         (0.0148)        (0.0139)        (0.0148)        (0.0147)   
lny1 # lny2               -0.0499***      -0.0478***      -0.0515***      -0.0516***
                         (0.0112)        (0.0107)        (0.0112)        (0.0112)   
lny1 # lnw1                0.0002          0.0014          0.0012          0.0008   
                         (0.0053)        (0.0050)        (0.0053)        (0.0053)   
lny2 # lnw1                0.0136          0.0052          0.0098          0.0113   
                         (0.0140)        (0.0135)        (0.0141)        (0.0141)   
trend # lny1               0.0045***       0.0046***       0.0043**        0.0046***
                         (0.0017)        (0.0016)        (0.0017)        (0.0017)   
trend # lny2               0.0094***       0.0099***       0.0098***       0.0093***
                         (0.0030)        (0.0028)        (0.0030)        (0.0030)   
trend # lnw1              -0.0086***      -0.0087***      -0.0084***      -0.0087***
                         (0.0025)        (0.0024)        (0.0025)        (0.0025)   
trend # trend # half       0.0766***       0.8308***       0.0810***       0.0764***
                         (0.0017)        (0.3148)        (0.0026)        (0.0017)   
Constant                  12.5398***      10.7597***      12.8159***      12.3189***
                         (3.1876)        (3.4275)        (2.6238)        (2.6156)   
------------------------------------------------------------------------------------
ln[var(vit)]                                                                        
Constant                  -4.5121***      -4.6378***      -4.5188***      -4.5163***
                         (0.0314)        (0.0314)        (0.0314)        (0.0314)   
------------------------------------------------------------------------------------
ln[var(ui)]                                                                         
Constant                  -3.6457***      -2.7813***      -3.7279***      -3.6967***
                         (0.0703)        (0.1215)        (0.0768)        (0.0747)   
------------------------------------------------------------------------------------
E[ui|z]                                                                             
Constant                   0.7969         17.0536*         0.6192***       0.6738***
                         (1.8311)       (10.1015)        (0.0839)        (0.1107)   
------------------------------------------------------------------------------------
beta[t]                                                                             
gamma                                     -0.5944***      -0.0380***       0.0138** 
                                         (0.1290)        (0.0136)        (0.0069)   
delta                                      0.0962***      -0.0048**                 
                                         (0.0205)        (0.0023)                   
------------------------------------------------------------------------------------
AIC                       -2.3421          1.1384         -2.3338         -2.3343   
BIC                       -2.3054          1.1752         -2.2971         -2.2976   
$\hat\sigma$               0.1869          1.0652          0.1877          0.1876   
RSS                       1710.04       317650.01         1158.39         1317.51   
log-likelihood            1496.47         1623.02         1501.31         1498.55   
$N$                           500             500             500             500   
$\sum T_{i}$                 2546            2546            2546            2546   
------------------------------------------------------------------------------------
* p<0.10, ** p<0.05, *** p<0.01
```

## Execution Speed (in seconds), Models 5-8

The estimation is extremely fast. Time in seconds:

``` stata
. timer list 5
   5:      0.56 /        1 =       0.5610

. timer list 6
   6:      3.05 /        1 =       3.0480

. timer list 7
   7:      0.68 /        1 =       0.6810

. timer list 8
   8:      0.57 /        1 =       0.5690
```

# Determinants of inefficiency and heterskedasticity function

In this section time-invariant determinants of inefficiency are used.

# Models 9-12: with determinants

The following specifications are used for the next models. The conditional mean of the truncated distribution can be made observation specific by imposing linear structure on it as in [Kumbhakar et al (1991)](https://doi.org/10.1080/07350015.1991.10509853),

$$
\mu_{i} = E[u_{i}|\boldsymbol{z}_{u_{2i}}] = \boldsymbol{z}_{u_{2i}} \boldsymbol{\delta}_{u},
$$

so that $u_{i} \sim N^{+}( \boldsymbol{z}_{u_{2i}} \boldsymbol{\delta}_{u}, \sigma_{u}^{2} )$ and only $\boldsymbol{\delta}_{u}$ parameters are estimated instead of $\mu_{i}, i=1,\ldots, N$. 

With half-normal or truncated normal distribution, the inefficiency term is still *i.i.d.*, which can introduce bias. If the heteroskedasticity function is known, bias can be eliminated by making $\sigma_{u}^{2}$ observation specific (see [Caudill, Ford and Gropper, 1995](https://doi.org/10.1080/07350015.1995.10524583))

$$
\ln{ \sigma_{u_{i}}^{2}} = \boldsymbol{z}_{u_{1i}} \boldsymbol{\gamma}_{u}
$$

Besides, since $E(u_{i})= (2/\pi)^{1/2}\sigma_{u_{i}} = (2/\pi)^{1/2} \exp{0.5 \boldsymbol{z}_{ui} \boldsymbol{\gamma}_{u}}$, the $\boldsymbol{z}_{ui}$ variables can be viewed as determinants of persistent inefficiency.

Similarly, the heteroskedasticity function of noise can be specified
$$
\ln{ \sigma_{v_{it}}^{2}} = \boldsymbol{z}_{v_{it}} \boldsymbol{\gamma}_{v}.
$$

In the following, $\mu_i$ is defined by `uimean`, $\ln{ \sigma_{u_{i}}^{2}}$ is defined by `uilnvariance`, and $\ln{ \sigma_{v_{it}}^{2}}$ is defined by `vitlnvariance`. 

## Model 9

``` stata
timer clear 9
timer on 9
xtsf1g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnvariance(er_ave) vitlnvariance(lnta) $decimalopt
timer off 9
predict M9_te, te_jlms_mean
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M9
quietly predictnl M9_mW1 = $monotonW1 if e(sample)
quietly predictnl M9_mW2 = 1-M9_mW1 if e(sample)
quietly predictnl M9_eY1 = $monotonY1 if e(sample)
quietly predictnl M9_eY2 = $monotonY2 if e(sample)
generate double M9_rts = 1/(M9_eY1 + M9_eY2)
tabstat M9_*, stat(min p25 mean p75 max sd) format(%9.4f)

. timer clear 9

. timer on 9

. xtsf1g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnvariance(er_ave) vitlnvariance(lnta) $decimalopt


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =  530.64988  (not concave)
iter 1:  Log-likelihood =  1274.0126  (not concave)
iter 2:  Log-likelihood =  1423.8577  
iter 3:  Log-likelihood =  1502.7847  
iter 4:  Log-likelihood =  1521.1266  
iter 5:  Log-likelihood =  1523.9925  
iter 6:  Log-likelihood =  1524.0591  
iter 7:  Log-likelihood =  1524.0591  

Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.9976
 Adj R-squared    = 0.9970
 AIC              = -2.2954
 BIC              = -2.2587
 Root MSE         = 0.1913
-----------------------------

 The first-generation estimator of
 the cost stochastic frontier model for panel data,
 where effciency is time-invariant

--------------------------------------------------------------------------------------
                 lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
---------------------+----------------------------------------------------------------
                lny1 |      0.302      0.141     2.15    0.03        0.026       0.578
                lny2 |     -1.249      0.392    -3.18    0.00       -2.017      -0.480
                lnw1 |      0.167      0.161     1.04    0.30       -0.148       0.482
               trend |     -0.543      0.039   -14.11    0.00       -0.619      -0.468
                     |
c.half#c.lny1#c.lny1 |      0.053      0.005    11.09    0.00        0.044       0.062
                     |
c.half#c.lny2#c.lny2 |      0.219      0.032     6.91    0.00        0.157       0.281
                     |
c.half#c.lnw1#c.lnw1 |     -0.044      0.014    -3.15    0.00       -0.071      -0.017
                     |
       c.lny1#c.lny2 |     -0.060      0.011    -5.42    0.00       -0.082      -0.038
                     |
       c.lny1#c.lnw1 |     -0.000      0.005    -0.02    0.98       -0.010       0.010
                     |
       c.lny2#c.lnw1 |      0.002      0.013     0.15    0.88       -0.023       0.027
                     |
      c.trend#c.lny1 |      0.004      0.002     2.21    0.03        0.000       0.007
                     |
      c.trend#c.lny2 |      0.013      0.003     4.36    0.00        0.007       0.019
                     |
      c.trend#c.lnw1 |     -0.008      0.003    -3.12    0.00       -0.013      -0.003
                     |
     c.trend#c.trend#|
              c.half |      0.076      0.002    44.14    0.00        0.073       0.079
                     |
               _cons |      9.823      2.582     3.80    0.00        4.763      14.883
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
                lnta |      0.430      0.085     5.04    0.00        0.263       0.597
               _cons |     -9.502      0.996    -9.54    0.00      -11.453      -7.550
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
              er_ave |     -3.430      2.589    -1.32    0.19       -8.504       1.644
               _cons |     -3.399      0.282   -12.07    0.00       -3.951      -2.847
---------------------+----------------------------------------------------------------
E[ui|z]              |
             llp_ave |      0.000      0.000    11.81    0.00        0.000       0.000
              la_ave |      0.400      0.048     8.33    0.00        0.306       0.495
--------------------------------------------------------------------------------------

. timer off 9

. predict M9_te, te_jlms_mean
 (872 missing values generated)

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -2.2954302

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -2.2587152

. eststo M9

. quietly predictnl M9_mW1 = $monotonW1 if e(sample)

. quietly predictnl M9_mW2 = 1-M9_mW1 if e(sample)

. quietly predictnl M9_eY1 = $monotonY1 if e(sample)

. quietly predictnl M9_eY2 = $monotonY2 if e(sample)

. generate double M9_rts = 1/(M9_eY1 + M9_eY2)
(872 missing values generated)

. tabstat M9_*, stat(min p25 mean p75 max sd) format(%9.4f)

   Stats |     M9_te    M9_mW1    M9_mW2    M9_eY1    M9_eY2    M9_rts
---------+------------------------------------------------------------
     Min |    1.0096   -0.1121    0.9162   -0.1612    0.2442    0.9010
     p25 |    1.2933   -0.0075    0.9717    0.1498    0.5997    1.0902
    Mean |    1.4539    0.0086    0.9914    0.1768    0.6805    1.1781
     p75 |    1.5453    0.0283    1.0075    0.2101    0.7585    1.2498
     Max |    3.0491    0.0838    1.1121    0.3337    1.2298    1.9093
      SD |    0.2417    0.0283    0.0283    0.0545    0.1206    0.1206
----------------------------------------------------------------------
```

## Model 10

``` stata
timer clear 10
timer on 10
xtsf2g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnvariance(er_ave) vitlnvariance(lnta) model("K1990") $decimalopt
timer off 10
predict M10_te, te_jlms_mean
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M10
quietly predictnl M10_mW1 = $monotonW1 if e(sample)
quietly predictnl M10_mW2 = 1-M10_mW1 if e(sample)
quietly predictnl M10_eY1 = $monotonY1 if e(sample)
quietly predictnl M10_eY2 = $monotonY2 if e(sample)
generate double M10_rts = 1/(M10_eY1 + M10_eY2)
tabstat M10_*, stat(min p25 mean p75 max sd) format(%9.4f)

. timer clear 10

. timer on 10

. xtsf2g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnva
> riance(er_ave) vitlnvariance(lnta) model("K1990") $decimalopt


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =  685.59629  (not concave)
iter 1:  Log-likelihood =  1271.6597  (not concave)
iter 2:  Log-likelihood =  1317.3539  (not concave)
iter 3:  Log-likelihood =  1407.2609  
iter 4:  Log-likelihood =   1517.782  (not concave)
...
iter 998: Log-likelihood =  1546.8553  (not concave)
iter 999: Log-likelihood =  1546.8778  (not concave)
iter 1000: Log-likelihood =  1546.9004  (not concave)
convergence not achieved


Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.9979
 Adj R-squared    = 0.9974
 AIC              = -2.2840
 BIC              = -2.2473
 Root MSE         = 0.1924
-----------------------------

 The second-generation estimator of
 the cost stochastic frontier model for panel data,
 where effciency is time-varying

--------------------------------------------------------------------------------------
                 lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
---------------------+----------------------------------------------------------------
                lny1 |      0.659      0.142     4.65    0.00        0.381       0.937
                lny2 |     -0.527      0.391    -1.35    0.18       -1.293       0.239
                lnw1 |      0.427      0.159     2.68    0.01        0.115       0.739
               trend |     -0.673      0.040   -16.86    0.00       -0.751      -0.595
                     |
c.half#c.lny1#c.lny1 |      0.048      0.005    10.15    0.00        0.038       0.057
                     |
c.half#c.lny2#c.lny2 |      0.186      0.031     5.98    0.00        0.125       0.247
                     |
c.half#c.lnw1#c.lnw1 |     -0.042      0.014    -3.06    0.00       -0.069      -0.015
                     |
       c.lny1#c.lny2 |     -0.087      0.011    -7.81    0.00       -0.109      -0.065
                     |
       c.lny1#c.lnw1 |     -0.002      0.005    -0.35    0.72       -0.012       0.008
                     |
       c.lny2#c.lnw1 |     -0.020      0.013    -1.55    0.12       -0.045       0.005
                     |
      c.trend#c.lny1 |      0.004      0.002     2.19    0.03        0.000       0.007
                     |
      c.trend#c.lny2 |      0.012      0.003     4.18    0.00        0.007       0.018
                     |
      c.trend#c.lnw1 |     -0.008      0.002    -3.35    0.00       -0.013      -0.003
                     |
     c.trend#c.trend#|
              c.half |      0.103      0.003    33.44    0.00        0.097       0.109
                     |
               _cons |      3.826      2.605     1.47    0.14       -1.280       8.932
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
                lnta |      0.443      0.086     5.16    0.00        0.275       0.612
               _cons |     -9.689      1.005    -9.64    0.00      -11.658      -7.720
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
              er_ave |     -4.248      2.663    -1.60    0.11       -9.468       0.971
               _cons |     -2.550      0.295    -8.63    0.00       -3.129      -1.970
---------------------+----------------------------------------------------------------
E[ui|z]              |
             llp_ave |      0.000      0.000    11.13    0.00        0.000       0.001
              la_ave |      0.550      0.073     7.50    0.00        0.406       0.694
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |     -0.818          .        .       .            .           .
               delta |      0.133      0.003    44.49    0.00        0.127       0.139
--------------------------------------------------------------------------------------
Note, the function of inefficiency change over time is
beta[t] = ( 1 + exp(gamma*t + delta*t^2) )^-1: 
 the larger is the beta[t] the larger is the inefficiency

. timer off 10

. predict M10_te, te_jlms_mean
 (872 missing values generated)

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -2.2839873

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -2.2472723

. eststo M10

. quietly predictnl M10_mW1 = $monotonW1 if e(sample)

. quietly predictnl M10_mW2 = 1-M10_mW1 if e(sample)

. quietly predictnl M10_eY1 = $monotonY1 if e(sample)

. quietly predictnl M10_eY2 = $monotonY2 if e(sample)

. generate double M10_rts = 1/(M10_eY1 + M10_eY2)
(872 missing values generated)

. tabstat M10_*, stat(min p25 mean p75 max sd) format(%9.4f)

   Stats |    M10_te   M10_mW1   M10_mW2   M10_eY1   M10_eY2   M10_rts
---------+------------------------------------------------------------
     Min |    1.0068   -0.1345    0.8841   -0.1484    0.3370    0.8216
     p25 |    1.2592   -0.0084    0.9714    0.1368    0.6276    1.0857
    Mean |    1.4119    0.0088    0.9912    0.1696    0.7054    1.1504
     p75 |    1.5019    0.0286    1.0084    0.2058    0.7779    1.2075
     Max |    3.0785    0.1159    1.1345    0.3579    1.3408    1.4940
      SD |    0.2410    0.0304    0.0304    0.0579    0.1228    0.0938
----------------------------------------------------------------------
```

## Model 11

``` stata
timer clear 11
timer on 11
xtsf2g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnvariance(er_ave) vitlnvariance(lnta) model("K1990modified") $decimalopt
timer off 11
predict M11_te, te_jlms_mean
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M11
quietly predictnl M11_mW1 = $monotonW1 if e(sample)
quietly predictnl M11_mW2 = 1-M11_mW1 if e(sample)
quietly predictnl M11_eY1 = $monotonY1 if e(sample)
quietly predictnl M11_eY2 = $monotonY2 if e(sample)
generate double M11_rts = 1/(M11_eY1 + M11_eY2)
tabstat M11_*, stat(min p25 mean p75 max sd) format(%9.4f)

. timer clear 11

. timer on 11

. xtsf2g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnva
> riance(er_ave) vitlnvariance(lnta) model("K1990modified") $decimalopt


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =  530.64988  (not concave)
iter 1:  Log-likelihood =  1275.4967  (not concave)
iter 2:  Log-likelihood =   1418.682  
iter 3:  Log-likelihood =  1495.3091  
iter 4:  Log-likelihood =  1508.3066  
iter 5:  Log-likelihood =  1527.6979  
iter 6:  Log-likelihood =  1530.5061  
iter 7:  Log-likelihood =  1535.1122  
iter 8:  Log-likelihood =  1535.3355  
iter 9:  Log-likelihood =  1535.3362  
iter 10: Log-likelihood =  1535.3362  

Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.9977
 Adj R-squared    = 0.9971
 AIC              = -2.2876
 BIC              = -2.2509
 Root MSE         = 0.1921
-----------------------------

 The second-generation estimator of
 the cost stochastic frontier model for panel data,
 where effciency is time-varying

--------------------------------------------------------------------------------------
                 lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
---------------------+----------------------------------------------------------------
                lny1 |      0.361      0.141     2.57    0.01        0.086       0.637
                lny2 |     -1.244      0.393    -3.17    0.00       -2.013      -0.474
                lnw1 |      0.191      0.159     1.21    0.23       -0.120       0.502
               trend |     -0.577      0.039   -14.64    0.00       -0.655      -0.500
                     |
c.half#c.lny1#c.lny1 |      0.051      0.005    10.94    0.00        0.042       0.060
                     |
c.half#c.lny2#c.lny2 |      0.223      0.031     7.09    0.00        0.161       0.284
                     |
c.half#c.lnw1#c.lnw1 |     -0.042      0.014    -3.05    0.00       -0.070      -0.015
                     |
       c.lny1#c.lny2 |     -0.064      0.011    -5.81    0.00       -0.086      -0.043
                     |
       c.lny1#c.lnw1 |      0.002      0.005     0.35    0.73       -0.008       0.012
                     |
       c.lny2#c.lnw1 |     -0.002      0.013    -0.19    0.85       -0.028       0.023
                     |
      c.trend#c.lny1 |      0.003      0.002     1.71    0.09       -0.000       0.006
                     |
      c.trend#c.lny2 |      0.016      0.003     5.18    0.00        0.010       0.021
                     |
      c.trend#c.lnw1 |     -0.007      0.003    -2.95    0.00       -0.012      -0.003
                     |
     c.trend#c.trend#|
              c.half |      0.080      0.002    33.13    0.00        0.075       0.085
                     |
               _cons |      9.501      2.598     3.66    0.00        4.408      14.594
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
                lnta |      0.412      0.085     4.85    0.00        0.246       0.579
               _cons |     -9.308      0.995    -9.36    0.00      -11.258      -7.359
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
              er_ave |     -3.854      2.655    -1.45    0.15       -9.057       1.349
               _cons |     -3.540      0.290   -12.20    0.00       -4.109      -2.971
---------------------+----------------------------------------------------------------
E[ui|z]              |
             llp_ave |      0.000      0.000    12.12    0.00        0.000       0.000
              la_ave |      0.347      0.048     7.22    0.00        0.253       0.442
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |     -0.078      0.021    -3.81    0.00       -0.119      -0.038
               delta |     -0.009      0.004    -2.40    0.02       -0.017      -0.002
--------------------------------------------------------------------------------------
Note, the function of inefficiency change over time is
beta[t] = 1 + gamma*(t-T_i) + delta*(t-T_i)^2: 
 the larger is the beta[t] the larger is the inefficiency

. timer off 11

. predict M11_te, te_jlms_mean
 (872 missing values generated)

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -2.2875846

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -2.2508695

. eststo M11

. quietly predictnl M11_mW1 = $monotonW1 if e(sample)

. quietly predictnl M11_mW2 = 1-M11_mW1 if e(sample)

. quietly predictnl M11_eY1 = $monotonY1 if e(sample)

. quietly predictnl M11_eY2 = $monotonY2 if e(sample)

. generate double M11_rts = 1/(M11_eY1 + M11_eY2)
(872 missing values generated)

. 
. tabstat M11_*, stat(min p25 mean p75 max sd) format(%9.4f)

   Stats |    M11_te   M11_mW1   M11_mW2   M11_eY1   M11_eY2   M11_rts
---------+------------------------------------------------------------
     Min |    1.0083   -0.1097    0.9189   -0.1551    0.2430    0.8731
     p25 |    1.2865   -0.0059    0.9721    0.1460    0.6038    1.0845
    Mean |    1.4430    0.0096    0.9904    0.1739    0.6879    1.1723
     p75 |    1.5357    0.0279    1.0059    0.2073    0.7688    1.2440
     Max |    3.1478    0.0811    1.1097    0.3327    1.2664    1.8977
      SD |    0.2435    0.0271    0.0271    0.0542    0.1252    0.1219
----------------------------------------------------------------------
```

## Model 12

``` stata
timer clear 12
timer on 12
xtsf2g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnvariance(er_ave) vitlnvariance(lnta) model("BC1992") $decimalopt
timer off 12
predict M12_te, te_jlms_mean
estadd scalar AIC = e(aic)
estadd scalar BIC = e(bic)
eststo M12
quietly predictnl M12_mW1 = $monotonW1 if e(sample)
quietly predictnl M12_mW2 = 1-M12_mW1 if e(sample)
quietly predictnl M12_eY1 = $monotonY1 if e(sample)
quietly predictnl M12_eY2 = $monotonY2 if e(sample)
generate double M12_rts = 1/(M12_eY1 + M12_eY2)
tabstat M12_*, stat(min p25 mean p75 max sd) format(%9.4f)

. timer clear 12

. timer on 12

. xtsf2g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnva
> riance(er_ave) vitlnvariance(lnta) model("BC1992") $decimalopt


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =  530.64988  (not concave)
iter 1:  Log-likelihood =  1274.1847  (not concave)
iter 2:  Log-likelihood =  1420.8998  
iter 3:  Log-likelihood =  1505.2472  
iter 4:  Log-likelihood =  1526.9576  
iter 5:  Log-likelihood =  1531.4269  
iter 6:  Log-likelihood =  1531.6695  
iter 7:  Log-likelihood =  1531.6716  
iter 8:  Log-likelihood =  1531.6716  

Sample:----------------------
 Number of obs    = 2546
 Number of groups = 500
Diagnostics:-----------------
 R-squared        = 0.9976
 Adj R-squared    = 0.9970
 AIC              = -2.2825
 BIC              = -2.2458
 Root MSE         = 0.1926
-----------------------------

 The second-generation estimator of
 the cost stochastic frontier model for panel data,
 where effciency is time-varying

--------------------------------------------------------------------------------------
                 lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
---------------------+----------------------------------------------------------------
                lny1 |      0.352      0.140     2.51    0.01        0.077       0.627
                lny2 |     -1.219      0.391    -3.12    0.00       -1.985      -0.453
                lnw1 |      0.196      0.159     1.23    0.22       -0.116       0.507
               trend |     -0.549      0.038   -14.40    0.00       -0.624      -0.474
                     |
c.half#c.lny1#c.lny1 |      0.051      0.005    10.98    0.00        0.042       0.061
                     |
c.half#c.lny2#c.lny2 |      0.219      0.031     7.00    0.00        0.158       0.281
                     |
c.half#c.lnw1#c.lnw1 |     -0.041      0.014    -2.96    0.00       -0.068      -0.014
                     |
       c.lny1#c.lny2 |     -0.063      0.011    -5.74    0.00       -0.085      -0.042
                     |
       c.lny1#c.lnw1 |      0.001      0.005     0.17    0.87       -0.009       0.011
                     |
       c.lny2#c.lnw1 |     -0.002      0.013    -0.17    0.87       -0.027       0.023
                     |
      c.trend#c.lny1 |      0.003      0.002     1.82    0.07       -0.000       0.006
                     |
      c.trend#c.lny2 |      0.015      0.003     5.07    0.00        0.009       0.021
                     |
      c.trend#c.lnw1 |     -0.008      0.003    -3.12    0.00       -0.013      -0.003
                     |
     c.trend#c.trend#|
              c.half |      0.075      0.002    43.53    0.00        0.072       0.079
                     |
               _cons |      9.342      2.585     3.61    0.00        4.276      14.408
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
                lnta |      0.423      0.085     4.97    0.00        0.256       0.589
               _cons |     -9.429      0.993    -9.49    0.00      -11.376      -7.482
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
              er_ave |     -3.388      2.571    -1.32    0.19       -8.427       1.651
               _cons |     -3.542      0.280   -12.63    0.00       -4.092      -2.993
---------------------+----------------------------------------------------------------
E[ui|z]              |
             llp_ave |      0.000      0.000    12.32    0.00        0.000       0.000
              la_ave |      0.374      0.045     8.23    0.00        0.285       0.463
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |      0.032      0.008     3.93    0.00        0.016       0.048
--------------------------------------------------------------------------------------
Note, the function of inefficiency change over time is
beta[t] = exp( -gamma*(t-T_i) ): 
 the larger is the beta[t] the larger is the inefficiency

. timer off 12

. predict M12_te, te_jlms_mean
 (872 missing values generated)

. estadd scalar AIC = e(aic)

added scalar:
                e(AIC) =  -2.2824787

. estadd scalar BIC = e(bic)

added scalar:
                e(BIC) =  -2.2457637

. eststo M12

. quietly predictnl M12_mW1 = $monotonW1 if e(sample)

. quietly predictnl M12_mW2 = 1-M12_mW1 if e(sample)

. quietly predictnl M12_eY1 = $monotonY1 if e(sample)

. quietly predictnl M12_eY2 = $monotonY2 if e(sample)

. generate double M12_rts = 1/(M12_eY1 + M12_eY2)
(872 missing values generated)

. tabstat M12_*, stat(min p25 mean p75 max sd) format(%9.4f)

   Stats |    M12_te   M12_mW1   M12_mW2   M12_eY1   M12_eY2   M12_rts
---------+------------------------------------------------------------
     Min |    1.0090   -0.1067    0.9180   -0.1562    0.2445    0.8826
     p25 |    1.3004   -0.0055    0.9720    0.1471    0.5997    1.0904
    Mean |    1.4585    0.0100    0.9900    0.1749    0.6828    1.1778
     p75 |    1.5551    0.0280    1.0055    0.2080    0.7624    1.2494
     Max |    3.1618    0.0820    1.1067    0.3320    1.2529    1.8928
      SD |    0.2464    0.0266    0.0266    0.0541    0.1233    0.1211
----------------------------------------------------------------------
```

## Results of Models 9-12

Use the *estout* command for this. Note that parameters $\delta$ and $\gamma$ have different interpretation in models 2-4.


``` stata
estout M9 M10 M11 M12, ///
 cells(b(star fmt(%9.4f)) se(par)) ///
 stats(AIC BIC shat RSS ll N sumTi, ///
 labels("AIC" "BIC" "\$\hat\sigma\$" "RSS" "log-likelihood" "\$N\$" "\$\sum T_{i}\$") ///
 fmt(%9.4f %9.4f %9.4f %9.2f %9.2f %9.0f %9.0f)) ///
 starlevels(* 0.10 ** 0.05 *** 0.01) ///
 varlabels(_cons Constant ) ///
 substitute("_ " "Frontier") ///
 legend label collabels(none) mlabels(none) replace
 
------------------------------------------------------------------------------------
Frontier                                                                                  
lny1                       0.3022**        0.6570***       0.3611**        0.3518** 
                         (0.1407)        (0.1432)        (0.1406)        (0.1404)   
lny2                      -1.2486***      -0.5262         -1.2435***      -1.2189***
                         (0.3923)        (0.3977)        (0.3927)        (0.3909)   
lnw1                       0.1669          0.4350***       0.1913          0.1955   
                         (0.1606)        (0.1602)        (0.1587)        (0.1589)   
trend                     -0.5434***      -0.6384***      -0.5774***      -0.5488***
                         (0.0385)        (0.0389)        (0.0394)        (0.0381)   
half # lny1 # lny1         0.0531***       0.0463***       0.0512***       0.0515***
                         (0.0048)        (0.0047)        (0.0047)        (0.0047)   
half # lny2 # lny2         0.2191***       0.1878***       0.2228***       0.2195***
                         (0.0317)        (0.0315)        (0.0314)        (0.0313)   
half # lnw1 # lnw1        -0.0438***      -0.0470***      -0.0424***      -0.0408***
                         (0.0139)        (0.0139)        (0.0139)        (0.0138)   
lny1 # lny2               -0.0602***      -0.0872***      -0.0642***      -0.0634***
                         (0.0111)        (0.0112)        (0.0111)        (0.0110)   
lny1 # lnw1               -0.0001          0.0007          0.0018          0.0009   
                         (0.0053)        (0.0052)        (0.0052)        (0.0052)   
lny2 # lnw1                0.0020         -0.0220*        -0.0024         -0.0021   
                         (0.0130)        (0.0130)        (0.0128)        (0.0129)   
trend # lny1               0.0038**        0.0048***       0.0029*         0.0031*  
                         (0.0017)        (0.0017)        (0.0017)        (0.0017)   
trend # lny2               0.0128***       0.0124***       0.0155***       0.0152***
                         (0.0029)        (0.0030)        (0.0030)        (0.0030)   
trend # lnw1              -0.0079***      -0.0073***      -0.0075***      -0.0079***
                         (0.0025)        (0.0025)        (0.0025)        (0.0025)   
trend # trend # half       0.0761***       0.0940***       0.0800***       0.0754***
                         (0.0017)        (0.0025)        (0.0024)        (0.0017)   
Constant                   9.8233***       3.7296          9.5008***       9.3421***
                         (2.5817)        (2.6552)        (2.5984)        (2.5848)   
------------------------------------------------------------------------------------
ln[var(vit)]                                                                        
lnta                       0.4298***       0.4423***       0.4123***       0.4228***
                         (0.0852)        (0.0861)        (0.0851)        (0.0850)   
Constant                  -9.5016***      -9.6649***      -9.3084***      -9.4288***
                         (0.9958)        (1.0073)        (0.9945)        (0.9934)   
------------------------------------------------------------------------------------
ln[var(ui)]                                                                         
er_ave                    -3.4300         -5.3783*        -3.8540         -3.3880   
                         (2.5888)        (2.8642)        (2.6546)        (2.5708)   
Constant                  -3.3989***      -2.0839***      -3.5403***      -3.5422***
                         (0.2815)        (0.3280)        (0.2903)        (0.2805)   
------------------------------------------------------------------------------------
E[ui|z]                                                                             
llp_ave                    0.0003***       0.0005***       0.0003***       0.0003***
                         (0.0000)        (0.0000)        (0.0000)        (0.0000)   
la_ave                     0.4004***       0.5800***       0.3475***       0.3743***
                         (0.0481)        (0.0947)        (0.0481)        (0.0455)   
------------------------------------------------------------------------------------
beta[t]                                                                             
gamma                                     -0.4181         -0.0784***       0.0319***
                                              (.)        (0.0206)        (0.0081)   
delta                                      0.0722***      -0.0092**                 
                                         (0.0027)        (0.0038)                   
------------------------------------------------------------------------------------
AIC                       -2.2954         -2.2840         -2.2876         -2.2825   
BIC                       -2.2587         -2.2473         -2.2509         -2.2458   
$\hat\sigma$               0.1913          0.1924          0.1921          0.1926   
RSS                        428.34          376.56          414.74          434.83   
log-likelihood            1524.06         1539.31         1535.34         1531.67   
$N$                           500             500             500             500   
$\sum T_{i}$                 2546            2546            2546            2546   
------------------------------------------------------------------------------------
* p<0.10, ** p<0.05, *** p<0.01
```

## Execution Speed (in seconds), Models 9-12

The estimation is extremely fast. Time in seconds:


``` stata
. timer list 9
   9:      0.51 /        1 =       0.5070

. timer list 10
  10:      4.62 /        1 =       4.6180

. timer list 11
  11:      0.67 /        1 =       0.6730

. timer list 12
  12:      0.60 /        1 =       0.6020
```
