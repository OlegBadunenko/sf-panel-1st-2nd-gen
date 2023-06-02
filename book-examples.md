# Replication

In this article, the functionality of the commands **xtsf1g** and **xtsf2g** is
showcased. At the same time, the code provided below can be used to replicate the results in the chapter:

Stochastic frontier analysis in Stata: using existing and coding new
commands in “Handbook of Research Methods and Applications in Empirical
Microeconomics” (edited by Nigar Hashimzade and Michael A. Thornton),
2021, Chapter 17, ***Edward Elgar Publishing***, [DOI
<img src="man/figures/doi.png"  width="12" height="12">](https://doi.org/10.4337/9781788976480.00027)

Throughout the article, the subset of the banking data will be used:

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
```

# No determinants of inefficiency

In the first 8 models the specification models don't include determinants of inefficiency.

The inefficiency follow the $u_{it} = \beta(t) u_{i}$. In the 4 models below, $\mathcal{A}(t)$ is defined as follows:

1. $\beta(t) = 1$ (time-invariant inefficiency)
2. $\beta(t) = ( 1 + \exp(\gamma t + \delta t^2) )^-1$ (time-variant inefficiency)
3. $\beta(t) = 1 + \gamma*(t-T_i) + \delta*(t-T_i)^2$ (time-variant inefficiency)
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
xtsf1g lnc $spetech if year > $year_c, cost iterate($itermax)
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

. xtsf1g lnc $spetech if year > $year_c, cost iterate($itermax)


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

----------------------------------------------------------------------------------------
                   lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-----------------------+----------------------------------------------------------------
                  lny1 |   .3965867   .1421524     2.79   0.005     .1179731    .6752004
                  lny2 |  -.6190153   .3989393    -1.55   0.121    -1.400922    .1628912
                  lnw1 |   .3851557   .1663171     2.32   0.021     .0591802    .7111312
                 trend |   -.561285    .040083   -14.00   0.000    -.6398463   -.4827237
                       |
  c.half#c.lny1#c.lny1 |   .0494929    .004811    10.29   0.000     .0400636    .0589222
                       |
  c.half#c.lny2#c.lny2 |   .1782781   .0322045     5.54   0.000     .1151585    .2413978
                       |
  c.half#c.lnw1#c.lnw1 |  -.0686835   .0138597    -4.96   0.000     -.095848    -.041519
                       |
         c.lny1#c.lny2 |  -.0669347   .0114179    -5.86   0.000    -.0893133   -.0445561
                       |
         c.lny1#c.lnw1 |  -.0025056   .0055485    -0.45   0.652    -.0133805    .0083694
                       |
         c.lny2#c.lnw1 |  -.0090899   .0134014    -0.68   0.498    -.0353561    .0171762
                       |
        c.trend#c.lny1 |   .0054134   .0017426     3.11   0.002      .001998    .0088288
                       |
        c.trend#c.lny2 |   .0121805    .003081     3.95   0.000     .0061418    .0182192
                       |
        c.trend#c.lnw1 |  -.0067636    .002597    -2.60   0.009    -.0118537   -.0016736
                       |
c.trend#c.trend#c.half |   .0764848   .0017912    42.70   0.000     .0729741    .0799955
                       |
                 _cons |   5.383199   2.641784     2.04   0.042     .2053973      10.561
-----------------------+----------------------------------------------------------------
ln[var(vit)]           |
                 _cons |  -4.422005    .031957  -138.37   0.000    -4.484639    -4.35937
-----------------------+----------------------------------------------------------------
ln[var(ui)]            |
                 _cons |  -2.521644   .0837447   -30.11   0.000     -2.68578   -2.357507
----------------------------------------------------------------------------------------

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

## Model 2 (2nd generation, time-varying inefficiency using the Kumbhakar specification)

``` stata
timer clear 2
timer on 2
xtsf2g lnc $spetech if year > $year_c, cost iterate($itermax)
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

. xtsf2g lnc $spetech if year > $year_c, cost iterate($itermax)


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

----------------------------------------------------------------------------------------
                   lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-----------------------+----------------------------------------------------------------
                  lny1 |   .4286034   .1402502     3.06   0.002     .1537181    .7034886
                  lny2 |  -.4833135   .3932531    -1.23   0.219    -1.254075    .2874484
                  lnw1 |   .4970021   .1624923     3.06   0.002     .1785231    .8154812
                 trend |  -.6561307   .0403878   -16.25   0.000    -.7352894   -.5769719
                       |
  c.half#c.lny1#c.lny1 |   .0467688   .0045977    10.17   0.000     .0377574    .0557802
                       |
  c.half#c.lny2#c.lny2 |   .1712858   .0314325     5.45   0.000     .1096792    .2328925
                       |
  c.half#c.lnw1#c.lnw1 |  -.0616038    .013398    -4.60   0.000    -.0878635   -.0353441
                       |
         c.lny1#c.lny2 |  -.0685929   .0111672    -6.14   0.000    -.0904802   -.0467056
                       |
         c.lny1#c.lnw1 |   .0001119   .0053354     0.02   0.983    -.0103454    .0105691
                       |
         c.lny2#c.lnw1 |  -.0236568   .0131298    -1.80   0.072    -.0493907    .0020771
                       |
        c.trend#c.lny1 |   .0054366   .0016856     3.23   0.001     .0021328    .0087404
                       |
        c.trend#c.lny2 |    .013075   .0030006     4.36   0.000      .007194     .018956
                       |
        c.trend#c.lnw1 |  -.0057521   .0025103    -2.29   0.022    -.0106723    -.000832
                       |
c.trend#c.trend#c.half |   .0936949   .0029256    32.03   0.000     .0879609     .099429
                       |
                 _cons |   4.496978   2.617188     1.72   0.086    -.6326162    9.626571
-----------------------+----------------------------------------------------------------
ln[var(vit)]           |
                 _cons |   -4.49755   .0323567  -139.00   0.000    -4.560968   -4.434132
-----------------------+----------------------------------------------------------------
ln[var(ui)]            |
                 _cons |  -2.314936   .0888831   -26.04   0.000    -2.489143   -2.140728
-----------------------+----------------------------------------------------------------
beta[t]                |
                 gamma |  -5.000087   1.898288    -2.63   0.008    -8.720663   -1.279512
                 delta |   .8193021   .3145215     2.60   0.009     .2028513    1.435753
----------------------------------------------------------------------------------------
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

## Model 3 (2nd generation, time-varying inefficiency using the modified Kumbhakar specification)


``` stata
timer clear 3
timer on 3
xtsf2g lnc $spetech if year > $year_c, model(K1990modified) cost iterate($itermax)
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

. xtsf2g lnc $spetech if year > $year_c, model(K1990modified) cost iterate($itermax)


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

----------------------------------------------------------------------------------------
                   lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-----------------------+----------------------------------------------------------------
                  lny1 |   .4747436   .1418124     3.35   0.001     .1967963    .7526909
                  lny2 |  -.5632274   .3944716    -1.43   0.153    -1.336377    .2099228
                  lnw1 |   .3932696   .1646438     2.39   0.017     .0705738    .7159655
                 trend |  -.5924608   .0405122   -14.62   0.000    -.6718631   -.5130584
                       |
  c.half#c.lny1#c.lny1 |   .0474519    .004691    10.12   0.000     .0382577     .056646
                       |
  c.half#c.lny2#c.lny2 |   .1796071   .0316048     5.68   0.000     .1176629    .2415514
                       |
  c.half#c.lnw1#c.lnw1 |  -.0663705   .0137783    -4.82   0.000    -.0933754   -.0393656
                       |
         c.lny1#c.lny2 |  -.0726608   .0113244    -6.42   0.000    -.0948563   -.0504653
                       |
         c.lny1#c.lnw1 |   .0000371   .0054986     0.01   0.995      -.01074    .0108141
                       |
         c.lny2#c.lnw1 |   -.012975    .013288    -0.98   0.329     -.039019     .013069
                       |
        c.trend#c.lny1 |   .0047703   .0017335     2.75   0.006     .0013727    .0081679
                       |
        c.trend#c.lny2 |   .0130351   .0030609     4.26   0.000     .0070357    .0190344
                       |
        c.trend#c.lnw1 |  -.0060754   .0025876    -2.35   0.019    -.0111471   -.0010038
                       |
c.trend#c.trend#c.half |      .0823   .0023315    35.30   0.000     .0777304    .0868696
                       |
                 _cons |   4.705405   2.630508     1.79   0.074    -.4502958    9.861105
-----------------------+----------------------------------------------------------------
ln[var(vit)]           |
                 _cons |  -4.435232   .0319327  -138.89   0.000    -4.497819   -4.372645
-----------------------+----------------------------------------------------------------
ln[var(ui)]            |
                 _cons |  -2.756972   .0997669   -27.63   0.000    -2.952511   -2.561432
-----------------------+----------------------------------------------------------------
beta[t]                |
                 gamma |  -.1251908   .0275921    -4.54   0.000    -.1792704   -.0711113
                 delta |  -.0197381   .0052411    -3.77   0.000    -.0300105   -.0094657
----------------------------------------------------------------------------------------
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
xtsf2g lnc $spetech if year > $year_c, model(BC1992) cost iterate($itermax)
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

. xtsf2g lnc $spetech if year > $year_c, model(BC1992) cost iterate($itermax)


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

----------------------------------------------------------------------------------------
                   lnc | Coefficient  Std. err.      z    P>|z|     [95% conf. interval]
-----------------------+----------------------------------------------------------------
                  lny1 |   .4592441   .1426825     3.22   0.001     .1795916    .7388967
                  lny2 |  -.4769618   .3985192    -1.20   0.231    -1.258045    .3041215
                  lnw1 |   .4166255   .1654059     2.52   0.012      .092436    .7408151
                 trend |  -.5561746   .0397535   -13.99   0.000    -.6340901   -.4782592
                       |
  c.half#c.lny1#c.lny1 |   .0481591   .0047315    10.18   0.000     .0388855    .0574327
                       |
  c.half#c.lny2#c.lny2 |   .1712548   .0318872     5.37   0.000      .108757    .2337526
                       |
  c.half#c.lnw1#c.lnw1 |  -.0658764   .0137573    -4.79   0.000    -.0928402   -.0389126
                       |
         c.lny1#c.lny2 |  -.0717502   .0113881    -6.30   0.000    -.0940705   -.0494299
                       |
         c.lny1#c.lnw1 |  -.0012537   .0054734    -0.23   0.819    -.0119815     .009474
                       |
         c.lny2#c.lnw1 |  -.0138682   .0133356    -1.04   0.298    -.0400056    .0122692
                       |
        c.trend#c.lny1 |   .0052777   .0017375     3.04   0.002     .0018724    .0086831
                       |
        c.trend#c.lny2 |   .0124644   .0030718     4.06   0.000     .0064438     .018485
                       |
        c.trend#c.lnw1 |  -.0063633   .0025944    -2.45   0.014    -.0114483   -.0012783
                       |
c.trend#c.trend#c.half |   .0761811   .0017917    42.52   0.000     .0726694    .0796928
                       |
                 _cons |   4.189569   2.657111     1.58   0.115    -1.018272     9.39741
-----------------------+----------------------------------------------------------------
ln[var(vit)]           |
                 _cons |   -4.42903     .03198  -138.49   0.000     -4.49171   -4.366351
-----------------------+----------------------------------------------------------------
ln[var(ui)]            |
                 _cons |  -2.638054    .092451   -28.53   0.000    -2.819255   -2.456854
-----------------------+----------------------------------------------------------------
beta[t]                |
                 gamma |   .0294572   .0093124     3.16   0.002     .0112052    .0477092
----------------------------------------------------------------------------------------
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
xtsf1g lnc $spetech if year > $year_c, cost iterate($itermax) distribution(t)
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

. xtsf1g lnc $spetech if year > $year_c, cost iterate($itermax) distribution(t)


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
                lny1 |   .2077101   .1414406     1.47   0.142    -.0695084    .4849285
                lny2 |  -1.735215   .3981144    -4.36   0.000    -2.515505   -.9549253
                lnw1 |   .0085818   .1781377     0.05   0.962    -.3405616    .3577252
               trend |   -.515125    .038677   -13.32   0.000    -.5909304   -.4393195
                     |
c.half#c.lny1#c.lny1 |   .0475832   .0047279    10.06   0.000     .0383166    .0568498
                     |
c.half#c.lny2#c.lny2 |   .2571491   .0323624     7.95   0.000     .1937199    .3205782
                     |
c.half#c.lnw1#c.lnw1 |  -.0298446   .0147587    -2.02   0.043    -.0587711   -.0009181
                     |
       c.lny1#c.lny2 |  -.0498619   .0111854    -4.46   0.000    -.0717848    -.027939
                     |
       c.lny1#c.lnw1 |   .0002169   .0053375     0.04   0.968    -.0102445    .0106782
                     |
       c.lny2#c.lnw1 |   .0136242   .0140488     0.97   0.332    -.0139109    .0411592
                     |
      c.trend#c.lny1 |   .0044511   .0016869     2.64   0.008     .0011448    .0077575
                     |
      c.trend#c.lny2 |   .0093749   .0029727     3.15   0.002     .0035485    .0152012
                     |
      c.trend#c.lnw1 |  -.0085535   .0025347    -3.37   0.001    -.0135214   -.0035856
                     |
     c.trend#c.trend#|
              c.half |   .0765679   .0017143    44.66   0.000     .0732079    .0799278
                     |
               _cons |   12.53978   3.187609     3.93   0.000     6.292178    18.78738
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
               _cons |  -4.512133   .0314075  -143.66   0.000    -4.573691   -4.450576
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
               _cons |  -3.645746   .0703237   -51.84   0.000    -3.783578   -3.507914
---------------------+----------------------------------------------------------------
E[ui|z]              |
               _cons |   .7969201   1.831099     0.44   0.663    -2.791968    4.385808
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

## Model 6 (2nd generation, time-varying inefficiency using the Kumbhakar specification)

``` stata
timer clear 6
timer on 6
xtsf2g lnc $spetech if year > $year_c, cost iterate($itermax) distribution(t)
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

. xtsf2g lnc $spetech if year > $year_c, cost iterate($itermax) distribution(t)


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
                lny1 |   .1778271   .1359689     1.31   0.191    -.0886671    .4443213
                lny2 |  -1.711872   .3821581    -4.48   0.000    -2.460888   -.9628562
                lnw1 |   .0734777   .1692998     0.43   0.664    -.2583439    .4052992
               trend |  -4.351412   1.593488    -2.73   0.006    -7.474591   -1.228233
                     |
c.half#c.lny1#c.lny1 |    .047979   .0044354    10.82   0.000     .0392858    .0566721
                     |
c.half#c.lny2#c.lny2 |   .2553466   .0308495     8.28   0.000     .1948827    .3158104
                     |
c.half#c.lnw1#c.lnw1 |  -.0233714    .013882    -1.68   0.092    -.0505797    .0038368
                     |
       c.lny1#c.lny2 |   -.047834   .0107315    -4.46   0.000    -.0688673   -.0268007
                     |
       c.lny1#c.lnw1 |   .0013558    .005017     0.27   0.787    -.0084774     .011189
                     |
       c.lny2#c.lnw1 |   .0052038   .0134594     0.39   0.699    -.0211761    .0315837
                     |
      c.trend#c.lny1 |   .0046392   .0015754     2.94   0.003     .0015514     .007727
                     |
      c.trend#c.lny2 |    .009932   .0028001     3.55   0.000     .0044439    .0154201
                     |
      c.trend#c.lnw1 |  -.0087413   .0023697    -3.69   0.000    -.0133859   -.0040968
                     |
     c.trend#c.trend#|
              c.half |   .8308269   .3147648     2.64   0.008     .2138992    1.447755
                     |
               _cons |   10.75969   3.427534     3.14   0.002     4.041847    17.47753
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
               _cons |  -4.637849   .0313975  -147.71   0.000    -4.699387   -4.576311
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
               _cons |  -2.781339   .1215107   -22.89   0.000    -3.019496   -2.543183
---------------------+----------------------------------------------------------------
E[ui|z]              |
               _cons |   17.05362   10.10155     1.69   0.091    -2.745051    36.85229
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |  -.5943576   .1290358    -4.61   0.000    -.8472632    -.341452
               delta |   .0962067   .0205007     4.69   0.000     .0560261    .1363874
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

## Model 7 (2nd generation, time-varying inefficiency using the modified Kumbhakar specification)


``` stata
timer clear 7
timer on 7
xtsf2g lnc $spetech if year > $year_c, model(K1990modified) cost iterate($itermax) distribution(t)
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

. xtsf2g lnc $spetech if year > $year_c, model(K1990modified) cost iterate($itermax) di
> stribution(t)


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
                lny1 |   .2311342    .141483     1.63   0.102    -.0461673    .5084357
                lny2 |  -1.787172   .3984823    -4.48   0.000    -2.568183   -1.006161
                lnw1 |   .0417005   .1789295     0.23   0.816    -.3089949    .3923959
               trend |  -.5345139   .0406138   -13.16   0.000    -.6141154   -.4549123
                     |
c.half#c.lny1#c.lny1 |   .0467602   .0046792     9.99   0.000     .0375893    .0559312
                     |
c.half#c.lny2#c.lny2 |   .2645153   .0323008     8.19   0.000     .2012068    .3278238
                     |
c.half#c.lnw1#c.lnw1 |    -.02998   .0147926    -2.03   0.043     -.058973    -.000987
                     |
       c.lny1#c.lny2 |  -.0515267     .01117    -4.61   0.000    -.0734195   -.0296339
                     |
       c.lny1#c.lnw1 |    .001227   .0053087     0.23   0.817    -.0091779    .0116319
                     |
       c.lny2#c.lnw1 |   .0097781   .0141188     0.69   0.489    -.0178942    .0374504
                     |
      c.trend#c.lny1 |   .0042695   .0016877     2.53   0.011     .0009617    .0075773
                     |
      c.trend#c.lny2 |   .0098409   .0029705     3.31   0.001     .0040188    .0156631
                     |
      c.trend#c.lnw1 |  -.0084499   .0025346    -3.33   0.001    -.0134177   -.0034821
                     |
     c.trend#c.trend#|
              c.half |    .081037   .0026074    31.08   0.000     .0759266    .0861474
                     |
               _cons |   12.81594   2.623768     4.88   0.000     7.673454    17.95843
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
               _cons |  -4.518785   .0314407  -143.72   0.000    -4.580408   -4.457163
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
               _cons |  -3.727864   .0767533   -48.57   0.000    -3.878298   -3.577431
---------------------+----------------------------------------------------------------
E[ui|z]              |
               _cons |   .6192252   .0839178     7.38   0.000     .4547493    .7837012
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |  -.0379511   .0135623    -2.80   0.005    -.0645328   -.0113694
               delta |  -.0048324   .0023034    -2.10   0.036    -.0093471   -.0003177
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
xtsf2g lnc $spetech if year > $year_c, model(BC1992) cost iterate($itermax) distribution(t)
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

. xtsf2g lnc $spetech if year > $year_c, model(BC1992) cost iterate($itermax) distribut
> ion(t)


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
                lny1 |   .2311192   .1415719     1.63   0.103    -.0463567    .5085951
                lny2 |  -1.714386   .3973282    -4.31   0.000    -2.493135    -.935637
                lnw1 |   .0254174   .1782129     0.14   0.887    -.3238735    .3747083
               trend |   -.505287   .0387943   -13.02   0.000    -.5813224   -.4292516
                     |
c.half#c.lny1#c.lny1 |   .0468804   .0046913     9.99   0.000     .0376857    .0560751
                     |
c.half#c.lny2#c.lny2 |   .2578118   .0321879     8.01   0.000     .1947248    .3208988
                     |
c.half#c.lnw1#c.lnw1 |  -.0282167   .0147394    -1.91   0.056    -.0571055     .000672
                     |
       c.lny1#c.lny2 |  -.0516376   .0111736    -4.62   0.000    -.0735375   -.0297377
                     |
       c.lny1#c.lnw1 |   .0007885   .0053044     0.15   0.882    -.0096079    .0111849
                     |
       c.lny2#c.lnw1 |   .0113034   .0140648     0.80   0.422    -.0162632    .0388699
                     |
      c.trend#c.lny1 |   .0045902   .0016842     2.73   0.006     .0012892    .0078911
                     |
      c.trend#c.lny2 |   .0093314   .0029666     3.15   0.002      .003517    .0151457
                     |
      c.trend#c.lnw1 |  -.0087169   .0025341    -3.44   0.001    -.0136837     -.00375
                     |
     c.trend#c.trend#|
              c.half |   .0763892   .0017169    44.49   0.000      .073024    .0797543
                     |
               _cons |   12.31889   2.615637     4.71   0.000      7.19234    17.44545
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
               _cons |  -4.516273   .0314477  -143.61   0.000    -4.577909   -4.454636
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
               _cons |  -3.696693   .0746589   -49.51   0.000    -3.843021   -3.550364
---------------------+----------------------------------------------------------------
E[ui|z]              |
               _cons |    .673796   .1107047     6.09   0.000     .4568189    .8907732
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |   .0138492   .0068658     2.02   0.044     .0003925    .0273059
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

The following specifications are used for the next models. The conditional mean of the truncated distribution can be made observation specific by imposing linear structure on it as in Kumbhakar et al (1991),

$$
\mu_{i} = E[u_{i}|\boldsymbol{z}_{u_{2i}}] = \boldsymbol{z}_{u_{2i}} \boldsymbol{\delta}_{u},
$$

so that $u_{i} \sim N^{+}( \boldsymbol{z}_{u_{2i}} \boldsymbol{\delta}_{u}, \sigma_{u}^{2} )$ and only $\boldsymbol{\delta}_{u}$ parameters are estimated instead of $\mu_{i}, i=1,\ldots, N$. 

With half-normal or truncated normal distribution, the inefficiency term is still *i.i.d.*, which can introduce bias. If the heteroskedasticity function is known, bias can be eliminated by making $\sigma_{u}^{2}$ observation specific (see Caudill et al 95))

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
xtsf1g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnvariance(er_ave) vitlnvariance(lnta)
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

. xtsf1g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnva
> riance(er_ave) vitlnvariance(lnta)


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
                lny1 |   .3021508   .1406718     2.15   0.032     .0264393    .5778624
                lny2 |  -1.248608     .39229    -3.18   0.001    -2.017483   -.4797342
                lnw1 |   .1668581   .1605561     1.04   0.299    -.1478261    .4815423
               trend |   -.543417   .0385162   -14.11   0.000    -.6189073   -.4679267
                     |
c.half#c.lny1#c.lny1 |    .053077   .0047843    11.09   0.000        .0437     .062454
                     |
c.half#c.lny2#c.lny2 |   .2190942      .0317     6.91   0.000     .1569633     .281225
                     |
c.half#c.lnw1#c.lnw1 |  -.0438382   .0139156    -3.15   0.002    -.0711123   -.0165641
                     |
       c.lny1#c.lny2 |  -.0602073   .0111059    -5.42   0.000    -.0819744   -.0384402
                     |
       c.lny1#c.lnw1 |   -.000108   .0052721    -0.02   0.984    -.0104412    .0102251
                     |
       c.lny2#c.lnw1 |   .0020095   .0130036     0.15   0.877    -.0234772    .0274961
                     |
      c.trend#c.lny1 |   .0037667   .0017034     2.21   0.027     .0004282    .0071053
                     |
      c.trend#c.lny2 |   .0128497   .0029479     4.36   0.000     .0070719    .0186275
                     |
      c.trend#c.lnw1 |  -.0078657   .0025196    -3.12   0.002    -.0128041   -.0029274
                     |
     c.trend#c.trend#|
              c.half |   .0760539   .0017229    44.14   0.000     .0726771    .0794307
                     |
               _cons |   9.823269   2.581683     3.80   0.000     4.763264    14.88327
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
                lnta |    .429784   .0852055     5.04   0.000     .2627843    .5967837
               _cons |  -9.501589   .9957925    -9.54   0.000    -11.45331   -7.549871
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
              er_ave |  -3.429966   2.588808    -1.32   0.185    -8.503938    1.644005
               _cons |  -3.398929   .2815485   -12.07   0.000    -3.950754   -2.847104
---------------------+----------------------------------------------------------------
E[ui|z]              |
             llp_ave |   .0003162   .0000268    11.81   0.000     .0002637    .0003687
              la_ave |   .4004424   .0480506     8.33   0.000      .306265    .4946198
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
xtsf2g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnvariance(er_ave) vitlnvariance(lnta) model("K1990")
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
> riance(er_ave) vitlnvariance(lnta) model("K1990")


Log-likelihood maximization using optimize() 

iter 0:  Log-likelihood =  685.59629  (not concave)
iter 1:  Log-likelihood =  1271.6597  (not concave)
iter 2:  Log-likelihood =  1317.3539  (not concave)
iter 3:  Log-likelihood =  1407.2609  
iter 4:  Log-likelihood =   1517.782  (not concave)
...
iter 101: Log-likelihood =  1539.2934  (not concave)
iter 102: Log-likelihood =  1539.3056  (not concave)
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
                lny1 |   .6569742   .1432205     4.59   0.000     .3762672    .9376812
                lny2 |  -.5262071    .397703    -1.32   0.186    -1.305691    .2532763
                lnw1 |   .4349992   .1601837     2.72   0.007     .1210448    .7489536
               trend |    -.63838   .0389399   -16.39   0.000    -.7147008   -.5620591
                     |
c.half#c.lny1#c.lny1 |    .046269   .0047232     9.80   0.000     .0370118    .0555262
                     |
c.half#c.lny2#c.lny2 |   .1878206   .0315442     5.95   0.000     .1259952     .249646
                     |
c.half#c.lnw1#c.lnw1 |  -.0469864   .0138988    -3.38   0.001    -.0742275   -.0197452
                     |
       c.lny1#c.lny2 |  -.0871997   .0112404    -7.76   0.000    -.1092305   -.0651689
                     |
       c.lny1#c.lnw1 |   .0007312    .005196     0.14   0.888    -.0094528    .0109151
                     |
       c.lny2#c.lnw1 |  -.0220436   .0130154    -1.69   0.090    -.0475533     .003466
                     |
      c.trend#c.lny1 |   .0047959   .0016961     2.83   0.005     .0014716    .0081203
                     |
      c.trend#c.lny2 |   .0123547   .0029877     4.14   0.000      .006499    .0182104
                     |
      c.trend#c.lnw1 |  -.0073209    .002492    -2.94   0.003    -.0122051   -.0024368
                     |
     c.trend#c.trend#|
              c.half |   .0939812   .0025075    37.48   0.000     .0890666    .0988958
                     |
               _cons |   3.729634   2.655195     1.40   0.160    -1.474453    8.933722
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
                lnta |   .4423142   .0861097     5.14   0.000     .2735424     .611086
               _cons |  -9.664923   1.007308    -9.59   0.000    -11.63921   -7.690637
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
              er_ave |  -5.378346   2.864169    -1.88   0.060    -10.99201    .2353224
               _cons |  -2.083921   .3280216    -6.35   0.000    -2.726832   -1.441011
---------------------+----------------------------------------------------------------
E[ui|z]              |
             llp_ave |   .0005282   .0000487    10.85   0.000     .0004327    .0006236
              la_ave |   .5800198    .094659     6.13   0.000     .3944915    .7655481
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |  -.4181429          .        .       .            .           .
               delta |   .0722293   .0026967    26.78   0.000      .066944    .0775147
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
xtsf2g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnvariance(er_ave) vitlnvariance(lnta) model("K1990modified")
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
> riance(er_ave) vitlnvariance(lnta) model("K1990modified")


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
                lny1 |    .361094   .1406042     2.57   0.010     .0855148    .6366732
                lny2 |  -1.243547   .3926977    -3.17   0.002     -2.01322   -.4738737
                lnw1 |   .1913169   .1587078     1.21   0.228    -.1197447    .5023785
               trend |  -.5774155   .0394469   -14.64   0.000    -.6547299   -.5001011
                     |
c.half#c.lny1#c.lny1 |   .0512063   .0046788    10.94   0.000     .0420361    .0603766
                     |
c.half#c.lny2#c.lny2 |   .2228153   .0314442     7.09   0.000     .1611859    .2844447
                     |
c.half#c.lnw1#c.lnw1 |  -.0424087   .0138853    -3.05   0.002    -.0696234    -.015194
                     |
       c.lny1#c.lny2 |  -.0641994   .0110566    -5.81   0.000      -.08587   -.0425288
                     |
       c.lny1#c.lnw1 |   .0018054   .0052178     0.35   0.729    -.0084213    .0120321
                     |
       c.lny2#c.lnw1 |  -.0023851   .0128492    -0.19   0.853     -.027569    .0227988
                     |
      c.trend#c.lny1 |   .0029335   .0017122     1.71   0.087    -.0004224    .0062894
                     |
      c.trend#c.lny2 |   .0155042   .0029943     5.18   0.000     .0096355    .0213728
                     |
      c.trend#c.lnw1 |  -.0074652   .0025292    -2.95   0.003    -.0124223   -.0025081
                     |
     c.trend#c.trend#|
              c.half |    .080009   .0024149    33.13   0.000     .0752758    .0847422
                     |
               _cons |    9.50084   2.598404     3.66   0.000     4.408063    14.59362
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
                lnta |   .4123455   .0850729     4.85   0.000     .2456056    .5790853
               _cons |  -9.308407   .9945235    -9.36   0.000    -11.25764   -7.359177
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
              er_ave |  -3.854034   2.654558    -1.45   0.147    -9.056873    1.348804
               _cons |   -3.54028   .2902774   -12.20   0.000    -4.109213   -2.971346
---------------------+----------------------------------------------------------------
E[ui|z]              |
             llp_ave |   .0002999   .0000247    12.12   0.000     .0002514    .0003484
              la_ave |   .3474996   .0481268     7.22   0.000     .2531728    .4418265
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |  -.0783596   .0205507    -3.81   0.000    -.1186381    -.038081
               delta |  -.0091516   .0038063    -2.40   0.016    -.0166117   -.0016915
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
xtsf2g lnc $spetech if year > $year_c, cost uimean(llp_ave la_ave, noconstant) uilnvariance(er_ave) vitlnvariance(lnta) model("BC1992")
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
> riance(er_ave) vitlnvariance(lnta) model("BC1992")


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
                lny1 |   .3517897   .1403674     2.51   0.012     .0766746    .6269048
                lny2 |  -1.218866   .3909032    -3.12   0.002    -1.985022   -.4527098
                lnw1 |   .1955138   .1589117     1.23   0.219    -.1159473     .506975
               trend |  -.5488284   .0381213   -14.40   0.000    -.6235448   -.4741121
                     |
c.half#c.lny1#c.lny1 |   .0514759   .0046869    10.98   0.000     .0422897     .060662
                     |
c.half#c.lny2#c.lny2 |   .2194645   .0313478     7.00   0.000     .1580239     .280905
                     |
c.half#c.lnw1#c.lnw1 |  -.0408463   .0137916    -2.96   0.003    -.0678774   -.0138152
                     |
       c.lny1#c.lny2 |  -.0633507   .0110306    -5.74   0.000    -.0849703   -.0417311
                     |
       c.lny1#c.lnw1 |   .0008703   .0051906     0.17   0.867    -.0093031    .0110438
                     |
       c.lny2#c.lnw1 |  -.0021251   .0128614    -0.17   0.869    -.0273329    .0230827
                     |
      c.trend#c.lny1 |   .0031144   .0017108     1.82   0.069    -.0002388    .0064675
                     |
      c.trend#c.lny2 |   .0152092   .0029994     5.07   0.000     .0093306    .0210879
                     |
      c.trend#c.lnw1 |  -.0078867   .0025243    -3.12   0.002    -.0128343   -.0029392
                     |
     c.trend#c.trend#|
              c.half |   .0754254   .0017327    43.53   0.000     .0720293    .0788214
                     |
               _cons |   9.342057   2.584798     3.61   0.000     4.275946    14.40817
---------------------+----------------------------------------------------------------
ln[var(vit)]         |
                lnta |   .4227814       .085     4.97   0.000     .2561845    .5893782
               _cons |  -9.428763   .9934331    -9.49   0.000    -11.37586    -7.48167
---------------------+----------------------------------------------------------------
ln[var(ui)]          |
              er_ave |  -3.387995   2.570759    -1.32   0.188     -8.42659    1.650599
               _cons |  -3.542199   .2804514   -12.63   0.000    -4.091874   -2.992525
---------------------+----------------------------------------------------------------
E[ui|z]              |
             llp_ave |   .0003079    .000025    12.32   0.000      .000259    .0003569
              la_ave |   .3743385   .0454699     8.23   0.000     .2852192    .4634578
---------------------+----------------------------------------------------------------
beta[t]              |
               gamma |   .0318869   .0081228     3.93   0.000     .0159665    .0478074
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
