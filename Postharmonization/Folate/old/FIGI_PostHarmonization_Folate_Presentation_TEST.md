FIGI Folate Post Harmonization
========================================================
author: Andre Kim
date: 11/26/2018
autosize: true

Folate Harmonization
========================================================
<br>

FFQ multi-step data harmonization

reference time = enrollment (cohort)

dietary folate equivalent (DFE) - mcg folate occurring in foods

Folate Harmonization
========================================================
### Dietary Folate
- prior fortification: mcg natural food folate
- after fortification: mcg natural food folate + 1.7 * mcg folic acid from fortified foods

### Supplemental Folate
- folic acid from supplements (single or multivitamins)
- for studies with binary variables (regular user yes/no), assume 400mcg/day or 400mcg/tablets

### Total Folate
- Total = mcg dietary folate + 1.7 folic acid from supplements



Slide With Code
========================================================


```r
summary(cars)
```

```
     speed           dist       
 Min.   : 4.0   Min.   :  2.00  
 1st Qu.:12.0   1st Qu.: 26.00  
 Median :15.0   Median : 36.00  
 Mean   :15.4   Mean   : 42.98  
 3rd Qu.:19.0   3rd Qu.: 56.00  
 Max.   :25.0   Max.   :120.00  
```

Slide With Plot
========================================================

![plot of chunk unnamed-chunk-2](FIGI_PostHarmonization_Folate_Presentation_TEST-figure/unnamed-chunk-2-1.png)
