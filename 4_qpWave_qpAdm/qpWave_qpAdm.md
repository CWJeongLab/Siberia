# 4. qpWave and qpAdm analysis

## Introduction

We used qpwave and qpadm functions from the R library ADMIXTOOLS2 v2.0.0. for admixture modeling analysis.  
We used the following populations as a base outgroup set for both qpWave and qpAdm analysis:  
1. Present-day Central African hunter-gatherers Mbuti (1240K: n=5; HumanOrigins: n=10)  
2. Taiwanese Aborigines Ami (1240K: n=2; HumanOrigins: n=10)
3. Native Americans Mixe (1240K: n=5; HumanOrigins: n=10)
4. Indigenous Andamanese islander Onge (1240K: n=2; HumanOrigins: n=11)
5. Early Neolithic Iranians from the Ganj Dareh site Iran_N (n=8) (Lazaridis et al., 2016; Narasimhan et al., 2019)
6. Epipaleolithic European Villabruna (n=1) (Fu et al., 2016)
7. Early Neolithic farmers from western Anatolia Anatolia_N  (n=23) (Mathieson et al., 2015)
8. Early Neolithic northern East Asian Yumin from Inner Mongolia (n=1) (Yang et al., 2020)
9. Neolithic southern Russia West_Siberia_N (n=3) (Narasimhan et al., 2019)
  
In addition, when multiple admixture models were feasible, qpAdm rotating approach, which systematically shifts candidates from source to outgroup, was used to find the best proximal source.  
  
## qpWave analysis in R
    # Run below code in R
    # The code below gives an example of qpWave analysis on UKY, Kolyma_M, and irk030
    library(admixtools)
    my_plink_dir <- "realpath_of_plink_data"
    left <- c("UKY","Kolyma_M","irk030")
    right <- c("Mbuti.DG","Ami.DG","Mixe.DG","Onge.DG","Iran_N","Villabruna","Anatolia_N","Yumin","West_Siberia_N")
    qpwave(my_plink_dir,left,right)
## qpAdm analysis in R
    # Run below code in R
    # The code below gives an example of qpAdm analysis testing WestBaikal_LNBA as an admixture between EastBaikal_N and irk030
    library(admixtools)
    my_plink_dir <- "realpath_of_plink_data"
    left <- c("UKY","Kolyma_M","irk030")
    right <- c("Mbuti.DG","Ami.DG","Mixe.DG","Onge.DG","Iran_N","Villabruna","Anatolia_N","Yumin","West_Siberia_N")
    qpwave(my_plink_dir,left,right)
    left <- c("EastBaikal_N", "irk030")
    right <- c("Mbuti.DG","Ami.DG","Mixe.DG","Onge.DG","Iran_N","Villabruna","Anatolia_N","Yumin","West_Siberia_N")
    target <- c("WestBaikal_LNBA")
    qpadm(my_plink_dir, left, right, target)
