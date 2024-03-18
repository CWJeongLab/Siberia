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
