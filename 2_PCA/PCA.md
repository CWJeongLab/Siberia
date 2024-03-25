# 2. Principal Component Analysis

## Introduction

We conducted principal components analysis (PCA) with present-day individuals genotyped on the autosomal part of the HumanOrigins array (n=593,124) using smartpca v18140 from EIGENSOFT v8.0.0 (Patterson et al., 2006).  
We used two population sets, the first including present-day Eurasian and American (2,270) and the second including Eurasian only (2,077). 
We projected ancient individuals not included in the PC calculation using the ‘lsqrproject: YES’ option. 
Samples used in PCA are listed in Table S2.


## Run smartpca
    # Prepare paramter file
    pt1=($(pwd)"/")
    fn1="realpath_of_EIGENSTRAT_file"
    of1="PCA"
    tn1="temp1_"${of1}
    for K in $(seq 1 2); do
        of2=${of1}"_"${K}; parf=${of2}".par"
        echo "genotypename: "${fn1}".geno" > ${parf}
        echo "snpname: "${fn1}".snp" >> ${parf}
        echo "indivname: "${fn1}".ind" >> ${parf}
        echo "evecoutname: "${pt1}${of2}".evec" >> ${parf}
        echo "evaloutname: "${pt1}${of2}".eval" >> ${parf}
        echo "poplistname: "${pt1}${of2}".pops" >> ${parf}
        echo -e "altnormstype: NO\nnumoutevec: 20\nnumoutlieriter: 0\nnumoutlierevec: 0" >> ${parf}
        echo -e "outliersigmathresh: 6.0\nnumthreads: 8\nqtmode: 0\nlsqproject: YES" >> ${parf}
    done
    
    # Run smartpca
    for K in $(seq 1 2); do
        of2=${of1}"_"${K}
        smartpca -p ${of2}.par > ${of2}.log
    done      
    
