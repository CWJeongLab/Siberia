# 5. qpGraph analysis

## Introduction
In order to test component-wise admixture models, graph-based analysis was implemented using the qpgraph function from the R library ADMIXTOOLS2 v2.0.0.
Before graph fitting, f₂ statistics between all pairs of targets were calculated by the extract_f2 function in ADMIXTOOLS2 with the ‘max_miss=0’ option, the same as the ‘allsnps: NO’ option from the previous version.
The number of SNPs remaining by applying this option was 182,628.
  
Mbuti population was also used as an outgroup in this analysis, and the following populations were used for distal representatives: MA1 for ANE; WHG for Mesolithic hunter-gatherers from Europe; USR1 for Native Americans; EastBaikal_N for ANA.  
Then, Middle Holocene populations were systematically added by following orders: irk030, Dzhylinda-1, WestBaikal_EN, Yakutia_MN, Saqqaq, WestBaikal_LNBA, and Yakutia_LN.
The estimated branch length and admixture proportions were converted to dot file by following code, and we plotted admixture graph using Graphviz 6.0.1.
