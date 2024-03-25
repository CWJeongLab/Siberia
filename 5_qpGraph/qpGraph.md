# 5. qpGraph analysis

## Introduction
In order to test component-wise admixture models, graph-based analysis was implemented using the qpgraph function from the R library ADMIXTOOLS2 v2.0.0.
Before graph fitting, f₂ statistics between all pairs of targets were calculated by the extract_f2 function in ADMIXTOOLS2 with the ‘max_miss=0’ option, the same as the ‘allsnps: NO’ option from the previous version.
The number of SNPs remaining by applying this option was 182,628.
  
Mbuti population was also used as an outgroup in this analysis, and the following populations were used for distal representatives: MA1 for ANE; WHG for Mesolithic hunter-gatherers from Europe; USR1 for Native Americans; EastBaikal_N for ANA.  
Then, Middle Holocene populations were systematically added by following orders: irk030, Dzhylinda-1, WestBaikal_EN, Yakutia_MN, Saqqaq, WestBaikal_LNBA, and Yakutia_LN.
The estimated branch length and admixture proportions were converted to dot file by following code, and we plotted admixture graph using Graphviz 6.0.1.

## Run qpGraph in R
    # Run below code in R  
  
    library(admixtools)  
    library(tidyverse)  
    library(magrittr)  
    library(plotly)  
    library(igraph)  
    my_plink_dir <- "realpath of PACKEDPED genotype data"  
    my_pops <- c("Mbuti.DG","WHG","EastBaikal_N","MA1","USR1","irk030","Dzhylinda-1","WestBaikal_LNBA","WestBaikal_EN","Yakutia_MN","Yakutia_LN","Saqqaq.SG")  
    f2 <- f2_from_geno(my_plink_dir,pops=my_pops,maxmiss=0,poly_only=F)  
  
    file <- "edge.csv"  # edge file  
    edge <- read.csv(paste(file,".csv",sep=""),header=F)  
    qpg_results = qpgraph(f2,edge,return_fstats = T)  
    Z = round(qpg_results$f4['z'] %>% slice_min(z,n=1,with_ties=F),2)  
    score = round(qpg_results$score,2)  
    plot_graph(qpg_results$edges,title=paste("worst Z: ",Z,"\nscore: ",score,sep=""))  

    # Save the results as edge file
    edges <- as.data.frame(qpg_results$edges)  
    write.csv(edges,"result_edge.csv", row.names=F,quote=F)  

## Convert edge file into dot file in python
    # Run below code in python
    import pandas as pd
    INPUT = "result_edge.csv"
    df = pd.read_csv(INPUT)

    node_list = list(set(df['from'].to_list() + df['to'].to_list()))
    pop_list = ["R","Mbuti","WHG","EastBaikal_N","MA1","USR1","irk030","Dzhylinda_1","WestBaikal_LNBA",
             "WestBaikal_EN","Yakutia_MN","Yakutia_LN","Saqqaq"] # target populations
    outfile = open(INPUT.replace('_edge.csv','dot'),'w')
    print("digraph G {", file=outfile, end="\n")
    print('size = \"7.5,10\" ;', file=outfile, end="\n")
    print('labelloc = \"t\" ;', file=outfile, end="\n")
    print('label = \"worst Z = -2.77\" ;', file=outfile, end="\n")
    for node in node_list:
    if node not in pop_list:
        print(node, ' [label = \"\" shape=\"point\" size=1 color=\"white\"];',file=outfile,end="\n")

    for i in range(len(df)):
    f, t = df.loc[i,['from','to']]
    if df.loc[i,'type'] == 'edge':
        label = int(round(df.loc[i,'weight']*1000))
        print(f,"->",t, "[ label = ",'\"',label, '\"];',sep="",file=outfile,end="\n")
    else:
        label = int(round(df.loc[i,'weight']*100))
        print(f,"->",t, "[ style=dotted, label = ",'\"',label, '%\"];',sep="",file=outfile, end="\n")
    print("}",file=outfile)
    outfile.close()
    
## Plot the results using Graphviz
    dot -Tpdf result.dot -o result.pdf

If you want to change the color of edge or node, you can find the manual of Graphviz and dot language in here: https://graphviz.org/
