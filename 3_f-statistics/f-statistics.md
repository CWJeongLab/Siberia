# 3. f-statistics calculation

## Introduction

We calculated the f-statistics by qp3pop and f₄ functions from the R library ADMIXTOOLS2 v2.0.0. (https://github.com/uqrmaie1/admixtools, publication pending).  
We calculated outgroup-fз using the central African population Mbuti as an outgroup to measure shared genetic drift between target populations.  
Likewise, Mbuti was used as an outgroup to calculate f₄ statistics in the form of f₄(Mbuti, X; target1, target2) for testing symmetricity between targets or searching additional admixture sources.  
Populations used in f-statistics are listed in Table S2, and the results of f-statistics are summarized in Table S7-S9.  

## Calculate outgroup-f3 statistics
    # Run below code in R
    # The code below gives an example of outgroup f3 calculation for Dzhylinda-1
    library(admixtools)
    library(tidyverse)
    library(magrittr)
    library(plotly)
    library(ggplot2)
    library(forcats)
    library(dplyr)
    library(patchwork)
    setwd("/Users/hcgill/projects/Samoyedic/paper/outgroupf3")
    my_plink_dir = "realpath_of_plink_data"
    pop_list <- as.vector(read.table("Fstat_1240K_pop_list.txt") # list of worldwide population 
    target <- "Dzhylnida-1"
    df_3pop <- qp3pop(my_plink_dir, "Mbuti.DG", pop_list$V1,target)
    
    df_3pop <- df_3pop[!(df_3pop$pop2 == target),]
    df_3pop['upper3'] <- df_3pop['est'] + 3*df_3pop['se']
    df_3pop['lower3'] <- df_3pop['est'] - 3*df_3pop['se']
    df_3pop['upper2'] <- df_3pop['est'] + 2*df_3pop['se']
    df_3pop['lower2'] <- df_3pop['est'] - 2*df_3pop['se']
    df_3pop['upper1'] <- df_3pop['est'] + 1*df_3pop['se']
    df_3pop['lower1'] <- df_3pop['est'] - 1*df_3pop['se']
    
    par(mar=c(2,4,2,4))
    df_3pop_ordered <- df_3pop[order(-df_3pop['est']),]
    p <- df_3pop_ordered[c(1:20),] %>%
      mutate(pop2 = fct_reorder(pop2,est)) %>%
      ggplot(aes(x=est, y= pop2)) + geom_pointrange(aes(xmin=lower3, xmax=upper3),size=.5) + 
      geom_pointrange(aes(xmin=lower2, xmax=upper2),size=1.2) + geom_pointrange(aes(xmin=lower1, xmax=upper1),fatten=1.5,size=2) +
      ggtitle(paste('f3(Mbuti;worldwide,',target,')',sep="")) + ylab("") +
      theme(panel.border = element_rect(color = "black", fill = NA, size = 0.7),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(size = 0.15, linetype = "solid", color = "lightgray"),
            panel.grid.minor = element_line(size = 0.15,
                                            linetype = "solid", color = "lightgray"),
            plot.title = element_text(hjust = 0.5),
            axis.text.y = element_text(size = 10,face='bold'),
            legend.position="none")
    
    pdf(paste(target,"_outgroupf3.pdf",sep=""),8,8)
    print(p)
    dev.off()  

## Calculate f4 statistics
    # Run below code in R
    # The code below gives an example of f4 symmetry between Yakutia_MN and WestBaikal_LNBA
    library(admixtools)
    library(tidyverse)
    library(magrittr)
    library(plotly)
    library(ggplot2)
    library(forcats)
    library(dplyr)
    library(patchwork)
  
    my_plink_dir = "realpath_of_plink_data"
    pop_list <- as.vector(read.table("Fstat_1240K_pop_list.txt"))
    
    target1 <- "Yakutia_MN"
    target2 <- "WestBaikal_LNBA"
    df_f4 <- f4(my_plink_dir,"Mbuti.DG",pop_list$V1,target1,target2)
    
    df_f4 <- df_f4[!(df_f4$pop2 == target1),]
    df_f4 <- df_f4[!(df_f4$pop2 == target2),]
    title <- paste("f4symmetry_",target1,"_",target2,".1240K.pdf",sep="")
    plottitle <- paste("f4(Mbuti,worldwide;",target1,",",target2,")",sep="")
    
    df_f4_order <- df_f4[order(-df_f4['est']),]
    df_f4_summarized <- df_f4_order[c(1:15,(nrow(df_f4)-14):nrow(df_f4)),]
    df_f4_summarized$head_tail <- factor(c(rep('upper 15',15),rep('lower 15',15)),levels=c('upper 15','lower 15'))
    p <- df_f4_summarized %>%
      mutate(pop2 = fct_reorder(pop2, est)) %>%
      ggplot(aes(x=est, y= pop2, color=head_tail)) + ylab("") + 
      geom_pointrange(aes(xmin=est-3*se, xmax=est+3*se),size=.5) + geom_pointrange(aes(xmin=est-2*se, xmax=est+2*se),size=1,fatten=.5) +
      geom_pointrange(aes(xmin=est-se, xmax=est+se),size=1.5,fatten=.5) + ggtitle(plottitle) +
      geom_vline(xintercept=0, linetype='dotted', color='black') +
      theme(panel.border = element_rect(color = "black", fill = NA, size = 0.7),
            panel.background = element_rect(fill = "white"),
            panel.grid.major = element_line(size = 0.15, linetype = "solid", color = "lightgray"),
            panel.grid.minor = element_line(size = 0.15,
                                            linetype = "solid", color = "lightgray"),
            axis.text.y = element_text(size = 7, face='bold'),
            plot.title = element_text(hjust = 0.5))
    
    pdf(title,7,7)
    print(p)
    dev.off()

