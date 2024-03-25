# 1. Genotype Data Preparation  

## Introduction  
We prepared the genotype data of modern and ancient target individuals whose genotype data are not available but fastq or BAM files are available.  
We used two SNP panel, 1240K(Mathieson et al., 2015; Fu et al., 2016) and Human Origins(Patterson et al., 2012; Lazaridis et al., 2016; Flegontov et al., 2017; Jeong et al., 2019).  
The genotype data of ancient individuals for the 1240K panel, excluding those for whom we took publicly available 1240K panel genotype calls, have been deposited in the Edmond Data Repository of the Max Planck Society [https://edmond.mpg.de/dataset.xhtml?persistentId=doi:10.17617/3.QZBM1X].

You can run the whole process with the two below wrapper scripts  
The first wrapper script runs following steps:  
1. AdapterRemoval
2. Bwa mapping and duplicate removal
3. MapDamage
4. Coverage calculation
5. X-based contamination rate estimation (ANGSD)
6. MT-based contamination rate estimation (schmutzi)

The second wrapper script runs following steps:
1. Trimming Bam files
2. PileupCaller
3. MT haplogroup assignment
4. Y haplogroup assignment

Before you run the wrapper scripts, you have to prepare the list file containing each sample's sequencing information.
You can check the example list file (Please check Sample_example.txt).
Also, you need the reference genome list file (Please check reflist.txt).
If you already have BAM files, then you can skip the AdapterRemoval step.

## Wrapper_1
    pt1=($(pwd)"/")
    
    listf="Sample_example.txt"   ## A list of samples to be processed (LID, run, FastQs, type)
    reff="reflist.txt"    ## the reference genome list file
    inum=1                ## the number of sample to be processed
    
    iid=($(awk -v inum="$inum" '{if (inum == NR) print $1}' ${listf}))   ## Individual ID
    lid=($(awk -v inum="$inum" '{if (inum == NR) print $2}' ${listf}))   ## library ID
    rid=($(awk -v inum="$inum" '{if (inum == NR) print $3}' ${listf}))   ## run ID
    tfqs=($(awk -v inum="$inum" '{if (inum == NR) print $4}' ${listf}))  ## comma separated FastQ files
    type=($(awk -v inum="$inum" '{if (inum == NR) print $5}' ${listf}))  ## library type
    spp=($(awk -v inum="$inum" '{if (inum == NR) print $6}' ${listf}))   ## comma separated species list
    
    
    seedval="32"; if [[ "$type" == "ssLib" ]] || [[ "$type" == "non-UDG" ]]; then seedval="9999"; fi
    nrmax=100000
    
    ## Create directories
    echo -e "0. Create individual directory\n\n"
    mkdir -p ./${iid}/FastQ/ ./${iid}/BAM/ ./${iid}/qualimap/ ./${iid}/mapDamage/
    mkdir -p ./${iid}/coverage/ ./${iid}/1240Kbam ./${iid}/Xcont ./${iid}/schmutzi
    
    ##############################################################
    ## Part 1. Run AdapterRemoval (keep 35 bp or longer reads)  ##
    
    ## Common options for AdapterRemoval; keep minlength 35 bps
    optstr="--gzip --threads 8 --trimns --trimqualities "
    optstr+="--adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC "
    optstr+="--adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA "
    optstr+="--minlength 35 --minquality 20 --minadapteroverlap 1"
    
    cd ${pt1}${iid}/FastQ/
    mkdir -p ${rid}; cd ${rid}
    
    echo -e "1. Adapter Removal\n"
    
    ## Take a list of FastQ files and decide if they are SE or PE
    stval=($(echo ${tfqs} | sed s/","/"\n"/g | awk '$1 ~ /_R2/ || $1 ~ /_2.fastq.gz/ || $1 ~ /_2.fq.gz/' | wc -l | awk '{if ($1 > 0) print "PE"; else print "SE"}'))
    
    echo -e "This library was "${stval}" sequenced.\n"
    
    fq1s=""; while read fq; do
        if [[ "$fq1s" == "" ]]; then fq1s+=${fq}; else fq1s+=" "${fq}; fi
    done < <(echo ${tfqs} | sed s/","/"\n"/g | awk '$1 !~ /_R2/ && $1 !~ /_2.fastq.gz/ && $1 !~ /_2.fq.gz/' )
    if [[ "$stval" == "PE" ]]; then
        fq2s=""; while read fq; do
            if [[ "$fq2s" == "" ]]; then fq2s+=${fq}; else fq2s+=" "${fq}; fi
        done < <(echo ${tfqs} | sed s/","/"\n"/g | awk '$1 ~ /_R2/ || $1 ~ /_2.fastq.gz/ || $1 ~ /_2.fq.gz/')
    fi
    
    echo -e "Fastq1s:  "${fq1s}", Fastq2s: "${fq2s}"\n"
    
    fqout=${rid}".L35.fq"
    
    ## Run AdapterRemoval
    if [[ "$stval" == "SE" ]]; then
        AdapterRemoval --file1 ${fq1s} --basename ${fqout} ${optstr}
    else
        AdapterRemoval --file1 ${fq1s} --file2 ${fq2s} --basename ${fqout} ${optstr} --collapse
        zcat ${fqout}.collapsed.gz ${fqout}.collapsed.truncated.gz ${fqout}.singleton.truncated.gz \
        | gzip > ${fqout}.combined.gz
    fi
    
    echo -e "Adapter Removal is complete\n\n\n\n"
    
    #######################################################################
    ## Part 2. Map adapter-removed FastQ files to the reference genome   ##
    ##         and remove duplicates and filter reads                    ##
    
    ## DeDup was changed from 0.12.5 to 0.12.8 in the new setup
    dedup="java -Xmx8192m -jar /opt/ohpc/pub/apps/dedup/0.12.8/DeDup-0.12.8.jar"
    
    ## Move into the target directory
    cd ${pt1}${iid}/BAM/
    
    echo -e "2. Read Mapping\n"
    
    ## Input FastQ (after running AdapterRemoval)
    if [[ "$stval" == "SE" ]]; then
        fqin=${pt1}${iid}"/FastQ/"${rid}"/"${fqout}".truncated.gz"
    else
        fqin=${pt1}${iid}"/FastQ/"${rid}"/"${fqout}".combined.gz"
    fi
    
    echo -e ${fqin}" is used for mapping\n"
    
    ## Repeat mapping steps for each reference species requested
    echo -e "Mapping on each reference\n"
    snum=0
    for sp in $(echo ${spp} | sed s/","/" "/g); do
        let snum+=1
        echo -e "2."${snum}" Read Mapping on "${sp}"\n"
    
        ref=($(awk -v sp="$sp" '{if ($1 == sp) print $2}' ${reff}))  ## Take the reference genome .fa file
        sid=${rid}"."${sp}
    
        ## Make directory for each run id + reference
        mkdir -p ${sid}; cd ${sid}
    
        ## Setup output BAM file prefix
        of1=${sid}".L35.mapped"
    
        ## Define read group
        RG="@RG\tID:"${rid}"\tSM:"${lid}"\tLB:"${rid}"\tPL:illumina"
    
        ## Run BWA aln (-l 9999 for non-UDG, -l 32 for UDG-half)
        bwa aln -t 8 -n 0.01 -l ${seedval} -f ${of1}.sai ${ref} ${fqin}
    
        ## Run BWA samse and filter out unmapped reads
        bwa samse -r ${RG} ${ref} ${of1}.sai ${fqin} | samtools view -h -F 0x0004 -o ${of1}.0.bam -
    
        ## Sort the output file and index
        samtools sort -@ 8 -m 15G ${of1}.0.bam -o ${of1}.bam
        samtools index ${of1}.bam
    
        ## Remove temporary files
        rm ${of1}.0.bam
        rm ${of1}.sai
    
        ## Remove duplicates using dedup
        ${dedup} -i ${of1}.bam -m -o ./
        mv ${of1}_rmdup.bam ${of1}.rmdup.bam
        samtools index ${of1}.rmdup.bam
        rm ${of1}.dedup.json
    
        mv ${of1}.hist ${of1}.rmdup.hist
        mv ${of1}.log ${of1}.rmdup.log
    
        ## Apply quality filter (-q30)
        samtools view -bh -q30 -o ${of1}.rmdup.q30.bam ${of1}.rmdup.bam
        samtools index ${of1}.rmdup.q30.bam
    
        echo -e "\nMapping on "${sp}" is finished\n\n\n"
    
        ## Move to BAM directory
        cd ..
    
    done
    
    echo -e "Mapping is totally complete\n\n\n\n"
    
    ############################
    ## Part 3. Run mapDamage  ##
    
    cd ${pt1}${iid}/mapDamage/
    
    echo -e "3. Check mapDamage\n"
    
    for sp in $(echo ${spp} | sed s/","/" "/g); do
        ref=($(awk -v sp="$sp" '{if ($1 == sp) print $2}' ${reff}))  ## Take the reference genome .fa file
        sid=${rid}"."${sp}
    
        ibam=($(realpath ${pt1}${iid}/BAM/${sid}/${sid}.L35.mapped.rmdup.q30.bam))
    
        ## Run mapDamage
        mapDamage -i ${ibam} -r ${ref} -d ./${sid} --merge-reference-sequences --no-stat -t ${sid} -n ${nrmax}
    done
    
    echo -e "mapDamage is complete \n\n\n\n"
      
    ###########################################################################
    ## Part 4. Calculate coverage on 1240K sites & mtDNA using unmasked data ##
    
    echo -e "4. Calculate coverage on 1240K sites and mtDNA"
    
    cd ${pt1}${iid}/coverage/
    
    for sp in $(echo ${spp} | sed s/","/" "/g); do
        ref=($(awk -v sp="$sp" '{if ($1 == sp) print $2}' ${reff}))  ## Take the reference genome .fa file
        sid=${rid}"."${sp}
        mkdir ${sid}; cd ${sid}
        ibam=($(realpath ${pt1}${iid}/BAM/${sid}/${sid}.L35.mapped.rmdup.q30.bam))
        of1=${sid}".1240K.coverage.txt"
        tn1="temp1_"${sid}".1240K.coverage.txt"
    
        samtools depth -q30 -Q30 -aa -b ${bedf} ${ibam} >> ${tn1}_1  ## 1240K sites
        samtools depth -q30 -Q30 -aa -r MT ${ibam} >> ${tn1}_1       ## MT
        awk 'BEGIN { xR = 0; yR = 0; aR = 0; mr = 0; xS = 0; yS = 0; aS = 0; mS = 0 }
         { chr = $1; pos = $2; cov = $3;
           if (chr == "X") { xR += cov; xS += 1 }
           else if (chr == "Y") { yR += cov; yS += 1 }
           else if (chr == "MT") { mR += cov; mS += 1 }
           else { aR += cov; aS += 1 }
         }
         END { OFS="\t"
           print("xCoverage", xS > 0 ? xR / xS : 0)
           print("yCoverage", yS > 0 ? yR / yS : 0)
           print("autCoverage", aS > 0 ? aR / aS : 0)
           print("mtCoverage", mS > 0 ? mR / mS : 0)
         }' ${tn1}_1 > ${of1}
        rm ${tn1}_1
    
        cd ..
    done
    
    echo -e "Coverage calculation is complete \n\n\n\n"
    
    ##########################################################
    ## Part 5. X-based contamination estimation using ANGSD ##
    
    echo -e "5. X-based contamination estimation"
    
    cd ${pt1}${iid}/Xcont/
    
    for sp in $(echo ${spp} | sed s/","/" "/g); do
        ref=($(awk -v sp="$sp" '{if ($1 == sp) print $2}' ${reff}))  ## Take the reference genome .fa file
        sid=${rid}"."${sp}
        mkdir ${sid}; cd ${sid}
    
        ibam=($(realpath ${pt1}${iid}/1240Kbam/${sid}/${sid}.L35.1240K.mapped.rmdup.q30.bam))
        of1="./"${sid}".angsdCounts"
        angsd -i ${ibam} -r X:5000000-154900000 -doCounts 1 -iCounts 1 -minMapQ 30 -minQ 30 -out ${of1}
    
        /opt/ohpc/pub/apps/angsd/misc/contamination -a ${of1}.icnts.gz \
        -h /opt/ohpc/pub/apps/angsd/RES/HapMapChrX.gz 2> ${of1}
        cd ..
    done
    
    echo -e "Xcont estimation is complete \n\n\n\n"
    
    ##########################################################
    ## Part 6. MT-based contamination estimation using ANGSD ##
    
    cd ${pt1}${iid}/schmutzi
    
    lenDeam="2"; if [[ "$type" == "ssLib" ]] || [[ "$type" == "non-UDG" ]]; then lenDeam="20"; fi
    lty="double"; if [[ "$type" == "ssLib" ]]; then lty="single"; fi
    MTref1="/home/References/Human/hg19_MT.fasta"
    MTref2="/home/References/Human/hg19_MT_500.fasta"
    ctgn="NC_012920.1"
    nrmax=50000  ## Use max 50,000 reads
    
    for sp in $(echo ${spp} | sed s/","/" "/g); do
        ref=($(awk -v sp="$sp" '{if ($1 == sp) print $2}' ${reff}))  ## Take the reference genome .fa file
        sid=${rid}"."${sp}
        mkdir ${sid}; cd ${sid}
    
        ibam=($(realpath ${pt1}${iid}/1240Kbam/${sid}/${sid}.L35.1240K.mapped.rmdup.q30.bam))
        rg="@RG\tID:"${rid}"\tSM:"${lid}"\tLB:"${rid}"\tPL:illumina"
        of1=${sid}".circmapper"
        bam1=${of1}".rmdup.q30.MD.small.bam"
        tn1="temp1_"${of1}
    
        mkdir -p ./circularMT; mkdir -p ./schmutzi; cd ./circularMT/
        samtools view -h ${ibam} MT | awk -v ctgn="$ctgn" '{OFS="\t"} {if ($1 == "@SQ" && $2 ~ /MT$/) print $1,"SN:"ctgn,$3;
        else if ($1 == "@SQ"); else if ($1 ~ /^@/) print $0; else { $3=ctgn; print $0 } }' | samtools view -bh -o ${of1}.rmdup.q30.bam -
        samtools index ${of1}.rmdup.q30.bam
    
        samtools fillmd -b ${of1}.rmdup.q30.bam ${MTref1} > ${of1}.rmdup.q30.MD.bam
        samtools index ${of1}.rmdup.q30.MD.bam
    
        ## Downsample reads to run Schmutzi easily (capped to 30,000 reads)
        nr1=($(samtools view -c ${of1}.rmdup.q30.MD.bam))
        pr1=($(echo ${nr1} | awk -v nrmax="$nrmax" '{printf "%.3f\n", nrmax/$1}'))
        if [ "$nr1" -le "$nrmax" ]; then
            cp ${of1}.rmdup.q30.MD.bam ${of1}.rmdup.q30.MD.small.bam
            samtools index ${of1}.rmdup.q30.MD.small.bam
        else
            samtools view -bh -s ${pr1} -o ${of1}.rmdup.q30.MD.small.bam ${of1}.rmdup.q30.MD.bam
            samtools index ${of1}.rmdup.q30.MD.small.bam
        fi
    
        ## Move into Schmutzi directory
        cd ../schmutzi/
    
        ## Run contDeam (Schmutzi step 1)
        contDeam.pl --library ${lty} --lengthDeam ${lenDeam} --out ./${of1}_nolen ${MTref1} ../circularMT/${bam1}
    
        ## Run schmutzi (Schmutzi step 2)
        schmutzi.pl -t 8 --notusepredC --lengthDeam ${lenDeam} --ref ${MTref1} ./${of1}_nolen /opt/ohpc/pub/apps/schmutzi/share/schmutzi/alleleFreqMT/197/freqs/ ../circularMT/${bam1}
    
        cd ..
    done
    
    echo -e "schmutzi is complete \n\n\n\n"
    
## Wrapper_2
    pt1=($(pwd)"/")
    
    listf="Sample_example.txt"   ## A list of samples to be processed (LID, run, FastQs, type)
    reff="reflist.txt"    ## the reference genome list file
    inum=1                ## the number of sample to be processed
    
    iid=($(awk -v inum="$inum" '{if (inum == NR) print $1}' ${listf}))   ## Individual ID
    lid=($(awk -v inum="$inum" '{if (inum == NR) print $2}' ${listf}))   ## library ID
    rid=($(awk -v inum="$inum" '{if (inum == NR) print $3}' ${listf}))   ## run ID
    tfqs=($(awk -v inum="$inum" '{if (inum == NR) print $4}' ${listf}))  ## comma separated FastQ files
    type=($(awk -v inum="$inum" '{if (inum == NR) print $5}' ${listf}))  ## library type
    spp=($(awk -v inum="$inum" '{if (inum == NR) print $6}' ${listf}))   ## comma separated species list
    mask=($(awk -v inum="$inum" '{if (inum == NR) print $7}' ${listf}))  ## the base pair to be trimmed
    
    ## 0. Make directories
    
    mkdir -p ./${iid}/Yhap
    
    ## 1. Trim the bam files 
    cd ${pt1}${iid}
    if [[ "$type" != "ssLib" ]]; then
        for sp in $(echo ${spp} | sed s/","/" "/g); do
            sid=${rid}"."${sp}
            ibam1="./BAM/${sid}/${sid}.L35.mapped.rmdup.q30.bam"
            obam1="./BAM/${sid}/${sid}.L35.mapped.rmdup.q30.m"${mask}".bam"
            bam trimBAM ${ibam1} ${obam1} ${mask}
            samtools index ${obam1}
        done
    fi
    
    
    ## 2. Run pileupCaller
    
    cd ${pt1}${iid}/genotypes
    ref="/home/References/Human/hs37d5.fa"
    posf="/home/References/Human/SNPCapBEDs/1240K.pos.list.txt"
    snpf1="/home/References/Human/SNPCapBEDs/1240K.snp"
    snpf2="/home/References/Human/SNPCapBEDs/HumanOrigins.snp"
    for sp in $(echo ${spp} | sed s/","/" "/g); do
        sid=${rid}"."${sp}
        mkdir ${sid}; cd ${sid}
        covf="${pt1}${iid}/coverage/${sid}/${sid}.1240K.coverage.txt"
        of1=${sid}".1240K"
        of2=${sid}".HO"
        tn1="temp1_"${of1}
        tn2="temp1_"${of2}
        ibam1="${pt1}${iid}/BAM/${sid}/${sid}.L35.mapped.rmdup.q30.bam"
        if [[ "$type" != "ssLib" ]]; then
            ## Run pileupCaller for raw bam file
            samtools mpileup -B -R -q30 -Q30 -l ${posf} -f ${ref} ${ibam1} > ${tn1}_1
            cut -f 1 ${tn1}_1 | sed s/'X'/'23'/g | sed s/'Y'/'24'/g | paste - ${tn1}_1 | cut -f 1,3- > ${of1}.pileup
            rm ${tn1}_1
            pileupCaller --randomHaploid --sampleNames ${rid} -f ${snpf1} -e ${tn1}_1 < ${of1}.pileup
    
            ## Run pileupCaller for trimmed bam file
            ibam2="${pt1}${iid}/BAM/${sid}/${sid}.L35.mapped.rmdup.q30.m${mask}.bam"
            samtools mpileup -B -R -q30 -Q30 -l ${posf} -f ${ref} ${ibam2} > ${tn1}_2
            cut -f 1 ${tn1}_2 | sed s/'X'/'23'/g | sed s/'Y'/'24'/g | paste - ${tn1}_2 | cut -f 1,3- > ${of1}.masked.pileup
            rm ${tn1}_2
            pileupCaller --randomHaploid --sampleNames ${rid} -f ${snpf1} -e ${tn1}_2 < ${of1}.masked.pileup
    
            ## Choose calls from masked data for C/T and G/A SNPs
            paste ${snpf1} ${tn1}_1.geno ${tn1}_2.geno | \
            awk '{if ($5$6 == "CT" || $5$6 == "TC" || $5$6 == "GA" || $5$6 == "AG") print $8; else print $7}' > ${of1}.geno        
        else
            ## Run pileupCaller for ibam1 only but with SingleStrandMode option
            samtools mpileup -B -R -q30 -Q30 -l ${posf} -f ${ref} ${ibam1} > ${tn1}_1
            cut -f 1 ${tn1}_1 | sed s/'X'/'23'/g | sed s/'Y'/'24'/g | paste - ${tn1}_1 | cut -f 1,3- > ${of1}.pileup
            rm ${tn1}_1
            pileupCaller --randomHaploid --singleStrandMode --sampleNames ${rid} -f ${snpf1} -e ${tn1}_1 < ${of1}.pileup
            mv ${tn1}_1.geno ${of1}.geno
        fi
        
        ## Update Sex information
        Ycov=($(cut -f 2 ${covf} | head -n 2 | tail -n 1))
        Acov=($(cut -f 2 ${covf} | head -n 3 | tail -n 1))
        echo $Ycov $Acov $sid | awk '{if ($1/($2 + 0.000001) > 0.3) print $3"\tM\t"$3;
        else if ($1/($2 + 0.000001) < 0.1) print $3"\tF\t"$3; else print $3"\tU\t"$3}' > ${of1}.ind    
    
        ## Extract HO snps
        awk '{print $1}' ${snpf2} > ${tn2}_1
        paste ${snpf1} ${of1}.geno | fgrep -wf ${tn2}_1 - | cut -f 7 > ${of2}.geno
        paste ${snpf1} ${of1}.geno | fgrep -wf ${tn2}_1 - | awk '{if ($2 <= 22) print $7}' > ${of2}.auto.geno
    
        ## Calculate the number of covered SNPs
        echo -e 'IID\tn1240K\tp1240K\tnHO\tpHO' > ${sid}.SNPcov.txt
        tnum1=($(awk '$1 != 9' ${of1}.geno | wc -l))
        tnum2=($(awk '$1 != 9' ${of2}.geno | wc -l))
        echo ${iid}" "${tnum1}" "${tnum2} | awk '{print $1,$2,$2/1233013,$3,$3/597573}' | sed s/" "/"\t"/g >> ${sid}.SNPcov.txt
    
        ## Remove temporary files
        rm ${tn1}*; rm ${tn2}*; rm *.pileup
    
        cd ..
    done
    
    
    ## 3. Run haplogrep for assigning MT haplogroups
    
    cd ${pt1}${iid}/schmutzi
    
    for sp in $(echo ${spp} | sed s/","/" "/g); do
        sid=${rid}"."${sp}
        cd ${sid}
        ## Retrieve consensus sequences
        faout1=${sid}".MTendo.schmutzi.nolen.fasta"
        faout2=${sid}".MTcont.schmutzi.nolen.fasta"
        logf1="./schmutzi/${sid}.circmapper_nolen_final_endo.log" ## log file for endogenous seq
        logf2="./schmutzi/${sid}.circmapper_nolen_final_cont.log" ## log file for contaminant seq
    
        for J in 10 20 30; do hv=">"${sid}".q"${J}
            if [ -f "$logf1" ]; then
                echo ${hv} >> ${faout1}; log2fasta -q ${J} ${logf1} | tail -n +2 >> ${faout1}
            fi
            if [ -f "$logf2" ]; then
                echo ${hv} >> ${faout2}; log2fasta -q ${J} ${logf2} | tail -n +2 >> ${faout2}
            fi
        done
        ## Run haplogrep in command line (works well!!)
        for fa in ${faout1} ${faout2}; do
            java -jar /opt/ohpc/pub/apps/haplogrep/haplogrep-2.1.20.jar --in ${fa} --format fasta --out ${fa}.txt
        done
        
        echo -e "Endo_q10\tEndo_q20\tEndo_q30\tCont_q10\tCont_q20\tCont_q30" > ${sid}.MThap.txt
        cat ${faout1}.txt ${faout2}.txt | sed s/'"'//g | awk '{if (NR % 4 != 1) print $3}' | tr '\n' '\t' | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6}' >> ${sid}.MThap.txt
        
        cd ..
    done
    
    
    ## 4. Run Yhaplo for assigning Y haplogroups (if the sample is male)
    
    
    cd ${pt1}${iid}/Yhap
    ref="/home/References/Human/hs37d5.fa"
    posf="/opt/ohpc/pub/apps/yhaplo/input/ISOGG_chrY_cleanup_170802_hs37d5.list"
    snpf="/opt/ohpc/pub/apps/yhaplo/input/ISOGG_chrY_cleanup_170802.snp"
    dataf="/opt/ohpc/pub/apps/yhaplo/input/ISOGG_chrY_cleanup_170802.txt"
    py="/opt/ohpc/pub/apps/yhaplo/callHaplogroups.py"
     
    for sp in $(echo ${spp} | sed s/","/" "/g); do
        sid=${rid}"."${sp}
        indf="${pt1}${iid}/genotypes/${sid}/${sid}.1240K.ind"
        sex=($(cut -f 2 ${indf}))
        mkdir -p ${sid}; cd ${sid}
        of1=${sid}".yhaplo"
        tn1="temp1_"${of1}
        if [[ "$sex" == "M" ]]; then
            ibam1="${pt1}${iid}/BAM/${sid}/${sid}.L35.mapped.rmdup.q30.bam"
            if [[ "$type" != "ssLib" ]]; then
                ## Run pileupCaller for raw bam file
                samtools mpileup -B -R -q30 -Q30 -l ${posf} -f ${ref} ${ibam1} > ${tn1}_1
                cut -f 1 ${tn1}_1 | sed s/'X'/'23'/g | sed s/'Y'/'24'/g | paste - ${tn1}_1 | cut -f 1,3- > ${of1}.pileup
                rm ${tn1}_1
                pileupCaller --majorityCall --sampleNames ${rid} -f ${snpf} -e ${tn1}_1 < ${of1}.pileup
    
                ## Run pileupCaller for trimmed bam file
                ibam2="${pt1}${iid}/BAM/${sid}/${sid}.L35.mapped.rmdup.q30.m${mask}.bam"
                samtools mpileup -B -R -q30 -Q30 -l ${posf} -f ${ref} ${ibam2} > ${tn1}_2
                cut -f 1 ${tn1}_2 | sed s/'X'/'23'/g | sed s/'Y'/'24'/g | paste - ${tn1}_2 | cut -f 1,3- > ${of1}.masked.pileup
                rm ${tn1}_2
                pileupCaller --majorityCall --sampleNames ${rid} -f ${snpf} -e ${tn1}_2 < ${of1}.masked.pileup
    
                ## Choose calls from masked data for C/T and G/A SNPs
                paste ${snpf} ${tn1}_1.geno ${tn1}_2.geno | \
                awk '{if ($5$6 == "CT" || $5$6 == "TC" || $5$6 == "GA" || $5$6 == "AG") print $8; else print $7}' > ${of1}.geno
                        
            else
                ## Run pileupCaller for ibam1 only but with SingleStrandMode option
                samtools mpileup -B -R -q30 -Q30 -l ${posf} -f ${ref} ${ibam1} > ${tn1}_1
                cut -f 1 ${tn1}_1 | sed s/'X'/'23'/g | sed s/'Y'/'24'/g | paste - ${tn1}_1 | cut -f 1,3- > ${of1}.pileup
                rm ${tn1}_1
                pileupCaller --majorityCall --singleStrandMode --sampleNames ${rid} -f ${snpf} -e ${tn1}_1 < ${of1}.pileup
                mv ${tn1}_1.geno ${of1}.geno
            fi
            
            cp ${snpf} ${of1}.snp
            cp ${indf} ${of1}.ind
    
            ## Change EIGENSTRAT format data into yHaplo input format
            /home/choongwon_jeong/bin/EIGENSTRAT_to_yhaplo_170629.py ${of1}
    
            /usr/bin/python ${py} --ancStopThresh 10 -i ${of1}.genos.txt -o ./output -c -hp -ds -dsd -as -asd
            echo -e "SID\tYhap" > ${of1}.txt
            cat ./*/haplogroups.${of1}.txt | sort -k1,1 | awk '{OFS="\t"} {if ($2 == $3) print $1,$4"("$2")"; else print $1,$4"("$2","$3")"}' >> ${of1}.txt
            rm ${tn1}*
            cd ..
        else
            echo -e "SID\tYhap" > ${of1}.txt
            echo -e ${sid}"\tNA" >> ${of1}.txt 
            cd ..
        fi
    done


