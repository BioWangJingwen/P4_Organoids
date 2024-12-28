#!/bin/bash

### call TF peak use narrow peak
# macs3 callpeak -t ChIP.bam -c Control.bam -f BAM -g hs -n test -B -q 0.01

### call Histone peak use broad peak
# macs3 callpeak -t ChIP.bam -c Control.bam --broad -g hs --broad-cutoff 0.1

# function
function CheckBam()
{
    bam=$1
    sam=$2
    if [ -s $bam ]
    then
        samtools quickcheck $bam
         if [ $? -eq 0 ]
            then
                touch ${bam}.done
            else
                echo "${bam} is truncked, ${bam} is deleted"
                rm $bam 
                rm $sam
            fi
        fi
}

# Set resources path
IndexHg19="/work/chenfeilab/Avocado/MyData/bowtieGenome/bowtie_hg19/bowtieHg19Index"
IndexHg38="/work/chenfeilab/Avocado/MyData/bowtieGenome/bowtie_GRCh38/bowtieHg38Index"
IndexGRCm39="/work/chenfeilab/Avocado/MyData/bowtieGenome/bowtie_GRCm39/bowtieGRCm39Index"
# blacklist="/share/home/Coconut/Project/Genome/hg19-blacklist.v2.bed"
blacklist="/work/chenfeilab/Avocado/MyData/hg19-blacklist.v2.bed"

# Set software path
Sdir="/work/chenfeilab/Avocado/Softwares/miniconda3/envs/Jingwen/bin"
seacr="/work/chenfeilab/MyPipline/SEACR-master/SEACR_1.3.sh"
gatk=/share/home/Coconut/Software/gatk-4.2.0.0/gatk
BASEDIR=/chenfeilab/Avocado/MyPipline/tobiasrauschATAC/ATACseq/src
# callpeak="callpeaks.py"

# setup tools
THREADS=30

# Setup paths
file_path=$1   # Run_Path
prefix=$2
SampleID=$3
Species=$4
# file_path=/chenfeilab/Avocado/P4_Organoids/P4_9CutNTag/P4_9_1Preprocess/P4_9_1_1_RunAlign/ 
# prefix=P4_9_1_1_
# SampleID=CUTTag_EBF_P003T_Rep1_240503_016_QZB
# Species=hg19
MT="chrM"
echo $Species
if [[ $Species == "hg19" ]]; then 
  Index=$IndexHg19
  ATYPE="hg19"
elif [[ $Species == "Hg38" ]]; then
  Index=${IndexHg38}
  ATYPE="hg38"
elif [[ $Species == "GRCm39" ]]; then
  Index=$IndexGRCm39
  ATYPE="mm39"
fi
echo $Index

raw_path=${file_path}/${prefix}00_rawdata
logs_path=${file_path}/${prefix}logs
fastqc_path=${file_path}/${prefix}01_fastqc

trimmedFastq_path=${file_path}/${prefix}02_trimmeddata
trimmedFastq_log_path=${file_path}/${prefix}logs/trimmeddata

align_exp_path=${file_path}/${prefix}03_bam
alignexp_log_path=${file_path}/${prefix}logs/align_exp
exp_bam_rmdup=${file_path}/${prefix}03_bam_rmdup
rmdup_exp_log=${file_path}/${prefix}logs/rmdup_state_exp

bw_path=${file_path}/${prefix}04_bigwig
bw_fulllength_dir=${bw_path}/${prefix}04_fulllength
bw_singlebase_dir=${bw_path}/${prefix}04_singlebase

peak_path=${file_path}/${prefix}05_peaks

motif_path=${file_path}/${prefix}06_motifs

#step 1.2
#### fastqc of raw data ####
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` fastqc of raw data at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir -p ${fastqc_path}

in_path=${raw_path}
out_path=${fastqc_path}

nohup_number=0
for FILE in `ls ${in_path}/${SampleID}*.gz`
do
    # FILE=${in_path}/${SampleID}*.gz
    sample=$(basename ${FILE/.fastq.gz/})
    if [ ! -s ${out_path}/"$(basename ${FILE/.fastq.gz/_fastqc.zip})" ]
    then
        echo "Generating file: ${out_path}/"$(basename ${FILE/.fastq.gz/_fastqc.zip})"..."
        fastqc $FILE -t 4 -o ${out_path}/ > ${out_path}/${sample}_fastqc.log 2>&1 &
    else
        echo "Fastqc finished ${SampleID} `printf "%0.s~" {1..16}`"
    fi
    
    nohup_number=`echo $nohup_number+4 | bc`
    if [[ $nohup_number -eq ${THREADS} ]]
    then
        echo "waiting..."
        wait
        nohup_number=0
    fi
done

#step 1.3
### merge reports of fastqc
if [ ! -s ${fastqc_path}/${SampleID}/rawdata_multiqc_${SampleID}.html ]
then
  multiqc ${fastqc_path}/${SampleID}* -n rawdata_multiqc_${SampleID} -o ${fastqc_path}/${SampleID} -q
else
  echo "Multiqc finished "${SampleID}
fi

#step 2.1
### Trimming adapters  (trim_galore)
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Trimming adapters at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir -p ${trimmedFastq_path}
mkdir -p ${trimmedFastq_log_path}

nohup_number=0
fq1=${raw_path}/${SampleID}_R1.fastq.gz
fq2=${fq1/R1.fastq.gz/R2.fastq.gz}
if [ ! -s ${trimmedFastq_path}/$SampleID"_R1_val_1.fq.gz" ]
then
   echo "Generating file: ${trimmedFastq_path}/"$(basename ${fq1/.fastq.gz/_val_1.fq.gz})";"
   trim_galore -q 25 --phred33 --length 36 -e 0.1 --fastqc \
      --stringency 4 --paired -o ${trimmedFastq_path} $fq1 $fq2 \
      > ${trimmedFastq_log_path}/"$(basename ${fq1/_R1.fastq.gz/_trimmed.log})" 2>&1 &
else
   echo "Trimming adapters done: ${SampleID} `printf "%0.s~" {1..16}`"
fi
wait
# PID 180723

#step 3.1
###Aligning to experimental genome#####
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Aligning to experimental genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir -p ${align_exp_path}
mkdir -p ${alignexp_log_path}

Out_file_sam="${align_exp_path}/$SampleID.sam"
Out_file_bam="${align_exp_path}/$SampleID.bam"
Out_file_sort_bam="${align_exp_path}/$SampleID.sorted.bam" 

CheckBam ${Out_file_bam} ${Out_file_sam}

if [ ! -s $Out_file_sam ]
then
  echo "Alignment $SampleID"
  F1=${trimmedFastq_path}/$SampleID"_R1_val_1.fq.gz"
  F2=${trimmedFastq_path}/$SampleID"_R2_val_2.fq.gz"
  $Sdir/bowtie2 -p 6 \
    --end-to-end --very-sensitive \
    --no-mixed --no-discordant \
    --phred33 \
    -I 10 \
    -X 700 \
    -x $Index \
    -1 $F1 \
    -2 $F2 \
    -S $Out_file_sam &> ${alignexp_log_path}/$SampleID.txt
fi

if [ ! -s $Out_file_bam ]
then  
  ## Filtering Bam
  echo "Filter and keep mapped pairs ${SampleID}"
  samtools view -bS -F 0x04 $Out_file_sam > $Out_file_bam
fi

if [ ! -s $Out_file_sort_bam ]
then 
  echo "Sort samples by name ${SampleID}"
  samtools sort $Out_file_bam > $Out_file_sort_bam
else
  echo "Alignment done ${SampleID}"
fi

#step 3.2
####### deduplicating for experimental genome
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Removing duplicates of experimental genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir -p ${exp_bam_rmdup}
mkdir -p ${rmdup_exp_log}
Out_file_rmdup=${exp_bam_rmdup}/${SampleID}_rmDup.bam

# samtools view -h $Out_file_sort_bam | grep -v chrM |\
#     samtools sort -m 300m -O bam  -@ 58 -o - >$Out_file_rmChrM
# #samtools flagstat ${exp_bam_rmdup}/$(basename ${SampleID}_rmChrM.bam) >${rmdup_exp_log}/$(basename ${SampleID}_rmChrM.stat)
# ## ref:https://www.biostars.org/p/170294/
# ## Calculate %mtDNA:
# samtools index -@ 58 $Out_file_sort_bam  # Needed for "samtools idxstats"
# mtReads=$(samtools idxstats $Out_file_sort_bam | grep 'chrM' | awk '{SUM += $3} END {print SUM}')
# totalReads=$(samtools idxstats $Out_file_sort_bam | awk '{SUM += $3} END {print SUM}')
# echo "${sample} ==> mtDNA Content: $(bc <<< "scale=2;100*$mtReads/$totalReads")%" >> ${rmdup_exp_log}/${SampleID}_mt.log

if [ ! -s $Out_file_rmdup ]
then
  samtools view -h ${Out_file_sort_bam} | grep -v $MT |
            samtools sort -m 300m -O bam  -@ 58 -o - >${exp_bam_rmdup}/${SampleID}_rmChrM.bam
        #samtools flagstat ${exp_bam_rmdup}/$(basename ${FILE%.bam}_rmChrM.bam) >${rmdup_exp_log}/$(basename ${FILE%.bam}_rmChrM.stat)
        ## ref:https://www.biostars.org/p/170294/
        ## Calculate %mtDNA:
  samtools index ${Out_file_sort_bam}  # Needed for "samtools idxstats"
  mtReads=$(samtools idxstats ${Out_file_sort_bam} | grep $MT | awk '{SUM += $3} END {print SUM}')
  totalReads=$(samtools idxstats ${Out_file_sort_bam} | awk '{SUM += $3} END {print SUM}')
  echo "${SampleID} ==> mtDNA Content: $(bc <<< "scale=2;100*$mtReads/$totalReads")%" >> ${rmdup_exp_log}/${SampleID}_mt.log
  
  echo "${SampleID} ==> mtDNA Content: $(bc <<< "scale=2;100*$mtReads/$totalReads")%" 
  echo "MarkDuplicates $SampleID"
  picard MarkDuplicates REMOVE_DUPLICATES=True \
    INPUT=${exp_bam_rmdup}/${SampleID}_rmChrM.bam \
    OUTPUT=$Out_file_rmdup \
    METRICS_FILE=${exp_bam_rmdup}/${SampleID}_rmDup.metrics \
    2>${exp_bam_rmdup}/${SampleID}_rmDup.log
else
  echo "Remove Duplicates done ${SampleID}"
fi

Out_file_rmdup_sort=${exp_bam_rmdup}/${SampleID}_rmDup_sorted.bam
if [ ! -s ${Out_file_rmdup_sort} ]
then
  samtools sort ${Out_file_rmdup} > ${Out_file_rmdup_sort}
  samtools index ${Out_file_rmdup_sort}  # Needed for "samtools idxstats"
  samtools flagstat ${Out_file_rmdup_sort} > ${exp_bam_rmdup}/${SampleID}_rmDup_sorted.stat
fi

# Out_file_shift=${exp_bam_rmdup}/${SampleID}_rmdup.q10.bam
# if [ ! -s $Out_file_shift ]
# then
#   samtools view -f 2 -q 10 -h $Out_file_rmdup |
#     samtools sort -m 300m -O bam  -@ 30 -o - > ${exp_bam_rmdup}/${SampleID}_rmdup.q10.bam
#     
#   
#   samtools flagstat ${exp_bam_rmdup}/${SampleID}_rmdup.q10.bam > ${exp_bam_rmdup}/${SampleID}_rmdup.q10.stat
#   
#   samtools index -@ 30 ${exp_bam_rmdup}/${SampleID}_rmdup.q10.bam
# #   alignmentSieve --numberOfProcessors 25 \
# #     --ATACshift -b ${exp_bam_rmdup}/${SampleID}_rmdup.q10.bam \
# #     -o ${exp_bam_rmdup}/${SampleID}_shift.bam 2>${exp_bam_rmdup}/${SampleID}_shift.log
# #   samtools sort -m 300m -O bam  -@ 30 -o $Out_file_shift ${exp_bam_rmdup}/${SampleID}_shift.bam
# #   samtools index -@ 30 $Out_file_shift
# # # bedtools bamtobed -i ${p_align}/${i}.last.bam  > ${p_align}/${i}.bed
# 
#   echo "rmdupQ10 flagstat has done;file has generated in ${exp_bam_rmdup}/${SampleID}_rmdup.q10.stat"
#   
#   # ## size distribution
#   # $gatk CollectInsertSizeMetrics -H ${rmdup_exp_log}/${SampleID}_InsertSize.pdf \
#   #   -I ${exp_bam_rmdup}/${SampleID}_shift.bam -O ${rmdup_exp_log}/${SampleID}_InsertSize.txt \
#   #   2>${rmdup_exp_log}/${SampleID}_InsertSize.log
# else
#   echo "rmdupQ10 flagstat has done: ${SampleID}"
# fi 

mkdir -p ${align_exp_path}/fragmentLen
Out_file_fLen=${align_exp_path}/fragmentLen/$SampleID.txt          
if [ ! -s $Out_file_fLen ]
then          
  ## Extract fragment length
  echo extracted fragment length $SampleID
  samtools view -F 0x04 $Out_file_sam |\
    awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' |\
    sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > $Out_file_fLen &
else
  echo "fragmentLen has done: ${SampleID}"
fi

#step 3.3
### calculate alignment information ###
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Calculating alignment information at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

exp_path=${rmdup_exp_log}
align_path=${alignexp_log_dir}
alignmyc_path=${alignmyc_log_dir}

mkdir -p ${logs_path}/scalefactor
scalefactor_file=${logs_path}/scalefactor/scalefactor_${SampleID}.csv

if [ ! -s ${scalefactor_file} ]
then  
  # $Out_file_rmdup

  echo -e "sample,ALLREADS,exp_READS,exp_mapping_RATIO(of total),mt_mapping_RATIO(of exp),exp_rmdup_RATIO(of exp_rmMT),exp_qc_READS,exp_qc_RATIO(of total)" > ${scalefactor_file}

  ALLREADS=$(cat ${alignexp_log_path}/${SampleID}.txt | grep "were paired; of these:$" | cut -d "(" -f 1 | awk '{print $1*2}')
  exp_READS=$(grep "aligned" ${alignexp_log_path}/${SampleID}.txt | awk '{print $1}' | xargs | awk 'END{print ($2+$3)*2}')
  exp_mapping_RATIO=$(cat ${alignexp_log_path}/${SampleID}.txt|grep "overall alignment rate"|cut -d "%" -f 1)
  mt_mapping_RATIO=$(cat ${rmdup_exp_log}/${SampleID}_mt.log | awk '{print $5}')
  exp_rmdup_RATIO=$(sed -n 8p ${exp_bam_rmdup}/${SampleID}_rmDup.metrics | awk -F'\t' '{print $9}')
  exp_qc_READS=$(cat ${exp_bam_rmdup}/${SampleID}_rmDup_sorted.stat|grep "total (QC-passed reads"|cut -d " " -f 1)
  exp_qc_RATIO=$(echo "${exp_qc_READS}/${ALLREADS}"|bc -l)

  #spike_qc_READS=$(cat ${spike_path}/${sample}_${spike_info}_rmdup.stat|grep "total (QC-passed reads"|cut -d " " -f 1)
  #QC_reads=$(echo "${exp_qc_READS}+${spike_qc_READS}"|bc )
  #spike_qc_RATIO_intotal=$(echo "${spike_qc_READS}/${ALLREADS}"|bc -l)
  #spike_qc_RATIO_inqc=$(echo "${spike_qc_READS}/${QC_reads}"|bc -l)
  #SCALEFACTOR=$(echo "1000000/${spike_qc_READS}"|bc -l)
  echo -e ${SampleID}","$ALLREADS","$exp_READS","${exp_mapping_RATIO}%","$mt_mapping_RATIO","$exp_rmdup_RATIO","$exp_qc_READS","$exp_qc_RATIO >> ${scalefactor_file}
  # "\t"$spike_qc_READS"\t"$spike_qc_RATIO_intotal"\t"$spike_qc_RATIO_inqc"\t"$SCALEFACTOR
else
  echo "Calculating alignment information has done: ${SampleID}"
fi

#step 4.1
### Making CPM-normalized bigWig files with full-length reads adjusted by spike-in / CPM equal to RPM###
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Making CPM-normalized bigWig files with full-length reads at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir ${bw_path}

Out_file_bw=${bw_path}/${SampleID}_CPM.bw
if [ ! -s $Out_file_bw ]
  then
  bamCoverage -p 60 \
    --skipNonCoveredRegions \
    --normalizeUsing CPM \
    --binSize 1 \
    --blackListFileName $blacklist \
    -b ${Out_file_rmdup_sort} \
    -o $Out_file_bw &
else
  echo "CPM-normalized bigWig has done: ${SampleID}"
fi

#step 5.1
### Peak Calling ###
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Calling peak at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir -p ${peak_path}

# --cutoff-analysis:decide proper cutoff,more than 30time to cost --bdg
# keep chr only: awk -v FS="\t" -v OFS="\t" 'length($1)<=5' | sort -k1,1 -k2,2n > ${Macs2Dir}/WT_SortPeaks.bed

Out_file_nPeak=${peak_path}/${SampleID}_peaks.final.narrowPeak
if [ ! -s $Out_file_nPeak ]
then
  macs2 callpeak -f BAMPE \
    --min-length 100 --keep-dup all --nolambda --bdg -g hs -n ${SampleID} \
    -t ${Out_file_rmdup_sort} --outdir ${peak_path} 2>${peak_path}/${SampleID}_macs2.log
  # filter peaks against ENCODE blacklist -- abnormal region that enriched in almost all the second generation high-throughput sequencing
  bedtools intersect -a ${peak_path}/${SampleID}_peaks.narrowPeak -b $blacklist -f 0.25 -v >$Out_file_nPeak
else
  echo "Calling peak has done: ${SampleID}"
fi

echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Motif discovery at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir -p ${motif_path}

# Motif discovery
Out_file_motif=${motif_path}/${SampleID}/homerResults.html
if [ ! -s $Out_file_motif ]
then
  ${BASEDIR}/motif.sh ${ATYPE} $Out_file_nPeak ${motif_path}/${SampleID} &
fi
wait

echo -e "\033[36m`printf "%0.s-" {1..60}`\033[0m"
echo -e "\033[36mJob finished for ${SampleID}.\033[0m"
echo -e "\033[36m DESTROYING ALL! \033[0m"
echo -e "\033[36m    Bazinga! \033[0m"
echo -e "\033[36m`printf "%0.s-" {1..60}`\033[0m"

# #step 4 
# ###Creation of bedgraph files#####
# echo -e "\n`printf "%0.s-" {1..90}`\nFILE Creation of bedgraph files at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`"
# 
# Out_file_bed=${align_exp_path}/${SampleID}.fragments.bedgraph
# if [ ! -s $Out_file_bed ]
# then
#   echo "Convert to bed ${SampleID}"
#   bedtools bamtobed -i $Out_file_sort_bam -bedpe > ${align_exp_path}/${SampleID}.bed
# 
#   echo "Keep read pairs on same chromosomes and fragment length less than 1000bp ${SampleID}"
#   awk '$1==$4 && $6-$2 < 1000 {print $0}' ${align_exp_path}/${SampleID}.bed > ${align_exp_path}/${SampleID}.clean.bed
# 
#   echo "Extract fragment related columns only ${SampleID}"
#   cut -f 1,2,6 ${align_exp_path}/${SampleID}.clean.bed | sort -k1,1 -k2,2n -k3,3n  > ${align_exp_path}/${SampleID}.fragments.bed
#   
#   echo "Make bedgraph ${SampleID}"
#   bedtools genomecov -bg -i ${align_exp_path}/${SampleID}.fragments.bed -g $chromSize > $Out_file_bed
# fi
# 
# #step 6 
# ###call peaks (use IgG control as reference) #####
# echo -e "\n`printf "%0.s-" {1..90}`\ncall peaks (use IgG control as reference) at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`"
# 
# mkdir -p ${peak_path}SEACR-stringent
# mkdir -p ${peak_path}/SEACR-relaxed
# 
# # for i in "${!SAMP[@]}"; do
# 
# Out_file_peak=${peak_path}/SEACR-relaxed/${SampleID}_seacr_top0.01.peaks
# if [ ! -s $Out_file_bed ]
# then
# 
# # echo Call peaks against IgG stringent ${SampleID}
# # bash $seacr $Out_file_bed ./mapped/bedgraph/A7_S7.fragments.bedgraph norm stringent ${peak_path}/SEACR-stringent/${SampleID}_seacr_IgG.peaks
# # 
# # echo Call peaks against IgG relaxed ${SampleID}
# # bash $seacr $Out_file_bed ./mapped/bedgraph/A7_S7.fragments.bedgraph norm relaxed ${peak_path}/SEACR-relaxed/${SampleID}_seacr_IgG.peaks
# 
# echo Call top 1% peaks stringent ${SampleID}
# bash $seacr $Out_file_bed 0.01 norm stringent ${peak_path}/SEACR-stringent/${SampleID}_seacr_top0.01.peaks
# 
# echo Call top 1% peaks relaxed ${SampleID}
# bash $seacr $Out_file_bed 0.01 norm relaxed $Out_file_peak
# 
# fi
# 
# # done
# 
# # mkdir -p ${peak_path}
# # 
# # Out_file_peak=${peak_path}/${SampleID}.peaks
# # 
# # $callpeak -b $Out_file_sort_bam -o ${Out_file_peak}
# 
# # =================================================================
# # ## PART 7 make bigwig files to check on igv
# # # ===============================================================
# 
# 
# ## Make merged bigwig files for WT and I-SceI
# 
# mkdir -p ${bw_path}
# 
# echo Index bam ${SAMP[i]}                                                   
# samtools index ./mapped/bam/${SAMP[i]}.sorted.bam
# 
# 
# 
# ## need to be in python virtual environment to use deeptools
# # source /home/sum879/python3.7.4/bin/activate
# 
# # echo Generate bigWig ${SAMP[i]}
# # bamCoverage -b ./mapped/bam/${SAMP[i]}.sorted.bam -o ./mapped/bigwig/${SAMP[i]}_raw.bw
# # echo Generate bigwig A3_S1
# # bamCoverage -b ./mapped/bam/A3_S1.sorted.bam -o ./mapped/bigwig/A3_S1_raw.bw
# # echo Generate bigwig I-SceI merged
# # bamCoverage -b ./mapped/bam/I-SceI.merged.sorted.bam -o ./mapped/bigwig/I-SceI_merge_raw.bw
# 	
# ## after done with python environment
# # deactivate
# 
# # done













