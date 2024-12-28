#!/bin/bash

# function
function CheckBam()
{
    bam=$1
    if [ -s $bam ]
    then
        samtools quickcheck $bam
         if [ $? -eq 0 ]
            then
                touch ${bam}.done
            else
                echo "${bam} is truncked, ${bam} is deleted"
                rm $bam 
            fi
        fi
}

# Set software path
Sdir="/work/chenfeilab/Avocado/Softwares/miniconda3/envs/Jingwen/bin"
GATK_PATH=/share/home/Coconut/Software/gatk-4.2.6.1/
BASEDIR=/chenfeilab/Avocado/MyPipline/tobiasrauschATAC/ATACseq/src

#reference genome
IndexHg19="/work/chenfeilab/Avocado/MyData/bowtieGenome/bowtie_hg19/bowtieHg19Index"
IndexHg38="/work/chenfeilab/Avocado/MyData/bowtieGenome/bowtie_GRCh38/bowtieHg38Index"
IndexGRCm39="/work/chenfeilab/Avocado/MyData/bowtieGenome/bowtie_GRCm39/bowtieGRCm39Index"
# GENOME_mycoplasma="/share/home/Coconut/Project/reference/genome/mycoplasma_index/mycoplasma"

RDNA_hg19="/work/chenfeilab/Avocado/MyData/bowtieGenome/bowtie_rDNA/Hg19rRNA"
RDNA_Hg38="/work/chenfeilab/Avocado/MyData/bowtieGenome/bowtie_rDNA/Hg38rRNA"
RDNA_GRCm39="/work/chenfeilab/Avocado/MyData/bowtieGenome/bowtie_rDNA/GRCm39rRNA"

blacklist_hg19="/work/chenfeilab/Avocado/MyData/hg19-blacklist.v2.bed"
blacklist_hg38="/work/chenfeilab/Avocado/MyData/hg38-blacklist.v2.bed"
blacklist_GRCm39="/work/chenfeilab/Avocado/MyData/mm10-blacklist.v2.Liftover.mm39.bed"
# download from https://github.com/Boyle-Lab/Blacklist/issues/25 mm10-blacklist.v2.Liftover.mm39.bed.txt

###Input
file_path=$1   # Run_Path
prefix=$2
SampleID=$3
Species=$4
# file_path=/work/chenfeilab/Avocado/P7_XXseq/P7_3ATAC/P7_3_1Preprocess/P7_3_1_1_RunAlign/ 
# prefix=P7_3_1_1_
# SampleID=ATAC_p300i_SMG_8Gy_Rep2_240720_185_RXX_GRCm39
# Species=GRCm39

echo $Species
if [[ $Species == "hg19" ]]; then 
  GENOME_EXP=${IndexHg19}
  GENOME_RDNA=${RDNA_hg19}
  blacklist=${blacklist_hg19}
  exp_info="hg19"
  ATYPE="hg19"
elif [[ $Species == "Hg38" ]]; then
  GENOME_EXP=${IndexHg38}
  GENOME_RDNA=${RDNA_Hg38}
  blacklist=${blacklist_hg38}
  exp_info="Hg38"
  ATYPE="hg38"
elif [[ $Species == "GRCm39" ]]; then
  GENOME_EXP=${IndexGRCm39}
  GENOME_RDNA=${RDNA_GRCm39}
  blacklist=${blacklist_GRCm39}
  exp_info="GRCm39"
  ATYPE="mm39"
fi
MT="chrM"
echo $GENOME_EXP

###file_path
# input_path=${all_path}/Raw_data/${datainfo}
# sampleinfo=${input_path}/sampleinfo_${datainfo}.txt
# file_path=${all_path}/Alignment/${datainfo}   # Run_Path

raw_dir=${file_path}/${prefix}00_rawdata
logs_dir=${file_path}/${prefix}logs
fastqc_dir=${file_path}/${prefix}01_fastqc
trimmedFastq_dir=${file_path}/${prefix}02.1_trimmeddata
trimmedFastq_log_dir=${file_path}/${prefix}logs/trimmeddata
trimmedFastq_fastqc_dir=${file_path}/${prefix}02.2_trimmeddata_fastqc

align_exp_dir=${file_path}/${prefix}03_bam
alignexp_log_dir=${file_path}/${prefix}logs/align_exp
align_myc_dir=${file_path}/${prefix}03_mycoplasmabam
alignmyc_log_dir=${file_path}/${prefix}logs/align_mycoplasma
exp_bam_rmdup=${file_path}/${prefix}03_bam_rmdup
rmdup_exp_log=${file_path}/${prefix}logs/rmdup_state_exp

bw_track_dir=${file_path}/${prefix}04_bigwig_track
peak_calling_dir=${file_path}/${prefix}04_peak

motif_path=${file_path}/${prefix}05_motifs

echo -e "*RNA-seq parameters:
          Run_path: ${file_path}
          prefix: ${prefix}
          SampleID: ${SampleID}
          Species: ${Species}"

# #================== rawdir ==================#
# ### rename if exists sampleinfo.txt
# #step 1.1 change file name#####
# echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Renaming files at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"
# 
# mkdir -p ${raw_dir}
# 
# cat $sampleinfo| while read id;
# do
#     arr=($id)
#     sample1=${arr[0]}
#     sample2=${arr[1]}
#     if [ ! -s ${raw_dir}/${sample2}_R1.fastq.gz ]
#     then
#         fq1=$(ls ${input_path}/*/Cleandata/*/*R1.f*q.gz|grep "$sample1")
#         fq2=$(ls ${input_path}/*/Cleandata/*/*R2.f*q.gz|grep "$sample1")
#         ln -s $fq1 ${raw_dir}/${sample2}_R1.fastq.gz
#         ln -s $fq2 ${raw_dir}/${sample2}_R2.fastq.gz
#     fi
# done


#step 1.2
####fastqc of raw data ####
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` fastqc of raw data at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir -p ${fastqc_dir}

in_path=${raw_dir}
out_path=${fastqc_dir}

nohup_number=0
for FILE in `ls ${in_path}/${SampleID}*.gz`
do
    sample=$(basename ${FILE/.fastq.gz/})
    if [ ! -s ${out_path}/"$(basename ${FILE/.fastq.gz/_fastqc.zip})" ]
    then
        echo "Generating file: ${out_path}/"$(basename ${FILE/.fastq.gz/_fastqc.zip})"..."
        fastqc $FILE -t 1 -o ${out_path}/ > ${out_path}/${sample}_fastqc.log 2>&1 &
        nohup_number=`echo $nohup_number+1 | bc`
    fi
    if [[ $nohup_number -eq 60 ]]
    then
        echo "waiting..."
        wait
        nohup_number=0
    fi
done

wait

#step 1.3
### merge reports of fastqc
if [ ! -s ${fastqc_dir}/${SampleID}/rawdata_multiqc_${SampleID}*.html ]
then
  multiqc ${fastqc_dir}/${SampleID}* -n rawdata_multiqc_${SampleID} -o ${fastqc_dir}/${SampleID} -q
else
  echo "Multiqc finished "${SampleID}
fi

#================== trim ==================#
#step 2.1
### Trimming adapters  (trim_galore)

echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Trimming adapters at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir -p ${trimmedFastq_dir}
mkdir -p ${trimmedFastq_log_dir}
mkdir -p ${trimmedFastq_fastqc_dir}

fq1=${in_path}/${SampleID}_R1.fastq.gz
fq2=${fq1/R1.fastq.gz/R2.fastq.gz}

if [ ! -s ${trimmedFastq_dir}/$SampleID"_R1_val_1.fq.gz" ]
then
    echo "Generating file: ${trimmedFastq_dir}/"$(basename ${fq1/.fastq.gz/_val_1.fq.gz})";"
    trim_galore -q 25 --phred33 --length 36 -e 0.1 --stringency 4 --paired \
    -o ${trimmedFastq_dir} $fq1 $fq2 \
    --fastqc_args "--outdir ${trimmedFastq_fastqc_dir}" \
    > ${trimmedFastq_log_dir}/"$(basename ${fq1/_R1.fastq.gz/_trimmed.log})" 2>&1 &
else
   echo "Trimming adapters done: ${SampleID} `printf "%0.s~" {1..16}`"
fi
wait

#================== align ==================#
#step 3.1
###Aligning to experimental genome#####
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Aligning to experimental genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir -p ${align_exp_dir}
mkdir -p ${alignexp_log_dir}

Out_bam=${align_exp_dir}/${SampleID}_${exp_info}.bam

CheckBam ${Out_bam}
if [ ! -s ${Out_bam} ]
then
  echo "aligning ${SampleID} to experimental genome"

  # -t:time -q:fastq(default)
  # -L:seedlength
  # no-mixed&no-discordant:for paired reads
  # Actually, the unpaired data has filitered here!
  ($Sdir/bowtie2 -p 58 -q -N 1 -L 25 -X 2000 \
      --no-mixed --no-discordant \
      -x ${GENOME_EXP} \
      -1 ${trimmedFastq_dir}/${SampleID}_R1_val_1.fq.gz \
      -2 ${trimmedFastq_dir}/${SampleID}_R2_val_2.fq.gz \
      2> ${alignexp_log_dir}/${SampleID}_align.log) |
      samtools sort -m 300M -@ 58 -O BAM -o ${Out_bam}
else
  echo "Aligning to experimental genome done: ${SampleID} `printf "%0.s~" {1..16}`"
fi

# #step 3.2
# ###Aligning to mycoplasmabam genome#####
# echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Aligning to mycoplasma genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"
# 
# mkdir -p ${align_myc_dir}
# mkdir -p ${alignmyc_log_dir}
# 
# for PAIR in $(ls ${trimmedFastq_dir} | sed 's/_R[1-2].*//' | uniq )
# do
#     CheckBam ${align_myc_dir}/${PAIR}_mycoplasma.bam
#     if [ ! -s ${align_myc_dir}/${PAIR}_mycoplasma.bam ]
#     then
#         echo "aligning ${PAIR} to mycoplasma genome"
# 
#         (bowtie2 -p 58 -q -N 1 -L 25 -X 2000 \
#             --no-mixed --no-discordant \
#             -x ${GENOME_mycoplasma} \
#             -1 ${trimmedFastq_dir}/${PAIR}_R1_val_1.fq.gz \
#             -2 ${trimmedFastq_dir}/${PAIR}_R2_val_2.fq.gz \
#             2> ${alignmyc_log_dir}/${PAIR}_mycoplasma_align.log) |
#             samtools sort -m 300M -@ 58 -O BAM -o ${align_myc_dir}/${PAIR}_mycoplasma.bam
#     fi
# done

#step 3.3
####### deduplicating for experimental genome
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Removing duplicates of experimental genome at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir -p ${exp_bam_rmdup}
mkdir -p ${rmdup_exp_log}

Out_rmdup_bam=${exp_bam_rmdup}/${SampleID}_rmdup.bam
if [ ! -s $Out_rmdup_bam ]
then
  samtools view -h ${Out_bam} | grep -v chrM |
    samtools sort -m 300m -O bam  -@ 58 -o - >${exp_bam_rmdup}/${SampleID}_rmChrM.bam
  #samtools flagstat ${exp_bam_rmdup}/${SampleID}_rmChrM.bam) >${rmdup_exp_log}/${SampleID}_rmChrM.stat)
  ## ref:https://www.biostars.org/p/170294/
  ## Calculate %mtDNA:
  samtools index -@ 58 ${Out_bam}  # Needed for "samtools idxstats"
  mtReads=$(samtools idxstats ${Out_bam} | grep 'chrM' | awk '{SUM += $3} END {print SUM}')
  totalReads=$(samtools idxstats ${Out_bam} | awk '{SUM += $3} END {print SUM}')
  echo "${SampleID} ==> mtDNA Content: $(bc <<< "scale=2;100*$mtReads/$totalReads")%" >> ${rmdup_exp_log}/${SampleID}_mt.log
  echo "${SampleID} ==> mtDNA Content: $(bc <<< "scale=2;100*$mtReads/$totalReads")%"

  picard MarkDuplicates REMOVE_DUPLICATES=True \
    INPUT=${exp_bam_rmdup}/${SampleID}_rmChrM.bam \
    OUTPUT=${Out_rmdup_bam} \
    METRICS_FILE=${exp_bam_rmdup}/${SampleID}_rmdup.metrics \
    2>${rmdup_exp_log}/${SampleID}_rmdup.log
else
  echo "Removing duplicates done: ${SampleID} `printf "%0.s~" {1..16}`"
fi

Out_shift_bam=${exp_bam_rmdup}/${SampleID}_shift.sorted.bam
if [ ! -s ${Out_shift_bam} ]
then
    samtools view -f 2 -q 10 -h ${Out_rmdup_bam} |
        samtools sort -m 300m -O bam  -@ 30 -o - >${exp_bam_rmdup}/${SampleID}_rmChrM_rmdup.q10.bam
    samtools flagstat ${exp_bam_rmdup}/${SampleID}_rmChrM_rmdup.q10.bam >${rmdup_exp_log}/${SampleID}_rmChrM_rmdup.q10.stat

    samtools index -@ 30 ${exp_bam_rmdup}/${SampleID}_rmChrM_rmdup.q10.bam
    alignmentSieve --numberOfProcessors 25 --ATACshift -b ${exp_bam_rmdup}/${SampleID}_rmChrM_rmdup.q10.bam -o ${exp_bam_rmdup}/${SampleID}_shift.bam 2>${rmdup_exp_log}/${SampleID}_shift.log
    samtools sort -m 300m -O bam  -@ 30 -o ${Out_shift_bam} ${exp_bam_rmdup}/${SampleID}_shift.bam
    samtools index -@ 30 ${Out_shift_bam}
    # bedtools bamtobed -i ${p_align}/${i}.last.bam  > ${p_align}/${i}.bed

    echo "rmdupQ10 flagstat has done;file has generated in ${rmdup_exp_log}/${SampleID}_rmChrM_rmdup.q10.stat"

    # ## size distribution
    # gatk=/share/home/Coconut/Software/gatk-4.2.0.0/gatk
    $GATK_PATH/gatk CollectInsertSizeMetrics \
        -H ${rmdup_exp_log}/${SampleID}_InsertSize.pdf \
        -I ${exp_bam_rmdup}/${SampleID}_shift.bam \
        -O ${rmdup_exp_log}/${SampleID}_InsertSize.txt \
        2>${rmdup_exp_log}/${SampleID}_InsertSize.log &
else
  echo "Removing duplicates done: ${SampleID} `printf "%0.s~" {1..16}`"
fi

#bedtools intersect -nonamecheck -v -a ${BamDir}/${Sample}_MD.bam -b ${Blacklist} > ${BamDir}/${Sample}_RMBL.bam
#samtools view -h -F 1804 -q 30



#step 3.3
### calculate alignment information ###
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Calculating alignment information at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

exp_path=${rmdup_exp_log}
align_path=${alignexp_log_dir}
alignmyc_path=${alignmyc_log_dir}

mkdir -p ${logs_dir}/scalefactor
scalefactor_file=${logs_dir}/scalefactor/scalefactor_${SampleID}.csv

echo -e "sample,ALLREADS,exp_READS,exp_mapping_RATIO(of total),mt_mapping_RATIO(of exp),exp_rmdup_RATIO(of exp_rmMT),exp_qc_READS,exp_qc_RATIO(of total)" > ${scalefactor_file}

ALLREADS=$(cat ${alignexp_log_dir}/${SampleID}_align.log|grep "were paired; of these:$"|cut -d "(" -f 1|awk '{print $1*2}')
exp_READS=$(grep "aligned" ${alignexp_log_dir}/${SampleID}_align.log | awk '{print $1}' | xargs | awk 'END{print ($2+$3)*2}')
exp_mapping_RATIO=$(cat ${alignexp_log_dir}/${SampleID}_align.log|grep "overall alignment rate"|cut -d "%" -f 1)%
# myc_mapping_RATIO=$(cat ${alignmyc_path}/${SampleID}_mycoplasma_align.log | grep "overall alignment rate"|cut -d "%" -f 1)
mt_mapping_RATIO=$(cat  ${rmdup_exp_log}/${SampleID}_mt.log | awk '{print $5}')
exp_rmdup_RATIO=$(sed -n 8p ${exp_bam_rmdup}/${SampleID}_rmdup.metrics | awk -F'\t' '{print $9}')
exp_qc_READS=$(cat ${rmdup_exp_log}/${SampleID}_rmChrM_rmdup.q10.stat|grep "total (QC-passed reads"|cut -d " " -f 1)
exp_qc_RATIO=$(echo "${exp_qc_READS}/${ALLREADS}"|bc -l)

#spike_qc_READS=$(cat ${spike_path}/${SampleID}_${spike_info}_rmdup.stat|grep "total (QC-passed reads"|cut -d " " -f 1)
#QC_reads=$(echo "${exp_qc_READS}+${spike_qc_READS}"|bc )
#spike_qc_RATIO_intotal=$(echo "${spike_qc_READS}/${ALLREADS}"|bc -l)
#spike_qc_RATIO_inqc=$(echo "${spike_qc_READS}/${QC_reads}"|bc -l)
#SCALEFACTOR=$(echo "1000000/${spike_qc_READS}"|bc -l)
echo -e $SampleID","$ALLREADS","$exp_READS","$exp_mapping_RATIO","$mt_mapping_RATIO","$exp_rmdup_RATIO","$exp_qc_READS","$exp_qc_RATIO >> ${scalefactor_file}
# "\t"$spike_qc_READS"\t"$spike_qc_RATIO_intotal"\t"$spike_qc_RATIO_inqc"\t"$SCALEFACTOR



#step 4.1
### Making CPM-normalized bigWig files with full-length reads adjusted by spike-in / CPM equal to RPM###
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Making CPM-normalized bigWig files with full-length reads at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir ${bw_track_dir}

if [ ! -s ${bw_track_dir}/${SampleID}_CPM.bw ]
then
    bamCoverage -p 60 \
    --skipNonCoveredRegions \
    --normalizeUsing CPM \
    --binSize 1 \
    --blackListFileName $blacklist \
    -b ${Out_shift_bam} \
    -o ${bw_track_dir}/${SampleID}_CPM.bw &
fi

#step 5.1
### Peak Calling ###
echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Calling peak at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir -p ${peak_calling_dir}

# --cutoff-analysis:decide proper cutoff,more than 30time to cost --bdg
# keep chr only: awk -v FS="\t" -v OFS="\t" 'length($1)<=5' | sort -k1,1 -k2,2n > ${Macs2Dir}/WT_SortPeaks.bed
Out_peak=${peak_calling_dir}/${SampleID}_peaks.final.narrowPeak
if [ ! -s ${Out_peak} ]
then
  macs2 callpeak -f BAMPE \
    --min-length 100 --keep-dup all --nolambda --bdg -g hs -n ${SampleID} \
    -t ${Out_shift_bam} --outdir ${peak_calling_dir} 2>${peak_calling_dir}/${SampleID}_macs2.log
    
  bedtools intersect -a ${peak_calling_dir}/${SampleID}_peaks.narrowPeak -b $blacklist -f 0.25 -v > ${Out_peak} &
fi

echo -e "\033[31m`printf "%0.s-" {1..90}`\n `echo $SampleID` Motif discovery at $(date +%Y"-"%m"-"%d" "%H":"%M":"%S)\n`printf "%0.s-" {1..90}`\033[0m"

mkdir -p ${motif_path}

# Motif discovery
Out_file_motif=${motif_path}/${SampleID}/homerResults.html
if [ ! -s $Out_file_motif ]
then
  ${BASEDIR}/motif.sh ${ATYPE} ${Out_peak} ${motif_path}/${SampleID} &
fi
wait

echo -e "\033[36m`printf "%0.s-" {1..60}`\033[0m"
echo -e "\033[36mJob finished for ${SampleID}.\033[0m"
echo -e "\033[36m DESTROYING ALL! \033[0m"
echo -e "\033[36m    Bazinga! \033[0m"
echo -e "\033[36m`printf "%0.s-" {1..60}`\033[0m"
