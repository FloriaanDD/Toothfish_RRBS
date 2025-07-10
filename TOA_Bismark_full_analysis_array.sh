#!/bin/bash
#SBATCH --time=24:0:0
#SBATCH --mem=10G
#SBATCH --cpus-per-task=16
#SBATCH --account=OD-220670
#SBATCH --array=1-144%70
##SBATCH --array=1-168%70

echo "task_id =" $SLURM_ARRAY_TASK_ID
echo "number of tasks/cpus" $SLURM_NTASKS

echo $PATH
start=$(date)
echo $start
SECONDS=0

cd /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/


module load bowtie/2.4.4
module load samtools/1.12
module load bcftools/1.15.1
module load bismark/0.23.0
module load fastqc/0.12.1
module load trimgalore/0.6.10
module load gcc/10.3.0
module load perl/5.32.1
module load zlib/1.2.11
module load R/4.0.5
module load gatk/4.2.6.1
module load bcftools
module load vcftools
module load miniconda3
module load bbmap/39.06
module load cgmaptools/0.1.3


echo ">>>>>>>>>>>>>>>>>>>>setup files<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
SAMPLETABLE=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/Toothfish_files.tsv # Path to a fastq table where the 1st column is the prefix of the raw fastq files. The 4th column is the sample ID.
# SAMPLETABLE=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/Toothfish_files_2024.tsv # Path to a fastq table where the 1st column array nr. The 2nd column is the sample ID.


echo ">>>>>SET scratch3 DIR<<<<<<<<"
SCRATCH_PATH=/scratch3/dev093/1.Toothfish/
cd ${SCRATCH_PATH}

# rsync -a ${GENOME_PATH} ${IN_PATH} ${SCRATCH_PATH}

# PREFIX=Toothfish_RRBS_array2024_SE_mitogenome
# GENOME_PATH=/scratch3/dev093/1.Toothfish/Bismark_mtgenome_index/
# GENOME_FILE=/scratch3/dev093/1.Toothfish/Bismark_mtgenome_index/D.mawsoni_mtDNA.fa
# FQ_TRIM=/scratch3/dev093/1.Toothfish/Toothfish_RRBS_array2024_SE/fastq_trim
# BAM_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark


# PREFIX=Toothfish_RRBS_array2024_SE_lambda
# GENOME_PATH=/scratch3/dev093/1.Toothfish/Lambda_genome_index/
# GENOME_FILE=/scratch3/dev093/1.Toothfish/Lambda_genome_index/Lambda_phage.fasta
# FQ_TRIM=/scratch3/dev093/1.Toothfish/Toothfish_RRBS_array2024_SE/fastq_trim
# BAM_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark


# PREFIX=Toothfish_RRBS_array2024_SE
# GENOME_PATH=/scratch3/dev093/1.Toothfish/Bismark_genome_index/
# GENOME_FILE=/scratch3/dev093/1.Toothfish/Bismark_genome_index/D.mawsoni.fasta
# IN_PATH=/scratch3/dev093/1.Toothfish/Merged_files2/
# FASTQC=/scratch3/dev093/1.Toothfish/${PREFIX}/fastqc
# FQ_TRIM=/scratch3/dev093/1.Toothfish/${PREFIX}/fastq_trim
# BAM_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark


PREFIX=Toothfish_RRBS_array2018_SE
GENOME_PATH=/scratch3/dev093/1.Toothfish/Bismark_genome_index/
GENOME_FILE=/scratch3/dev093/1.Toothfish/Bismark_genome_index/D.mawsoni.fasta
IN_PATH0=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/1.Raw_fastq/Merged_files2
IN_PATH=/scratch3/dev093/1.Toothfish/Merged_files2/
FASTQC=/scratch3/dev093/1.Toothfish/${PREFIX}/fastqc
FQ_TRIM=/scratch3/dev093/1.Toothfish/${PREFIX}/fastq_trim
BAM_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark

# PREFIX=Toothfish_RRBS_array2018_SE_mitogenome
# GENOME_PATH=/scratch3/dev093/1.Toothfish/Bismark_mtgenome_index/
# GENOME_FILE=/scratch3/dev093/1.Toothfish/Bismark_mtgenome_index/D.mawsoni_mtDNA.fa
# FQ_TRIM=/scratch3/dev093/1.Toothfish/Toothfish_RRBS_array2018_SE/fastq_trim
# BAM_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark

# PREFIX=Toothfish_RRBS_array2018_SE_lambda
# GENOME_PATH=/scratch3/dev093/1.Toothfish/Lambda_genome_index/
# GENOME_FILE=/scratch3/dev093/1.Toothfish/Lambda_genome_index/D.Lambda_phage.fa
# FQ_TRIM=/scratch3/dev093/1.Toothfish/Toothfish_RRBS_array2018_SE/fastq_trim
# BAM_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark


echo ">>Set DIR<<"
rm -r ./${PREFIX}
mkdir ${PREFIX}
mkdir ${PREFIX}/fastqc
mkdir ${PREFIX}/fastqc/initial
mkdir ${PREFIX}/fastqc/final
mkdir ${PREFIX}/fastqc/final_100
mkdir ${PREFIX}/fastq_trim
mkdir ${PREFIX}/Bismark
mkdir ${PREFIX}/Bismark/BAM
mkdir ${PREFIX}/CGMapTools

cd ./${PREFIX}



echo -e "\n\n\n\n>>>>>>>>>>>>>>>>>>>>START ANALYSIS<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n\n\n"

Bismark - genome index ## ONLY USED MAX 6GB
echo -e "\n\n>>>>>>>>>>>>>>>>>>>>Genome prep<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"
bismark_genome_preparation --bowtie2 --parallel 8 --verbose ${GENOME_PATH}
samtools faidx ${GENOME_FILE}




SAMPLE_ID=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $7}' $SAMPLETABLE)
echo $SAMPLE_ID

# SAMPLE_ID=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $SAMPLETABLE)
# echo $SAMPLE_ID

FQ_FILE=${IN_PATH}${SAMPLE_ID}_R1.fastq.gz
BAM_FILE=${BAM_PATH}/${SAMPLE_ID}_R1_trimmed.fq_trimmed_bismark_bt2.bam
# BAM_FILE=${BAM_PATH}/${SAMPLE_ID}_R1_trimmed.fq_trimmed_100_bismark_bt2.bam
BAM_SORT=${BAM_PATH}/${SAMPLE_ID}_R1_trimmed_sort_bismark_bt2.bam
BSMK_REPORT=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/${SAMPLE_ID}_R1_trimmed_sort_bismark_bt2.CpG_report.txt.gz
ATCGmap=/scratch3/dev093/1.Toothfish/${PREFIX}/CGMapTools/${SAMPLE_ID}.ATCGmap.gz
CGmap=/scratch3/dev093/1.Toothfish/${PREFIX}/CGMapTools/${SAMPLE_ID}.CGmap.gz
VCF=/scratch3/dev093/1.Toothfish/${PREFIX}/CGMapTools/${SAMPLE_ID}.vcf
SNV=/scratch3/dev093/1.Toothfish/${PREFIX}/CGMapTools/${SAMPLE_ID}.snv
ASM=/scratch3/dev093/1.Toothfish/${PREFIX}/CGMapTools/${SAMPLE_ID}.asm
NUGEN_PATH=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/NuMetRRBS/





#################################################################################################################################################################################################################
## USES APPROX 60GB - 15h per sample on 16 cpus
# echo $FQ_FILE

echo -e "\n\n>>>>>>>>>>>>>>>>>>>>Running fastqc<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"
fastqc -t 16 ${FQ_FILE} --outdir ${FASTQC}/initial

echo -e "\n\n>>>>>>>>>>>>>>>>>>>>Running trim_galore<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"
trim_galore -j 16 -a AGATCGGAAGAGC --path_to_cutadapt /apps/cutadapt/4.3/bin/cutadapt --length 35 --nextseq 30 -o ${FQ_TRIM} ${FQ_FILE}

echo -e "\n\n>>>>>>>>>>>>>>>>>>>>Remove reads without MspI restriction site<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n"
		##USES APPROX 5GB - 2h per sample
		## If a restriction enzyme (like MSP or Taq1) recognition site (YGG or similar) is found at the beginning of the sequence, it trims off the bases from the beginning up to the YGG position.
		## it trims 5 bases from the 3' end for single-end reads, or 6 bases for paired-end reads. This is intended to remove low-quality bases or adapter sequences from the 3' end.
cd ${FQ_TRIM}
python ${NUGEN_PATH}trimRRBSdiversityAdaptCustomers_chatgpt.py3 -1 ${FQ_TRIM}/${SAMPLE_ID}_R1_trimmed.fq.gz


echo -e "\n\n>>>>>>>>>>>>>>>>>>>>Running fastqc<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n"
fastqc -t 16 ${FQ_TRIM}/${SAMPLE_ID}_R1_trimmed.fq_trimmed.fq.gz --outdir ${FASTQC}/final

echo -e "\n\n>>>>>>>>>>>>>>>>>>>>cutting to 100 bp <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"
	###forcetrimright=94 bececause NUGEN also cuts 5bp from 5` end for the 2018 100bp fastq files, 
bbduk.sh in=${FQ_TRIM}/${SAMPLE_ID}_R1_trimmed.fq_trimmed.fq.gz out=${FQ_TRIM}/${SAMPLE_ID}_R1_trimmed.fq_trimmed_100.fq.gz forcetrimright=94
echo -e "\n\n>>>>>>>>>>>>>>>>>>>>Running fastqc<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n"
fastqc -t 16 ${FQ_TRIM}/${SAMPLE_ID}_R1_trimmed.fq_trimmed_100.fq.gz --outdir ${FASTQC}/final_100





#################################################################################################################################################################################################################
echo -e "\n\n\n\n>>>>>>>>>>>>>>>>>>>>BISMARK - Mapping reads to genome<<<<<<<<<<<<<<<<<<<<<\n\n"
cd ${BAM_PATH}

bismark --bowtie2 --multicore 16 -p 2 --score_min L,0,-0.4 -N 1 --rg_tag --rg_id ${SAMPLE_ID} --rg_sample ${SAMPLE_ID} -o ${BAM_PATH} --genome ${GENOME_PATH} ${FQ_TRIM}/${SAMPLE_ID}_R1_trimmed.fq_trimmed.fq.gz
# bismark --bowtie2 --multicore 16 -p 2 --score_min L,0,-0.4 -N 1 --rg_tag --rg_id ${SAMPLE_ID} --rg_sample ${SAMPLE_ID} -o ${BAM_PATH} --genome ${GENOME_PATH} ${FQ_TRIM}/${SAMPLE_ID}_R1_trimmed.fq_trimmed_100.fq.gz

echo -e "\n\n>>>>>>>>>>>>>>>>>>>>SAMTOOLS - Sorting mapped reads<<<<<<<<<<<<<<<<<<<<<\n"
samtools sort ${BAM_FILE} -o ${BAM_SORT}
samtools index ${BAM_SORT}

mv ${BAM_FILE} ./BAM/
# BSMK_REPORT1=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/${SAMPLE_ID}_R1_trimmed.fq_trimmed_100_bismark_bt2_SE_report.txt
BSMK_REPORT1=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/${SAMPLE_ID}_R1_trimmed.fq_trimmed_bismark_bt2_SE_report.txt
BSMK_REPORT2=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/${SAMPLE_ID}_R1_trimmed_sort_bismark_bt2_SE_report.txt
mv ${BSMK_REPORT1} ${BSMK_REPORT2}

echo -e "\n\n>>>>>>>>>>>>>>>>>>>>BISMARK - Extracting methylation data<<<<<<<<<<<<<<<<<\n\n"
bismark_methylation_extractor --multicore 12 --single-end --gzip --cytosine_report --comprehensive --cutoff 4 --ignore 3 --bedGraph --genome_folder ${GENOME_PATH} ${BAM_SORT}
cd ../




#################################################################################################################################################################################################################
echo -e "\n\n>>>>>>>>>>>>>>>>>>>>SNP DATA<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"

# echo -e "\n\n>>>>>>>>>>>>>>>>>>>>double-masking BAM<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"

#### ONLY REQUIRES 20GB RAM

cd ${BAM_PATH}
source activate /datasets/work/ncmi-toa-rrbs/work/2.School_shark/revelio_env

REVELIO=/datasets/work/ncmi-toa-rrbs/work/2.School_shark/revelio/revelio.py
REV_TMP=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/revelio_temp/
BAM_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/
BAM_PATH2=/scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/masked/
GATK_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}/gatk/

rm -r ${REV_TMP}
mkdir ${REV_TMP}
mkdir ${BAM_PATH2}
mkdir ${GATK_PATH}

gatk CreateSequenceDictionary -R ${GENOME_FILE}


cd ${BAM_PATH2}

echo $SAMPLE_ID
GVCF=${GATK_PATH}${SAMPLE_ID}.g.vcf.gz
BAM_MASK=${BAM_PATH2}${SAMPLE_ID}_R1_trimmed_sort_bismark_bt2
samtools calmd --threads 12 -b ${BAM_SORT} ${GENOME_FILE} 1> ${BAM_MASK}_A.bam 2> /dev/null
samtools sort -@ 12 ${BAM_MASK}_A.bam -o ${BAM_MASK}_A.bam
samtools index -@ 12 ${BAM_MASK}_A.bam
$REVELIO -T 12 --quality --temp ${REV_TMP} ${BAM_MASK}_A.bam ${BAM_MASK}_B.bam
	
samtools sort -@ 12 ${BAM_MASK}_B.bam -o ${BAM_MASK}_D.bam
samtools index -@ 12 ${BAM_MASK}_D.bam

bam trimBam ${BAM_MASK}_D.bam ${BAM_MASK}_E.bam -L 3 -R 0 --clip
samtools sort -@ 12 ${BAM_MASK}_E.bam -o ${BAM_MASK}_E.bam


mv ${BAM_MASK}_E.bam ${BAM_PATH2}${SAMPLE_ID}.bam
samtools index -@ 12 ${BAM_PATH2}${SAMPLE_ID}.bam

rm ${BAM_MASK}_A.bam ${BAM_MASK}_A.bam.bai ${BAM_MASK}_B.bam ${BAM_MASK}_B.bam.bai ${BAM_MASK}_D.bam ${BAM_MASK}_D.bam.bai


gatk --java-options "-Xmx100g -XX:ParallelGCThreads=16" HaplotypeCaller \
	 -R ${GENOME_FILE} \
	 -I ${BAM_PATH2}${SAMPLE_ID}.bam \
	 -O ${GVCF} \
	 --sample-name $SAMPLE_ID \
	 --min-dangling-branch-length 1 --min-pruning 1 \
	 --add-output-vcf-command-line false \
	 --emit-ref-confidence GVCF \
	 --native-pair-hmm-threads 3 \
	 --heterozygosity 0.01 \
	 --sample-ploidy 2


gatk IndexFeatureFile -I ${GVCF} -O ${GATK_PATH}${SAMPLE_ID}.tbi ## vcf files need to be indexed

conda deactivate



#################################################################################################################################################################################################################
echo -e "\n\n>>>>>>>>>>>>>>>>>>>>CGmap - Extracting SNP data<<<<<<<<<<<<<<<<<<<<<<<<<<\n"
echo ">>Convert<<"
cd ./CGMapTools
cgmaptools convert bam2cgmap -b ${BAM_SORT} -g ${GENOME_FILE} -o ${SAMPLE_ID}
echo ">>SNV<<"
echo ">BAYES<"
zcat ${ATCGmap} | cgmaptools snv -m bayes --bayes-dynamicP -v ${VCF} -o ${SNV}
cgmaptools asm -r ${GENOME_FILE} -b ${BAM_SORT} -l ${VCF} -o ${ASM}
echo ">>Stats<<"
cgmaptools mbin -i ${CGmap} -B 500 -c 4 -f png -t ${SAMPLE_ID} -p ${SAMPLE_ID} > ${SAMPLE_ID}_mbin.data
cgmaptools mstat -i ${CGmap} -c 4 -f png -p ${SAMPLE_ID} -t ${SAMPLE_ID} > ${SAMPLE_ID}_mstat.data
cgmaptools oac bin -i ${ATCGmap} -B 1000 -f png -p ${SAMPLE_ID} -t ${SAMPLE_ID} > ${SAMPLE_ID}_oac_bin.data
cgmaptools oac stat -i ${ATCGmap} -p ${SAMPLE_ID} -f png > ${SAMPLE_ID}_oac_stat.data
cgmaptools mec bin -i ${CGmap} -B 1000 -f png -p ${SAMPLE_ID} -t ${SAMPLE_ID} > ${SAMPLE_ID}_mec_bin.data
cgmaptools mec stat -i ${CGmap} -p ${SAMPLE_ID} -f png > ${SAMPLE_ID}_mec_stat.data



echo -e "\n\n>>RE-FORMAT VCF<<\n\n"

bcftools reheader --fai ${GENOME_FILE}.fai ${SAMPLE_ID}.vcf > ${SAMPLE_ID}_temp.vcf
sed -e "5 i ##reference=file://${GENOME_FILE}" ${SAMPLE_ID}_temp.vcf > ${SAMPLE_ID}_temp2.vcf ## add reference tag to header
awk ' BEGIN { FS=OFS="\t" } {gsub(":",";",$8)}1' ${SAMPLE_ID}_temp2.vcf > ${SAMPLE_ID}_temp.vcf ## change colon to semicolin for INFO
awk ' BEGIN { FS=OFS="\t" } {gsub("0","99",$6)}1' ${SAMPLE_ID}_temp.vcf > ${SAMPLE_ID}_temp2.vcf ## adjust Quality scores
awk ' BEGIN { FS=OFS="\t" } {gsub(":0:",":99:",$10)}1' ${SAMPLE_ID}_temp2.vcf > ${SAMPLE_ID}_temp.vcf
	
sed -e "s/NA00001/${SAMPLE_ID}/g" ${SAMPLE_ID}_temp.vcf > ${SAMPLE_ID}.vcf ## add sample name
cat ${SAMPLE_ID}.vcf | grep -P '^#' > ${SAMPLE_ID}_temp.vcf
cat ${SAMPLE_ID}.vcf | awk -v OFS='\t' '{{ if(($4!="N" && ($5=="A"||$5=="C"||$5=="G"||$5=="T"))) {{ print }} }}' >> ${SAMPLE_ID}_temp.vcf	## remove ambigous ALT alleles

	
gatk SelectVariants -R ${GENOME_FILE} -V ${SAMPLE_ID}_temp.vcf \
	--select-type-to-include SNP -restrict-alleles-to BIALLELIC \
	--add-output-vcf-command-line false \
	-O ${SAMPLE_ID}_NEW.vcf ## double check to keep bialleleic SNPs

bcftools view ${SAMPLE_ID}_NEW.vcf -O z -o ./${SAMPLE_ID}.vcf.gz
bcftools index -t ./${SAMPLE_ID}.vcf.gz
vcf-validator ./${SAMPLE_ID}.vcf.gz
rm ${SAMPLE_ID}_temp.vcf
rm ${SAMPLE_ID}_NEW.vcf


echo -e "\n\n>>>>>>>>>>>>>>>>>>>>Move files<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"

OUT_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}
RESULT_PATH=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/
rsync -a ${OUT_PATH} ${RESULT_PATH}


duration=$SECONDS
echo "$(($duration / 3600)) hours, $((($duration / 60) % 60)) minutes and $(($duration % 60)) seconds elapsed."
