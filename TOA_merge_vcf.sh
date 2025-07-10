#!/bin/bash
#SBATCH --time=2:0:0
#SBATCH --nodes=1
#SBATCH --mem=10GB
#SBATCH --ntasks-per-node=8
#SBATCH --account=OD-220670
##SBATCH --output=/dev/null

echo $PATH
start=$(date)
echo $start
SECONDS=0



module load bowtie/2.4.4
module load samtools/1.18.0
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
module load muscle
module load mafft/7.490

module load multiqc/1.15


cd /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/


echo ">>>>>>>>>>>>>>>>>>>>setup files<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"


# rm -r /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/Toothfish_RRBS_array2024_SE_mitogenome
# rm -r /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/Toothfish_RRBS_array2024_SE_lambda
# rm -r /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/Toothfish_RRBS_array2024_SE
# rm -r /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/Toothfish_RRBS_array2018_SE_mitogenome
# rm -r /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/Toothfish_RRBS_array2018_SE_lambda
# rm -r /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/Toothfish_RRBS_array2018_SE
# rm -r /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/MITOGENOME
# rm -r /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/LAMBDA
# rm -r /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/GENOME

#############################################################################################################################
# echo -e "\n\n>>>>>>>>>>>>>>>>>>>MITOGENOME<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

COMB_PATH=/scratch3/dev093/1.Toothfish/MITOGENOME
COV_PATH=/scratch3/dev093/1.Toothfish/MITOGENOME/Combined_methylation
SNP_PATH=/scratch3/dev093/1.Toothfish/MITOGENOME/Combined_SNP

rm -r ${COMB_PATH}
mkdir ${COMB_PATH}
mkdir ${COV_PATH}
mkdir ${SNP_PATH}

PREFIX=Toothfish_RRBS_array2018_SE_mitogenome
rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.cov.gz ${COV_PATH}
rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.bam ${COV_PATH}
rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*_SE_report.txt ${COV_PATH}
rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.M-bias.txt ${COV_PATH}
rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/masked/*.bam ${SNP_PATH}
rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/gatk/*vcf.gz ${SNP_PATH}

rm -r /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/revelio_temp
OUT_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}
RESULT_PATH=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/ 
rsync -av ${OUT_PATH} ${RESULT_PATH}



PREFIX=Toothfish_RRBS_array2024_SE_mitogenome
GENOME_PATH=/scratch3/dev093/1.Toothfish/Bismark_mtgenome_index/
GENOME_FILE=/scratch3/dev093/1.Toothfish/Bismark_mtgenome_index/D.mawsoni_mtDNA.fa
cd ${SNP_PATH}
rename -n "_R1_trimmed_sort_mask_bismark_bt2_E" "" *.bam

# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.cov.gz ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.cov.gz ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.bam ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*_SE_report.txt ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.M-bias.txt ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/masked/*.bam ${SNP_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/gatk/*vcf.gz ${SNP_PATH}

cd /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/
rm *_R1_trimmed_sort_bismark_bt2_SE_report.txt
rename "_R1_trimmed.fq_trimmed_100_bismark_bt2_SE_report.txt" "_R1_trimmed_sort_bismark_bt2_SE_report.txt" *.txt


rm -r /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/revelio_temp
OUT_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}
RESULT_PATH=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/ 
rsync -av ${OUT_PATH} ${RESULT_PATH}


echo ">>>>>>>>>>>>>>>>>>>>BISMARK SUMMARY<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
cd ${COV_PATH}
rm *_R1_trimmed_sort_bismark_bt2_SE_report.txt
rename "_R1_trimmed.fq_trimmed_100_bismark_bt2_SE_report.txt" "_R1_trimmed_sort_bismark_bt2_SE_report.txt" *.txt
bismark2summary -o Toothfish_RRBS_array_ALL_SE_mitogenome


echo ">>>>>>>>>>>>>>>>>>>>BAM CONSENSUS<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
module unload bismark/0.23.0
module load samtools/1.18.0

cd ${SNP_PATH}
rm *_new.bam
rm *.fa
for file in *.bam; do \
	echo $file ; \
	samtools index ${file}; \
	samtools view -bo ${file%%.*}_new.bam ${file} "LC138011" ;\
	samtools index ${file%%.*}_new.bam ; \
    samtools consensus --mode "simple" --no-use-MQ --min-depth 3 --call-fract 0.5 --het-fract 0.4 --ambig --format fasta ${file%%.*}_new.bam --output ${file%%.*}.fa ;\
	sed -e "s/LC138011/${file%%.*}/g" ${file%%.*}.fa > ${file%%.*}_2.fa
	rm ${file%%.*}.fa
	mv ${file%%.*}_2.fa ${file%%.*}.fa
done ## calling consencsus from double-masked bam files

rm seqs.fa
cat $GENOME_FILE `ls *.fa` > seqs.fa
set OMP_NUM_THREADS=16
export OMP_NUM_THREADS=16
muscle -super5 seqs.fa -output 0_Toothfish_mtgenome_aligned.fa
muscle -super5 0_TOA_input_Seq_geneious.fasta -output 0_Toothfish_mtgenome_aligned_MUSCLE.fasta
rsync -av /scratch3/dev093/1.Toothfish/Toothfish_RRBS_array2024_SE_mitogenome/Combined_SNP/0_Toothfish_mtgenome_aligned_MUSCLE.fasta /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/Toothfish_RRBS_array2024_SE_mitogenome/Combined_SNP/
	### G-INS-i (suitable for sequences of similar lengths; recommended for <200 sequences; iterative refinement method incorporating global pairwise alignment information)
cd /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/MITOGENOME/Combined_SNP/
mafft --globalpair --retree 3 --maxiterate 100 --thread 16 --reorder 0_TOA_input_Seq_geneious.fasta > 0_Toothfish_mtgenome_aligned_MAFFT_geneious.fasta
mafft --globalpair --retree 3 --maxiterate 100 --thread 16 --reorder seqs.fa > 0_Toothfish_mtgenome_aligned_MAFFT.fasta
rsync -av /scratch3/dev093/1.Toothfish/MITOGENOME/Combined_SNP/0_Toothfish_mtgenome_aligned_MAFFT.fasta /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/MITOGENOME/Combined_SNP/
rsync -av /scratch3/dev093/1.Toothfish/Toothfish_RRBS_array2024_SE_mitogenome/Combined_SNP/*.fa /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/Toothfish_RRBS_array2024_SE_mitogenome/Combined_SNP/


RESULT_PATH=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/ 
rsync -a ${COMB_PATH} ${RESULT_PATH}

module unload samtools/1.18.0
module load bismark/0.23.0






#############################################################################################################################
echo -e "\n\n>>>>>>>>>>>>>>>>>>>LAMBDA<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

COMB_PATH=/scratch3/dev093/1.Toothfish/LAMBDA
COV_PATH=/scratch3/dev093/1.Toothfish/LAMBDA/Combined_methylation
SNP_PATH=/scratch3/dev093/1.Toothfish/LAMBDA/Combined_SNP


rm -r ${COMB_PATH}
mkdir ${COMB_PATH}
mkdir ${COV_PATH}
mkdir ${SNP_PATH}


PREFIX=Toothfish_RRBS_array2018_SE_lambda
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.cov.gz ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.bam ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*_SE_report.txt ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.M-bias.txt ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/masked/*.bam ${SNP_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/gatk/*vcf.gz ${SNP_PATH}

rm -r /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/revelio_temp
OUT_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}
RESULT_PATH=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/ 
# rsync -av ${OUT_PATH} ${RESULT_PATH}


PREFIX=Toothfish_RRBS_array2024_SE_lambda
GENOME_PATH=/scratch3/dev093/1.Toothfish/Lambda_genome_index/
GENOME_FILE=/scratch3/dev093/1.Toothfish/Lambda_genome_index/Lambda_phage.fasta

# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.cov.gz ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.bam ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*_SE_report.txt ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.M-bias.txt ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/masked/*.bam ${SNP_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/gatk/*vcf.gz ${SNP_PATH}


cd /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/
rm *_R1_trimmed_sort_bismark_bt2_SE_report.txt
rename "_R1_trimmed.fq_trimmed_100_bismark_bt2_SE_report.txt" "_R1_trimmed_sort_bismark_bt2_SE_report.txt" *.txt




echo ">>>>>>>>>>>>>>>>>>>>BISMARK SUMMARY<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
cd ${COV_PATH}
rm *_R1_trimmed_sort_bismark_bt2_SE_report.txt
rename "_R1_trimmed.fq_trimmed_100_bismark_bt2_SE_report.txt" "_R1_trimmed_sort_bismark_bt2_SE_report.txt" *.txt
bismark2summary -o Toothfish_RRBS_array_ALL_SE_lambda


RESULT_PATH=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/
rsync -a ${COMB_PATH} ${RESULT_PATH}










#############################################################################################################################
echo -e "\n\n>>>>>>>>>>>>>>>>>>>GENOME<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

COMB_PATH=/scratch3/dev093/1.Toothfish/GENOME
COV_PATH=/scratch3/dev093/1.Toothfish/GENOME/Combined_methylation
# SNP_PATH=/scratch3/dev093/1.Toothfish/GENOME/Combined_SNP
# CGMAP_PATH=/scratch3/dev093/1.Toothfish/GENOME/Combined_SNP_CGMAP
SNP_PATH=/scratch3/dev093/1.Toothfish/Toothfish_RRBS_array2018_SE/gatk
CGMAP_PATH=/scratch3/dev093/1.Toothfish/Toothfish_RRBS_array2018_SE/CGMapTools


rm -r ${COMB_PATH}
mkdir ${COMB_PATH}
mkdir ${COV_PATH}
mkdir ${SNP_PATH}
mkdir ${CGMAP_PATH}


PREFIX=Toothfish_RRBS_array2018_SE
GENOME_PATH=/scratch3/dev093/1.Toothfish/Bismark_genome_index/
GENOME_FILE=/scratch3/dev093/1.Toothfish/Bismark_genome_index/D.mawsoni.fasta
FASTQC0=/scratch3/dev093/1.Toothfish/${PREFIX}/fastqc/initial/
FASTQC=/scratch3/dev093/1.Toothfish/${PREFIX}/fastqc/final/



# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.cov.gz ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.bam ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*_SE_report.txt ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.M-bias.txt ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/masked/*.bam ${SNP_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/gatk/*vcf.gz ${SNP_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/CGMapTools/*vcf.gz ${CGMAP_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/CGMapTools/*.asm ${CGMAP_PATH}

echo -e "\n\n>>>>>>>>>>>>>>>>>>>>FASTQC SUMMARY<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n"
cd ${FASTQC0}
multiqc --filename ${PREFIX}_fastqc_initial ./
cd ${FASTQC}
multiqc --filename ${PREFIX}_fastqc_final ./


rm -r /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/revelio_temp
OUT_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}
RESULT_PATH=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/ 
rsync -av ${OUT_PATH} ${RESULT_PATH}



PREFIX=Toothfish_RRBS_array2024_SE
GENOME_PATH=/scratch3/dev093/1.Toothfish/Bismark_genome_index/
GENOME_FILE=/scratch3/dev093/1.Toothfish/Bismark_genome_index/D.mawsoni.fasta
FASTQC0=/scratch3/dev093/1.Toothfish/${PREFIX}/fastqc/initial/
FASTQC=/scratch3/dev093/1.Toothfish/${PREFIX}/fastqc/final/
FASTQC100=/scratch3/dev093/1.Toothfish/${PREFIX}/fastqc/final_100/
FASTQTRIM=/scratch3/dev093/1.Toothfish/${PREFIX}/fastq_trim/


cd ${SNP_PATH}
rename -n '_R1_trimmed_sort_bismark_bt2_E' '' *

# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.cov.gz ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.bam ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*_SE_report.txt ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/*.M-bias.txt ${COV_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/masked/*.bam ${SNP_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/gatk/*vcf.gz ${SNP_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/CGMapTools/*vcf.gz ${CGMAP_PATH}
# rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/CGMapTools/*.asm ${CGMAP_PATH}


echo -e "\n\n>>>>>>>>>>>>>>>>>>>>FASTQC SUMMARY<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n"
cd ${FASTQC0}
multiqc --filename ${PREFIX}_fastqc_initial ./
cd ${FASTQC}
multiqc --filename ${PREFIX}_fastqc_final ./
cd ${FASTQTRIM}
multiqc --filename ${PREFIX}_fastqc_intermed ./
cd ${FASTQC100}
multiqc --filename ${PREFIX}_fastqc_final_100 ./

cd /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/
rm *_R1_trimmed_sort_bismark_bt2_SE_report.txt
rename "_R1_trimmed.fq_trimmed_100_bismark_bt2_SE_report.txt" "_R1_trimmed_sort_bismark_bt2_SE_report.txt" *.txt


rm -r /scratch3/dev093/1.Toothfish/${PREFIX}/Bismark/revelio_temp
OUT_PATH=/scratch3/dev093/1.Toothfish/${PREFIX}
RESULT_PATH=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/ 
rsync -av ${OUT_PATH} ${RESULT_PATH}
rsync -av /scratch3/dev093/1.Toothfish/${PREFIX}/fastqc/final_100 /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/${PREFIX}/fastqc/

echo ">>>>>>>>>>>>>>>>>>>>BISMARK SUMMARY<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"
cd ${COV_PATH}
rm *_R1_trimmed_sort_bismark_bt2_SE_report.txt
rename "_R1_trimmed.fq_trimmed_100_bismark_bt2_SE_report.txt" "_R1_trimmed_sort_bismark_bt2_SE_report.txt" *.txt
bismark2summary -o Toothfish_RRBS_array_ALL_SE_genome



RESULT_PATH=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/
rsync -a ${COMB_PATH} ${RESULT_PATH}



##########################################################################################################################################################
echo -e "\n\n>>>>>>>>>>>>>>>>>>>>COMBINE VCF<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n\n"

# OUTPUT=Toothfish_Nuclear_SNP
OUTPUT=Toothfish_Nuclear_SNP_2018

###################################################################################
echo -e "\n\n>>>>gVCF from REVELIO<<<\n\n"
cd ${SNP_PATH}

echo -e "\n\n>>>>>gatk - merge vcf<<<\n\n"

echo -e "\n\n>>gatk - merge vcf<<\n\n"
gatk CreateSequenceDictionary -R ${GENOME_FILE}


find ./ -type f -name "TOA*.vcf.gz" > input.list
INTERVAL=${GENOME_PATH}/interval.list #scaffold names

for file in *.vcf.gz; do gatk IndexFeatureFile -I ${file} -O ${file}.tbi; done ## vcf files need to be indexed

rm -r ./TOA_revelio_database

gatk --java-options "-Xmx990G -XX:ParallelGCThreads=16" GenomicsDBImport \
  --variant input.list \
  --genomicsdb-workspace-path TOA_revelio_database \
  --intervals ${INTERVAL} \
  --genomicsdb-shared-posixfs-optimizations true \
  --bypass-feature-reader \
  --merge-contigs-into-num-partitions 100 \
  --add-output-vcf-command-line false \
  --reader-threads 1 ### TAKES 26h on 1TB RAM 16 cpu ### for 2018: ### 2018 - TAKES 14h on 1TB RAM 16 cpu (MAX 30GB)

 # --genomicsdb-update-workspace-path TOA_revelio_database ###THIS MAKES LARGE LOG FILES for GENDB://, MAKE SURE TO SET #SBTACH --output=/dev/null
  
gatk --java-options "-Xmx1000G -XX:ParallelGCThreads=16" GenotypeGVCFs \
	-R ${GENOME_FILE} \
	-V gendb://TOA_revelio_database \
	--standard-min-confidence-threshold-for-calling 0.0 \
	--intervals ${INTERVAL} \
	--add-output-vcf-command-line false \
	--sample-ploidy 2 \
	-O ${SNP_PATH}/${OUTPUT}_TOA_revelio.vcf.gz ### TAKES 3.2 days on 1TB RAM 16 cpu for 420,000,000 variants?? ### 2018 - TAKES 1day:14h on 1TB RAM 16 cpu for 226,520,449 variants (MAX 10GB)


du -ah /scratch3/dev093/1.Toothfish/Toothfish_RRBS_array2024_SE/Combined_SNP/Toothfish_Nuclear_SNP.vcf.gz

gatk SelectVariants -R ${GENOME_FILE} -V ${SNP_PATH}/${OUTPUT}_TOA_revelio.vcf.gz --select-type-to-include SNP -restrict-alleles-to BIALLELIC \
	--add-output-vcf-command-line false \
	-O ${SNP_PATH}/${OUTPUT}_TOA_revelio_biallelic.vcf.gz
	
bcftools +setGT ${SNP_PATH}/${OUTPUT}_TOA_revelio_biallelic.vcf.gz -- -t q -n . -i 'FORMAT/DP<10 || FORMAT/GQ < 30' > ${SNP_PATH}/${OUTPUT}_TOA_revelio_biallelic_filtered.vcf.gz
gatk ValidateVariants -V ${SNP_PATH}/${OUTPUT}_TOA_revelio_biallelic_filtered.vcf.gz -R ${GENOME_FILE}
rsync ${SNP_PATH}/${OUTPUT}* /datasets/work/ncmi-toa-rrbs/work/1.Toothfish/Toothfish_RRBS_array2018_SE/gatk/

###################################################################################
echo -e "\n\n>>bcftools - merge vcf<<\n\n"
bcftools merge --threads 16 --output-type z SS*.vcf.gz > ${OUTPUT}_3.vcf.gz
bcftools index -t ${OUTPUT}_3.vcf.gz
gatk SelectVariants -R ${GENOME_FILE} -V ${OUTPUT}_3.vcf.gz --select-type-to-include SNP -restrict-alleles-to BIALLELIC \
	--add-output-vcf-command-line false \
	-O ${OUTPUT}_joint3.vcf.gz
source activate /datasets/work/ncmi-toa-rrbs/work/3.Epaulette_Shark/pacbio
vcf-merge -R 0/0 SS*.vcf.gz | bgzip -c > ${OUTPUT}_4.vcf.gz
conda deactivate


###################################################################################
echo -e "\n\n>>GLnexus - merge vcf<<\n\n"

	### REQUIRES 1766.51 GB RAM & 2h
cd /datasets/work/ncmi-toa-rrbs/work/3.Epaulette_Shark/GLnexus
chmod +x glnexus_cli
/datasets/work/ncmi-toa-rrbs/work/3.Epaulette_Shark/GLnexus/glnexus_cli \
	--dir ${SNP_PATH}/TOA_GLnexus_gatk.DB \
	--config gatk_unfiltered \
	--threads 16 \
    ${SNP_PATH}/TOA*.vcf.gz > ${SNP_PATH}/all_data.bcf
bcftools view ${SNP_PATH}/all_data.bcf | bgzip -@ 4 -c > ${SNP_PATH}/${OUTPUT}_TOA_GLnexus_gatk.vcf.gz

# 	--config gatk_unfiltered_rrbs \
# 	--config gatk_unfiltered \


	####### double check to keep bialleleic SNPs
gatk IndexFeatureFile -I ${SNP_PATH}/${OUTPUT}_TOA_GLnexus_gatk.vcf.gz -O ${SNP_PATH}/${OUTPUT}_TOA_GLnexus_gatk.vcf.gz.tbi

gatk SelectVariants -R ${GENOME_FILE} -V ${SNP_PATH}/${OUTPUT}_TOA_GLnexus_gatk.vcf.gz --select-type-to-include SNP -restrict-alleles-to BIALLELIC \
	--add-output-vcf-command-line false \
	-O ${SNP_PATH}/${OUTPUT}_TOA_GLnexus_gatk_biallelic.vcf.gz
	
bcftools +setGT ${SNP_PATH}/${OUTPUT}_TOA_GLnexus_gatk_biallelic.vcf.gz -- -t q -n . -i 'FORMAT/DP<10 || FORMAT/GQ < 30' -o ${SNP_PATH}/${OUTPUT}_TOA_GLnexus_gatk_biallelic_filtered.vcf.gz

gatk ValidateVariants -V ${SNP_PATH}/${OUTPUT}_TOA_GLnexus_gatk_biallelic_filtered.vcf.gz -R ${GENOME_FILE}









###################################################################################
echo -e "\n\n>>>>gVCF from CGMAPTOOLS<<<\n\n"
cd ${CGMAP_PATH}

echo -e "\n\n>>gatk - merge vcf<<\n\n"

### gatk CreateSequenceDictionary -R ${GENOME_FILE}

find ./ -type f -name "TOA*.vcf.gz" > input.list
INTERVAL=${GENOME_PATH}/interval.list #scaffold names

for file in *.vcf.gz; do gatk IndexFeatureFile -I ${file} -O ${file}.tbi; done ## vcf files need to be indexed

gatk --java-options "-Xmx990G -XX:ParallelGCThreads=16" GenomicsDBImport \
  --variant input.list \
  --genomicsdb-workspace-path TOA_cgmaptools_database \
  --intervals ${INTERVAL} \
  --genomicsdb-shared-posixfs-optimizations true \
  --bypass-feature-reader \
  --merge-contigs-into-num-partitions 100 \
  --add-output-vcf-command-line false \
  --reader-threads 1 ### TAKES 26h on 1TB RAM 16 cpu

  
gatk --java-options "-Xmx1000G -XX:ParallelGCThreads=16" GenotypeGVCFs \
	-R ${GENOME_FILE} \
	-V gendb://TOA_cgmaptools_database \
	--standard-min-confidence-threshold-for-calling 0.0 \
	--intervals ${INTERVAL} \
	--add-output-vcf-command-line false \
	--sample-ploidy 2 \
	-O ${CGMAP_PATH}/${OUTPUT}_TOA_cgmaptools.vcf.gz ### TAKES 3.2 days on 1TB RAM 16 cpu for 420,000,000 variants????


du -ah /scratch3/dev093/1.Toothfish/Toothfish_RRBS_array2024_SE/Combined_SNP/Toothfish_Nuclear_SNP.vcf.gz

gatk SelectVariants -R ${GENOME_FILE} -V ${CGMAP_PATH}/${OUTPUT}_TOA_cgmaptools.vcf.gz --select-type-to-include SNP -restrict-alleles-to BIALLELIC \
	--add-output-vcf-command-line false \
	-O ${CGMAP_PATH}/${OUTPUT}_TOA_cgmaptools_biallelic.vcf.gz
	
bcftools +setGT ${CGMAP_PATH}/${OUTPUT}_TOA_cgmaptools_biallelic.vcf.gz -- -t q -n . -i 'FORMAT/DP<10 || FORMAT/GQ < 30' -o ${CGMAP_PATH}/${OUTPUT}_TOA_cgmaptools_biallelic_filtered.vcf.gz

gatk ValidateVariants -V ${CGMAP_PATH}/${OUTPUT}_TOA_cgmaptools_biallelic_filtered.vcf.gz -R ${GENOME_FILE}



###################################################################################
echo -e "\n\n>>bcftools - merge vcf<<\n\n"
bcftools merge --threads 16 --output-type z SS*.vcf.gz > ${OUTPUT}_3.vcf.gz
bcftools index -t ${OUTPUT}_3.vcf.gz
gatk SelectVariants -R ${GENOME_FILE} -V ${OUTPUT}_3.vcf.gz --select-type-to-include SNP -restrict-alleles-to BIALLELIC \
	--add-output-vcf-command-line false \
	-O ${OUTPUT}_joint3.vcf.gz
source activate /datasets/work/ncmi-toa-rrbs/work/3.Epaulette_Shark/pacbio
vcf-merge -R 0/0 SS*.vcf.gz | bgzip -c > ${OUTPUT}_4.vcf.gz
conda deactivate


###################################################################################
echo -e "\n\n>>GLnexus - merge vcf<<\n\n"

	### REQUIRES 1766.51 GB RAM & 2h
cd /datasets/work/ncmi-toa-rrbs/work/3.Epaulette_Shark/GLnexus
chmod +x glnexus_cli
/datasets/work/ncmi-toa-rrbs/work/3.Epaulette_Shark/GLnexus/glnexus_cli \
	--dir ${CGMAP_PATH}/TOA_GLnexus_cgmaptools.DB \
	--config gatk_unfiltered \
	--threads 16 \
    ${CGMAP_PATH}/TOA*.vcf.gz > ${CGMAP_PATH}/all_data.bcf
bcftools view ${CGMAP_PATH}/all_data.bcf | bgzip -@ 4 -c > ${CGMAP_PATH}/${OUTPUT}_TOA_GLnexus_cgmaptools.vcf.gz

# 	--config gatk_unfiltered_rrbs \
# 	--config gatk_unfiltered \


	####### double check to keep bialleleic SNPs
gatk IndexFeatureFile -I ${CGMAP_PATH}/${OUTPUT}_TOA_GLnexus_cgmaptools.vcf.gz -O ${CGMAP_PATH}/${OUTPUT}_TOA_GLnexus_cgmaptools.vcf.gz.tbi

gatk SelectVariants -R ${GENOME_FILE} -V ${CGMAP_PATH}/${OUTPUT}_TOA_GLnexus_cgmaptools.vcf.gz --select-type-to-include SNP -restrict-alleles-to BIALLELIC \
	--add-output-vcf-command-line false \
	-O ${CGMAP_PATH}/${OUTPUT}_TOA_GLnexus_cgmaptools_biallelic.vcf.gz
	
bcftools +setGT ${CGMAP_PATH}/${OUTPUT}_TOA_GLnexus_cgmaptools_biallelic.vcf.gz -- -t q -n . -i 'FORMAT/DP<10 || FORMAT/GQ < 30' > ${CGMAP_PATH}/${OUTPUT}_TOA_GLnexus_cgmaptools_biallelic_filtered.vcf.gz

gatk ValidateVariants -V ${CGMAP_PATH}/${OUTPUT}_TOA_GLnexus_cgmaptools_biallelic_filtered.vcf.gz -R ${GENOME_FILE}








##########################################################################################################################################################
echo -e "\n\n>>>>>>>>>>>>>>>>>>>>COMBINE VCF FROM POOLED SAMPLES<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"

echo -e "\n>>>>gVCF from REVELIO<<<\n"

cd ${SNP_PATH}
find ./ -type f -name "*Pool*.vcf.gz" > input2.list
INTERVAL=${GENOME_PATH}/interval.list #scaffold names

gatk CreateSequenceDictionary -R ${GENOME_FILE}

### REQUIRES 106 GB RAM & 2h
gatk --java-options "-Xmx100G -XX:ParallelGCThreads=10" GenomicsDBImport \
  --variant input2.list \
  --genomicsdb-workspace-path TOA_pool_revelio_database \
  --intervals ${INTERVAL} \
  --merge-contigs-into-num-partitions 100 \
  --add-output-vcf-command-line false \
  --reader-threads 1
  
 # --genomicsdb-update-workspace-path SS_revelio_database ###THIS MAKES LARGE LOG FILES for GENDB://, MAKE SURE TO SET #SBTACH --output=/dev/null
  
gatk --java-options "-Xmx100G -XX:ParallelGCThreads=10" GenotypeGVCFs \
	-R ${GENOME_FILE} \
	-V gendb://TOA_pool_revelio_database \
	--standard-min-confidence-threshold-for-calling 0.0 \
	--intervals ${INTERVAL} \
	--add-output-vcf-command-line false \
	--sample-ploidy 2 \
	-O ${SNP_PATH}/${OUTPUT}_POOL_revelio.vcf.gz
	
	
# double check to keep bialleleic SNPs
gatk SelectVariants -R ${GENOME_FILE} -V ${SNP_PATH}/${OUTPUT}_POOL_revelio.vcf.gz --select-type-to-include SNP -restrict-alleles-to BIALLELIC \
	--add-output-vcf-command-line false \
	-O ${SNP_PATH}/${OUTPUT}_POOL_revelio_biallelic.vcf.gz
gatk ValidateVariants -V ${SNP_PATH}/${OUTPUT}_POOL_revelio_biallelic.vcf.gz -R ${GENOME_FILE}

mv ${SNP_PATH}/${OUTPUT}_POOL.vcf.gz ${SNP_PATH}/${OUTPUT}_POOL_revelio.vcf.gz
mv ${SNP_PATH}/${OUTPUT}_POOL_biallelic.vcf.gz ${SNP_PATH}/${OUTPUT}_POOL_revelio_biallelic.vcf.gz



###################################################################################
echo -e "\n\n>>GLnexus - merge vcf<<\n\n"

### REQUIRES 20 GB RAM & 1min
cd /datasets/work/ncmi-toa-rrbs/work/3.Epaulette_Shark/GLnexus
chmod +x glnexus_cli
/datasets/work/ncmi-toa-rrbs/work/3.Epaulette_Shark/GLnexus/glnexus_cli \
	--dir ${SNP_PATH}/TOA_pool_GLnexus_gatk.DB \
	--config gatk_unfiltered \
	--threads 16 \
    ${SNP_PATH}/*Pool*.vcf.gz > ${SNP_PATH}/all_data.bcf
bcftools view ${SNP_PATH}/all_data.bcf | bgzip -@ 4 -c > ${SNP_PATH}/${OUTPUT}_POOL.vcf.gz

# 	--config gatk_unfiltered_rrbs \
# 	--config gatk_unfiltered \


# double check to keep bialleleic SNPs

gatk IndexFeatureFile -I ${SNP_PATH}/${OUTPUT}_POOL_GLnexus_gatk.vcf.gz -O ${SNP_PATH}/${OUTPUT}_POOL_GLnexus_gatk.vcf.gz.tbi

gatk SelectVariants -R ${GENOME_FILE} -V ${SNP_PATH}/${OUTPUT}_POOL_GLnexus_gatk.vcf.gz --select-type-to-include SNP -restrict-alleles-to BIALLELIC \
	--add-output-vcf-command-line false \
	-O ${SNP_PATH}/${OUTPUT}_POOL_GLnexus_gatk_biallelic.vcf.gz
gatk ValidateVariants -V ${SNP_PATH}/${OUTPUT}_POOL_GLnexus_gatk_biallelic.vcf.gz -R ${GENOME_FILE}



###################################################################################
echo -e "\n>>>>POPOOLATION2<<<\n"
cd ${SNP_PATH}
POPOOLATION_PATH=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/popoolation2
POPOOLATION_RESULTS=/scratch3/dev093/1.Toothfish/${PREFIX}/Combined_SNP/popoolation
mkdir ${POPOOLATION_RESULTS}
rsync -av ${SNP_PATH}/*Pool*.bam ${POPOOLATION_RESULTS}/
cd ${POPOOLATION_RESULTS}

## REQUIRES 15G & 2h

samtools mpileup -f ${GENOME_FILE} *Pool*.bam > TOA_pool_all_indiv.mpileup

perl ${POPOOLATION_PATH}/mpileup2sync.pl --fastq-type sanger --min-qual 20 --input ${POPOOLATION_RESULTS}/TOA_pool_all_indiv.mpileup --output ${POPOOLATION_RESULTS}/TOA_pool_all_indiv.sync
java -ea -Xmx10g -jar ${POPOOLATION_PATH}/mpileup2sync.jar --input TOA_pool_all_indiv.mpileup --output TOA_pool_all_indiv.sync --fastq-type sanger --min-qual 20 --threads 8

Making mpileup files for each sample. Can be used for calculating nucleotide diversity for each sample.
for file in *Pool*.bam; do \
		echo $file ;\
		samtools coverage ${file} -o ${file}.coverage ;\
		samtools mpileup -f ${GENOME_FILE} ${file} > ${file}.mpileup ;\
		perl ${POPOOLATION_PATH}/Variance-sliding.pl --measure pi --input ${file}.mpileup --output ${file}.pi --min-count 2 --min-coverage 10 --max-coverage 200 --window-size 10000 --step-size 10000 --pool-size 80 --fastq-type sanger ;\
 done


perl ${POPOOLATION_PATH}/snp-frequency-diff.pl --input TOA_pool_all_indiv.sync --output-prefix TOA_pool_all_indiv --min-count 6 --min-coverage 50 --max-coverage 200
perl ${POPOOLATION_PATH}/fst-sliding.pl --input TOA_pool_all_indiv.sync  --output TOA_pool_all_indiv.fst --suppress-noninformative --min-count 6 --min-coverage 50 --max-coverage 200 --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size 80
perl ${POPOOLATION_PATH}/fst-sliding.pl --input TOA_pool_all_indiv.sync  --output TOA_pool_all_indiv_w5000.fst --min-count 6 --min-coverage 50 --max-coverage 200 --min-covered-fraction 1 --window-size 5000 --step-size 5000 --pool-size 80
perl ${POPOOLATION_PATH}/fisher-test.pl --input TOA_pool_all_indiv.sync  --output TOA_pool_all_indiv.fet --min-count 6 --min-coverage 50 --max-coverage 200 --suppress-noninformative
perl ${POPOOLATION_PATH}/cmh-test.pl --input p1_p2_p1_p2.sync --output p1_p2_p1_p2.cmh --min-count 12 --min-coverage 50 --max-coverage 200 --population 1-2,3-4
perl ${POPOOLATION_PATH}/create-genewise-sync.pl --input p1_p2.sync --gtf 2R_exons.gtf --output p1_p2_genes.sync
perl ${POPOOLATION_PATH}/fst-sliding.pl --min-count 6 --min-coverage 50 --max-coverage 200 --pool-size 80 --min-covered-fraction 0.0 --window-size 1000000 --step-size 1000000 --input p1_p2_genes.sync --output p1_p2_genewise.fst

echo ">>>>>>>>>>>>>>>>>>>>MOVE FILES<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

RESULT_PATH=/datasets/work/ncmi-toa-rrbs/work/1.Toothfish/
rsync -a ${COMB_PATH} ${RESULT_PATH}


duration=$SECONDS
echo "$(($duration / 3600)) hours, $((($duration / 60) % 60)) minutes and $(($duration % 60)) seconds elapsed."
