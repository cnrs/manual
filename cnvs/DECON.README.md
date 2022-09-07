# 利用DECoN从WES数据中检测CNV  

# R
```
if (!requireNamespace("BiocManager", quietly = TRUE))  
    install.packages("BiocManager")  
BiocManager::install(version = "3.15")  

BiocManager::install("ExomeDepth")  
BiocManager::install("optparse")  
BiocManager::install("getopt")  
```

# conda
```
conda create -n decon
conda activate decon
conda install -c bioconda gatk4 samtools bwa fastp trimmomatic  
conda install r-biocmanager r-exomedepth r-ggplot2 r-getopt r-optparse r-reshape bioconductor-genomicranges  
```


# Download UCSC hg19 genome
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz  
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz  
perl ext_ucsc_human_canonical_fa.pl hg19.fa > hg19.fa.1  
mv hg19.fa.1 hg19.fa  
```


# 下载变异注释文件
```
wget -c -r -nd -np -k -L -p ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19

或者：

wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_omni2.5.hg19.sites.vcf.idx.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.indels.hg19.sites.vcf.idx.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/1000G_phase1.snps.high_confidence.hg19.sites.vcf.idx.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/CEUTrio.HiSeq.WGS.b37.bestPractices.hg19.vcf.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/CEUTrio.HiSeq.WGS.b37.bestPractices.hg19.vcf.idx.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.excluding_sites_after_129.vcf.idx.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.idx.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3_hg19_pop_stratified_af.vcf.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3_hg19_pop_stratified_af.vcf.gz.tbi
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/hapmap_3.3.hg19.sites.vcf.idx.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.sites.vcf.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.sites.vcf.idx.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.vcf.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/NA12878.HiSeq.WGS.bwa.cleaned.raw.subset.hg19.vcf.idx.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/NA12878.knowledgebase.snapshot.20131119.hg19.vcf.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/NA12878.knowledgebase.snapshot.20131119.hg19.vcf.idx.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.dict.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.fai.gz
wget -c -t 0 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/ucsc.hg19.fasta.gz
```


# 1. Reads mapping
```
bwa index -a bwtsw hg19.fa  

ll *_1.fastq | sed -e 's/_1.fastq//g' | awk '{print "trimmomatic PE -threads 18 -phred33 " $9 "_1.fastq " $9 "_2.fastq " $9 "_1.fq " $9 "_1.unpaired.fq " $9 "_2.fq "  $9 "_2.unpaired.fq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40 AVGQUAL:20"}'  

# ID:   样品的ID号
# PU:   测序仪器
# PL:   测序平台
# LB:   文库名
# SM:   样品名 

#例如：@RG ID:ZX1_ID SM:ZX1 LB:PE400 PU:Illumina PL:Miseq

samples=`ll *_R1.fq.gz | sed -e 's/_R1.fq.gz//g' | awk '{print $9}'`
for SAMPLE in ${samples}
do
  echo ${SAMPLE}
  #ID=LANE1
  ID=${SAMPLE}
  PU="Illumina"
  PL="NovaSeq6000"
  LB="PE400"
  SM=${SAMPLE}
  
  nohup bwa mem -t 18 -M -R "@RG\tID:${ID}\tPU:${PU}\tPL:${PL}\tLB:${LB}\tSM:${SM}" -o ${SAMPLE}.sam /usr/local/db/ucsc/human/BWA_hg19_genome/hg19.fa ${SAMPLE}_R1.fq.gz ${SAMPLE}_R2.fq.gz > ${SAMPLE}.nohup 2>&1 &
  sleep 1s
done

# 注意要用全路径，应为索引文件都在这个目录下，不要用链接：
# ln -s /usr/local/db/ucsc/human/BWA_hg19_genome/hg19.fa

#ll *_1.fq | sed -e 's/_1.fq//g'| awk '{print "bwa mem -t 18 -M -R "@RG\tID:${ID}\tPU:${PU}\tPL:${PL}\tLB:${LB}\tSM:${SM}\" -o \${SAMPLE}\.sam /usr/local/db/ucsc/human/BWA_hg19_genome/hg19.fa " $9 "_1.fq " $9 "_2.fq & sleep 1s"}'
ll *.sam | sed -e 's/.sam//g' | awk '{print "samtools view -bS -o " $9 ".tmp " $9 ".sam & sleep 1s"}'
ll *.tmp | sed -e 's/.tmp//g' | awk '{print "samtools sort -o " $9 ".bam " $9 ".tmp & sleep 1s"}'
ll *.bam | awk '{print "samtools index " $9 " & sleep 1s"}'


标记PCR重复序列并建立索引
samples=`ll *.bam | sed -e 's/.bam//g' | awk '{print $9}'`
for SAMPLE in ${samples}
do
  #echo ${SAMPLE}
  #nohup gatk MarkDuplicates -I ${SAMPLE}.bam -O ${SAMPLE}.M.bam -M ${SAMPLE}.metrics > ${SAMPLE}.nohup 2>&1 &
  #nohup picard MarkDuplicates -I ${SAMPLE}.bam -O ${SAMPLE}.M.bam -M ${SAMPLE}.metrics > ${SAMPLE}.nohup 2>&1 &
  echo "nohup picard MarkDuplicates -I ${SAMPLE}.bam -O ${SAMPLE}.M.bam -M ${SAMPLE}.metrics > ${SAMPLE}.nohup 2>&1 &"
  #sleep 1s
done

ll *.bam | awk '{print "samtools index " $9 " & sleep 1s"}'


${GATK} MarkDuplicates -I ${WORK_DIR}/mapping/${SAMPLE}.sorted.bam -O ${WORK_DIR}/mapping/${SAMPLE}.sorted.MarkDuplicates.bam -M ${WORK_DIR}/mapping/${SAMPLE}.sorted.bam.metrics
${SAMTOOLS} index ${WORK_DIR}/mapping/${SAMPLE}.sorted.MarkDuplicates.bam

重新校正碱基质量值（BQSR）
#下载变异注释文件 
wget  -c -r -nd -np -k -L -p ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19

${GATK} BaseRecalibrator -R ${GENOME} -I ${WORK_DIR}/mapping/${SAMPLE}.sorted.MarkDuplicates.bam -O ${WORK_DIR}/mapping/${SAMPLE}.recal_data.table --known-sites ${known1000G_indels} --known-sites ${GoldStandard_indels} --known-sites ${dbSNP}
${GATK} ApplyBQSR -R ${GENOME} -I ${WORK_DIR}/mapping/${SAMPLE}.sorted.MarkDuplicates.bam -bqsr ${WORK_DIR}/mapping/${SAMPLE}.recal_data.table -O ${WORK_DIR}/mapping/${SAMPLE}.bam

rm ${WORK_DIR}/mapping/${SAMPLE}*.sam
rm ${WORK_DIR}/mapping/${SAMPLE}*.tmp
rm ${WORK_DIR}/mapping/${SAMPLE}*.sorted.bam
rm ${WORK_DIR}/mapping/${SAMPLE}*.sorted.MarkDuplicates.bam
```

# 2. ReadInBams
```
Rscript ReadInBams.R --bam bamlist.txt --bed skin.bed --fsa hg19.fa --out CNV_COUNTS
```

bamlist.txt
```
/home/wangk/lab/baiyun/FASTQ/SAMPLE01.bam
/home/wangk/lab/baiyun/FASTQ/SAMPLE02.bam
/home/wangk/lab/baiyun/FASTQ/SAMPLE03.bam
/home/wangk/lab/baiyun/FASTQ/SAMPLE04.bam
/home/wangk/lab/baiyun/FASTQ/SAMPLE05.bam
/home/wangk/lab/baiyun/FASTQ/SAMPLE06.bam
/home/wangk/lab/baiyun/FASTQ/SAMPLE07.bam
/home/wangk/lab/baiyun/FASTQ/SAMPLE08.bam
/home/wangk/lab/baiyun/FASTQ/SAMPLE09.bam
/home/wangk/lab/baiyun/FASTQ/SAMPLE10.bam
/home/wangk/lab/baiyun/FASTQ/SAMPLE11.bam
/home/wangk/lab/baiyun/FASTQ/SAMPLE12.bam
/home/wangk/lab/baiyun/FASTQ/SAMPLE13.bam
```

skin.bed

```
chr1    156084501       156085065       LMNA    NM_170707.1
chr1    156100407       156100564       LMNA    NM_170707.2
chr1    156104193       156104319       LMNA    NM_170707.3
chr1    156104595       156104766       LMNA    NM_170707.4
chr1    156104977       156105103       LMNA    NM_170707.5
chr1    156105691       156105912       LMNA    NM_170707.6
chr1    156106004       156106227       LMNA    NM_170707.7
chr1    156106711       156106819       LMNA    NM_170707.8
chr1    156106903       156107023       LMNA    NM_170707.9
chr1    156107444       156107534       LMNA    NM_170707.10
chr1    156108278       156108548       LMNA    NM_170707.11
chr1    156108870       156109872       LMNA    NM_170707.12
chr1    156095960       156096014       LMNA    NM_001257374.1
chr1    156100407       156100564       LMNA    NM_001257374.2
chr1    156104193       156104319       LMNA    NM_001257374.3
chr1    156104595       156104766       LMNA    NM_001257374.4
chr1    156104977       156105103       LMNA    NM_001257374.5
chr1    156105691       156105912       LMNA    NM_001257374.6
chr1    156106004       156106227       LMNA    NM_001257374.7
chr1    156106711       156106819       LMNA    NM_001257374.8
chr1    156106903       156107023       LMNA    NM_001257374.9
chr1    156107444       156107534       LMNA    NM_001257374.10
chr1    156108278       156108548       LMNA    NM_001257374.11
chr1    156108870       156108893       LMNA    NM_001257374.12
chr1    156109560       156109872       LMNA    NM_001257374.13
chr1    156084460       156085065       LMNA    NM_170708.1
chr1    156100407       156100564       LMNA    NM_170708.2
chr1    156104193       156104319       LMNA    NM_170708.3
chr1    156104595       156104766       LMNA    NM_170708.4
chr1    156104977       156105103       LMNA    NM_170708.5
chr1    156105691       156105912       LMNA    NM_170708.6
chr1    156106004       156106227       LMNA    NM_170708.7
chr1    156106711       156106819       LMNA    NM_170708.8
chr1    156106903       156107023       LMNA    NM_170708.9
chr1    156108278       156108548       LMNA    NM_170708.10
chr1    156108870       156109880       LMNA    NM_170708.11
chr1    156084460       156085065       LMNA    NM_005572.1
chr1    156100407       156100564       LMNA    NM_005572.2
```

# 3. IdentifyFailures
```
Rscript IdentifyFailures.R --rds CNV_COUNTS.RData --cor 0.97 --cov 50 --out CNV_COUNTS_OUT
```

产生文件：CNV_COUNTS_OUT_Failures.txt

```
Sample  Gene    Exon    Type    Info
All     INDEX:1|LMNA    INDEX:1|NM_170707.1     Whole exon      Low median read depth (FPKM):  0
All     INDEX:2|LMNA    INDEX:2|NM_170707.2     Whole exon      Low median read depth (FPKM):  0
All     INDEX:3|LMNA    INDEX:3|NM_170707.3     Whole exon      Low median read depth (FPKM):  0
All     INDEX:4|LMNA    INDEX:4|NM_170707.4     Whole exon      Low median read depth (FPKM):  0
All     INDEX:5|LMNA    INDEX:5|NM_170707.5     Whole exon      Low median read depth (FPKM):  0
All     INDEX:6|LMNA    INDEX:6|NM_170707.6     Whole exon      Low median read depth (FPKM):  0
All     INDEX:7|LMNA    INDEX:7|NM_170707.7     Whole exon      Low median read depth (FPKM):  0
All     INDEX:8|LMNA    INDEX:8|NM_170707.8     Whole exon      Low median read depth (FPKM):  0
All     INDEX:9|LMNA    INDEX:9|NM_170707.9     Whole exon      Low median read depth (FPKM):  0
All     INDEX:67|LMNA   INDEX:67|NM_001282626.8 Whole exon      Low median read depth (FPKM):  0
All     INDEX:139|ATM   INDEX:139|NM_000051.45  Whole exon      Low median read depth (FPKM):  0
All     INDEX:140|ATM   INDEX:140|NM_000051.46  Whole exon      Low median read depth (FPKM):  0
All     INDEX:141|ATM   INDEX:141|NM_000051.47  Whole exon      Low median read depth (FPKM):  0
All     INDEX:142|ATM   INDEX:142|NM_000051.48  Whole exon      Low median read depth (FPKM):  0
All     INDEX:143|ATM   INDEX:143|NM_000051.49  Whole exon      Low median read depth (FPKM):  0
All     INDEX:144|ATM   INDEX:144|NM_000051.50  Whole exon      Low median read depth (FPKM):  0
All     INDEX:145|ATM   INDEX:145|NM_000051.51  Whole exon      Low median read depth (FPKM):  0
All     INDEX:146|ATM   INDEX:146|NM_000051.52  Whole exon      Low median read depth (FPKM):  0
All     INDEX:147|ATM   INDEX:147|NM_000051.53  Whole exon      Low median read depth (FPKM):  0
All     INDEX:148|ATM   INDEX:148|NM_000051.54  Whole exon      Low median read depth (FPKM):  0
All     INDEX:149|ATM   INDEX:149|NM_000051.55  Whole exon      Low median read depth (FPKM):  0
All     INDEX:150|ATM   INDEX:150|NM_000051.56  Whole exon      Low median read depth (FPKM):  0
All     INDEX:151|ATM   INDEX:151|NM_000051.57  Whole exon      Low median read depth (FPKM):  0
All     INDEX:152|ATM   INDEX:152|NM_000051.58  Whole exon      Low median read depth (FPKM):  0
All     INDEX:153|ATM   INDEX:153|NM_000051.59  Whole exon      Low median read depth (FPKM):  0
All     INDEX:155|ATM   INDEX:155|NM_000051.61  Whole exon      Low median read depth (FPKM):  0
All     INDEX:157|ATM   INDEX:157|NM_000051.63  Whole exon      Low median read depth (FPKM):  0
All     INDEX:158|BRCA2 INDEX:158|NM_000059.1   Whole exon      Low median read depth (FPKM):  0
All     INDEX:159|BRCA2 INDEX:159|NM_000059.2   Whole exon      Low median read depth (FPKM):  0
All     INDEX:183|BRCA2 INDEX:183|NM_000059.26  Whole exon      Low median read depth (FPKM):  0
All     INDEX:184|BRCA2 INDEX:184|NM_000059.27  Whole exon      Low median read depth (FPKM):  0
All     INDEX:185|SERPING1      INDEX:185|NM_000062.1   Whole exon      Low median read depth (FPKM):  0
All     INDEX:186|SERPING1      INDEX:186|NM_000062.2   Whole exon      Low median read depth (FPKM):  0
All     INDEX:187|SERPING1      INDEX:187|NM_000062.3   Whole exon      Low median read depth (FPKM):  0
All     INDEX:188|SERPING1      INDEX:188|NM_000062.4   Whole exon      Low median read depth (FPKM):  0
All     INDEX:189|SERPING1      INDEX:189|NM_000062.5   Whole exon      Low median read depth (FPKM):  0
All     INDEX:190|SERPING1      INDEX:190|NM_000062.6   Whole exon      Low median read depth (FPKM):  0
```

去掉未通过的exons：
```
perl exclude_failures.pl CNV_COUNTS_OUT_Failures.txt skin.bed > skin.pass.bed
```
或者使用生成的CNV_COUNTS_OUT.exon.bed重新作分析，从ReadInBams开始。

# 4. cnvcalls/makeCNVcalls
```
#Rscript makeCNVcalls.R --rds CNV_COUNTS.RData --out CNV_result --plot All --plotFolder CNV_Plots
ll *.rds | sed -e 's/.rds//g' | awk '{print "Rscript cnvcalls.R --rds " $9 ".RData --out " $9 ".CNV & sleep 1s"}'
Rscript cnvcalls.R --rds CNV_COUNTS.RData --out CNV_COUNTS.CNV & sleep 1s
or
Rscript cnvcalls.R --bam bamlist.txt --bed skin.pass.bed --fsa hg19.fa --out CNV_results
```

