# 利用lumpy从WGS数据中检测CNV

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
conda create -n lumpy
conda activate lumpy
conda install -c bioconda lumpy-sv samblaster svtyper scipy gatk4 samtools bwa trimmomatic picard r-biocmanager r-ggplot2 r-getopt r-optparse r-reshape bioconductor-genomicranges
#conda install r-biocmanager r-ggplot2 r-getopt r-optparse r-reshape
#conda install -c bioconda lumpy-sv samblaster svtyper scipy gatk4 samtools bwa trimmomatic picard
```


# Download UCSC hg19 genome
```
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/genes/hg19.refGene.gtf.gz
perl ext_ucsc_human_canonical_fa.pl hg19.fa > hg19.fa.1  
mv hg19.fa.1 hg19.fa
```


# 1. Reads mapping
```
bwa index -a bwtsw hg19.fa  

ln -s rawdata/Sample_R22035863-RNDVA2946Q-VA2946Q/R22035863-RNDVA2946Q-VA2946Q_combined_R1.fastq.gz VA2946Q_R1.fq.gz
ln -s rawdata/Sample_R22035863-RNDVA2946Q-VA2946Q/R22035863-RNDVA2946Q-VA2946Q_combined_R2.fastq.gz VA2946Q_R2.fq.gz
ln -s rawdata/Sample_R22035864-RNDCA2519B-CA2519B/R22035864-RNDCA2519B-CA2519B_combined_R1.fastq.gz CA2519B_R1.fq.gz
ln -s rawdata/Sample_R22035864-RNDCA2519B-CA2519B/R22035864-RNDCA2519B-CA2519B_combined_R2.fastq.gz CA2519B_R2.fq.gz
ln -s rawdata/Sample_R22035865-RNDCA2594B-CA2594B/R22035865-RNDCA2594B-CA2594B_combined_R1.fastq.gz CA2594B_R1.fq.gz
ln -s rawdata/Sample_R22035865-RNDCA2594B-CA2594B/R22035865-RNDCA2594B-CA2594B_combined_R2.fastq.gz CA2594B_R2.fq.gz
ln -s rawdata/Sample_R22035866-RNDCA2608B-CA2608B/R22035866-RNDCA2608B-CA2608B_combined_R1.fastq.gz CA2608B_R1.fq.gz
ln -s rawdata/Sample_R22035866-RNDCA2608B-CA2608B/R22035866-RNDCA2608B-CA2608B_combined_R2.fastq.gz CA2608B_R2.fq.gz
ln -s rawdata/Sample_R22035867-RNDCT3-4-CT3-4/R22035867-RNDCT3-4-CT3-4_combined_R1.fastq.gz CT3_4_R1.fq.gz
ln -s rawdata/Sample_R22035867-RNDCT3-4-CT3-4/R22035867-RNDCT3-4-CT3-4_combined_R2.fastq.gz CT3_4_R2.fq.gz

# ll *_1.fastq | sed -e 's/_1.fastq//g' | awk '{print "trimmomatic PE -threads 18 -phred33 " $9 "_1.fastq " $9 "_2.fastq " $9 "_1.fq " $9 "_1.unpaired.fq " $9 "_2.fq "  $9 "_2.unpaired.fq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40 AVGQUAL:20"}'  

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
  
  #speedseq align -R "@RG\tID:${ID}\tPU:${PU}\tPL:${PL}\tLB:${LB}\tSM:${SM}" hg19.fa ${SAMPLE}_R1.fq.gz ${SAMPLE}_R2.fq.gz
  nohup bwa mem -t 18 -M -R "@RG\tID:${ID}\tPU:${PU}\tPL:${PL}\tLB:${LB}\tSM:${SM}" -o ${SAMPLE}.sam /usr/local/db/ucsc/human/BWA_hg19_genome/hg19.fa ${SAMPLE}_R1.fq.gz ${SAMPLE}_R2.fq.gz > ${SAMPLE}.nohup 2>&1 &
  sleep 1s
done

# 注意参考基因组要用全路径，因为索引文件都在这个目录下，不要用链接：
# ln -s /usr/local/db/ucsc/human/BWA_hg19_genome/hg19.fa

如果忘记添加文库信息，可以用：
picard AddOrReplaceReadGroups -I ${SAMPLE}.sam -O ${SAMPLE}.H.sam -R hg19.fa -ID ${ID} -PU ${PU} -PL ${PL} -LB ${LB} -SM ${SM}

#ll *_1.fq | sed -e 's/_1.fq//g'| awk '{print "bwa mem -t 18 -M -R "@RG\tID:${ID}\tPU:${PU}\tPL:${PL}\tLB:${LB}\tSM:${SM}\" -o \${SAMPLE}\.sam /usr/local/db/ucsc/human/BWA_hg19_genome/hg19.fa " $9 "_1.fq " $9 "_2.fq & sleep 1s"}'
#ll *.sam | sed -e 's/.sam//g' | awk '{print "samtools view -bS -o " $9 ".tmp " $9 ".sam & sleep 1s"}'
#ll *.tmp | sed -e 's/.tmp//g' | awk '{print "samtools sort -o " $9 ".bam " $9 ".tmp & sleep 1s"}'
#ll *.bam | awk '{print "samtools index " $9 " & sleep 1s"}'
```


# 2. 用samblaster进行markduplicate
```
samples=`ll *.sam | sed -e 's/.sam//g' | awk '{print $9}'`
for SAMPLE in ${samples}
do
  #echo ${SAMPLE}
  nohup samblaster -i ${SAMPLE}.sam -o ${SAMPLE}.RM.sam --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 > ${SAMPLE}.nohup 2>&1 & sleep 1s
  #echo "nohup picard MarkDuplicates -I ${SAMPLE}.bam -O ${SAMPLE}.M.bam -M ${SAMPLE}.metrics > ${SAMPLE}.nohup 2>&1 &"
  #sleep 1s
done

或者：

ll *.sam | sed -e 's/.sam//g' | awk '{print "samblaster -i " $9 ".sam -o " $9 ".RM.sam --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 & sleep 1s"}'
ll *.RM.sam | sed -e 's/.RM.sam//g' | awk '{print "samtools view -Sb -o " $9 ".tmp " $9 ".RM.sam & sleep 1s"}' 
ll *.tmp | sed -e 's/.tmp//g' | awk '{print "samtools sort -o " $9 ".bam " $9 ".tmp & sleep 1s"}'
ll *.bam | awk '{print "samtools index " $9 " & sleep 1s"}'

samblaster 主要参数：
-i --input 输入sam文件（必须包含header且按reads id排序）
-o --output 输出sam文件
-d --discordantFile 输出discordant read pairs
-s --splitterFile 输出split reads
-u --unmappedFile 输出unmapped/clipped reads

其他参数:
-a --acceptDupMarks 不去重
-e --excludeDups 去掉discordant, splitter, and/or unmapped等重复（具体定义详见samblaster主页）
-r --removeDups 去掉重复(-e --excludeDups类似)
--addMateTags 添加MC and MQ tags

-M 与bwa mem -M 类似
```

# 3. Extract the discordant paired-end alignments
```
ll *.bam | sed -e 's/.bam//g' | awk '{print "samtools view -b -F 1294 -o " $9 ".discordants.unsorted.bam " $9 ".bam & sleep 1s"}'

samtools view -b -F 1294 -o CA2519B.discordants.unsorted.bam CA2519B.bam & sleep 1s
samtools view -b -F 1294 -o CA2594B.discordants.unsorted.bam CA2594B.bam & sleep 1s
samtools view -b -F 1294 -o CA2608B.discordants.unsorted.bam CA2608B.bam & sleep 1s
samtools view -b -F 1294 -o CT3_4.discordants.unsorted.bam CT3_4.bam & sleep 1s
samtools view -b -F 1294 -o VA2946Q.discordants.unsorted.bam VA2946Q.bam & sleep 1s
```


# 4. Extract the split-read alignments
```
ll *.discordants.unsorted.bam | sed -e 's/.discordants.unsorted.bam//g' | awk '{print "samtools view -h " $9 ".bam | extractSplitReads_BwaMem -i stdin | samtools view -Sb -o " $9 ".splitters.unsorted.bam - & sleep 1s"}'

samtools view -h CA2519B.bam | extractSplitReads_BwaMem -i stdin | samtools view -Sb -o CA2519B.splitters.unsorted.bam - & sleep 1s
samtools view -h CA2594B.bam | extractSplitReads_BwaMem -i stdin | samtools view -Sb -o CA2594B.splitters.unsorted.bam - & sleep 1s
samtools view -h CA2608B.bam | extractSplitReads_BwaMem -i stdin | samtools view -Sb -o CA2608B.splitters.unsorted.bam - & sleep 1s
samtools view -h CT3_4.bam | extractSplitReads_BwaMem -i stdin | samtools view -Sb -o CT3_4.splitters.unsorted.bam - & sleep 1s
samtools view -h VA2946Q.bam | extractSplitReads_BwaMem -i stdin | samtools view -Sb -o VA2946Q.splitters.unsorted.bam - & sleep 1s
```


# 5. Sort discordant & split alignments
```
ll *.discordants.unsorted.bam | sed -e 's/.discordants.unsorted.bam//g' | awk '{print "samtools sort -o " $9 ".discordants.bam " $9 ".discordants.unsorted.bam & sleep 1s"}'
ll *.splitters.unsorted.bam | sed -e 's/.splitters.unsorted.bam//g' | awk '{print "samtools sort -o " $9 ".splitters.bam " $9 ".splitters.unsorted.bam & sleep 1s"}'

samtools sort -o CA2519B.discordants.bam CA2519B.discordants.unsorted.bam & sleep 1s
samtools sort -o CA2594B.discordants.bam CA2594B.discordants.unsorted.bam & sleep 1s
samtools sort -o CA2608B.discordants.bam CA2608B.discordants.unsorted.bam & sleep 1s
samtools sort -o CT3_4.discordants.bam CT3_4.discordants.unsorted.bam & sleep 1s
samtools sort -o VA2946Q.discordants.bam VA2946Q.discordants.unsorted.bam & sleep 1s

samtools sort -o CA2519B.splitters.bam CA2519B.splitters.unsorted.bam & sleep 1s
samtools sort -o CA2594B.splitters.bam CA2594B.splitters.unsorted.bam & sleep 1s
samtools sort -o CA2608B.splitters.bam CA2608B.splitters.unsorted.bam & sleep 1s
samtools sort -o CT3_4.splitters.bam CT3_4.splitters.unsorted.bam & sleep 1s
samtools sort -o VA2946Q.splitters.bam VA2946Q.splitters.unsorted.bam & sleep 1s
```


# 6. Running lumpyexpress
```
samples=`ll *.discordants.bam | sed -e 's/.discordants.bam//g' | awk '{print $9}'`
for SAMPLE in ${samples}
do
  nohup lumpyexpress -B ${SAMPLE}.bam -S ${SAMPLE}.splitters.bam -D ${SAMPLE}.discordants.bam -o ${SAMPLE}.vcf > ${SAMPLE}.nohup 2>&1 & sleep 1s
done

ll *.discordants.bam | sed -e 's/.discordants.bam//g' | awk '{print "lumpyexpress -B " $9 ".bam -S " $9 ".splitters.bam -D " $9 ".discordants.bam -o " $9 ".vcf"}'

# single sample
lumpyexpress -B NA12877.bam -S NA12877.splitters.bam -D NA12877.discordants.bam -o NA12877.vcf

# multiple samples
lumpyexpress -B NA12877.bam,NA12878.bam,NA12879.bam -S NA12877.splitters.bam,NA12878.splitters.bam,NA12879.splitters.bam -D NA12877.discordants.bam,NA12878.discordants.bam,NA12879.discordants.bam -o multi_sample.vcf

# tumor-normal pair
lumpyexpress -B tumor.bam,normal.bam -S tumor.splitters.bam,normal.splitters.bam -D tumor.discordants.bam,normal.discordants.bam -o tumor_normal.vcf
```


# 7. Call genotypes using svtyper
```
samples=`ll *.vcf | sed -e 's/.vcf//g' | awk '{print $9}'`
for SAMPLE in ${samples}
do
  nohup svtyper -B ${SAMPLE}.bam -i ${SAMPLE}.vcf -o ${SAMPLE}.genotypes.vcf > ${SAMPLE}.nohup 2>&1 & sleep 1s
done
```

# 8. AnnotSV 注释
```
ll *.vcf | sed -e 's/.vcf//g' |awk '{print "AnnotSV -genomeBuild GRCh37 -SVinputFile " $9 ".vcf -outputFile " $9 ".annotated.tsv -svtBEDcol 4"}'
```
AnnotSV -genomeBuild GRCh37 -SVinputFile NA12877.vcf -outputFile NA12877.annotated.tsv -svtBEDcol 4


# 9. ClassifyCNV致病性注释

https://github.com/Genotek/ClassifyCNV

The numeric pathogenicity score, calculated by ClassifyCNV, is converted to pathogenicity classification using the following cutoffs:
```
≤ 0.99: benign variant
0.90 .. 0.98: likely benign variant
0.89 .. 0.89: variant of uncertain significance
0.90 .. 0.98: likely pathogenic variant
≥ 0.99: pathogenic variant
```

运行实例：
```
# python3 ClassifyCNV.py --infile YourCNVFile.bed --GenomeBuild {hg19,hg38}
ll *.bed |sed -e 's/.bed//g' | awk '{print "ClassifyCNV.py --infile " $9 ".bed --outdir " $9 " --GenomeBuild hg19"}'

ClassifyCNV.py --infile NA12877.bed --outdir NA12877 --GenomeBuild hg19
```

Input

ClassifyCNV accepts a BED file as input. The file must include the following columns in this order:

chromosome
CNV start position
CNV end position
CNV type (DEL or DUP)
CNVs on alternative contigs are not evaluated. Both hg19 and hg38 coordinates are supported.


NA12877.bed数据格式参考ClassifyCNV网站 

https://github.com/Genotek/ClassifyCNV/blob/master/Examples/
hg19.bed：
```
chr1  668630  850204  DUP
chr1  738570  742020  DEL
chr1  766600  769112  DEL
chr1  775292  791968  DEL
chr1  794496  799549  DEL
chr1  873391  874042  DEL
chr1  939918  939968  DEL
chr1  947113  948003  DEL
chr1  963826  974172  DUP
chr1  988572  988623  DEL
chr1  1019492  1019819  DEL
```


