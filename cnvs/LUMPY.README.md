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

遇到的问题：

```
问题1：
libcrypto.so.1.0.0: cannot open shared object file: No such file or directory

cd /usr/local/prog/anaconda/envs/lumpy/lib
ln -s libcrypto.so.3 libcrypto.so.1.0.0

问题2：
Checking for required python modules (/usr/local/prog/anaconda/envs/lumpy/bin/python)...
/usr/local/prog/anaconda/envs/lumpy/bin/lumpyexpress:行15: -n：未找到命令
/usr/local/prog/anaconda/envs/lumpy/bin/lumpyexpress: 第 16 行：[: ==：需要一元表达式

在/usr/local/prog/anaconda/envs/lumpy/bin/lumpyexpress.config中添加 
HEXDUMP=`which hexdump || true`
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


# 统计对照

https://pubmed.ncbi.nlm.nih.gov/24970577/  
LUMPY: a probabilistic framework for structural variant discovery  

To estimate sensitivity and FDR, we compared predictions made by each tool to two truth sets: 1) 3,376 validated, non-overlapping deletions from the 1000 Genomes Project [12] (Additional file 3); and 2) 4,095 deletions that were detected by at least one tool in the 50X dataset, or that were reported by Mills et al. [12] (which used numerous SV detection tools), and that were validated by split-read mapping analysis of independent long-read sequencing data from PacBio or Illumina Moleculo platforms (Additional file 4).



13059_2013_3363_MOESM1_ESM.zip  
Additional file 1: This file contains the breakpoints used for the homozygous variant simulation. The format is BEDPE. Each line has a ‘TYPE:’ field that indicates DELETION, DUPLICATION, INVERSION, or TRANSLOCATION. (ZIP 315 KB)

13059_2013_3363_MOESM2_ESM.zip  
Additional file 2: This file contains the breakpoints used for the heterogeneous tumor simulation. The format is BEDPE. These deletions are based on the variants released by the 1000 Genomes Project in [29]. We selected non-overlapping deletions that were at least 50 bases long and successfully lifted over from build 36 of the human reference genome to build 37. (ZIP 170 KB)

13059_2013_3363_MOESM3_ESM.zip  
Additional file 3: This file contains the breakpoints used for the Mills et al . truth set. The format is BEDPE, which is described by [30]. The breakpoints are the non-overlapping validated deletions observed in NA12878 and are based on the variants given in [31]. (ZIP 83 KB)

13059_2013_3363_MOESM4_ESM.zip  
Additional file 4: This file contains the breakpoint intervals for the deletion predictions that were made by LUMPY (pe?+?sr, trio, prior, pe?+?sr&rd), GASVPro, DELLY, Pindel, or the 1000 Genomes Project [12] and that were validated by long read alignments from PacBio and/or Illumina Moleculo sequencing. The format is BEDPE, and the score field indicates the number of overlapping predictions. Note: in some cases one algorithm made two predictions that contributed to a single call. For example, there are two calls with a score of 8. In both cases Pindel contributed two very similar calls. (ZIP 77 KB)

13059_2013_3363_MOESM5_ESM.zip  
Additional file 5: This file contains the calls made by LUMPY for NA12878 with paired-end and split-read evidence that were also validated with PacBio/Moleculo data. The format is BEDPE and the score field is the total amount of supporting evidence. This file contains extra fields that are described at [16]. (ZIP 162 KB)

https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2014-15-6-r84/MediaObjects/13059_2013_3363_MOESM1_ESM.zip
https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2014-15-6-r84/MediaObjects/13059_2013_3363_MOESM2_ESM.zip
https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2014-15-6-r84/MediaObjects/13059_2013_3363_MOESM3_ESM.zip
https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2014-15-6-r84/MediaObjects/13059_2013_3363_MOESM4_ESM.zip
https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2014-15-6-r84/MediaObjects/13059_2013_3363_MOESM5_ESM.zip

```
perl lumpy_cnv_lens_distribution.pl LUMPY.VCF 50 > lens_distribution.txt
可以看到lumpy预测的cnv窦在2500 bp 以下

wc -l  3717462611446476_add3.bedpe
3376 3717462611446476_add3.bedpe
cp 3717462611446476_add3.bedpe 3376_deletions.bedpe

wc -l 3717462611446476_add4.bedpe
4095 3717462611446476_add4.bedpe
cp 3717462611446476_add4.bedpe 4095_deletions.bedpe

samples=`ll *.vcf | sed -e 's/.vcf//g' | awk '{print $9}'`
for SAMPLE in ${samples}
do
  echo "perl ext_lumpy_cnv.pl 3376_deletions.bedpe ${SAMPLE}.vcf DEL > overlapping_lumpy.${SAMPLE}.txt"
done

Validated variants from the NA12878 genome. 
http://www.nature.com/nature/journal/v470/n7332/extref/nature09708-s5.zip

perl ext_lumpy_cnv.pl 3376_deletions.bedpe CT3_4.vcf DEL > CT3_4.txt
2540 / 3376 (75.24 %) of DEL CNV deletected!


[wangk@cs angpu]$ awk '$5 == "<DEL>" {print}' CT3_4.vcf | wc -l
4281
[wangk@cs angpu]$ sort -u CT3_4.txt | wc -l 
2449

该覆盖度下性能：
sensibility = 2540 / 3376 (75.24 %)
specificity = 2449 / 4281 (57.21 %)

```


# 优化指标：
```
测试WGS CNV检出软件，优化检出流程，需达到如下指标：
1、本地测序的4个CNV阳性样本，阳性位点全检出。
2、下载5例30X左右的已知CNV结果的WGS数据检测流程的灵敏度和特异性，灵敏度希望能达到95%以上，特异性越高越好，希望能达到85%以上。
3、本地测了NA12878样本，编号CT3-4。可用于检测流程的灵敏度和特异性，NA12878 CNV结果可从公开数据下载，阳性位点可参考各检测平台检出的可信CNV区段。

本地测试的4个阳性样本阳性位点如下：

CNV
VA2946Q PRKN (PARK2) 基因 3 号外显子杂合
CA2519B PRKN 基因 3 号外显子和 5-7 号外显子区域存在杂合缺失
CA2594B RB1
CA2608B PMP22，TEKT3
```



```
samples="PRKN RB1 PMP22 TEKT3"
for SAMPLE in ${samples}
do
  #echo ${SAMPLE}
  echo "awk '\$4 == \"${SAMPLE}\" {print}' hg19.exon.bed >> T.bed"
  #awk '$4 == $SAMPLE {print}' hg19.exon.bed >> T.bed
done

# T.bed
chr6    161768448       161771243       PRKN    NM_004562.1
chr6    161781119       161781237       PRKN    NM_004562.2
chr6    161807825       161807909       PRKN    NM_004562.3
chr6    161969885       161970035       PRKN    NM_004562.4
chr6    161990386       161990448       PRKN    NM_004562.5
chr6    162206803       162206940       PRKN    NM_004562.6
chr6    162394333       162394449       PRKN    NM_004562.7
chr6    162475122       162475206       PRKN    NM_004562.8
chr6    162622162       162622284       PRKN    NM_004562.9
chr6    162683556       162683797       PRKN    NM_004562.10
chr6    162864341       162864505       PRKN    NM_004562.11
chr6    163148693       163148798       PRKN    NM_004562.12
chr6    161768448       161771243       PRKN    NM_013988.1
chr6    161781119       161781237       PRKN    NM_013988.2
chr6    161807825       161807909       PRKN    NM_013988.3
chr6    161969885       161970035       PRKN    NM_013988.4
chr6    161990386       161990448       PRKN    NM_013988.5
chr6    162206803       162206940       PRKN    NM_013988.6
chr6    162394333       162394449       PRKN    NM_013988.7
chr6    162864341       162864505       PRKN    NM_013988.8
chr6    163148693       163148798       PRKN    NM_013988.9
chr6    161768448       161771243       PRKN    NM_013987.1
chr6    161781119       161781237       PRKN    NM_013987.2
chr6    161807825       161807909       PRKN    NM_013987.3
chr6    161969885       161970035       PRKN    NM_013987.4
chr6    161990386       161990448       PRKN    NM_013987.5
chr6    162206803       162206940       PRKN    NM_013987.6
chr6    162394333       162394449       PRKN    NM_013987.7
chr6    162622162       162622284       PRKN    NM_013987.8
chr6    162683556       162683797       PRKN    NM_013987.9
chr6    162864341       162864505       PRKN    NM_013987.10
chr6    163148693       163148798       PRKN    NM_013987.11
chr13   48877882        48878185        RB1     NM_000321.1
chr13   48881415        48881542        RB1     NM_000321.2
chr13   48916734        48916850        RB1     NM_000321.3
chr13   48919215        48919335        RB1     NM_000321.4
chr13   48921960        48921999        RB1     NM_000321.5
chr13   48923091        48923159        RB1     NM_000321.6
chr13   48934152        48934263        RB1     NM_000321.7
chr13   48936950        48937093        RB1     NM_000321.8
chr13   48939029        48939107        RB1     NM_000321.9
chr13   48941629        48941739        RB1     NM_000321.10
chr13   48942662        48942740        RB1     NM_000321.11
chr13   48947540        48947628        RB1     NM_000321.12
chr13   48951053        48951170        RB1     NM_000321.13
chr13   48953729        48953786        RB1     NM_000321.14
chr13   48954188        48954220        RB1     NM_000321.15
chr13   48954300        48954377        RB1     NM_000321.16
chr13   48955382        48955579        RB1     NM_000321.17
chr13   49027128        49027247        RB1     NM_000321.18
chr13   49030339        49030485        RB1     NM_000321.19
chr13   49033823        49033969        RB1     NM_000321.20
chr13   49037866        49037971        RB1     NM_000321.21
chr13   49039133        49039247        RB1     NM_000321.22
chr13   49039340        49039504        RB1     NM_000321.23
chr13   49047495        49047526        RB1     NM_000321.24
chr13   49050836        49050979        RB1     NM_000321.25
chr13   49051490        49051540        RB1     NM_000321.26
chr13   49054133        49056026        RB1     NM_000321.27
chr17   15133095        15134397        PMP22   NM_000304.1
chr17   15142787        15142928        PMP22   NM_000304.2
chr17   15162410        15162510        PMP22   NM_000304.3
chr17   15163966        15164078        PMP22   NM_000304.4
chr17   15168470        15168643        PMP22   NM_000304.5
chr17   15133095        15134397        PMP22   NM_001281456.1
chr17   15142787        15142928        PMP22   NM_001281456.2
chr17   15162410        15162510        PMP22   NM_001281456.3
chr17   15163966        15164074        PMP22   NM_001281456.4
chr17   15168470        15168643        PMP22   NM_001281456.5
chr17   15133093        15134397        PMP22   NM_001281455.1
chr17   15142787        15142928        PMP22   NM_001281455.2
chr17   15162410        15162510        PMP22   NM_001281455.3
chr17   15163966        15164078        PMP22   NM_001281455.4
chr17   15165741        15165906        PMP22   NM_001281455.5
chr17   15138535        15138601        PMP22   NM_001330143.1
chr17   15142787        15142928        PMP22   NM_001330143.2
chr17   15162410        15162510        PMP22   NM_001330143.3
chr17   15163966        15164078        PMP22   NM_001330143.4
chr17   15168470        15168643        PMP22   NM_001330143.5
chr17   15133095        15134397        PMP22   NM_153322.1
chr17   15142787        15142928        PMP22   NM_153322.2
chr17   15162410        15162510        PMP22   NM_153322.3
chr17   15163966        15164137        PMP22   NM_153322.4
chr17   15133095        15134397        PMP22   NM_153321.1
chr17   15142787        15142928        PMP22   NM_153321.2
chr17   15162410        15162510        PMP22   NM_153321.3
chr17   15163966        15164078        PMP22   NM_153321.4
chr17   15165745        15165869        PMP22   NM_153321.5
chr17   15207128        15207469        TEKT3   NM_031898.1
chr17   15211980        15212135        TEKT3   NM_031898.2
chr17   15215575        15215798        TEKT3   NM_031898.3
chr17   15217403        15217547        TEKT3   NM_031898.4
chr17   15222393        15222464        TEKT3   NM_031898.5
chr17   15231308        15231392        TEKT3   NM_031898.6
chr17   15234323        15234931        TEKT3   NM_031898.7
chr17   15243344        15243378        TEKT3   NM_031898.8
chr17   15244834        15244914        TEKT3   NM_031898.9
```

查看是否预测到：
```
perl ext_lumpy_cnv_gene_overlapping.pl T.bed VA2946Q.vcf PRKN > VA2946Q.PRKN.txt

3号外显子检测到一个DUP
chr6    67315800        2378    N       <DUP>   .       .       SVTYPE=DUP;STRANDS=-+:10;SVLEN=100162734;END=167478534;CIPOS=-10,9;CIEND=-10,9;CIPOS95=0,0;CIEND95=0,0;SU=10;PE=0;SR=10        GT:SU:PE:SR     ./.:10:0:10     |       chr6    161807825       161807909     PRKN     NM_013987.3

其他外显子也多检测到DUP


perl ext_lumpy_cnv_gene_overlapping.pl T.bed CA2519B.vcf PRKN > CA2519B.PRKN.txt

3号外显子检测到一个DEL，一个DUP
chr6    5162756 2364    N       <DEL>   .       .       SVTYPE=DEL;STRANDS=+-:6;SVLEN=-156727823;END=161890579;CIPOS=-6,5;CIEND=-6,5;CIPOS95=0,0;CIEND95=0,0;SU=6;PE=0;SR=6    GT:SU:PE:SR     ./.:6:0:6       |       chr6    161807825       161807909       PRKN    NM_013987.3
chr6    67315800        2375    N       <DUP>   .       .       SVTYPE=DUP;STRANDS=-+:7;SVLEN=100162734;END=167478534;CIPOS=-5,4;CIEND=-5,4;CIPOS95=0,0;CIEND95=0,0;SU=7;PE=0;SR=7     GT:SU:PE:SR     ./.:7:0:7       |       chr6    161807825       161807909       PRKN  NM_013987.3

5号外显子检测到一个DUP
chr6    67315800        2375    N       <DUP>   .       .       SVTYPE=DUP;STRANDS=-+:7;SVLEN=100162734;END=167478534;CIPOS=-5,4;CIEND=-5,4;CIPOS95=0,0;CIEND95=0,0;SU=7;PE=0;SR=7     GT:SU:PE:SR     ./.:7:0:7       |       chr6    161990386       161990448       PRKN  NM_013987.5

6号外显子检测到一个DEL，一个DUP
chr6    162018387       2365    N       <DEL>   .       .       SVTYPE=DEL;STRANDS=+-:27;SVLEN=-534841;END=162553228;CIPOS=-8,9;CIEND=-10,9;CIPOS95=0,1;CIEND95=0,0;IMPRECISE;SU=27;PE=21;SR=6 GT:SU:PE:SR     ./.:27:21:6     |       chr6    162206803       162206940     PRKN     NM_013987.6
chr6    67315800        2375    N       <DUP>   .       .       SVTYPE=DUP;STRANDS=-+:7;SVLEN=100162734;END=167478534;CIPOS=-5,4;CIEND=-5,4;CIPOS95=0,0;CIEND95=0,0;SU=7;PE=0;SR=7     GT:SU:PE:SR     ./.:7:0:7       |       chr6    162206803       162206940       PRKN  NM_013987.6

7号外显子检测到一个DEL，一个DUP
chr6    162018387       2365    N       <DEL>   .       .       SVTYPE=DEL;STRANDS=+-:27;SVLEN=-534841;END=162553228;CIPOS=-8,9;CIEND=-10,9;CIPOS95=0,1;CIEND95=0,0;IMPRECISE;SU=27;PE=21;SR=6 GT:SU:PE:SR     ./.:27:21:6     |       chr6    162394333       162394449     PRKN     NM_013987.7
chr6    67315800        2375    N       <DUP>   .       .       SVTYPE=DUP;STRANDS=-+:7;SVLEN=100162734;END=167478534;CIPOS=-5,4;CIEND=-5,4;CIPOS95=0,0;CIEND95=0,0;SU=7;PE=0;SR=7     GT:SU:PE:SR     ./.:7:0:7       |       chr6    162394333       162394449       PRKN  NM_013987.7

其他外显子也多测到DUP


perl ext_lumpy_cnv_gene_overlapping.pl T.bed CA2594B.vcf RB1 > CA2594B.RB1.txt
1-2号外显子检测到DEL
chr13   48573292        4237    N       <DEL>   .       .       SVTYPE=DEL;STRANDS=+-:20;SVLEN=-319902;END=48893194;CIPOS=-10,9;CIEND=-10,9;CIPOS95=0,0;CIEND95=0,0;SU=20;PE=15;SR=5   GT:SU:PE:SR     ./.:20:15:5     |       chr13   48877882        48878185        RB1   NM_000321.1
chr13   48573292        4237    N       <DEL>   .       .       SVTYPE=DEL;STRANDS=+-:20;SVLEN=-319902;END=48893194;CIPOS=-10,9;CIEND=-10,9;CIPOS95=0,0;CIEND95=0,0;SU=20;PE=15;SR=5   GT:SU:PE:SR     ./.:20:15:5     |       chr13   48881415        48881542        RB1   NM_000321.2


perl ext_lumpy_cnv_gene_overlapping.pl T.bed CA2608B.vcf PMP22 > CA2608B.PMP22.txt
未检测到PMP22的CNV变异


perl ext_lumpy_cnv_gene_overlapping.pl T.bed CA2608B.vcf TEKT3 > CA2608B.TEKT3.txt
未检测到TEKT3的CNV变异



```


```
#==============================================
方法1 conda安装python 2.7环境
conda create -n lumpy
WARNING: A conda environment already exists at '/public/home/users/fdu010/.prog/anaconda/envs/lumpy'
Remove existing environment (y/[n])? y
conda activate lumpy

conda install -c bioconda lumpy-sv
#conda install -c bioconda svtyper
https://github.com/arq5x/lumpy-sv/archive/refs/heads/master.zip
cd /public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit


在/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/lumpyexpress.config中添加 
HEXDUMP=`which hexdump || true`


#!/bin/bash -e

# Find original directory of file, resovling symlinks
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in/246128#246128
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

# general
LUMPY_HOME=$DIR

LUMPY=`which lumpy || true`
SAMBLASTER=`which samblaster || true`
# either sambamba or samtools is required
SAMBAMBA=`which sambamba || true`
SAMTOOLS=`which samtools || true`

# python 2.7 or newer, must have pysam, numpy installed
PYTHON=`which python || true`

# python scripts
PAIREND_DISTRO=$LUMPY_HOME/scripts/pairend_distro.py
BAMGROUPREADS=$LUMPY_HOME/scripts/bamkit/bamgroupreads.py
BAMFILTERRG=$LUMPY_HOME/scripts/bamkit/bamfilterrg.py
BAMLIBS=$LUMPY_HOME/scripts/bamkit/bamlibs.py


方法2 下载源码安装：
python 3.10环境：
https://blog.csdn.net/kx453653102/article/details/107686297
https://www.cnblogs.com/crazychris/p/4213029.html
GCC-4.9.4:
http://ftp.tsukuba.wide.ad.jp/software/gcc/releases/gcc-4.9.4/gcc-4.9.4.tar.gz

export http_proxy="192.168.155.150:443"
export https_proxy="192.168.155.150:443"

sed -e 's/ftp/--no-check-certificate https/g' ./contrib/download_prerequisites 
sed -e 's/wget/wget --no-check-certificate/g' download_prerequisites

sh ./contrib/download_prerequisites
mkdir build
cd build

../configure --prefix=/public/home/users/fdu010/.prog/gcc --disable-multilib --enable-bootstrap --enable-threads=posix
#The following languages will be built: c,c++,fortran,java,lto,objc
make
make install


checking for texi2dvi... no
checking for C compiler vendor... gnu
checking CFLAGS_WARN for maximum warnings... -Wall
checking for an ANSI C-conforming const... make[2]: *** [configure-stage1-cloog] Killed
make[2]: Leaving directory `/public/home/users/fdu010/.prog/local/pkgs/gcc-4.9.4/build'
make[1]: *** [stage1-bubble] Error 2
make[1]: Leaving directory `/public/home/users/fdu010/.prog/local/pkgs/gcc-4.9.4/build'
make: *** [all] Error 2

cd /public/home/users/fdu010/.prog/local/pkgs/gcc-4.9.4/cloog
../configure --prefix=/public/home/users/fdu010/.prog/gcc
make
make install


----------------------------------------------------------------------
Libraries have been installed in:
   /public/home/users/fdu010/.prog/gcc/lib/../lib64

If you ever happen to want to link against installed libraries
in a given directory, LIBDIR, you must either use libtool, and
specify the full pathname of the library, or use the `-LLIBDIR'
flag during linking and do at least one of the following:
   - add LIBDIR to the `LD_LIBRARY_PATH' environment variable
     during execution
   - add LIBDIR to the `LD_RUN_PATH' environment variable
     during linking
   - use the `-Wl,-rpath -Wl,LIBDIR' linker flag
   - have your system administrator add LIBDIR to `/etc/ld.so.conf'

See any operating system documentation about shared libraries for
more information, such as the ld(1) and ld.so(8) manual pages.
----------------------------------------------------------------------


[fdu010@a110 ~]$ gcc --version
gcc (GCC) 4.9.4
Copyright (C) 2014 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


先安装：
yum install libffi-devel -y

openssl:
https://www.openssl.org/
https://www.openssl.org/source/openssl-1.1.1l.tar.gz


./config --prefix=/public/home/users/fdu010/.prog/anaconda/envs/lumpy
make && make install


最先安装：
python=3.10.0 
https://www.python.org/downloads/source/
https://www.python.org/ftp/python/3.10.0/Python-3.10.0.tgz
./configure --prefix=/public/home/users/fdu010/.prog/anaconda/envs/lumpy
make &
make install

ln -s python3 python
ln -s pip3 pip

actv lumpy 使用新的环境
(lumpy) [fdu010@a110 pysam-0.17.0]$ python setup.py install --prefix=/public/home/users/fdu010/.prog/anaconda/envs/lumpy

查看已经安装的模块：
>>> help('modules')


libffi:
http://www.sourceware.org/libffi/
https://github.com/libffi/libffi/releases/download/v3.4.2/libffi-3.4.2.tar.gz
./configure --prefix=/public/home/users/fdu010/.prog/anaconda/envs/lumpy
make
make install

setuptools==58.4.0  openssl-python 0.1.1 pysam NumPy 
python setup.py install --prefix=/public/home/users/fdu010/.prog/anaconda/envs/lumpy
ln -s python3 python

pip install pysam
https://github.com/pysam-developers/pysam

https://files.pythonhosted.org/packages/bf/a6/b1802ac26bb501ce8dada6f41dd83cdf591ecc834f2c208b07c6d97eba84/pysam-0.17.0.tar.gz
https://files.pythonhosted.org/packages/5f/d6/ad58ded26556eaeaa8c971e08b6466f17c4ac4d786cd3d800e26ce59cc01/numpy-1.21.3.zip

export http_proxy="http://192.168.155.150:443"
export https_proxy="https://192.168.155.150:443"




lumpy：
https://github.com/arq5x/lumpy-sv
https://github.com/hall-lab/speedseq

/public/home/users/fdu010/.prog/anaconda/envs/lumpy 
cp -rf scripts /public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/

bwa:
https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2
make
cp -rf * /public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin


samtools, htslib：
http://www.htslib.org/
https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2
./configure --prefix=/public/home/users/fdu010/.prog/anaconda/envs/lumpy
make
make install

samblaster：
https://github.com/GregoryFaust/samblaster
https://github.com/GregoryFaust/samblaster/archive/refs/heads/master.zip
cd samblaster
make
cp samblaster /public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/

sambamba：
https://github.com/lomereiter/sambamba
https://github.com/biod/sambamba/archive/refs/tags/v0.8.1.tar.gz
https://github.com/biod/sambamba/releases/download/v0.8.1/sambamba-0.8.1-linux-amd64-static.gz
mv sambamba-0.8.1-linux-amd64-static sambamba
cp sambamba /public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/

gawk:
https://ftp.gnu.org/gnu/gawk/
https://ftp.gnu.org/gnu/gawk/gawk-5.1.1.tar.gz
./configure --prefix=/public/home/users/fdu010/.prog/anaconda/envs/lumpy
make
make install

lumpy-sv：
https://github.com/arq5x/lumpy-sv
wget --no-check-certificate https://github.com/arq5x/lumpy-sv/archive/refs/heads/master.zip -O lumpy-sv.zip
https://github.com/samtools/htslib/archive/d2d9c76ade2df2b63b9cf79ae8decda1dfadc042.zip
下载到lib/htslib/


/usr/lib64/libz.so
export ZLIB_PATH="/usr/lib64"; 


#conda install -c bioconda lumpy-sv


actv lumpy
make
cp bin/* /public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/
cp -rf scripts /public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/
ln -s scripts/*

bamkit:
https://github.com/hall-lab/bamkit/archive/refs/heads/master.zip
cd /public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit/
```




# 安装运行出错解决方案：

```
错误0：
安装时报错ModuleNotFoundError: No module named '_ctypes'的解决办法
1、执行如下命令：
yum install libffi-devel 
2、从"./configure ..."重新安装

错误1：
<string>:1: DeprecationWarning: the imp module is deprecated in favour of importlib; see the module's documentation for alternative uses

(lumpy) [fdu010@a110 lumpy]$ grep imp ~/.prog/anaconda/envs/lumpy/bin/lumpyexpress
    $PYTHON_TEST -c "import imp; imp.find_module('pysam')"
    $PYTHON_TEST -c "import imp; imp.find_module('numpy')"
修改成
    $PYTHON_TEST -c "import importlib; importlib.import_module('pysam')"
    $PYTHON_TEST -c "import importlib; importlib.import_module('numpy')"

AttributeError: module 'importlib' has no attribute 'find_module'
find_module
换成：
import_module


错误2：
File "/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit/bamlibs.py", line 61
print ','.join(lib_rg[lib])
改成：
print （','.join(lib_rg[lib])）
需要加括号


错误3：
File "/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit/bamlibs.py", line 80/128
    except IOError, e:
                  ^
SyntaxError: invalid syntax


    except IOError, e:
改成
    except IOError as e:

(lumpy) [fdu010@a110 lumpy]$ grep IOError /public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit/*py
/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit/bamcleanheader.py:    except IOError, e:
/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit/bamfilterrg.py:    except IOError, e:
/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit/bamfixflags.py:    except IOError, e:
/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit/bamgroupreads.py:    except IOError, e:
/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit/bamheadrg.py:    except IOError, e:
/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit/bamlibs.py:    except IOError as e:
/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit/bamtofastq.py:    except IOError, e:


错误4：
########
Checking for required python modules (/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/python)...
/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/lumpyexpress: line 15: -n: command not found
/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/lumpyexpress: line 16: [: ==: unary operator expected
Calculating insert distributions...

(lumpy) [fdu010@a110 lumpy]$ Library read groups: LIBRARY
Library read length: 150
Removed 251 outliers with isize >= 689
done
0
Running LUMPY...
########

在/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/lumpyexpress.config中添加 
HEXDUMP=`which hexdump || true`

===========================================================================
#!/bin/bash -e

# Find original directory of file, resovling symlinks
# http://stackoverflow.com/questions/59895/can-a-bash-script-tell-what-directory-its-stored-in/246128#246128
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
    DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
    SOURCE="$(readlink "$SOURCE")"
    [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"

# general
LUMPY_HOME=$DIR
HEXDUMP=`which hexdump || true`
LUMPY=`which lumpy || true`
SAMBLASTER=`which samblaster || true`
# either sambamba or samtools is required
SAMBAMBA=`which sambamba || true`
SAMTOOLS=`which samtools || true`

# python 2.7 or newer, must have pysam, numpy installed
PYTHON=`which python || true`

# python scripts
PAIREND_DISTRO=$LUMPY_HOME/scripts/pairend_distro.py
BAMGROUPREADS=$LUMPY_HOME/scripts/bamkit/bamgroupreads.py
BAMFILTERRG=$LUMPY_HOME/scripts/bamkit/bamfilterrg.py
BAMLIBS=$LUMPY_HOME/scripts/bamkit/bamlibs.py
===========================================================================



错误5：
(lumpy) [fdu010@a110 lumpy]$ Sourcing executables from /public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/lumpyexpress.config ...
Checking for required python modules (/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/python)...
[E::idx_find_and_load] Could not retrieve index file for 'ERR194147.bam'
Calculating insert distributions...
Library read groups: LIBRARY
Library read length: 150
[E::idx_find_and_load] Could not retrieve index file for '-'
Removed 64 outliers with isize >= 677
done
0
Running LUMPY...

samtools sort  ERR194147.bam  -o ERR194147.S.bam &
samtools index  ERR194147.S.bam &
lumpyexpress -B  ERR194147.S.bam -S ERR194147.splitters.bam -D ERR194147.discordants.bam -o ERR194147.vcf


测试：
#!/usr/bin/python
# -*- coding: UTF-8 -*-
 
str = "-";
seq = ("a", "b", "c");
print str.join( seq );

(lumpy) [fdu010@a110 lumpy]$ python 1.py
File "/public/home/users/fdu010/wangkai/lumpy/1.py", line 3
print str.join( seq );
          ^
SyntaxError: invalid syntax


错误6：
configure error C compiler cannot create executables错误解决
我们在编译软件的时候，是不是经常遇到下面的错误信息呢？


checking build system type... i686-pc-linux-gnu
checking host system type... i686-pc-linux-gnu
checking for gcc... gcc
checking for C compiler default output file name...
configure: error: C compiler cannot create executables
See `config.log' for more details.

这个错误产生的原因其实很简单： 由于我们在编译软件之前，进行了export操作，改变了CFLAGS和LIBS的值。这个时候只要讲这个值清空就可以了。
export LIBS=""
export DYLD_LIBRARY_PATH=""
export LD_LIBRARY_PATH=""

sh   export LIBS=""
sh   export CFLAGS=""

echo $LD_LIBRARY_PATH

*** LIBRARY_PATH shouldn't contain the current directory when

export LIBRARY_PATH="/public/software/compiler/intel/composer_xe_2015.2.164/compiler/lib/intel64:/public/software/compiler/intel/composer_xe_2015.2.164/mkl/lib/intel64"

错误7：
ModuleNotFoundError: No module named '_ctypes'



/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/python: can't open file '/public/home/users/fdu010/.prog/anaconda/envs/lumpy/bin/scripts/bamkit/bamlibs.py': [Errno 2] No such file or directory


实验验证数据：
SV detection in the NA12878 genome:
Although it is difficult to precisely measure the sensitivity and accuracy of SV predictions made from a real data set, it is also important to evaluate each tool’s performance when confronted with real data containing artifacts that are not easily captured by simulations (for example, PCR artifacts, chimeric molecules, reads from poorly assembled genomic regions, and so on). In this experiment we compared SV detection performance in the NA12878 individual by analyzing the Illumina Platinum Genomes dataset, which represents approximately 50X coverage of the NA12878 genome (European Nucleotide Archives; ERA172924). We additionally subsampled this dataset to approximately 5X coverage to assess SV detection in low coverage scenarios.

To estimate sensitivity and FDR, we compared predictions made by each tool to two truth sets: 
1) 3,376 validated, non-overlapping deletions from the 1000 Genomes Project [12] (Additional file 3); and 
2) 4,095 deletions that were detected by at least one tool in the 50X dataset, or that were reported by Mills et al. [12] (which used numerous SV detection tools), 
and that were validated by split-read mapping analysis of independent long-read sequencing data from PacBio or Illumina Moleculo platforms (Additional file 4).

Additional file 3:
https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r84#MOESM3
https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2014-15-6-r84/MediaObjects/13059_2013_3363_MOESM3_ESM.zip

Additional file 4:
https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-6-r84#MOESM4
https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2014-15-6-r84/MediaObjects/13059_2013_3363_MOESM4_ESM.zip




```

