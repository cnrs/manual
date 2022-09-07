https://github.com/Illumina/SMNCopyNumberCaller

适合whole-genome sequencing (WGS)和特定基因
SMNCopyNumberCaller is a tool to call the copy number of full-length SMN1, full-length SMN2, as well as SMN2Δ7–8 (SMN2 with a deletion of Exon7-8) from a whole-genome sequencing (WGS) BAM file.
This caller works with standard WGS sequencing depth (>=30X).


#/public/home/users/fdu010/anaconda3/envs/smncopynumbercaller

https://github.com/Illumina/SMNCopyNumberCaller/blob/master/requirements.txt

/public/home/users/fdu010/.prog/anaconda/envs/smncopynumbercaller
conda create -n smncopynumbercaller
conda activate smncopynumbercaller
#conda install pysam=0.15.3
conda install python=3.7 reportlab numpy=1.16 scipy=1.2 statsmodels=0.9 pysam=0.15.3
#conda install reportlab numpy=1.16 scipy=1.2 statsmodels=0.9 pysam=0.15.3

reportlab
numpy >=1.16
scipy >=1.2
pysam >=0.15.3
statsmodels >=0.9


https://github.com/Illumina/SMNCopyNumberCaller/archive/refs/heads/master.zip


示例：
smn_caller.py --manifest MANIFEST_FILE --genome [19/37/38] --prefix OUTPUT_FILE_PREFIX --outDir OUTPUT_DIRECTORY --threads NUMBER_THREADS

usage: smn_caller.py [-h] --manifest MANIFEST --genome GENOME --outDir OUTDIR --prefix PREFIX [--threads THREADS] [--reference REFERENCE]
                     [--countFilePath COUNTFILEPATH]

Call Copy number of full-length SMN1, full-length SMN2 and SMN* (Exon7-8 deletion) from a WGS bam file.

optional arguments:
  -h, --help            show this help message and exit
  --manifest MANIFEST   Manifest listing absolute paths to input BAM/CRAM files
  --genome GENOME       Reference genome, select from 19, 37, or 38
  --outDir OUTDIR       Output directory
  --prefix PREFIX       Prefix to output file
  --threads THREADS     Number of threads to use. Default is 1
  --reference REFERENCE
                        Optional path to reference fasta file for CRAM
  --countFilePath COUNTFILEPATH
                        Optional path to count files




/public/home/users/fdu010/wangkai/CNV/smncopynumbercaller

数据路径已经上传至上海超算
/public/home/users/fdu010/wangkai/WGS_data
样本说明及报告见report_otherinfor.zip及附件；主要测试动态突变及CNV分析模块



ll *_R1_001.fastq.gz | sed -e 's/_R1_001.fastq.gz//g'| awk '{print "bwa mem -t 24 /public/home/users/fdu010/.db/ucsc/hg19.fa " $9 "_R1_001.fastq.gz " $9 "_R2_001.fastq.gz > " $9 ".sam"}'

bwa mem -t 24 /public/home/users/fdu010/.db/ucsc/hg19.fa KYNGQX2101608801_S220_L004_R1_001.fastq.gz KYNGQX2101608801_S220_L004_R2_001.fastq.gz > KYNGQX2101608801_S220_L004.sam
bwa mem -t 24 /public/home/users/fdu010/.db/ucsc/hg19.fa KYNGQX2101608901_S254_L004_R1_001.fastq.gz KYNGQX2101608901_S254_L004_R2_001.fastq.gz > KYNGQX2101608901_S254_L004.sam
bwa mem -t 24 /public/home/users/fdu010/.db/ucsc/hg19.fa KYNGQX2101609001_S255_L004_R1_001.fastq.gz KYNGQX2101609001_S255_L004_R2_001.fastq.gz > KYNGQX2101609001_S255_L004.sam
bwa mem -t 24 /public/home/users/fdu010/.db/ucsc/hg19.fa KYNGQX2101609101_S136_L002_R1_001.fastq.gz KYNGQX2101609101_S136_L002_R2_001.fastq.gz > KYNGQX2101609101_S136_L002.sam
bwa mem -t 24 /public/home/users/fdu010/.db/ucsc/hg19.fa KYNGQX2101609101_S65_L001_R1_001.fastq.gz KYNGQX2101609101_S65_L001_R2_001.fastq.gz > KYNGQX2101609101_S65_L001.sam
bwa mem -t 24 /public/home/users/fdu010/.db/ucsc/hg19.fa KYNGQX2101609201_S137_L002_R1_001.fastq.gz KYNGQX2101609201_S137_L002_R2_001.fastq.gz > KYNGQX2101609201_S137_L002.sam
bwa mem -t 24 /public/home/users/fdu010/.db/ucsc/hg19.fa KYNGQX2101609301_S138_L002_R1_001.fastq.gz KYNGQX2101609301_S138_L002_R2_001.fastq.gz > KYNGQX2101609301_S138_L002.sam
bwa mem -t 24 /public/home/users/fdu010/.db/ucsc/hg19.fa KYNGQX2101609301_S66_L001_R1_001.fastq.gz KYNGQX2101609301_S66_L001_R2_001.fastq.gz > KYNGQX2101609301_S66_L001.sam


ll *.sam | awk '{print "nohup samtools view -bS " $9 " -o " $9 ".tmp > " $9 ".t.nohup 2>&1 & sleep 1s"}'
ll *.tmp | awk '{print "nohup samtools sort -o " $9 ".bam " $9  " > " $9 ".s.tmp.nohup 2>&1 & sleep 1s"}' | sed -e 's/sam.tmp.bam/bam/g'
ll *.bam | awk '{print "samtools index " $9 }'



(smncopynumbercaller) [wangk@cs data]$ 
cd  /public/home/users/fdu010/wangkai/CNV/smncopynumbercaller
ll *.bam | sed -e 's/.bam//g' | awk '{print "echo /public/home/users/fdu010/wangkai/WGS_data/rawdata/" $9 ".bam > " $9 ".txt"}'
ll *.bam | sed -e 's/.bam//g' | awk '{print "echo /public/home/users/fdu010/wangkai/CNV/smncopynumbercaller/" $9 ".bam > " $9 ".txt"}'


cd  /public/home/users/fdu010/wangkai/CNV/smncopynumbercaller



ll *.bam | awk '{print "samtools index " $9 " & sleep 1s"}'
samtools index CA0933B.bam & sleep 1s
samtools index CA1185B.bam & sleep 1s



#/public/home/users/fdu010/.prog/anaconda/envs/smncopynumbercaller
(smncopynumbercaller) [wangk@cs data]$ 

ll *.bam | sed -e 's/.bam//g' | awk '{print "smn_caller.py --manifest " $9 ".txt --genome 19 --prefix " $9 " --outDir " $9 " --threads 24"}'


smn_caller.py --manifest KYNGQX2101608801_S220_L004.txt --genome 19 --prefix S220_L004 --outDir S220_L004 --threads 24
smn_caller.py --manifest KYNGQX2101608901_S254_L004.txt --genome 19 --prefix S254_L004 --outDir S254_L004 --threads 24
smn_caller.py --manifest KYNGQX2101609001_S255_L004.txt --genome 19 --prefix S255_L004 --outDir S255_L004 --threads 24
smn_caller.py --manifest KYNGQX2101609101_S136_L002.txt --genome 19 --prefix S136_L002 --outDir S136_L002 --threads 24
smn_caller.py --manifest KYNGQX2101609101_S65_L001.txt --genome 19 --prefix S65_L001 --outDir S65_L001 --threads 24
smn_caller.py --manifest KYNGQX2101609201_S137_L002.txt --genome 19 --prefix S137_L002 --outDir S137_L002 --threads 24
smn_caller.py --manifest KYNGQX2101609301_S138_L002.txt --genome 19 --prefix S138_L002 --outDir S138_L002 --threads 24
smn_caller.py --manifest KYNGQX2101609301_S66_L001.txt --genome 19 --prefix S66_L001 --outDir S66_L001 --threads 24

  


可视化：
smn_charts.py -s SMN_JSON_FILE -o OUTPUT_DIRECTORY

(smncopynumbercaller) [fdu010@a110 smncopynumbercaller]$  ll |grep " S" | awk '{print "./smn_charts.py -s " $9 "/" $9 ".json -o " $9}'
./smn_charts.py -s S136_L002/S136_L002.json -o S136_L002
./smn_charts.py -s S137_L002/S137_L002.json -o S137_L002
./smn_charts.py -s S138_L002/S138_L002.json -o S138_L002
./smn_charts.py -s S220_L004/S220_L004.json -o S220_L004
./smn_charts.py -s S254_L004/S254_L004.json -o S254_L004
./smn_charts.py -s S255_L004/S255_L004.json -o S255_L004
./smn_charts.py -s S65_L001/S65_L001.json -o S65_L001
./smn_charts.py -s S66_L001/S66_L001.json -o S66_L001




需要图形支持，需要在本地服务器上测试下。


出现这个问题是由于没有加解释器（加上就可以了）
#!/usr/bin/env python3
  


例子：
A .pdf file is produced for each sample in OUTPUT_DIRECTORY, containing four plots.?
https://github.com/Illumina/SMNCopyNumberCaller/blob/master/charts/data/smn_HG03458.pdf

