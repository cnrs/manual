使用范围：
全基因组

CNVnator是一个基于read-depth方法检测CNV的软件，该软件
1. 无需健康人对照样本，以自身为对照进行拷贝数变异检测；
2. 该软件比较适合高深度全基因组的CNV检测；
3. 可以单个样本检测；

相关文档：
#https://mp.weixin.qq.com/s?__biz=MzU1ODU4OTkzNA==&mid=2247484342&idx=1&sn=ccc000c980df3ba8ab83678473eaf07b&chksm=fc257137cb52f821f27d04431422e04a3388a0b7d089c0a6f74a85806a2ebd07d345d5802468&scene=21#wechat_redirect
http://www.360doc.com/content/19/1224/13/68068867_881785146.shtml
https://cloud.tencent.com/developer/news/271069
https://www.jianshu.com/p/7eeae8a19038
https://www.codetd.com/article/6237849

下载安装：
https://github.com/abyzovlab/CNVnator
https://github.com/abyzovlab/CNVnator/blob/master/INSTALL
https://www.jianshu.com/p/7eeae8a19038

conda create -n cnvnator
conda activate cnvnator
conda install -c bioconda cnvnator
conda install -c conda-forge samtools

git clone https://github.com/abyzovlab/CNVnator.git

cp *.py /public/home/users/fdu010/.prog/anaconda/envs/cnvnator/bin
cp *.pl /public/home/users/fdu010/.prog/anaconda/envs/cnvnator/bin
cp -rf pytools /public/home/users/fdu010/.prog/anaconda/envs/cnvnator/bin

CNVNATOR=/public/home/users/fdu010/.prog/anaconda/envs/cnvnator/bin
GENOME=/public/home/users/fdu010/.db/ucsc/hg19.fa
WORK_DIR=/public/home/users/fdu010/wangkai
BINSIZE=30


流程演示：
1. Extracting read mapping from bam/sam files
# If option -chrom is not used all chromosomes from bam file will be extracted.
cnvnator -root out.root -tree CA1248B.sorted.MarkDuplicates.BQSR.bam
sleep 1s;
cnvnator -root out.root  -merge file1.root ...

ll *.bam | awk '{print "cnvnator -root " $9 ".root -tree " $9 " & sleep 1s"}' | sed -e 's/bam.root/root/g'
cnvnator -root CA2106B.root -tree CA2106B.bam & sleep 1s
cnvnator -root CA2115B.root -tree CA2115B.bam & sleep 1s


2. Generating a read depth histogram
#Files with individual chromosome sequences (.fa) are required and should reside in the current directory or in the directory specified by the -d option. Files should be named as: chr1.fa, chr2.fa, etc.
这一步不消耗内存，因此可以同时对所有染色体进行。
它也可以用于染色体的一个子集。需要具有单个染色体序列（.fa）的文件，这些文件应位于当前目录或-d选项指定的目录中。文件应命名为：chr1.fa、chr2.fa等。

#https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
#cp chromFa.tar.gz /home/wangk/test/wes/CNV/genome
cnvnator -root ${SAMPLE}.root -his ${BINSIZE} -d genome

#不能使用整个基因组文件，或者将多个单个的chro文件打包：
#tar -f hg19.tar.gz -cz *.fa &
#cnvnator -root ${SAMPLE}.root -his ${BINSIZE} -fasta hg19.tar.gz    
sleep 1s;

ll *.root | awk '{print "cnvnator -root " $9 " -his 30 -d genome & sleep 1s"}'
#ll *.root | awk '{print "cnvnator -root " $9 " -his 30 -fasta hg19.tar.gz & sleep 1s"}'

3. calculating statistics
cnvnator -root ${SAMPLE}.${BINSIZE}.root -stat ${BINSIZE}
sleep 1s;
cnvnator -root ${SAMPLE}.${BINSIZE}.root -eval ${BINSIZE}
sleep 1s;

ll *.root | awk '{print "cnvnator -root " $9 " -stat 30 & sleep 1s"}'
cnvnator -root CA2106B.1.root -stat 30 & sleep 1s
cnvnator -root CA2106B.root -stat 30 & sleep 1s
cnvnator -root CA2115B.root -stat 30 & sleep 1s
cnvnator -root CA2116B.root -stat 30 & sleep 1s
cnvnator -root DYDFH09117.root -stat 30 & sleep 1s
cnvnator -root L18D00192.root -stat 30 & sleep 1s
cnvnator -root L19D0017956.root -stat 30 & sleep 1s
cnvnator -root L19D0017964.root -stat 30 & sleep 1s
cnvnator -root PD0595B.root -stat 30 & sleep 1s


ll *.root | awk '{print "cnvnator -root " $9 " -eval 30 > " $9 ".eval & sleep 1s"}'
cnvnator -root CA2106B.root -eval 30 > CA2106B.root.eval & sleep 1s
cnvnator -root CA2115B.root -eval 30 > CA2115B.root.eval & sleep 1s
cnvnator -root CA2116B.root -eval 30 > CA2116B.root.eval & sleep 1s
cnvnator -root DYDFH09117.root -eval 30 > DYDFH09117.root.eval & sleep 1s
cnvnator -root L18D00192.root -eval 30 > L18D00192.root.eval & sleep 1s
cnvnator -root L19D0017956.root -eval 30 > L19D0017956.root.eval & sleep 1s
cnvnator -root L19D0017964.root -eval 30 > L19D0017964.root.eval & sleep 1s
cnvnator -root PD0595B.root -eval 30 > PD0595B.root.eval & sleep 1s



#（his是经过调试后的，需要根据你测序深度来决定，否则一直死循环找不到）
#指的是bin size，可以通过-eval参数进行筛选，也可以根据经验值进行确定，
#一般测序深度2-3x选取500，20-30x选取bin size大小100，100x选取30
#检验mean to segma的参数为-eval，就是cnvnator 你的root文件 -eval 你设置的bin值
#这个root文件里可以存储多个bin值，当达到4-5的时候，再往下进行第四步 partition 第五步call 这样得到的结果才是最准确的
#要保证mean to sigma的值在4-5之间

4. rd signal partitioning
耗时较长
cnvnator -root ${SAMPLE}.${BINSIZE}.root -partition ${BINSIZE}
ll *.root | awk '{print "cnvnator -root " $9 " -partition 30 & sleep 1s"}'
cnvnator -root CA2106B.root -partition 30 & sleep 1s
cnvnator -root CA2115B.root -partition 30 & sleep 1s
cnvnator -root CA2116B.root -partition 30 & sleep 1s
cnvnator -root DYDFH09117.root -partition 30 & sleep 1s
cnvnator -root L18D00192.root -partition 30 & sleep 1s
cnvnator -root L19D0017956.root -partition 30 & sleep 1s
cnvnator -root L19D0017964.root -partition 30 & sleep 1s
cnvnator -root PD0595B.root -partition 30 & sleep 1s


5. cnv calling
cnvnator -root ${SAMPLE}.root -call ${BINSIZE} > ${WORK_DIR}/${SAMPLE}.CNVNATOR.CNV
ll *.root | awk '{print "cnvnator -root " $9 " -call 30 > " $9 ".CNVNATOR.CNV & sleep 1s"}' | sed -e 's/root.CNVNATOR.CNV/CNVNATOR.CNV/g'
cnvnator -root CA2106B.root -call 30 > CA2106B.CNVNATOR.CNV & sleep 1s
cnvnator -root CA2115B.root -call 30 > CA2115B.CNVNATOR.CNV & sleep 1s
cnvnator -root CA2116B.root -call 30 > CA2116B.CNVNATOR.CNV & sleep 1s
cnvnator -root DYDFH09117.root -call 30 > DYDFH09117.CNVNATOR.CNV & sleep 1s
cnvnator -root L18D00192.root -call 30 > L18D00192.CNVNATOR.CNV & sleep 1s
cnvnator -root L19D0017956.root -call 30 > L19D0017956.CNVNATOR.CNV & sleep 1s
cnvnator -root L19D0017964.root -call 30 > L19D0017964.CNVNATOR.CNV & sleep 1s
cnvnator -root PD0595B.root -call 30 > PD0595B.CNVNATOR.CNV & sleep 1s


#对于原始的cnv call的输出，还可以通过软件自带的脚本转换为VCF格式，代码如下
cnvnator2VCF.pl ${WORK_DIR}/${SAMPLE}.CNVNATOR.CNV > ${WORK_DIR}/${SAMPLE}.CNVNATOR.CNV.vcf
#cnvnator2VCF.pl -prefix study1 -reference GRCh37 sample1.cnvnator.out /path/to/individual/fasta_files


#The output columns are as follows:
#CNV_type coordinates CNV_size normalized_RD e-val1 e-val2 e-val3 e-val4 q0

https://github.com/abyzovlab/CNVnator/issues/210

导入VCF数据
# Import SNP data
cnvnator -root out.root -vcf file.vcf.gz

cnvnator -root out.root -vcf file.vcf.gz -addchr
  # Options -addchr or -rmchr can be used to add or remove the "chr" prefix from 
  # chromosome names in vcf file to match chromosom names from bam file. 

# Import mask data
cnvnator -root out.root -mask mask.fa.gz
  OR
cnvnator -root out.root -mask mask.fa.gz -addchr
  
# Generate SNP histograms
cnvnator -root out.root -baf ${BINSIZE}

# Ploting
cnvnator -root out.root -view ${BINSIZE}
>1:1M-50M
>1:1M-50M baf

# List root file content
cnvnator -root out.root -ls

# Copy RD and SNP data to new root file
cnvnator -root out.root -cptrees new_file.root

# Ploting RD and BAF whole genome circular plots using python tool:
plotcircular.py out.root

