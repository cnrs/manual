适用范围：
1. 全基因组/全外显子组数据/panel；
2. 有实验样本，最好有对照，但也可以自己生成参考对照；
3. 其他要求最好和DECoN一样；


相关文档：
https://cnvkit.readthedocs.io/en/stable/
https://cnvkit.readthedocs.io/en/stable/quickstart.html
https://cnvkit.readthedocs.io/en/stable/pipeline.html
https://www.yuque.com/biotrainee/wes/bmgldi
https://www.jianshu.com/p/140f885a36cf
https://www.sohu.com/a/231637430_99971433
https://blog.csdn.net/weixin_43569478/article/details/108079649
http://www.360doc.com/content/19/1224/13/68068867_881784738.shtml

#https://mp.weixin.qq.com/s?__biz=MzIwODA1MzI4Mg==&mid=2456012609&idx=1&sn=dd4a2d6525aee5f196e89470da609d6f&chksm=809f920cb7e81b1a5fa1532f070edb0574543f73d6a5cc614affb050d63c4447f1e6dc62bd64&scene=21#wechat_redirect

#https://cnvkit.readthedocs.io/en/stable/pipeline.html#access
#Other known unmappable, variable, or poorly sequenced regions can be excluded with the -x/--exclude option. This option can be used more than once to exclude several BED files listing different sets of regions. For example, regions of poor mappability have been precalculated by others and are available from the UCSC FTP Server (see here for hg19).
#An “access” file precomputed for the UCSC reference human genome build hg19, with some know low-mappability regions excluded, is included in the CNVkit source distribution under the data/ directory (data/access-5kb-mappable.hg19.bed).


准备数据：
# path for reference 
http://hgdownload.soe.ucsc.edu/downloads.html
GENOME=/public/home/users/fdu010/.db/ucsc/hg19.fa
CNVKIT=/public/home/users/fdu010/.prog/anaconda/envs/cnvkit/bin/
SAMPLE=CA1248B
actv cnvkit

注释文件：
#The access command computes the locations of the accessible sequence regions for a given reference genome based on these masked-out sequences, treating long spans of ‘N’ characters as the inaccessible regions and outputting the coordinates of the regions between them.
#http://genome.ucsc.edu/goldenpath/gbdDescriptionsOld.html#RefFlat
#http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz


grep 'NM_' hg19.refFlat.txt | awk '{print $3 "\t" $5 "\t" $6 "\t" $1}' | sort -u > ucsc_hg19_baits.bed 


下载安装：
https://github.com/etal/cnvkit
https://github.com/etal/cnvkit-examples
conda create -n cnvkit -c bioconda cnvkit
#actv  cnvkit

一步法：
全外显子组：
cnvkit.py batch mapping/DYDFH09117.bam mapping/L1*.bam \
--normal mapping/CA21*.bam mapping/PD0595B.bam \
-p 0 --method amplicon \
--drop-low-coverage \
--targets hg19.exon.bed \
--fasta hg19.fa --access access-5k-mappable.hg19.bed \
--output-reference reference.cnn \
--output-dir cnvkit_results \
--diagram --scatter

#如果targets格式文件如下：则不需要注释文件，否则需要hg19.refFlat.txt (--annotate hg19.refFlat.txt )
chr1        1508981 1509154 SSU72
chr1        2407978 2408183 PLCH2
chr1        2409866 2410095 PLCH2

对于无对照样本的情况（-n）：
cnvkit.py batch *Tumor.bam -n -t my_baits.bed -f hg19.fasta  --access data/access-5kb-mappable.hg19.bed  --output-reference my_flat_reference.cnn -d cnvkit_results

对于已经生成了reference file的情况：
cnvkit.py batch *Tumor.bam -r my_reference.cnn -p 0 --scatter --diagram -d cnvkit_results

An “access” file precomputed for the UCSC reference human genome build hg19, with some know low-mappability regions excluded, is included in the CNVkit source distribution under the data/ directory (data/access-5kb-mappable.hg19.bed).
https://github.com/etal/cnvkit/tree/master/data


全基因组：
cnvkit.py batch 700_bwa.sam.bam --annotate ucsc-human-refflat2.txt 
--normal 699_bwa.sam.bam \
--method wgs \
--fasta hg38.fa \
--output-reference my_flat_reference.cnn \
-d  cnvkit_results \

#============================================================
cd cnvkit_results

cnvkit.py export seg *bqsr.cns -o gistic.segments
sed 's/_bqsr//' gistic.segments

awk '{print FILENAME"\t"$0}' *bqsr.cns  | grep -v chromosome |sed 's/_bqsr.cns//g' |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$8"\t"$6}' >final.seg
#============================================================

分步法：
第一步：生成基线 reference.cnn 文件
# Reusing targets and antitargets to build a new reference, but no analysis
cnvkit.py batch --normal *Normal.bam --fasta hg19.fa --access access-5k-mappable.hg19.bed --targets target.bed --antitargets antitarget.bed --output-reference normal.reference.cnn

--normal：基线样本BAM文件，空分割
--fasta：参考基因组文件
--access：可比对区域文件（bin：5kb)
--targets ：panel 内 区域文件，可用 target 命令获得
--antitargets：off-panel 区域文件，可用 antitarget 命令获得

第二步：肿瘤样本 call cnv
# Reusing a reference for additional samples
cnvkit.py batch -p 0 -m amplicon --drop-low-coverage -r ${SAMPLE}.Reference.cnn -d  cnvkit_results ${SAMPLE}.bam

--reference ：上一步生成的 baseline 的 depth、log2 等文件
--output-dir ： 输出目录，默认是当前目录

${CNVKIT}/cnvkit.py batch ${SAMPLE}.bam -m amplicon -c --drop-low-coverage -p 0 -t hg19.exon.bed -f hg19.fa --access access-5k-mappable.hg19.bed --output-reference ${SAMPLE}_cnv_reference.cnn --output-dir cnvkit_results --diagram --scatter


=================================================================================
详解：
#这里可以输入多个tumor.bam和normal.bam
#batch --method amplicon

# target子命令用于处理目的区域的bed文件，可以添加对应的基因注释等信息，用法如下
cnvkit.py target hg19.exon.bed --annotate hg19.refFlat.txt -o my_targets.bed

#除了目的区域in-target，还需要计算off-set区域，也称之为antitarget,  in-target和off-target区域加起来就是基因组上所有可覆盖的区域。
#二测测序并不能达到100%的覆盖度，基因组上的高重复区域，端粒，着丝粒等区域就无法覆盖，
#所以cnvkit通过access子命令来计算基因组上可以覆盖到的区域，命令如下
cnvkit.py access hg19.fa -x excludes.bed -o access.hg19.bed

#access-5kb-mappable.hg19.bed
#access-5k-mappable.hg19.bed

#计算出可覆盖的区域之后，减去in-target区域， 就可以得到off-target区域，通过antitarget子命令来实现，代码如下
#cnvkit.py antitarget my_targets.bed -g access.hg19.bed -o my_antitargets.bed
cnvkit.py antitarget hg19.exon.bed -g access-5k-mappable.hg19.bed -o my_antitargets.bed

#coverage和autobin两个子命令都可以用来计算测序深度，以coverage为例，用法如下
#分别统计target和antitarget区域的测序深度信息，输出结果后缀为cnn,是cnvkit中定义的一种格式，专门用来存储测序深度信息。
cnvkit.py coverage ${SAMPLE}.bam hg19.exon.bed -o ${SAMPLE}.targetcoverage.cnn
cnvkit.py coverage ${SAMPLE}.bam my_antitargets.bed -o ${SAMPLE}.antitargetcoverage.cnn
#Provide the *.targetcoverage.cnn and *.antitargetcoverage.cnn files created by the coverage command:

#通过reference子命令来构建正常基因组的测序分布模型，采用对照样本的测序深度，校正GC含量等系统误差。有多个对照样本时，可以将所有的对照样本合并来创建，用法如下
#For target amplicon or whole-genome sequencing protocols, the “antitarget” BED and .cnn files can be omitted. 
#Otherwise, ensure the filename prefixes are the same for each pair of “.targetcoverage.cnn” and “.antitargetcoverage.cnn” files (as it’s done by default).
#cnvkit.py reference *coverage.cnn -f hg19.fa -o Reference.cnn
cnvkit.py reference ${SAMPLE}.targetcoverage.cnn -f hg19.fa -o ${SAMPLE}.Reference.cnn

#当没有对照样本时，软件可以模拟出一个正常的测序深度分布模型，用法如下
cnvkit.py reference -o ${SAMPLE}.Reference.cnn -f hg19.fa -t hg19.exon.bed -a my_antitargets.bed

#计算实验样本相对正常对照的log2 ratio, 输出结果后缀为cnr, 是cnvkit中定义的一种格式，专门用来存储log2ratio的信息。
cnvkit.py fix ${SAMPLE}.targetcoverage.cnn ${SAMPLE}.antitargetcoverage.cnn ${SAMPLE}.Reference.cnn -o ${SAMPLE}.cnr


#通过segment子命令进行segment的划分，输出结果后缀为cns, 是cnvkit中定义的一种格式，和SEG格式类似，用来存储CNV分析的结果。
cnvkit.py segment ${SAMPLE}.cnr -o ${SAMPLE}.cns

#接下来还可以通过call子命令，计算每个segment区域的绝对拷贝数，用法如下
cnvkit.py call ${SAMPLE}.cns -o ${SAMPLE}.call.cns


cnvkit.py diagram -s ${SAMPLE}.cns ${SAMPLE}.cnr
cnvkit.py scatter -s ${SAMPLE}.cns ${SAMPLE}.cnr
cnvkit.py heatmap *.cns
=================================================================================

