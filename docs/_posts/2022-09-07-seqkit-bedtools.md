---
layout: post
title: 使用seqkit和bedtools从基因组序列文件中提取片段
tags: [sequence, toolkits]
thumbnail-img: /assets/img/mds/DNA_sequence.png
---

fasta和fastq是两种存储核酸序列（DNA、RNA）或者蛋白质序列（AA）的文件格式，是众多生物信息分析的基础文件。在[NGS Analysis](https://learn.gencore.bio.nyu.edu/)中有对两种格式的详细说明。目前处理两种文件的工具有[seqkit](https://bioinf.shenwei.me/seqkit/)、[seqtk](https://github.com/lh3/seqtk)，以及一些python包，如[pyfastx](https://pyfastx.readthedocs.io/en/latest/)等。

根据分析需求，这里比较了使用seqkit和[bedtoos](https://bedtools.readthedocs.io/en/latest/index.html)从基因组提取特定bed文件对应序列的异同。使用的是[ENCODE](https://www.encodeproject.org/data-standards/reference-sequences/) 的小鼠参考基因组序列和注释信息（mm10_no_alt_analysis_set_ENCODE和NCFF871VGR）

### 定制gtf|bed文件

- 基因（protein_coding）

```shell
awk '{if($3=="gene" && $12=="\"protein_coding\";"){print $0}}' gencode.gtf
```

- 转录本（protein_coding)

```shell
awk '{if($3=="transcript" && $18=="\"protein_coding\";"){print $0}}' gencode.gtf
```

- 启动子区，注意正负链的延伸方向不一致

```shell
awk 'BEGIN{OFS=FS="\t"}{if($3=="gene") {if($7=="+") {start=$4-500; end=$4+100;} else {if($7=="-") start=$5-100; end=$5+500;} if(start<0) start=0; print $1,$2,$3,start,end,$6,$7,$8,$9;}}' | grep "protein_coding" > gencode.gtf
```

### 使用seqkit

- 提取序列

1、 基因区域（TSS-500bp to TES+500bp）

```shell
seqkit subseq --gtf mm10_ENCFF871VGR_gene.gtf  -u 500 -d 500 mm10_no_alt_analysis_set_ENCODE.fasta | seqkit sort -l --quiet > mm10_ENCFF871VGR_gene_u500_d500.fa
```

2、 启动子区域（TSS-500bp to TSS+100bp）

```shell
seqkit subseq --gtf mm10_ENCFF871VGR_promoter.gtf mm10_no_alt_analysis_set_ENCODE.fasta |seqkit seq -u > mm10_ENCFF871VGR_promoter.seqkit.fa
```

- 序列信息统计

```shell
seqkit stats mm10_ENCFF871VGR_gene_u500_d500.fa > mm10_ENCFF871VGR_gene_u500_d500.fa.stats.txt #序列统计信息
seqkit watch --fields ReadLen mm10_ENCFF871VGR_gene_u500_d500.fa -O mm10_ENCFF871VGR_gene_u500_d500.fa.lendis.png #直方图展示统计信息
seqkit fx2tab mm10_ENCFF871VGR_gene_promoter.fa -l -g -n -i -H > mm10_ENCFF871VGR_promoter.seqkit.fa.tab.txt #每条序列信息
```

### 使用bedtools

- 提取序列

```shell
bedtools getfasta -s -fi mm10_no_alt_analysis_set_ENCODE.fasta -bed mm10_ENCFF871VGR_promoter.gtf > mm10_ENCFF871VGR_promoter.bedtools.fa
```

- 随机背景序列

```shell
bedtools shuffle -seed 1234 -i mm10_ENCFF871VGR_promoter.gtf -g mm10_no_alt_analysis_set_ENCODE.fasta.fai > mm10_ENCFF871VGR_promoter_shuffle.gtf  #根据参考基因组坐标信息（fai|genome文件），在基因组上随机选取位置片段，片段长度与bed文件的长度一致
bedtools getfasta -s -fi mm10_no_alt_analysis_set_ENCODE.fasta -bed mm10_ENCFF871VGR_promoter_shuffle.gtf > mm10_ENCFF871VGR_promoter_shuffle.fa #提取序列
```

**总结一下：**  
seqkit-subseq默认对负链取反向互补，而bedtools-getfasta默认不会取反向互补，反向互补需要加参数-s；  
当选取区域长度大于参考序列长度时，seqkit会保留区域内参考序列，bedtools则会丢弃；  
bedtools的速度快于seqkit。

### 关于GTF、bed和fastaidx文件的说明

一些常见的生信数据格式可参考[UCSC-format](http://genome.ucsc.edu/FAQ/FAQformat.html#format4)

- GTF (gene transfer format)文件

[Gencode文档](https://www.gencodegenes.org/pages/data_format.html)

[WIKI基因结构说明](https://en.wikipedia.org/wiki/Gene_structure)

GTF格式文件用来对基因组进行注释，由以tab键分割的9列组成，分别是：

1. seq_id: 序列的编号，一般为chr或者scanfold编号；
2. source: 注释的来源，一般为数据库或者注释的机构，如果未知，则用点“.”代替；
3. type: 注释信息的类型，比如gene, transcript, CDS, exon, start_codon, stop_codon, UTR等；
4. start: 该基因或转录本在参考序列上的起始位置；
5. end: 该基因或转录本在参考序列上的终止位置；
6. score: 得分，数字，是注释信息可能性的说明，可以是序列相似性比对时的E-values值或者基因预测是的P-values值，“.”表示为空；
7. strand: 该基因或转录本位于参考序列的正链(+)或负链(-)上;
8. phase: 仅对注释类型为“CDS”有效，表示起始编码的位置，有效值为0、1、2(对于编码蛋白质的CDS来说，本列指定下一个密码子开始的位置。每3个核苷酸翻译一个氨基酸，从0开始，CDS的起始位置，除以3，余数就是这个值，表示到达下一个密码子需要跳过的碱基个数。该编码区第一个密码子的位置，取值0,1,2。0表示该编码框的第一个密码子第一个碱基位于其5'末端；1表示该编码框的第一个密码子的第一个碱基位于该编码区外；2表示该编码框的第一个密码子的第一、二个碱基位于该编码区外；如果Feature为CDS时，必须指明具体值；
9. attributes: 一个包含众多属性的列表，格式为“标签=值”（tag=value），标签与值之间以空格分开，且每个特征之后都要有分号（包括最后一个特征），其内容必须包括gene_id和transcript_id。

另一种基因注释文件为，GFF (general feature format)，与GTF定义基本一致。

- Fasta index文件

fai文件由tab键分割的5列组成，如下

1. NAME: 序列的名称，只保留“>”后，第一个空白之前的内容；
2. LENGTH: 序列的长度，单位为bp；
3. OFFSET: 第一个碱基的偏移量，从0开始计数，换行符也统计进行；
4. LINEBASES: 除了最后一行外，其他代表序列的行的碱基数，单位为bp；
5. LINEWIDTH: 行宽，除了最后一行外，其他代表序列的行的长度，包括换行符，在windows系统中换行符为\r\n, 要在序列长度的基础上加2。

- BED (Browser Extensible Data)文件

BED文件是一种为序列定义基本序列特征的简单方法。每行至少包括chrom、chromStart和chromEnd 3列，另外还可添加额为的9列，总共3-12列组成，顺序固定。

必须的3列为：

1. chrom: 染色体号，例如，chr1；
2. chromStart: feature在染色体上起始位置。从0开始算，染色体上第一个碱基位置标记为0，即0-base规则；
3. chromEnd: feature在染色体上终止位置。染色体上前100个碱基片段的位置位置标记为：chromStart=0, chromEnd=100；

可选的9列为：

1. name: BED行名，在基因组浏览器左边显示；
2. score: 在基因组浏览器中显示的灰度设定，值介于0-1000；
3. strand: 正负链标记，“.”表示不确定链性，“+” or “-”表示正链或负链；
4. thickStart: feature起始位置；
5. thickEnd: feature编码终止位置；
6. itemRgb: R,G,B (e.g. 255,0,0)值，当itemRgb 设置为 “On”，BED的行会显示颜色；
7. blockCount: blocks (exons)数目；
8. blockSizes: blocks (exons)大小列表，逗号分隔，对应于blockCount；
9. blockStarts: blocks (exons)起始位置列表，逗号分隔，对应于blockCount，与chromStart的一个相对位置。

---
2022年9月5日·四川泸定县6.8级地震 🙏🙏🙏
