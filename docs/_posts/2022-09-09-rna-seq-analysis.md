---
layout: post
title: 一次Hisat2+featureCounts+DESeq2的RNA-seq分析记录
tags: [RNA-seq, toolkits, pipline]
cover-img: /assets/img/covers/sanger_seq.png
thumbnail-img: /assets/img/mds/volcano_plot.png
---

在过去的十年中，RNA 测序 (RNA-seq) 成为了全转录组分析差异基因表达和mRNA差异剪接的重要工具。随着下一代测序技术的发展，RNA-seq方法已经被广泛用于研究RNA生物学的许多不同领域，包括单细胞基因表达、翻译和RNA结构等。另外，例如空间转录组学（spatialomics）等一些新兴的应用也正在探索之中。结合新的长读长和直接RNA-seq技术以及用于数据分析的更好的计算工具的开发，RNA-seq的创新将有助于更全面地了解RNA生物学，包括转录发生的时间和位置，以及RNA功能调控问题等。

RNA-seq的常规分析流程主要包括：

- 建库测序，获取reads；
- 对reads进行质控；
- 去除接头和低质量碱基；
- 将reads比对回参考基因组；
- 统计reads数表征基因表达水平；
- 基因差异表达分析或者联合多组学分析；
- 基因功能富集分析等。

常用工具如下：

1. 质控：
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)、[fastp](https://github.com/OpenGene/fastp)、[RSeQC](http://rseqc.sourceforge.net/)等
2. 预处理reads：
[Cutadapt](https://cutadapt.readthedocs.io/en/stable/)、[Trim Galore (Cutadapt+fastQC)](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)、[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)等
3. reads回贴：
[HISAT2](http://daehwankimlab.github.io/hisat2/)、[STAR](https://github.com/alexdobin/STAR)
4. 计算基因表达量：
[featureCounts](http://subread.sourceforge.net/featureCounts.html)、[StringTie](https://ccb.jhu.edu/software/stringtie/)、[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html)、[Cufflinks](https://github.com/cole-trapnell-lab/cufflinks)，其中Salmon可以选择跳过reads比对回参考基因组步骤。
5. 差异分析：
[Cufflinks-cuffdiff](https://github.com/cole-trapnell-lab/cufflinks)、[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)、[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)、[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
6. 功能富集分析：
[DAVID](https://david.ncifcrf.gov/home.jsp)、[GSEA](https://www.gsea-msigdb.org/gsea/index.jsp)、[Enrichr](https://maayanlab.cloud/Enrichr/)、[panther](http://pantherdb.org/)、[clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)等

[Snakemake](https://snakemake.readthedocs.io/en/stable/)是基于Python语言实现的一个工作流管理系统，用于创建可重复和可扩展的数据分析流程的工具。Snakemake 可以根据分析所需软件描述，自动部署到任何执行环境。此外，可以无缝扩展到服务器、集群、网格和云环境，无需修改工作流定义。

**Ref**：  
[A survey of best practices for RNA-seq data analysis. 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0881-8)  
[Exaggerated false positives by popular differential expression methods when analyzing human population samples. 2022](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02648-4)

### reads数据处理

如下为基于Snakemake定义的RNA-seq分析流程，使用了Hisat2+featureCounts+DESeq2的组合工具，对双端测序（pair-end）数据使用，命令行运行 *snakemake -s snakefile --cores 12 --keep-going*：

```python
###################
# Genome files #
genome_index="/home/user/genomes/hisat2Index/hg38/genome"  #参考基因组
genome_gtf="/home/user/genomes/hg38.ncbiRefSeq.gtf"        #参考基因组注释文件
###################
# input files #
species="Human/Homo sapiens"  #
reads=["_R1","_R2"]           #或者_1,_2
ext=".fastq.gz"               #或者.fq.gz
path_origin=""                #测序原始数据存放路径
samples=[]                    #样本名，不加后缀
###################
# pre-pipeline #
import os
import sys
from pathlib import Path


p = Path(path_origin)
ss=list(p.glob("**/*"+ext))   # 注意，**目录链接无法找到文件
l=len(ext)+len(reads[0])

Path("rawdata").mkdir(parents=True, exist_ok=True)
pwd=Path.cwd()
for s in ss:
    f=s.name
    if f[0:-l] in samples:
        s2=pwd / "rawdata"/ f
        #s2.unlink(missing_ok=True)  #python 3.8支持
        if s2.exists():
            os.remove(s2)
        s2.symlink_to(s)  #注意s2为生成的软链接文件
        #os.symlink(s,s2)

print("#"*10)
print("Species is:")
print("All samples:", samples)
print("#"*10)
###################
rule all:
    input:
        expand("fastQC/{sample}{read}_fastqc.html",sample=samples,read=reads),
        expand("alignment/{sample}.filter.sorted.bam.bai",sample=samples),
        expand("counts/{sample}.tab",sample=samples),
        expand("featurecounts/{sample}.count.txt",sample=samples),
        #"counts/gene_count_matrix.csv"

rule trimgalore:
    input:
        r1="rawdata/{sample}"+reads[0]+ext,
        r2="rawdata/{sample}"+reads[1]+ext
    output:
        r1="trimdata/{sample}"+reads[0]+ext,
        r2="trimdata/{sample}"+reads[1]+ext
    params:
        tmp1="trimdata/{sample}"+reads[0]+"_val_1.fq.gz",
        tmp2="trimdata/{sample}"+reads[1]+"_val_2.fq.gz"
    threads: 4
    log:
        "trimdata/logs/{sample}.log"
    shell:
        "trim_galore --paired --stringency 3 -q 20 --length 20 -o trimdata --cores {threads} {input.r1} {input.r2} 2> {log}; "
        "mv {params.tmp1} {output.r1}; "
        "mv {params.tmp2} {output.r2}"

rule fastqc:
    input:
        "trimdata/{sample}{read}"+ext
    output:
        "fastQC/{sample}{read}_fastqc.html"
    log:
        "fastQC/logs/{sample}{read}.log"
    shell:
        "fastqc -o fastQC -t {threads} -q {input} 2> {log}"

rule hisat2:
    input:
        r1="trimdata/{sample}"+reads[0]+ext,
        r2="trimdata/{sample}"+reads[1]+ext
    output:
        bam="alignment/{sample}.bam",
        summary="alignment/{sample}.summary"
    threads: 12
    params:
        index=genome_index,
        splice="alignment/{sample}.splice.txt",
        met="alignment/{sample}.matrix.txt"
    log:
        "alignment/logs/{sample}.align.log"
    shell:
        "hisat2 -q --rna-strandness RF --dta-cufflinks -x {params.index} -1 {input.r1} -2 {input.r2} -p {threads} --novel-splicesite-outfile {params.splice} --summary-file {output.summary} --met-file {params.met} | "
        "samtools view -Sbh - > {output.bam} 2> {log}"
#链特异性，如果使用dUTP: --rna-strandness RF，默认非链特异性

rule filter:
    input:
        "alignment/{sample}.bam"
    output:
        "alignment/{sample}.filter.sorted.bam"
    threads: 12
    log:
        "alignment/logs/{sample}.filter.log"
    shell:
        "samtools view -b -F 780 -f 2 -q 30 -@ {threads} {input} | "
        "samtools sort -@ {threads} -o {output} 2> {log}"

rule bamindex:
    input:
        "alignment/{sample}.filter.sorted.bam"
    output:
        "alignment/{sample}.filter.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule stringtie:
    input:
        "alignment/{sample}.filter.sorted.bam"
    output:
        gtf="counts/{sample}/{sample}.gtf",
        tab="counts/{sample}.tab"
    params:
        ann=genome_gtf
    threads: 12
    log:
        "counts/logs/{sample}.log"
    shell:
        "stringtie -o {output.gtf} -A {output.tab} -e -G {params.ann} --rf -p {threads} {input} 2> {log}"

rule featurecounts:
    input:
        "alignment/{sample}.filter.sorted.bam"
    output:
        "featurecounts/{sample}.count.txt"
    params:
        gtf=genome_gtf
    log:
        "featurecounts/logs/{sample}.count.log"
    shell:
        "featureCounts -s 2 -p -B -t exon -g gene_id -a {params.gtf} -o {output} {input} 2> {log}"   
```

关于**StringTie表达值的转换**

StringTie计算得到的基因表达量已经标准化为FPKM和TPM，在使用DESeq2做差异分析时需要提供Count值，因此官方提供了python脚本[prepDE.py3](https://github.com/gpertea/stringtie/blob/master/prepDE.py3)来完成转换，（计算方法为coverage\*transcript_len/read_len，而非FPKM转换，注意计算差异时删除生成文件的最后一行），使用方法为 *python prepDE.py3 -i PATH_TO_STRINGTIE_OUTDIR -g gene_count_matrix.csv -t transcript_count_matrix.csv -l READS_LENGTH*。

### 差异表达分析

- 使用Cufflinks-cuffdiff进行差异分析，从bam文件开始，不需要定量表达水平。

```shell
cuffdiff [options]* transcripts.gtf \
sample1_replicate1.sam[,…,sample1_replicateM.sam] \
sample2_replicate1.sam[,…,sample2_replicateM.sam] … \
sampleN.sam_replicate1.sam[,…,sample2_replicateM.sam]

#[options]*可设置为：
--library-type fr-firststrand --min-reps-for-js-test M --labels sample1,sample2,...,sampleN -p 12
```

需要注意，*差异倍数计算方向是输入文件的后者-前者，如sample2-sample1，可以对照实际表达值校对*

- 使用R包DESeq2进行差异分析，

```R
library(tidyverse)

setwd('featurecounts')

#合并featurecounts结果
files <- list.files(pattern = '.count.txt$')
files
counts_d <- NULL
for (i in files){
  d <- read_tsv(i,skip = 1)
  n <- str_sub(i,end=-11)
  n2 <- paste0(n,'_fpkm')
  d <- select(d,c(1,6,7))
  colnames(d) <- c('gene_id','length',n)
  N <- sum(d[,n])
  d[,n2] <- round(10^9*d[,n]/(d[,'length']*N),2)
  # d <- mutate(d,fpkm=round(10^9*d[,n]/(d[,'length']*N),2))#计算fpkm
  # colnames(d) <- c('gene_id','length',n,n2)
  #colnames 无法改名？？，
  if (length(counts_d)==0){
    counts_d <- d
  }else{
    counts_d <- merge(counts_d,d,by=c('gene_id','length'))
  }
}
head(counts_d)
write.table(counts_d,'expr_of_rnaseq.txt',sep='\t',row.names = F)

#相关性分析
corr_d <- counts_d[,c(4,6,8,10)]  #取fpkm值做相关性分析，注意样本数量
x <- rowSums(corr_d)
corr_d <- corr_d[x>0,]

# library(ggcorrplot)
# corr <- round(cor(log2(corr_d+1)),2)
# ggcorrplot(
#   corr,
#   hc.order = TRUE,
#   # type = "lower",
#   outline.color = "white",
#   ggtheme = ggplot2::theme_gray,
#   colors = c("#6D9EC1", "white", "#E46726"),
#   lab=TRUE
# )
# ggsave('corrplot.jpg',width = 14,height = 14,dpi = 300,units = 'cm')

library(GGally)

ggpairs(log2(corr_d+1),columnLabels=c('ctrl_#1','treat_#1','ctrl_#2','treat_#2'))+theme_bw()
ggsave('pairsplot.jpg',width = 14,height = 14,dpi = 300,units = 'cm')

## 基因差异表达分析
library(DESeq2)
cts <- counts_d%>%
  select(1,3,5,7,9)%>%
  column_to_rownames('gene_id') # 基因表达矩阵，基因名作为行名，样本名作为列名，至少有2个重复

sampleinfo <- data.frame(
  sample = colnames(cts),
  condition = c("ctrl", "treat", "ctrl", "treat"),
  batch = c("I", "I", "II", "II")
) #设置批次，去除批次效应
coldata <- data.frame(row.names = sampleinfo$sample, 
                      condition = factor(sampleinfo$condition, levels = c("treat", "ctrl")), 
                      batch = factor(sampleinfo$batch))
dds <- DESeqDataSetFromMatrix(cts, colData = coldata, design = ~ batch + condition) # 去除批次效应~batch+condition
dds$condition <- factor(dds$condition, levels = c("treat", "ctrl")) # 设置比较方向，如treat-ctrl
dds <- DESeq(dds)

# PCA 分析及可视化
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, "condition")
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$batch)
plotPCA(vsd, "condition")  #观察去除批次效应后的效果
# pcaData <- plotPCA(vsd, ntop = 2000, intgroup = "condition", returnData = TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# pcaData %>%
#   ggplot(aes(PC1, PC2, color = condition)) +
#   geom_point(size = 5) +
#   geom_text(aes(label = name)) +
#   scale_color_manual(values = c("darkred", "darkblue")) +
#   xlab(paste0("PC1: ", percentVar[1], "% variance")) +
#   ylab(paste0("PC2: ", percentVar[2], "% variance")) +
#   theme_base()

# 差异表达分析
res <- results(dds, contrast = c("condition", "treat", "ctrl"))
res <- res[order(res$pvalue), ]
res <- na.omit(as.data.frame(res))
res$gene_id <- rownames(res)

# 定义显著性阈值，筛选差异表达分析
p <- log10(0.05) * (-1)
fc <- log2(2)
degs <- res %>%
  select(gene_id,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj) %>%
  mutate(logp = log10(pvalue) * (-1)) %>%
  mutate(type = if_else(logp > p, if_else(log2FoldChange > fc, "Up", if_else(log2FoldChange < (-fc), "Down", "N.s")), "N.s")) %>%
  arrange(type)

table(degs$type)
write.table(degs, "degs_DESeq2.txt", row.names = F, sep='\t')

# 火山图可视化
degs %>%
  filter(log2FoldChange < 10) %>%
  filter(logp < 50) %>%
  ggplot(aes(log2FoldChange, logp)) +
  geom_point(aes(colour = type), alpha = 0.8, size = 0.8) +
  scale_colour_manual(values = c("#3C5488FF", "grey", "#DC0000FF")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  # geom_vline(xintercept = c(-fc, fc), linetype = 2) +  #添加指示线
  # geom_hline(yintercept = p * (-1), linetype = 2) +
  scale_x_continuous(breaks = seq(-4, 4, 2), labels = seq(-4, 4, 2), limits = c(-5, 5)) +
  labs(x = quote(log[2] ~ FoldChange), y = quote(-log[10] ~ pvalue), colour = "") +
  theme_light() +
  theme(
    panel.border = element_rect(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(size = 14),
    legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm")
  )
ggsave("vocalno.jpg", width = 10, height = 8, dpi = 1200, units = "cm")

```

### 功能富集分析

- DAVID结果气泡图展示

```R
library(tidyverse)  #数据输入和变换
library(stringr)
library(openxlsx)  #打开Excel
library(ggthemes)  #图形主题
library(hrbrthemes)
library(ggpubr)
library(RColorBrewer)  #图形色彩
library(ggsci)
library(scico)
library(viridis)
library(viridisLite)
library(scales)
library(ggbreak)  #打断坐标轴
library(ggpointdensity)  #点密度图
library(pheatmap)  #热图

mypal <- pal_npg("nrc",alpha = 0.6)(4)  #颜色选择
show_col(mypal)

go <- read_tsv(RESULT_OF_DAVID)

p <- go %>%
  slice(1:10) %>% #选取前10个term
  mutate(logp=log10(PValue) * (-1),Percentage=`%`)%>%   
  separate(col=Term,into=c('ID','Term'),sep=':')%>%  #提取term，KEGG为~，GO为:
  separate(Term,into=c('first','rest'),sep = 1)%>%
  mutate(first_=lapply(first, str_to_upper))%>%
  unite("Term",first_,rest,sep='')%>%  #term的首个字母转变为大写
  arrange(Count)%>%
  mutate(Term=str_wrap(Term,70))%>%  #长度过长的term换行，每行至多70个字符
  mutate(Term2=factor(Term,levels = Term))%>%
  ggplot(aes(x = Count, y =Term2 , size = Percentage, colour = logp)) +
  geom_point() +
  scale_colour_gradient(low = "#788BDA", high = "#2E629F",breaks=c(2,4,6)) + ## 渐变色方向
  scale_size_continuous(range = c(2, 5),breaks=c(1,2,3)) +
  scale_x_continuous(breaks = c(3,6,9), labels = c(3,6,9),limits = c(2,10)) +
  labs(y = "", colour = quote(-log[10] ~ pvalue), size = "Percentage") + ## 注意下标的标识方法
  theme_light() +
  theme(
    panel.border = element_rect(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    axis.text = element_text(colour = "black", size = 14),
    axis.title = element_text(size = 14),
    legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm")
  )
p
ggsave(paste0("RESULT_OF_DAVID", ".jpg"), egg::set_panel_size(p, width=unit(2.6, "cm"), height=unit(7, "cm")),height = 10, width = 21, units = "cm", dpi = 1200)  #egg::set_panel_size定义绘图区大小，同一批次图片保持一致
```

- 带误差线的柱状图

```R
data <- read_tsv(filename)
data %>%
  pivot_longer(cols = rep1:rep3, names_to = "rep") %>%  # 数据转形为长型数据
  group_by(group, time) %>%
  summarise(mean = mean(value), sd = sd(value), .groups = "drop") %>%
  mutate(time2 = factor(time, levels = c("G0/G1", "S", "G2/M")), group2 = factor(group, levels = c("NC", "KD"))) %>%
  ggplot(aes(x = time2, y = mean, fill = group2)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black", width = 0.8, show.legend = F) +
  geom_signif(map_signif_level = TRUE, y_position = c(54, 84), xmin = c(0.8, 1.8), xmax = c(1.2, 2.2), tip_length = 0.02, annotations = "***", textsize = 3, size = 0.4, vjust = 0.2) +
  scale_fill_manual(values = mycol, name = "") +
  geom_errorbar(aes(ymin = mean, ymax = mean + sd), width = 0.3, position = position_dodge(0.8)) +
  ylim(0, 100) +
  theme_pubr() +
  labs(x = "", y = "Percent of Cell Number (%)", title = "") +
  theme(axis.title = element_text(size = 10))
ggsave("cell_cycle.jpg", width = 6, height = 6, units = "cm", dpi = 1200)
```

**prepDE.py3代码如下：**

```python
#!/usr/bin/env python3
import re, csv, sys, os, glob, warnings, itertools
from math import ceil
from optparse import OptionParser
from operator import itemgetter
from collections import defaultdict

parser=OptionParser(description='Generates two CSV files containing the count matrices for genes and transcripts, using the coverage values found in the output of `stringtie -e`')
parser.add_option('-i', '--input', '--in', default='.', help="a folder containing all sample sub-directories, or a text file with sample ID and path to its GTF file on each line [default: %default/]")
parser.add_option('-g', default='gene_count_matrix.csv', help="where to output the gene count matrix [default: %default")
parser.add_option('-t', default='transcript_count_matrix.csv', help="where to output the transcript count matrix [default: %default]")
parser.add_option('-l', '--length', default=75, type='int', help="the average read length [default: %default]")
parser.add_option('-p', '--pattern', default=".", help="a regular expression that selects the sample subdirectories")
parser.add_option('-c', '--cluster', action="store_true", help="whether to cluster genes that overlap with different gene IDs, ignoring ones with geneID pattern (see below)")
parser.add_option('-s', '--string', default="MSTRG", help="if a different prefix is used for geneIDs assigned by StringTie [default: %default]")
parser.add_option('-k', '--key', default="prepG", help="if clustering, what prefix to use for geneIDs assigned by this script [default: %default]")
parser.add_option('-v', action="store_true", help="enable verbose processing")

parser.add_option('--legend', default="legend.csv", help="if clustering, where to output the legend file mapping transcripts to assigned geneIDs [default: %default]")
(opts, args)=parser.parse_args()

samples = [] # List of tuples. If sample list, (first column, path). Else, (subdirectory name, path to gtf file in subdirectory)
if (os.path.isfile(opts.input)):
    # gtfList = True
    try:
        fin = open(opts.input, 'r')
        for line in fin:
            if line[0] != '#':
                lineLst = tuple(line.strip().split(None,2))
                if (len(lineLst) != 2):
                    print("Error: line should have a sample ID and a file path:\n%s" % (line.strip()))
                    exit(1)
                if lineLst[0] in samples:
                    print("Error: non-unique sample ID (%s)" % (lineLst[0]))
                    exit(1)
                if not os.path.isfile(lineLst[1]):
                    print("Error: GTF file not found (%s)" % (lineLst[1]))
                    exit(1)
                samples.append(lineLst)
    except IOError:
        print("Error: List of .gtf files, %s, doesn't exist" % (opts.input))
        exit(1)
else:
    # gtfList = False
    ## Check that opts.input directory exists
    if not os.path.isdir(opts.input):
      parser.print_help()
      print(" ")
      print("Error: sub-directory '%s' not found!" % (opts.input))
      sys.exit(1)
    #####
    ## Collect all samples file paths and if empty print help message and quit
    #####
    samples = []
    for i in next(os.walk(opts.input))[1]:
        if re.search(opts.pattern,i):
         for f in glob.iglob(os.path.join(opts.input,i,"*.gtf")):
            samples.append((i,f)) 

if len(samples) == 0:
  parser.print_help()
  print(" ")
  print("Error: no GTF files found under base directory %s !" % (opts.input))
  sys.exit(1)

RE_GENE_ID=re.compile('gene_id "([^"]+)"')
RE_GENE_NAME=re.compile('gene_name "([^"]+)"')
RE_TRANSCRIPT_ID=re.compile('transcript_id "([^"]+)"')
RE_COVERAGE=re.compile('cov "([\-\+\d\.]+)"')
RE_STRING=re.compile(re.escape(opts.string))

RE_GFILE=re.compile('\-G\s*(\S+)') #assume filepath without spaces..


#####
## Sort the sample names by the sample ID
#####

samples.sort()

#if opts.v:
#  print "Sample GTFs found:"
#  for s in samples:
#     print s[1]

#####
## Checks whether a given row is a transcript 
## other options: ex. exon, transcript, mRNA, 5'UTR
#####
def is_transcript(x):
  return len(x)>2 and x[2]=="transcript"

def getGeneID(s, ctg, tid):
  r=RE_GENE_ID.search(s)
  #if r: return r.group(1)
  rn=RE_GENE_NAME.search(s)
  #if rn: return ctg+'|'+rn.group(1)
  if r:
    if rn: 
      return r.group(1)+'|'+rn.group(1)
    else:
      return r.group(1)
  return tid

def getCov(s):
  r=RE_COVERAGE.search(s)
  if r:
    v=float(r.group(1))
    if v<0.0: v=0.0
    return v
  return 0.0

def is_overlap(x,y): #NEEDS TO BE INTS!
  return x[0]<=y[1] and y[0]<=x[1]


def t_overlap(t1, t2): #from badGenes: chromosome, strand, cluster, start, end, (e1start, e1end)...
    if t1[0] != t2[0] or t1[1] != t2[1] or t1[5]<t2[4]: return False
    for i in range(6, len(t1)):
        for j in range(6, len(t2)):
            if is_overlap(t1[i], t2[j]): return True
    return False

## Average Readlength
read_len=opts.length

## Variables/Matrices to store t/g_counts
t_count_matrix, g_count_matrix=[],[]

##Get ready for clustering, stuff is once for all samples##
geneIDs=defaultdict(lambda: str) #key=transcript, value=cluster/gene_id


## For each of the sorted sample paths
for s in samples:
    badGenes=[] #list of bad genes (just ones that aren't MSTRG)
    try:
        ## opts.input = parent directory of sample subdirectories
        ## s = sample currently iterating through
        ## os.path.join(opts.input,s,"*.gtf") path to current sample's GTF
        ## split = list of lists: [[chromosome, ...],...]

        #with open(glob.iglob(os.path.join(opts.input,s,"*.gtf")).next()) as f:
        #    split=[l.split('\t') for l in f.readlines()]
#        if not gtfList:
#            f = open(glob.iglob(os.path.join(opts.input,s[1],"*.gtf")).next())
#        else:
#            f = open(s[1])
        with open(s[1]) as f:
            split=[l.split('\t') for l in f.readlines()]

        ## i = numLine; v = corresponding i-th GTF row
        for i,v in enumerate(split):
            if is_transcript(v):
                try:
                  t_id=RE_TRANSCRIPT_ID.search(v[8]).group(1)
                  g_id=getGeneID(v[8], v[0], t_id)
                except:
                  print("Problem parsing file %s at line:\n:%s\n" % (s[1], v))
                  sys.exit(1)
                geneIDs[t_id]=g_id
                if not RE_STRING.match(g_id):
                    badGenes.append([v[0],v[6], t_id, g_id, min(int(v[3]),int(v[4])), max(int(v[3]),int(v[4]))]) #chromosome, strand, cluster/transcript id, start, end
                    j=i+1
                    while j<len(split) and split[j][2]=="exon":
                        badGenes[len(badGenes)-1].append((min(int(split[j][3]), int(split[j][4])), max(int(split[j][3]), int(split[j][4]))))
                        j+=1

    except StopIteration:
        warnings.warn("Didn't get a GTF in that directory. Looking in another...")

    else: #we found the "bad" genes!
        break

##THE CLUSTERING BEGINS!##
if opts.cluster and len(badGenes)>0:
    clusters=[] #lists of lists (could be sets) or something of transcripts
    badGenes.sort(key=itemgetter(3)) #sort by start coord...?
    i=0
    while i<len(badGenes): #rather un-pythonic
        temp_cluster=[badGenes[i]]

        k=0
        while k<len(temp_cluster):
            j=i+1
            while j<len(badGenes):
                if t_overlap(temp_cluster[k], badGenes[j]):
                    temp_cluster.append(badGenes[j])
                    del badGenes[j]
                else:
                    j+=1
            k+=1
        if len(temp_cluster)>1:
            clusters.append([t[2] for t in temp_cluster])
        i+=1

    print(len(clusters))

    for c in clusters:
        c.sort()

    clusters.sort(key=itemgetter(0))
    legend=[]
    for u,c in enumerate(clusters):
        my_ID=opts.key+str((u+1))
        legend.append(list(itertools.chain.from_iterable([[my_ID],c]))) #my_ID, clustered transcript IDs
        for t in c:
            geneIDs[t]=my_ID
##            geneIDs[t]="|".join(c) #duct-tape transcript IDs together, disregarding ref_gene_names and things like that

    with open(opts.legend, 'w') as l_file:
        my_writer=csv.writer(l_file)
        my_writer.writerows(legend)

geneDict=defaultdict(lambda: defaultdict(lambda: 0)) #key=gene/cluster, value=dictionary with key=sample, value=summed counts
t_dict=defaultdict(lambda: defaultdict(lambda: 0))
guidesFile='' # file given with -G for the 1st sample
for q, s in enumerate(samples):
    if opts.v:
       print(">processing sample %s from file %s" % s)
    lno=0
    try:
        #with open(glob.iglob(os.path.join(opts.input,s,"*.gtf")).next()) as f: #grabs first .gtf file it finds inside the sample subdirectory
#        if not gtfList:
#            f = open(glob.iglob(os.path.join(opts.input,s[1],"*.gtf")).next())
#        else:
        f = open(s[1])
        transcript_len=0
        
        for l in f:
            lno+=1
            if l.startswith('#'):
                if lno==1:
                    ei=l.find('-e')
                    if ei<0:
                       print("Error: sample file %s was not generated with -e option!" % ( s[1] ))
                       sys.exit(1)
                    gf=RE_GFILE.search(l)
                    if gf:
                       gfile=gf.group(1)
                       if guidesFile:
                          if gfile != guidesFile:
                             print("Warning: sample file %s generated with a different -G file (%s) than the first sample (%s)" % ( s[1], gfile, guidesFile ))
                       else:
                          guidesFile=gfile
                    else:
                       print("Error: sample %s was not processed with -G option!" % ( s[1] ))
                       sys.exit(1)
                continue
            v=l.split('\t')
            if v[2]=="transcript":
                if transcript_len>0:
##                        transcriptList.append((g_id, t_id, int(ceil(coverage*transcript_len/read_len))))
                    t_dict[t_id][s[0]] = int(ceil(coverage*transcript_len/read_len))
                t_id=RE_TRANSCRIPT_ID.search(v[len(v)-1]).group(1)
                #g_id=RE_GENE_ID.search(v[len(v)-1]).group(1)
                g_id=getGeneID(v[8], v[0], t_id)
                #coverage=float(RE_COVERAGE.search(v[len(v)-1]).group(1))
                coverage=getCov(v[8])
                transcript_len=0
            if v[2]=="exon":
                transcript_len+=int(v[4])-int(v[3])+1 #because end coordinates are inclusive in GTF

##            transcriptList.append((g_id, t_id, int(ceil(coverage*transcript_len/read_len))))
        t_dict[t_id][s[0]]=int(ceil(coverage*transcript_len/read_len))

    except StopIteration:
#        if not gtfList:
#            warnings.warn("No GTF file found in " + os.path.join(opts.input,s[1]))
#        else:
        warnings.warn("No GTF file found in " + s[1])


##        transcriptList.sort(key=lambda bla: bla[1]) #gene_id
    
    for i,v in t_dict.items():
##        print i,v
       try:
          geneDict[geneIDs[i]][s[0]]+=v[s[0]]
       except KeyError:
          print("Error: could not locate transcript %s entry for sample %s" % ( i, s[0] ))
          raise

if opts.v:
   print("..writing %s " % ( opts.t ))
with open(opts.t, 'w') as csvfile:
   my_writer = csv.DictWriter(csvfile, fieldnames = ["transcript_id"] + [x for x,y in samples])
   my_writer.writerow(dict((fn,fn) for fn in my_writer.fieldnames))
   for i in t_dict:
        t_dict[i]["transcript_id"] = i
        my_writer.writerow(t_dict[i])
if opts.v:
   print("..writing %s " % ( opts.g ))
with open(opts.g, 'w') as csvfile:
   my_writer = csv.DictWriter(csvfile, fieldnames = ["gene_id"] + [x for x,y in samples])
##    my_writer.writerow([""]+samples)
##    my_writer.writerows(geneDict)
   my_writer.writerow(dict((fn,fn) for fn in my_writer.fieldnames))
   for i in geneDict:
        geneDict[i]["gene_id"] = i #add gene_id to row
        my_writer.writerow(geneDict[i])
if opts.v:
   print("All done.")
```

---
缅怀毛主席逝世46周年·1976年9月9日至2022年9月9日 💮🕯🕯🕯
