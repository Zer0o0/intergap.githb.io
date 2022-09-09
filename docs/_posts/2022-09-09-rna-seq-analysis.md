---
layout: post
title: 一次Hisat2+featureCounts+DESeq的RNA-seq分析记录
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

1. 质控
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)、[fastp](https://github.com/OpenGene/fastp)、[RSeQC](http://rseqc.sourceforge.net/)等
2. 预处理reads
[Cutadapt](https://cutadapt.readthedocs.io/en/stable/)、[Trim Galore (Cutadapt+fastQC)](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)、[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)等
3. reads回贴
[HISAT2](http://daehwankimlab.github.io/hisat2/)、[STAR](https://github.com/alexdobin/STAR)
4. 计算基因表达量
[featureCounts](http://subread.sourceforge.net/featureCounts.html)、[StringTie](https://ccb.jhu.edu/software/stringtie/)、[Salmon](https://salmon.readthedocs.io/en/latest/salmon.html)、[Cufflinks](https://github.com/cole-trapnell-lab/cufflinks)
Salmon可以选择跳过reads比对回参考基因组步骤。
5. 差异分析
[Cufflinks-cuffdiff](https://github.com/cole-trapnell-lab/cufflinks)、[DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)、[limma](https://bioconductor.org/packages/release/bioc/html/limma.html)、[edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
6. 功能富集分析
[DAVID](https://david.ncifcrf.gov/home.jsp)、[GSEA](https://www.gsea-msigdb.org/gsea/index.jsp)、[Enrichr](https://maayanlab.cloud/Enrichr/)、[panther](http://pantherdb.org/)、[GREAT](http://bejerano.stanford.edu/great/public/html/)、[clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html)等

[Snakemake](https://snakemake.readthedocs.io/en/stable/)是基于Python语言实现的一个工作流管理系统，用于创建可重复和可扩展的数据分析流程的工具。Snakemake 可以根据分析所需软件描述，自动部署到任何执行环境。此外，可以无缝扩展到服务器、集群、网格和云环境，无需修改工作流定义。

### reads数据处理

如下为基于Snakemake定义的RNA-seq分析流程，使用了Hisat2+featureCounts+DESeq2的组合工具，针对双端测序（pair-end）使用：
```
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
    threads: 24
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
    threads: 16
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
    log:
        "counts/logs/{sample}.log"
    shell:
        "stringtie -o {output.gtf} -A {output.tab} -e -G {params.ann} --rf {input} 2> {log}"

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

### 差异表达分析

使用R包DESeq2进行差异分析，
```
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
  colnames(d) <- c('geneid','length',n)
  N <- sum(d[,n])
  d[,n2] <- round(10^9*d[,n]/(d[,'length']*N),2)
  # d <- mutate(d,fpkm=round(10^9*d[,n]/(d[,'length']*N),2))#计算fpkm
  # colnames(d) <- c('geneid','length',n,n2)
  #colnames 无法改名？？，
  if (length(counts_d)==0){
    counts_d <- d
  }else{
    counts_d <- merge(counts_d,d,by=c('geneid','length'))
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
cts <- counts_d%>%select(1,3,5,7,9) # 基因表达矩阵，基因名作为行名，样本名作为列名，至少有2个重复
rownames(cts) <- cts$geneid
cts$geneid <- NULL
design <- data.frame(
  sample = colnames(cts),
  condition = c("ctrl", "treat", "ctrl", "treat"),
  batch = c("I", "I", "II", "II")
) #设置批次，去除批次效应
coldata <- data.frame(row.names = design$sample, condition = factor(design$condition), batch = factor(design$batch))
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
res$geneid <- rownames(res)

# 定义显著性阈值，筛选差异表达分析
p <- log10(0.05) * (-1)
fc <- log2(2)
degs <- res %>%
  select(geneid, log2FoldChange, pvalue, padj) %>%
  mutate(logp = log10(pvalue) * (-1)) %>%
  mutate(type = if_else(logp > p, if_else(log2FoldChange > fc, "Up", if_else(log2FoldChange < (-fc), "Down", "N.s")), "N.s"))
table(degs$type)
write.table(degs, "degs_of_rnaseq.txt", row.names = F, sep='\t')

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

---
缅怀毛主席逝世46周年·1976年9月9日至2022年9月9日 💮🕯🕯🕯
