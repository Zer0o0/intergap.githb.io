---
layout: post
title: 一次Hisat2+featureCounts+DESeq的RNA-seq分析流程
tags: [RNA-seq, toolkits, pipline]
cover-img: /assets/img/mds/sanger_seq.png
thumbnail-img: /assets/img/mds/volcano_plot.png
---

在过去的十年中，RNA 测序 (RNA-seq) 成为了全转录组分析差异基因表达和mRNA差异剪接的重要工具。随着下一代测序技术的发展，RNA-seq方法已经被广泛用于研究RNA生物学的许多不同领域，包括单细胞基因表达、翻译和RNA结构等。另外，例如空间转录组学（spatialomics）等一些新兴的应用也正在探索之中。结合新的长读长和直接RNA-seq技术以及用于数据分析的更好的计算工具的开发，RNA-seq的创新将有助于更全面地了解RNA生物学，包括转录发生的时间和位置，以及RNA功能调控问题等。

RNA-seq的常规分析流程主要包括：1）建库测序，获取reads；2）对reads进行质控；3）去除接头和低质量碱基；4）将reads比对回参考基因组；5）统计reads数表征基因表达水平；6）基因差异表达分析或者联合多组学分析；7）基因功能富集分析等。常用工具如下：

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
