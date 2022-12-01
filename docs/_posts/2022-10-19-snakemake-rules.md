---
layout: post
title: 二代测序应用分析时，一些有用的Snakemake Rule
tags: [snakemake, NGS]
---

```python
###################
# Genome files//以Ensembl数据库中小鼠GRCm38.p6参考基因组为例
genome_index_ht="/home/user/genomes/hisat2Index/GRCm38/genome"
genome_index_bt="/home/user/genomes/bowtie2Index/GRCm38/genome"
genome_index_st="/home/user/genomes/starIndex/GRCm38"

genome_gtf="/home/user/genomes/Mus_musculus.GRCm38.102.gtf"

blacklist="/home/user/genomes/mm10_blacklist_ENCFF547MET.bed"
effectiveGenomeSize="2730871774"  #human=3096649726
###################
# input files #
species="Mouse/Mus musculus"  # "Human/Homo sapiens" 或者其他物种
reads=["_R1","_R2"]           # 单端测序使用不设置此选项
ext=".fastq.gz"               # 或者.fq.gz
path_origin=""                # 原始数据存放路径
samples=[]                    # 样本名称，注意不带后缀名
###################
# pre-pipeline #
import os
import sys
from pathlib import Path

p = Path(path_origin)
ss=list(p.glob("**/*"+ext))
l=len(ext)+len(reads[0])

Path("rawdata").mkdir(parents=True, exist_ok=True)
pwd=Path.cwd()
for s in ss:
    f=s.name
    if f[0:-l] in samples:
        s2=pwd / "rawdata"/ f
        #s2.unlink(missing_ok=True)  # python 3.8支持
        if s2.exists():
            os.remove(s2)
        #s2.symlink_to(s)             # s2为生成的软链接文件
        os.symlink(s,s2)

print("#"*10)
print("Species is:")
print("All samples:", samples)
print("#"*10)

# 目录结构
# workdir
# |--rawdata
# |--trimdata
# |--fastQC
# |--chipQC          # ChIP-seq
# |--alignment
# |  |--filter
# |--deeptools       # ChIP-seq
# |  |--fragsize
# |  |--bigwig
# |  |--profile
# |--counts
# |  |--stringtie
# |  |--featurecounts

###################
# rules, 双端测序使用正常代码，单端测序使用注释代码
rule trimming:
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

# rule trimming:
#     input:
#         "rawdata/{sample}"+ext
#     output:
#         "trimdata/{sample}"+ext
#     params:
#         "trimdata/{sample}"+"_trimmed.fq.gz"
#     threads: 4
#     log:
#         "trimdata/logs/{sample}.log"
#     shell:
#         "trim_galore --stringency 3 -q 20 --length 20 -o trimdata --cores {threads} {input} 2> {log}; "
#         "mv {params} {output}"

rule fastQC:
    input:
        "trimdata/{sample}{read}"+ext
    output:
        "fastQC/{sample}{read}_fastqc.html"
    threads: 12
    log:
        "fastQC/logs/{sample}{read}.log"
    shell:
        "fastqc -o fastQC -t {threads} -q {input} 2> {log}"

# rule fastQC:
#     input:
#         "trimdata/{sample}"+ext
#     output:
#         "fastQC/{sample}_fastqc.html"
#     threads: 12
#     log:
#         "fastQC/logs/{sample}.log"
#     shell:
#         "fastqc -o fastQC -t {threads} -q {input} 2> {log}"

rule hisat2:
    input:
        r1="trimdata/{sample}"+reads[0]+ext,
        r2="trimdata/{sample}"+reads[1]+ext
    output:
        bam="alignment/{sample}.sorted.bam",
        summary="alignment/{sample}.summary"
    threads: 12
    params:
        index=genome_index_ht,
        splice="alignment/{sample}.splice.txt",
        met="alignment/{sample}.matrix.txt"
    log:
        "alignment/logs/{sample}.align.log"
    shell:
        "hisat2 -q --rna-strandness RF --dta-cufflinks -x {params.index} -1 {input.r1} -2 {input.r2} -p {threads} --novel-splicesite-outfile {params.splice} --summary-file {output.summary} --met-file {params.met} | "
        "samtools view -Sbh - | samtools sort -@ {threads} -o {output.bam} 2> {log}"

# rule star:
#     input:
#         r1="trimdata/{sample}"+reads[0]+ext,
#         r2="trimdata/{sample}"+reads[1]+ext
#     output:
#         "alignment/{sample}.sorted.bam"
#     threads: 12
#     params:
#         index=genome_index_st,
#         gtf=genome_gtf,
#         prefix="alignment/{sample}_"
#     log:
#         "alignment/logs/{sample}.align.log"
#     shell:
#     "STAR --runThreadN {threads} --genomeDir {params.index} --readFilesIn {input.r1} {input.r2} "
#     "--readFilesCommand zcat --sjdbGTFfile {params.gtf} --outFileNamePrefix {params.prefix} --outSAMtype BAM Unsorted --outStd BAM_Unsorted - | samtools sort -@ {threads} -o {output} 2> {log}"

# rule hisat2:
#     input:
#         "trimdata/{sample}"+ext
#     output:
#         bam="alignment/{sample}.sorted.bam",
#         summary="alignment/{sample}.summary"
#     threads: 12
#     params:
#         index=genome_index_ht,
#         splice="alignment/{sample}.splice.txt",
#         met="alignment/{sample}.matrix.txt"
#     log:
#         "alignment/logs/{sample}.align.log"
#     shell:
#         "hisat2 -q --rna-strandness RF --dta-cufflinks -x {params.index} -U {input} -p {threads} --novel-splicesite-outfile {params.splice} --summary-file {output.summary} --met-file {params.met} | "
#         "samtools view -Sbh - | samtools sort -@ {threads} -o {output.bam} 2> {log}"

rule bowtie2:
    input:
        r1="trimdata/{sample}"+reads[0]+ext,
        r2="trimdata/{sample}"+reads[1]+ext
    output:
        bam="alignment/{sample}.sorted.bam",
        summary="alignment/{sample}.summary"
    threads: 12
    params:
        index=genome_index_bt
    log:
        "alignment/logs/{sample}.align.log"
    shell:
        "bowtie2 -q -X 1000 -x {params.index} -1 {input.r1} -2 {input.r2} -p {threads} 2> {output.summary} "
        "| samtools view -Sbh - | samtools sort -@ {threads} -o {output.bam} 2> {log}"

# rule bowtie2:
#     input:
#         "trimdata/{sample}"+ext
#     output:
#         bam="alignment/{sample}.sorted.bam",
#         summary="alignment/{sample}.summary"
#     threads: 12
#     params:
#         index=genome_index_bt
#     log:
#         "alignment/logs/{sample}.align.log"
#     shell:
#         "bowtie2 -q -x {params.index} -U {input} -p {threads} 2> {output.summary} "
#         "| samtools view -Sbh - | samtools sort -@ {threads} -o {output.bam} 2> {log}"

rule filter4chip:
    input:
        "alignment/{sample}.sorted.bam"
    output:
        "alignment/filter/{sample}.filter.bam"
    threads: 12
    log:
        "alignment/filter/logs/{sample}.filter.log"
    shell:
        "samtools view -b -F 1804 -f 2 -q 30 -@ {threads} {input} > {output} 2> {log}"
    #可通过https://www.samformat.info/sam-format-flag查阅过滤规则

rule filter4rna:
    input:
        "alignment/{sample}.sorted.bam"
    output:
        "alignment/filter/{sample}.filter.bam"
    threads: 12
    log:
        "alignment/filter/logs/{sample}.filter.log"
    shell:
        "samtools view -b -F 780 -f 2 -q 30 -@ {threads} {input} > {output} 2> {log}"

rule bamindex:
    input:
        "alignment/{sample}.sorted.bam"
    output:
        "alignment/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule bamindex2:
    input:
        "alignment/filter/{sample}.filter.bam"
    output:
        "alignment/filter/{sample}.filter.bam.bai"
    shell:
        "samtools index {input}"

rule stringtie:
    input:
        "alignment/filter/{sample}.filter.bam"
    output:
        ogtf="counts/stringtie/{sample}/{sample}.gtf",
        tab="counts/stringtie/{sample}.tab"
    params:
        gtf=genome_gtf
    log:
        "counts/stringtie/logs/{sample}.log"
    shell:
        "stringtie -o {output.ogtf} -A {output.tab} -e -G {params.gtf} --rf {input} 2> {log}"

rule featurecounts:
    input:
        "alignment/filter/{sample}.filter.bam"
    output:
        "counts/featurecounts/{sample}.count.txt"
    params:
        gtf=genome_gtf
    log:
        "counts/featurecounts/logs/{sample}.count.log"
    shell:
        "featureCounts -s 2 -p -B -t exon -g gene_id -a {params.gtf} -o {output} {input} 2> {log}"

rule chipQC:
    input:
        "alignment/{sample}.sorted.bam"
    output:
        "chipQC/{sample}_qc.txt"
    log:
        "chipQC/logs/{sample}_qc.log"
    shell:{% raw %}"""
        bedtools bamtobed -bedpe -i {input} \
        | awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$6}}' \
        | grep -v 'chrM' | sort | uniq -c \
        | awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} END{{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}}' > {output}
        """{% endraw %}

rule fragmentszie:
    input:
        "alignment/{sample}.sorted.bam"
    output:
        png="deeptools/fragsize/{sample}_fragsize.png",
        hist="deeptools/fragsize/{sample}_fragsize.hist.txt"
    threads: 12
    log:
        "deeptools/logs/{sample}.fragsize.log"
    shell:
        "bamPEFragmentSize -b {input} -hist {output.png} --outRawFragmentLengths {output.hist} --maxFragmentLength 1000"

rule deduplicates:
    input:
        "alignment/filter/{sample}.filter.bam"
    output:
        bam="alignment/filter/{sample}.filter.dedup.bam",
        matrix="alignment/filter/{sample}.filter.dedup.txt"
    threads: 12
    log:
        "alignment/filter/logs/{sample}.dedup.log"
    shell:
        "picard MarkDuplicates INPUT={input} OUTPUT={output.bam} METRICS_FILE={output.matrix} "
        "REMOVE_DUPLICATES=true CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT 2> {log}"

rule bamcoverage:
    input:
        "alignment/filter/{sample}.filter.dedup.bam"
    output:
        "deeptools/bigwig/{sample}.filter.dedup.bw"
    params:
        black=blacklist
    threads: 12
    log:
        "deeptools/logs/{sample}.bamcoverage.log"
    shell:
        "bamCoverage -b {input} -o {output} -p {threads} "
        "--binSize 10 --smoothLength 30 --normalizeUsing RPKM --ignoreForNormalization chrX chrY chrM --extendReads 2> {log}"

rule profile:
    input:
        "deeptools/bigwig/{sample}.filter.dedup.bw"
    output:
        bed="deeptools/profile/{sample}_TSS.bed",
        matrix="deeptools/profile/{sample}_TSS.matrix.gz"
    params:
        ref=genome_gtf
    threads: 12
    log:
        "deeptools/logs/{sample}.profile.log"
    shell:
        "computeMatrix reference-point --referencePoint TSS -b 2500 -a 2500 -R {params.ref} "
        "-S {input} -o {output.matrix} --outFileSortedRegions {output.bed} --skipZeros -p {threads} 2> {log}"

rule plotprofile:
    input:
        "deeptools/profile/{sample}_TSS.matrix.gz"
    output:
        "deeptools/profile/{sample}_TSS.png"
    shell:
        "plotHeatmap -m {input} -out {output} --dpi 360"

rule plotfingerprint:
    input:
        expand("alignment/filter/{sample}.filter.dedup.bam",sample=samples)
    output:
        matrix="deeptools/plotfingerprint.txt",
        png="deeptools/plotfingerprint.png"
    threads: 12
    log:
        "deeptools/logs/plotfingerprint.log"
    shell:
        "plotFingerprint -b {input} -o {output.png} --outQualityMetrics {output.matrix} "
        "-p {threads} --ignoreDuplicates --extendReads 2> {log}"

rule correlation:
    input:
        expand("alignment/filter/{sample}.filter.dedup.bam",sample=samples)
    output:
        counts="deeptools/multibamsummary.counts.txt",
        npz="deeptools/multibamsummary.counts.npz"
    threads: 12
    params:
        black=blacklist
    log:
        "deeptools/logs/multibamsummary.log"
    shell:
        "multiBamSummary bins --bamfiles {input} --outRawCounts {output.counts} -o {output.npz} --extendReads --binSize 2000 -p {threads} 2> {log}"

rule plotcorrelation:
    input:
        "deeptools/multibamsummary.counts.npz"
    output:
        "deeptools/correlation.png"
    shell:
        "plotCorrelation -in {input} -o {output} --corMethod pearson --whatToPlot scatterplot --skipZeros --removeOutliers --log1p"

# rule all 模块
# RNA-seq使用
rule all:
    input:
        expand("fastQC/{sample}{read}_fastqc.html",sample=samples,read=reads),
        expand("alignment/filter/{sample}.filter.bam.bai",sample=samples),
        expand("counts/stringtie/{sample}/{sample}.gtf",sample=samples),
        expand("counts/featurecounts/{sample}.count.txt",sample=samples),

#ChIP-seq使用
rule all:
    input:
        expand("fastQC/{sample}{read}_fastqc.html",sample=samples,read=reads),
        expand("alignment/{sample}.sorted.bam.bai",sample=samples),
        expand("alignment/filter/{sample}.filter.bam.bai",sample=samples),
        #expand("chipQC/{sample}_qc.txt", sample=samples),
        expand("deeptools/bigwig/{sample}.filter.dedup.bw", sample=samples),
        expand("deeptools/fragsize/{sample}_fragsize.png", sample=samples),
        expand("deeptools/profile/{sample}_TSS.bed",sample=samples),
        expand("deeptools/profile/{sample}_TSS.png",sample=samples),
        "deeptools/plotfingerprint.png",
        "deeptools/multibamsummary.counts.txt",
        "deeptools/correlation.png",

```

其他的二代测序应用，如MeDIP-seq、RIP-seq等也可借鉴使用。
