---
layout: post
title: 一次Bowtie2+deepTools+MACS2的ChIP-seq分析记录
tags: [ChIP-seq, toolkits, pipline]
cover-img: /assets/img/covers/sanger_seq.png
thumbnail-img: /assets/img/mds/chip_profile.png
---

蛋白质-DNA相互作用和表观遗传标记的全基因组图谱对于全面理解转录调控过程，破译构成各种生物过程的基因调控网络至关重要。核小体定位与 DNA 和组蛋白的动态修饰相结合，在基因调控中发挥关键作用，并且指导发育和分化。染色质状态的改变，如通过允许或阻止 DNA 结合蛋白的接触，或者修饰核小体表面以增强或阻碍效应蛋白复合物的募集，能够直接影响转录。最近的研究进展表明，染色质状态和转录之间的这种相互作用是动态的，并且比以前认识到的更为复杂，因此需要对多种细胞类型和阶段的表观基因组进行系统分析，以了解发育过程和疾病状态。  
染色质免疫沉淀后测序 (ChIP-seq, Chromatin immunoprecipitation followed by sequencing) 是一种对 DNA 结合蛋白、组蛋白修饰或核小体定位进行全基因组分析的技术。由于新一代测序技术的巨大进步，ChIP-seq 提供了更高的分辨率、更低的噪音和更大的覆盖范围。ChIP-seq 实验产生大量数据，有效的计算分析将极大的促进生物学机制的揭示。

ChIP-seq的一般分析过程包括：

- 建库测序，获取reads；
- 对reads进行质控；
- 去除接头和低质量碱基；
- 将reads比对回参考基因组；
- ChIP质量评估；
- 去掉重复的reads；
- 检测信号富集区域（peakcalling）及注释；
- IGV可视化peak；
- Motif分析；
- 基因功能富集分析等；
- 差异分析或者联合多组学分析。

常用工具有：

1. 质控：
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)、[fastp](https://github.com/OpenGene/fastp)、[RSeQC](http://rseqc.sourceforge.net/)等
2. 预处理reads：
[Cutadapt](https://cutadapt.readthedocs.io/en/stable/)、[Trim Galore (Cutadapt+fastQC)](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)、[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)等
3. reads回贴：
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)、[BWA](http://bio-bwa.sourceforge.net/)
4. 过滤重复reads：
[Picard](http://broadinstitute.github.io/picard/)、[Samtools](https://www.htslib.org/doc/samtools.html)
5. peakcalling和注释:
[MACS](https://github.com/macs3-project/MACS)、[Homer](http://homer.ucsd.edu/homer/)、[ChIP](https://bioconductor.org/packages/release/bioc/vignettes/ChIPseeker/inst/doc/ChIPseeker.html)等
6. Motif分析：
[MEME](https://meme-suite.org/meme/)、[Homer](http://homer.ucsd.edu/homer/)
7. 功能富集分析：
[GREAT](http://bejerano.stanford.edu/great/public/html/)等
8. 差异分析：
[DiffBind](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)
9. 可视化：
[deepTools](https://deeptools.readthedocs.io/en/develop/)、[IGV](https://igv.org/)、[bedtools](https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)、[WWebLogo](https://weblogo.threeplusone.com/)等

### reads数据处理

基于Snakemake定义的RNA-seq分析流程，使用了Bowtie2+deepTools的组合工具，针对双端测序（pair-end）使用：
```
###################
# Genome files #
genome_index="/home/user/genomes/bowtie2Index/mm10/DNA/genome"
genome_gtf="/home/user/genomes/mm10.refGene.gtf"
genome_gtf_gene="/home/user/genomes/mm10.refGene.transcript.gtf"
blacklist="/home/user/genomes/mm10_blacklist_ENCFF547MET.bed"
###################
# input files #
species="Mouse/Mus musculus"  #"Human/Homo sapiens"
reads=["_R1","_R2"]           #_1,_2
ext=".fastq.gz"               #.fq.gz
path_origin=""                #测序原始数据存放路径
sample_input=[]               #input样本，不加后缀
sample_chip=[]                #chip样本
###################
# pre-pipeline #
import os
import sys
from pathlib import Path

samples=sample_input+sample_chip
p = Path(path_origin)
ss=list(p.glob("**/*"+ext))  ## 注意，**目录链接无法找到文件
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
        s2.symlink_to(s)  #s2为生成的软链接文件
        #os.symlink(s,s2)

print("#"*10)
print("Species is:")
print("All samples:", samples)
print("#"*10)
###################
rule all:
    input:
        expand("fastQC/{sample}{read}_fastqc.html",sample=samples,read=reads),
        expand("alignment/{sample}.sorted.bam",sample=samples),
        expand("alignment/{sample}.sorted.bam.bai",sample=samples),
        #expand("chipQC/{sample}_qc.txt", sample=samples),
        expand("deeptools/bigwig/{sample}.filter.dedup.bw", sample=samples),
        expand("deeptools/fragsize/{sample}_fragsize.png", sample=samples),
        expand("deeptools/profile/{sample}_gene.bed",sample=samples),
        expand("deeptools/profile/{sample}_gene.png",sample=samples),
        "deeptools/plotfingerprint.png",
        "deeptools/multibamsummary.counts.txt",
        "deeptools/correlation.png",

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

rule bowtie2:
    input:
        r1="trimdata/{sample}"+reads[0]+ext,
        r2="trimdata/{sample}"+reads[1]+ext
    output:
        bam="alignment/{sample}.sorted.bam",
        summary="alignment/{sample}.summary"
    threads: 12
    params:
        index=genome_index
    log:
        "alignment/logs/{sample}.align.log"
    shell:
        "bowtie2 -q -X 1000 -x {params.index} -1 {input.r1} -2 {input.r2} -p {threads} 2> {output.summary} "
        "| samtools view -Sbh - | samtools sort -@ {threads} -o {output.bam} 2> {log}"

rule samindex:
    input:
        "alignment/{sample}.sorted.bam"
    output:
        "alignment/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input}"

rule chipQC:
    input:
        "alignment/{sample}.sorted.bam"
    output:
        "chipQC/{sample}_qc.txt"
    log:
        "chipQC/logs/{sample}_qc.log"
    shell:"""
        bedtools bamtobed -bedpe -i {input} \
        | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$6}}' \
        | grep -v 'chrM' | sort | uniq -c \
        | awk 'BEGIN{{mt=0;m0=0;m1=0;m2=0}} ($1==1){{m1=m1+1}} ($1==2){{m2=m2+1}} {{m0=m0+1}} {{mt=mt+$1}} END{{m1_m2=-1.0; if(m2>0) m1_m2=m1/m2; printf "%d\\t%d\\t%d\\t%d\\t%f\\t%f\\t%f\\n",mt,m0,m1,m2,m0/mt,m1/m0,m1_m2}}' > {output}
        """
#echo -e "TotalReadPairs\\tDistinctReadPairs\\tOneReadPair\\tTwoReadPairs\\tNRF\\tPBC1\\tPBC2" > {output};  #首行插入表头
#sed -i '1i/TotalReadPairs\tDistinctReadPairs\tOneReadPair\tTwoReadPairs\tNRF\tPBC1\tPBC2' {output}  #效果同上

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

rule filter:
    input:
        "alignment/{sample}.sorted.bam"
    output:
        "alignment/filter/{sample}.filter.bam"
    threads: 12
    log:
        "alignment/filter/logs/{sample}.filter.log"
    shell:
        "samtools view -b -F 1804 -f 2 -q 30 -@ {threads} {input} > {output} 2> {log}"

rule deduplicates:
    input:
        "alignment/filter/{sample}.filter.bam"
    output:
        bam="alignment/filter/{sample}.filter.dedup.bam",
        matrix="alignment/filter/{sample}.filter.dedup.txt"
    threads: 16
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
# --blackListFileName {params.black}

rule profile:
    input:
        "deeptools/bigwig/{sample}.filter.dedup.bw"
    output:
        bed="deeptools/profile/{sample}_gene.bed",
        matrix="deeptools/profile/{sample}_gene.matrix.gz"
    params:
        ref=genome_gtf_gene
    threads: 12
    log:
        "deeptools/logs/{sample}.profile.log"
    shell:
        "computeMatrix reference-point --referencePoint TSS -b 2500 -a 2500 -R {params.ref} "
        "-S {input} -o {output.matrix} --outFileSortedRegions {output.bed} --skipZeros -p {threads} 2> {log}"

rule plotprofile:
    input:
        "deeptools/profile/{sample}_gene.matrix.gz"
    output:
        "deeptools/profile/{sample}_gene.png"
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
    threads:12
    params:
        black=blacklist
    log:
        "deeptools/logs/multibamsummary.log"
    shell:
        "multiBamSummary bins --bamfiles {input} --blackListFileName {params.black} --outRawCounts {output.counts} -o {output.npz} --extendReads --binSize 2000 -p {threads} 2> {log}; "

rule plotcorrelation:
    input:
        "deeptools/multibamsummary.counts.npz"
    output:
        "deeptools/correlation.png"
    shell:
        "plotCorrelation -in {input} -o {output} --corMethod pearson --whatToPlot scatterplot --skipZeros --removeOutliers --log1p"
```

### 使用MACS2 callpeak以及Homer注释peak

```
#callpeak
macs2 callpeak -f BAMPE -t CHIP.bam -c INPUT.bam --outdir PATH_TO_OUTDIR -n PRE_NAME -g mm|hs -p 0.01;  #sharp-peak，-p保留peak的阈值
macs2 callpeak -f BAMPE -t CHIP.bam -c INPUT.bam --outdir PATH_TO_OUTDIR -n PRE_NAME -g mm|hs -p 1e-5 --broad --broad-cutoff 0.01  #broad-peak，-p设置不同的值，将会影响broad-peak的pz值

#过滤peak
cat PRE_NAME.Peak |awk -F'\t' '{if ($8>5) print $0}' > PRE_NAME_p1e5.Peak

#合并peak
mergePeaks -venn <filename> <primary peak file> [additional peak/annotation files...] > outfile

#Homer注释peak
annotatePeaks.pl <peak/BED file> <genome> > <output file>

#Motif分析
# DNA
findMotifsGenome.pl peakfile.bed /mnt/d/ubuntu/genomes/mouse/mm10.fa motif_homer -size 200 -mask -dumpFasta
# RNA
findMotifsGenome.pl peakfile.bed /mnt/d/ubuntu/genomes/mouse/mm10.fa motif_homer -rna
```

### 饼图展示Homer注释的peak分布情况

这里要注意Homer注释的文件说明，  
1) The first determines the distance to the nearest TSS and assigns the peak to that gene.  
2) The second determines the genomic annotation of the region occupied by the center of the peak/region.  
也就是说，Distance to TSS和Annotation 两列的信息不完全一致，分别对应以上两步骤的对应基因。

```R
library(tidyverse)

# peak 靶基因的归属，包括Promoter <Distance to TSS (-2500~500 bp)>, Intron <Annotation>, Exon <Annotation>，Intergenic<Annotation>

## gene id转换
# BiocManager::install("rtracklayer")
library(rtracklayer)
library(ggsci)
library(scales)
# 读取gff文件
gff <- readGFF() # gff文件，如mm10.refGene.gtf
# head(gff)
id_map <- gff[, c("transcript_id", "gene_id", "gene_name")]
id_map <- unique(id_map) # gene id对照表
# id 转换
f <- "" # 注释文件
data <- read_tsv(f)
data2 <- data %>%
  select(1:4, 8, 10, 12:17, 19) %>%
  separate(col = Annotation, into = c("feature", "annotation"), sep = " \\(") %>%
  separate(col = annotation, into = c("transcript_id", "xxx"), sep = "[,\\)]")
data2 <- merge(data2, id_map, by = "transcript_id", all.x = T)
data2 <- select(data2, 2:6, 1, 17, 8:15)
# colnames(data2)
#  [1] "PeakID"
#  [2] "Chr"
#  [3] "Start"
#  [4] "End"
#  [5] "feature"
#  [6] "transcript_id"
#  [7] "gene_name"
#  [8] "Distance to TSS"
# [9] "Entrez ID"
# [10] "Nearest Unigene"
# [11] "Nearest Refseq"
# [12] "Nearest Ensembl"
# [13] "Gene Name"
# [14] "Gene Alias"
# [15] "Gene Type"
# 基因特征归类
feature2 <- as.character(data2$feature)
feature2[feature2 == "non-coding"] <- "Intergenic"
feature2[feature2 == "intron"] <- "Intron"
feature2[feature2 == "TTS"] <- "Exon"
feature2[feature2 == "3' UTR"] <- "Exon"
feature2[feature2 == "5' UTR"] <- "Exon"
feature2[feature2 == "exon"] <- "Exon"
feature2[feature2 == "promoter-TSS"] <- "Promoter"
feature2[feature2 == "promoter"] <- "Promoter"
# table(feature2)
data2$feature2 <- feature2
data2 <- data2[!is.na(feature2), ] # 去掉没有feature的peak
# 定义promoter区域大小，如 -2500bp ~ 500bp
feature2 <- as.character(data2$feature2)
feature3 <- c()
dist <- as.numeric(as.character(data2$`Distance to TSS`))
for (i in 1:length(dist)) {
  if (dist[i] > (-2500) & dist[i] < 500 & (feature2[i] != "Promoter")) {
    feature3 <- c(feature3, "Promoter")
  } else {
    feature3 <- c(feature3, feature2[i])
  }
}
# table(feature3)
data2$feature3 <- feature3
write.table(data2, paste0(f, "2"), row.names = F, sep = "\t")
# 饼图展示HOMER 注释peak的分布情况
mypal <- pal_npg("nrc", alpha = 0.6)(4)
show_col(mypal)
names(mypal) <- c("Intergenic", "Intron", "Exon", "Promoter")

data_pie <- data2 %>%
  group_by(feature3) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  mutate(prop = round(n / sum(n), 3) * 100) %>%
  mutate(ypos = cumsum(prop) - 0.5 * prop)
data_pie$feature3 <- factor(data_pie$feature3, levels = rev(data_pie$feature3))
data_pie$prop2 <- sprintf("%0.2f", data_pie$prop) ## 保留小数后的0
data_pie <- na.omit(data_pie)
mypalx <- mypal[as.character(data_pie$feature3)]
data_pie %>%
  ggplot(aes(x = 2, y = prop, fill = feature3)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  # scale_fill_npg(alpha=0.8,)+
  scale_fill_manual(values = mypalx) +
  # geom_text(aes(label = paste0(prop, "%")), position = position_stack(vjust=0.5),colour='white')+
  geom_text(aes(y = ypos, label = paste0(prop2, "%")), color = "black", size = 4, nudge_x = 0.8) +
  theme_void() +
  labs(fill = "Features") +
  theme(text = element_text(size = 14), legend.key.height = unit(0.4, "cm"), legend.key.width = unit(0.4, "cm"))
ggsave(paste0(f, ".png"), width = 10, height = 10, units = "cm", dpi = 1200)
```

### IGV 可视化peak

可自行生成参考基因组及相应注释文件的genome文件。

---

铭记九一八事变·1931年9月18日日本全面侵华 😠💢💢


