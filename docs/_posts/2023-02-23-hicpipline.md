---
layout: post
title: 三维基因组
tags: [HiC, Juicer, HiC-Pro, Docker]
---

### HiC 原理简单介绍

[A 3D Map of the Human Genome at Kilobase Resolution Reveals Principles of Chromatin Looping](https://www.cell.com/fulltext/S0092-8674(14)01497-4)

- Hi-C protocol

![](https://www.cell.com/cms/attachment/6599fde9-c119-4849-bb0c-e01eb117e87e/gr1.jpg)

- Chromatin architecture:

1. Megadomains: 5–20 Mb intervals  
2. Topologically associated domains (TADs): ∼1 Mb  
3. Contact domains are often preserved across cell types: 40 kb to 3 Mb (median size 185 kb)  
4. Loops are short (<2 Mb) and strongly conserved across cell types and between human and mouse  
5. Eigenvector: Compartment A is highly enriched for open chromatin; compartment B is enriched for closed chromatin  
6. Genes whose promoters are associated with a loop are much more highly expressed than genes whose promoters are not associated with a loop (6-fold)

![]({{ '/assets/img/hic/chromatin.png' | relative_url }})
[ref: Current Opinion in Genetics & Development, Izabela Harabula. 2021](https://doi.org/10.1016/j.gde.2020.12.008)

---

### 电脑配置信息

```sh
uname -a  查看内核/操作系统
#Linux chunyu-PowerEdge-R720 5.13.0-43-generic #48~20.04.1-Ubuntu SMP Thu May 12 12:59:19 UTC 2022 x86_64 x86_64 x86_64 GNU/Linux
df -h  #硬盘信息
#文件系统                                          容量  已用  可用 已用% 挂载点
#rpool/USERDATA/lxh_gscj6g                        46T   298G  46T  1%   /home/lxh
#rpool/ROOT/ubuntu_ktv56z/usr/local               46T   3.7G  46T  1%   /usr/local
free -h  #内存信息
#              总计         已用        空闲      共享    缓冲/缓存    可用
#内存：       503Gi         260Gi       241Gi     10Mi   1.3Gi        240Gi
#交换：       2.0Gi         0B          2.0Gi
cat /proc/cpuinfo  #查看CPU信息
nvidia-smi  #查看GPU信息和使用情况
cut -d: -f1 /etc/passwd  #查看系统所有用户
```

**CPU信息展示**:

![]({{ '/assets/img/hic/cpuinfo.png' | relative_url }})

**GPU信息展示**：

![]({{ '/assets/img/hic/gpuinfo.png' | relative_url }})

>表头释义：  
>Fan：显示风扇转速，数值在0到100%之间，是计算机的期望转速，如果计算机不是通过风扇冷却或者风扇坏了，显示出来就是N/A；  
>Temp：显卡内部的温度，单位是摄氏度；  
>Perf：表征性能状态，从P0到P12，P0表示最大性能，P12表示状态最小性能；  
>Pwr：能耗表示；  
>Bus-Id：涉及GPU总线的相关信息；  
>Disp.A：是Display Active的意思，表示GPU的显示是否初始化；  
>Memory Usage：显存的使用率；  
>Volatile GPU-Util：浮动的GPU利用率；  
>Compute M：计算模式；  

### 安装docker

参考教程：[Install Docker Engine on Ubuntu](https://docs.docker.com/engine/install/ubuntu/)

```sh
#Uninstall old versions
sudo apt-get remove docker docker-engine docker.io containerd runc
sudo rm -rf /var/lib/docker
sudo rm -rf /var/lib/containerd

##Install using the repository
#Set up the repository
sudo apt-get update
sudo apt-get install \
    ca-certificates \
    curl \
    gnupg \
    lsb-release
sudo mkdir -m 0755 -p /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

#Install Docker Engine
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
sudo docker run hello-world

#添加用户进docker组，支持非root使用
sudo gpasswd -a username docker
newgrp docker
```

### CUDA and NVIDIA GPU介绍

CUDA Toolkit 的[下载](https://developer.nvidia.com/cuda-downloads)、[安装](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html)和[使用](https://docs.nvidia.com/cuda/cuda-quick-start-guide/index.html)

### 创建hic_pipeline的docker环境

```sh
# hic-pipeline/
cd docker/hic-pipeline
docker build -t encodedcc/hic-pipeline:1.15.1 .
cd docker/delta
docker build -t encodedcc/hic-pipeline:1.15.1_delta .
cd docker/hiccups
docker build -t encodedcc/hic-pipeline:1.15.1_hiccups .

#或者从docker hub下载镜像
docker pull encodedcc/hic-pipeline:1.15.1
docker pull encodedcc/hic-pipeline:1.15.1_delta
docker pull encodedcc/hic-pipeline:1.15.1_hiccups

#列出安装的镜像
docker images
#REPOSITORY               TAG              IMAGE ID       CREATED         SIZE
#encodedcc/hic-pipeline   1.15.1_delta     41fc8e6dd680   2 hours ago     5.81GB
#encodedcc/hic-pipeline   1.15.1           95c8bf750530   3 months ago    2.21GB
#encodedcc/hic-pipeline   1.15.1_hiccups   363f736bd618   7 months ago    2.52GB
#hello-world              latest           feb5d9fea6a5   17 months ago   13.3kB
```

---

### Juicer & Juicerbox

[Juicer](https://github.com/aidenlab/juicer) is a one-click pipeline for processing terabase scale Hi-C datasets.
![](https://github.com/theaidenlab/juicer/wiki/images/graphic_juicer.png)

基于Juicer的 ENCODE HiC-pipline 分析流程搭建  
The [ENCODE pipeline for processing Hi-C data](https://github.com/ENCODE-DCC/hic-pipeline) based on Juicer.

#### 简单使用

```sh
git clone https://github.com/ENCODE-DCC/hic-pipeline.git
cd hic-pipeline
pip install caper  #依赖caper实现，首次运行时安装
caper run hic.wdl -i tests/functional/json/test_hic.json --docker
```

#### 使用指南

1. 在hic-pipeline的根目录下使用，每次使用将会在目录下生成日志文件cromwell.out和目录cromwell-workflow-logs，以及结果文件存放目录hic，每次运行将随机产生序列号
2. 配置文件（json格式）的设置，参考[Workflows说明](https://github.com/ENCODE-DCC/hic-pipeline/blob/dev/docs/reference.md)

3. 输入文件  
需要准备文件包括测序数据，参考基因组序列文件（生成酶切片段），参考基因组bwa索引文件，参考基因组染色体大小文件，参考基因组酶切片段坐标信息文件。  
hic.wdl的配置文件json内容如下：

```json
{
  "hic.assembly_name": "ce10",
  "hic.reference_index": "tests/data/ce10_selected.tar.gz",
  "hic.chrsz": "tests/data/ce10_selected.chrom.sizes.tsv",
  "hic.reference_fasta": "",

  "hic.restriction_enzymes": ["MboI"],  #MboI, HindIII, DpnII, and none
  "hic.restriction_sites": "tests/data/ce10_selected_MboI.txt.gz",
  "hic.ligation_site_regex": "",

  "hic.fastq": [
    [
      {
        "read_1": "tests/data/merged_read1.fastq.gz",
        "read_2": "tests/data/merged_read2.fastq.gz"
      },
    ],
  ],  #From fastq

  "hic.input_hic": "",  #From .hic
  "hic.normalization_methods": "",  #VC, VC_SQRT, KR, SCALE, GW_KR, GW_SCALE, GW_VC, INTER_KR and INTER_SCALE

  "hic.no_pairs": false,
  "hic.no_call_loops": false,  #requires GPUs
  "hic.no_call_tads": false,
  "hic.no_delta": false,
  "hic.no_slice": false,  #requires GPUs
  "hic.no_eigenvectors": false,

  "hic.align_num_cpus": 32,
  "hic.align_ram_gb_in_situ": 64,
  "hic.align_ram_gb_intact": 88,
  "hic.align_disk_size_gb_in_situ": 1000,
  "hic.align_disk_size_gb_intact": 1500,
  "hic.chimeric_sam_nonspecific_disk_size_gb": 7000,
  "hic.chimeric_sam_specific_disk_size_gb": 1500,
  "hic.dedup_ram_gb_in_situ": 32,
  "hic.dedup_ram_gb_intact": 48,
  "hic.dedup_disk_size_gb_in_situ": 5000,
  "hic.dedup_disk_size_gb_intact": 7500,
  "hic.create_hic_num_cpus": ,
  "hic.create_hic_ram_gb": ,
  "hic.create_hic_juicer_tools_heap_size_gb": ,
  "hic.create_hic_disk_size_gb": ,
  "hic.add_norm_num_cpus": ,
  "hic.add_norm_ram_gb": ,
  "hic.add_norm_disk_size_gb": ,
  "hic.create_accessibility_track_ram_gb": ,
  "hic.create_accessibility_track_disk_size_gb": ,
}
```

4. 输出文件

|task|outfiles|description|
|--|--|--|
|get_ligation_site_regex|ligation_site_regex.txt|酶切位点|
|normalize_assembly_name|normalized_assembly_name.txt, is_supported.txt|参考基因组标识，如GRCm38|
|align|aligned.bam, ligation_count.txt|数据比对结果文件|
|chimeric_sam_specific|chimeric_sam_specific.bam, result_norm.txt.res.txt|比对文件中的chimeric reads，保留的unambiguous部分？|
|merge|merged.bam|合并比对结果文件，包括R1和R2的正常比对以及部分的unambiguous chimeric read|
|dedup|merged_dedup.bam|去除PCR重复后的比对文件|
|bam_to_pre|merged_nodups_~{quality}.txt.gz, merged_nodups_~{quality}_index.txt.gz|生成.hic文件的中间文件，详见[Juicer-Pre](https://github.com/aidenlab/juicer/wiki/Pre#file-format)|
|pre_to_pairs|pairix.bsorted.pairs.gz|-|
|calculate_stats|stats_~{quality}.txt, stats_~{quality}.json, stats_~{quality}_hists.m|关于比对结果的数据统计信息，默认保留MAPQ>0和MAPQ>30两种结果|
|create_hic|inter_~{quality}_unnormalized.hic|生成原始.hic文件，及未经过标准化|
|add_norm|inter_~{quality}.hic|标准化的.hic文件|
|arrowhead|_~{quality}.bedpe.gz|Arrowhead Algorithm for Domain Annotation，检测TADs|
|hiccups|merged_loops_~{quality}.bedpe.gz|Hi-C Computational Unbiased Peak Search for Peak Calling，检测chromatin loops|
|hiccups_2|merged_loops_~{quality}.bedpe.gz|-|
|localizer|localized_loops_~{quality}.bedpe.gz|-|
|delta|predicted_loops_merged.bedpe.gz, predicted_domains_merged.bedpe.gz, predicted_stripes_merged.bedpe.gz, predicted_loop_domains_merged.bedpe.gz|检测TAD、loops，参考[Delta](https://delta.ngdc.cncb.ac.cn/)|
|create_eigenvector|eigenvector_~{resolution}.wig, eigenvector_~{resolution}.bw|Compartment A is highly enriched for open chromatin; compartment B is enriched for closed chromatin，分析区别compartment A和compartment B|
|slice|slice_subcompartment_clusters_~{resolution}.bed.gz|-|
|create_accessibility_track|inter_30.bw|-|
||||

**问题集合**:

1. docker: permission denied while trying to connect to the Docker daemon socket at unix:///var/run/docker.sock  
解决方法：将用户加入docker组，  
*sudo gpasswd -a username docker;*  
*newgrp docker*

2. Unable to find image 'encodedcc/hic-pipeline@ locally  
可能是因为网络原因，无法直接从[docker hub](https://hub.docker.com/) pull相关的image，  
解决方法：  
1）创建本地目录encodedcc，将docker目录中的镜像文件Dockerfile进行复制，  
*cp -r docker/ encodedcc/*  
这样的话好像首次次运行会根据Dockerfile创建相应的容器，会耗费一定时间。
2）也可自行创建容器，  
*cd docker/hic-pipeline*  
*docker build -t encodedcc/hic-pipeline:1.15.1 .*

#### HiC数据说明

参考[Data Production and Processing Standard of 
the Hi-C Mapping Center](https://www.encodeproject.org/documents/75926e4b-77aa-4959-8ca7-87efcba39d79/@@download/attachment/comp_doc_7july2018_final.pdf)

#### Quality Control

1. Hi-C Library Statistics  
2. Aggregate Peak Analysis (APA) for Peak Calling Quality Control  
3. Experimental validation

![]({{ '/assets/img/hic/tab2.png' | relative_url }})
![]({{ '/assets/img/hic/fig2.png' | relative_url }})
![]({{ '/assets/img/hic/fig3.png' | relative_url }})

#### Hi-C Data Processing

Juicer pipeline consists of two main parts:  

1. the transformation of raw reads into a highly compressed binary file that provides fast random access to contact matrices stored at multiple resolutions with multiple normalization schemes;  
2. the analysis portion, which operates on the binary file and automatically annotates contact domains, loops, anchors, and compartments.

![]({{ '/assets/img/hic/fig6.png' | relative_url }})
![]({{ '/assets/img/hic/tab3.png' | relative_url }})
![]({{ '/assets/img/hic/tab4.png' | relative_url }})

- Fastqs to contact matrices

1. Sequence alignment  
1.1 Filtering of abnormal alignments  
**normal**:  each read in a read pair will align to a single site in the genome.  
**chimeric**: at least one of the two reads comprises multiple subsequences, each of which align to different parts of the genome. For instance, the first 50 base pairs might map perfectly to one position, whereas the next 50 map perfectly to a second position several megabases away. Chimeric read pairs are classified as **"unambiguous"** or **"ambiguous"**, ambiguous chimeric reads are not included in Hi-C maps.  
**unalignable**: they have at least one end that cannot be successfully aligned.  
1.2 Filtering of duplicates  
1.3 Filtering of low-quality alignments  
1.4 Construction of contact matrices  
默认bin大小包括2.5 Mb, 1 Mb, 500 kb, 250 kb, 100 kb, 50 kb, 25 kb, 10 kb, and 5 kb，可通过参数设置  

2. Contact Matrix Normalization  
2.1 Vanilla coverage normalization (VC normalization)  
2.2 A role for matrix balancing in Hi-C  
2.3 New methods for matrix balancing

3. Diploid and Higher Order Maps  
3.1 Construction of Diploid Hi-C maps  
3.2 Higher-order Contact Analysis  

- Feature Annotation and Other Downstream Analysis Methods

1. HiCCUPS: Hi-C Computational Unbiased Peak Search for Peak Calling  
Much of the work on genome architecture so far has centered on the study of chromatin looping. **Chromatin loops manifest as local peaks in a proximity ligation dataset**, which occur between two points whenever they interact with each other significantly more than with random points in their neighborhood.

- Visualization and Data Integration

1. [Juicebox](https://github.com/aidenlab/Juicebox) for Data Visualization and Integration with ENCODE tracks

![]({{ '/assets/img/hic/fig22.png' | relative_url }})

---

### HiC-Pro 分析

[HiC-Pro](https://github.com/nservant/HiC-Pro) was designed to process Hi-C data, from raw fastq files (paired-end Illumina data) to normalized contact maps.
![](http://nservant.github.io/HiC-Pro/_images/hicpro_wkflow.png)

#### 简单使用

```sh
HiC-Pro -i FULL_PATH_TO_DATA_FOLDER -o FULL_PATH_TO_OUTPUTS -c MY_LOCAL_CONFIG_FILE

# FULL_PATH_TO_DATA_FOLDER structure
   + FULL_PATH_TO_DATA_FOLDER
     + sample1
       ++ file1_R1.fastq.gz
       ++ file1_R2.fastq.gz
       ++ ...
     + sample2
       ++ file1_R1.fastq.gz
       ++ file1_R2.fastq.gz
     *...
```

#### 使用指南

详细参考[HIC-PRO](http://nservant.github.io/HiC-Pro/MANUAL.html)

- 安装（基于conda环境）

```sh
git clone https://github.com/nservant/HiC-Pro.git
cd HiC-Pro
conda env create -f environment.yml  #将name改为hicenv
conda activate hicenv

make configure
make install
echo 'export PATH="/usr/local/bin/HiC-Pro_3.1.0/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

HiC-Pro --help
```

- 输入文件

测序数据，参考基因组bowtie2索引文件，参考基因组染色体大小文件及参考基因组酶切片段坐标信息文件（bed格式），  
config配置文件如下：

```python
# Please change the variable settings below if necessary

#########################################################################
## Paths and Settings  - Do not edit !
#########################################################################

TMP_DIR = tmp
LOGS_DIR = logs
BOWTIE2_OUTPUT_DIR = bowtie_results
MAPC_OUTPUT = hic_results
RAW_DIR = rawdata

#######################################################################
## SYSTEM AND SCHEDULER - Start Editing Here !!
#######################################################################
N_CPU = 12
SORT_RAM = 1000M
LOGFILE = hicpro.log

JOB_NAME = 
JOB_MEM = 
JOB_WALLTIME = 
JOB_QUEUE = 
JOB_MAIL = 

#########################################################################
## Data
#########################################################################

PAIR1_EXT = _R1
PAIR2_EXT = _R2

#######################################################################
## Alignment options
#######################################################################

MIN_MAPQ = 10

BOWTIE2_IDX_PATH = /home/lxh/HiC_test/refgenomes/indexs/GRCh38/bowtie2index
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder

#######################################################################
## Annotation files
#######################################################################

REFERENCE_GENOME = refgen
GENOME_SIZE = /home/lxh/HiC_test/refgenomes/chromsizes/GRCh38.chrom.sizes

#######################################################################
## Allele specific analysis
#######################################################################

ALLELE_SPECIFIC_SNP = 

#######################################################################
## Capture Hi-C analysis
#######################################################################

CAPTURE_TARGET =
REPORT_CAPTURE_REPORTER = 1

#######################################################################
## Digestion Hi-C
#######################################################################

GENOME_FRAGMENT = /home/lxh/HiC_test/refgenomes/enzymesites/hicpro/GRCh38_mboi_or_dpnii.bed
LIGATION_SITE = GATC
MIN_FRAG_SIZE = 
MAX_FRAG_SIZE =
MIN_INSERT_SIZE =
MAX_INSERT_SIZE =

#######################################################################
## Hi-C processing
#######################################################################

MIN_CIS_DIST =
GET_ALL_INTERACTION_CLASSES = 1
GET_PROCESS_SAM = 0
RM_SINGLETON = 1
RM_MULTI = 1
RM_DUP = 1

#######################################################################
## Contact Maps
#######################################################################

BIN_SIZE = 20000 40000 150000 500000 1000000
MATRIX_FORMAT = upper

#######################################################################
## Normalization
#######################################################################
MAX_ITER = 100
FILTER_LOW_COUNT_PERC = 0.02
FILTER_HIGH_COUNT_PERC = 0
EPS = 0.1
```

- 输出文件

---

### 参考基因组

```sh
##下载序列文件和注释文件
#对于人类和小鼠，从GENCODE (https://www.gencodegenes.org/) 下载参考基因组序列文件和注释文件
#Genome assembly GRCh38 /2023-02-08
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gff3.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz

#Genome assembly GRCm38 /2020-04-29
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gff3.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz

##构建index文件
#GRCh38
bowtie2-build --threads 24 /home/lxh/HiC_test/refgenomes/sequences/GRCh38/GRCh38.primary_assembly.genome.fa refgen
bwa index -p refgen -a bwtsw /home/lxh/HiC_test/refgenomes/sequences/GRCh38/GRCh38.primary_assembly.genome.fa
#GRCm38
bowtie2-build --threads 24 /home/lxh/HiC_test/refgenomes/sequences/GRCm38/GRCm38.primary_assembly.genome.fa refgen
bwa index -p refgen -a bwtsw /home/lxh/HiC_test/refgenomes/sequences/GRCm38/GRCm38.primary_assembly.genome.fa

#染色体大小
samtools faidx GRCh38.primary_assembly.genome.fa
cut -f1,2 GRCh38.primary_assembly.genome.fa.fai > GRCh38.chrom.sizes

##酶切片段文件
#HindIII
----------------
5' A^AGCTT 3'
3'  TTCGA^A 5'
----------------
#BglII
----------------
5' A^GATCT 3'
3'  TCTAG^A 5'
----------------
#DpnII
----------------
5' ^GATC 3'
3'  CTAG^ 5'
----------------
#MboI
----------------
5' ^GATC 3'
3'  CTAG^ 5'
----------------
##HiC-Pro
#GRCh38
/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py -r A^AGCTT -o GRCh38_hindiii.bed /home/lxh/HiC_test/refgenomes/sequences/GRCh38/GRCh38.primary_assembly.genome.fa  #HindIII
/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py -r hindiii -o GRCh38_hindiii.bed /home/lxh/HiC_test/refgenomes/sequences/GRCh38/GRCh38.primary_assembly.genome.fa  #效果同上
# /usr/local/bin/HiC-Pro_3.1.0/utils/digest_genome.py -r hindiii dpnii -o GRCh38_hindiii_dpnii.bed /home/lxh/HiC_test/refgenome/sequences/GRCh38/GRCh38.primary_assembly.genome.fa  #双酶切
/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py -r mboi -o GRCh38_mboi_or_dpnii.bed /home/lxh/HiC_test/refgenomes/sequences/GRCh38/GRCh38.primary_assembly.genome.fa
/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py -r bglii -o GRCh38_bglii.bed /home/lxh/HiC_test/refgenomes/sequences/GRCh38/GRCh38.primary_assembly.genome.fa

#GRCm38
/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py -r hindiii -o GRCm38_hindiii.bed /home/lxh/HiC_test/refgenomes/sequences/GRCm38/GRCm38.primary_assembly.genome.fa
/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py -r mboi -o GRCm38_mboi_or_dpnii.bed /home/lxh/HiC_test/refgenomes/sequences/GRCm38/GRCm38.primary_assembly.genome.fa
/usr/local/bin/HiC-Pro_3.1.0/bin/utils/digest_genome.py -r bglii -o GRCm38_bglii.bed /home/lxh/HiC_test/refgenomes/sequences/GRCm38/GRCm38.primary_assembly.genome.fa

##juicer
# 'HindIII'     : 'AAGCTT',
# 'DpnII'       : 'GATC',
# 'MboI'        : 'GATC',
# 'Sau3AI'      : 'GATC',
# 'Arima'       : [ 'GATC', 'GANTC' ]
#restriction_site_locations.json
{
  "make_restriction_site_locations.reference_fasta": "/home/lxh/HiC_test/refgenomes/sequences/GRCm38/GRCm38.primary_assembly.genome.fa",
  "make_restriction_site_locations.restriction_enzyme": "DpnII",
  "make_restriction_site_locations.assembly_name": "GRCm38"
}

caper run ../make_restriction_site_locations.wdl -i restriction_site_locations.json --docker
#也可直接下载，mm10的链接有访问权限，只能下载到GRCh38的文件
#DpnII, MboI  GRCh38
wget https://www.encodeproject.org/files/ENCFF132WAM/@@download/ENCFF132WAM.txt.gz
#HindIII  GRCh38
https://www.encodeproject.org/files/ENCFF984SUZ/@@download/ENCFF984SUZ.txt.gz

```

---
