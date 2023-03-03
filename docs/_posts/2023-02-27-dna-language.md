---
layout: post
title: DNA序列的数据表示方法
tags: [DNA, Computational biology]
---

### DNA相关知识

[Molecular Structure of Nucleic Acids: A Structure for Deoxyribose Nucleic Acid](https://www.nature.com/articles/171737a0)

|||
|--|--|
|![Animal cell diagram](https://upload.wikimedia.org/wikipedia/commons/1/11/Animal_Cell.svg)|Components of a typical animal cell: 1.Nucleolus, 2.Nucleus, 3.Ribosome, 4.Vesicle, 5.Rough endoplasmic reticulum, 6.Golgi apparatus, 7.Cytoskeleton, 8.Smooth endoplasmic reticulum, 9.Mitochondrion, 10.Vacuole, 11.Cytosol, 12.Lysosome, 13.Centrosome, 14.Cell membrane|

- Wiki

[DNA](https://en.wikipedia.org/wiki/DNA)  

![The structure of the DNA double helix (type B-DNA)](https://upload.wikimedia.org/wikipedia/commons/4/4c/DNA_Structure%2BKey%2BLabelled.pn_NoBB.png)

[RNA](https://en.wikipedia.org/wiki/RNA)  
[Protein](https://en.wikipedia.org/wiki/Protein)  

![A ribosome produces a protein using mRNA as template](https://upload.wikimedia.org/wikipedia/commons/b/b1/Ribosome_mRNA_translation_en.svg)
![3D structure of the protein myoglobin showing turquoise α-helices](https://upload.wikimedia.org/wikipedia/commons/6/60/Myoglobin.png)

[Gene](https://en.wikipedia.org/wiki/Gene)  
[Central dogma](https://en.wikipedia.org/wiki/Central_dogma_of_molecular_biology)  

![Information flow in biological systems](https://upload.wikimedia.org/wikipedia/commons/d/dd/Extended_Central_Dogma_with_Enzymes.jpg)

[DNA and RNA codon tables](https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables)  
[Transcription](https://en.wikipedia.org/wiki/Transcription_(biology))  

![Simplified diagram of mRNA synthesis and processing.](https://upload.wikimedia.org/wikipedia/commons/9/9b/MRNA.svg)

[Translation](https://en.wikipedia.org/wiki/Translation_(biology))  

![Overview of the translation of eukaryotic messenger RNA](https://upload.wikimedia.org/wikipedia/commons/4/44/Protein_synthesis.svg)
![Initiation and elongation stages of translation](https://upload.wikimedia.org/wikipedia/commons/0/01/Translation_-_Initiation_%26_Elongation.svg)

[Genome](https://en.wikipedia.org/wiki/Genome)  

![Human genome](https://upload.wikimedia.org/wikipedia/commons/b/b1/Human_karyotype_with_bands_and_sub-bands.png)

[Epigenetics](https://en.wikipedia.org/wiki/Epigenetics)

![Epigenetic mechanisms](https://upload.wikimedia.org/wikipedia/commons/f/fc/Epigenetic_mechanisms.png)

[Regulation of gene expression](https://en.wikipedia.org/wiki/Regulation_of_gene_expression)

![](https://www.science.org/cms/10.1126/science.aat8950/asset/d0a49575-0040-4f9e-8dc3-eb03cbe80970/assets/graphic/361_1332_f1.jpeg)

ref: [Chromatin plasticity: A versatile landscape that underlies cell fate and identity. 2018](https://www.science.org/doi/10.1126/science.aat8950)

![](https://www.science.org/cms/10.1126/science.aau0320/asset/3253216a-aff0-4b6a-9d8d-423fae8ace76/assets/graphic/361_1341_f2.jpeg)
![](https://www.science.org/cms/10.1126/science.aau0320/asset/412087a1-de48-433f-bc76-9ca9b70b9a1b/assets/graphic/361_1341_f4.jpeg)

ref: [Developmental enhancers and chromosome topology. 2018](https://www.science.org/doi/10.1126/science.aau0320)

![](https://ars.els-cdn.com/content/image/1-s2.0-S0959437X20301726-gr1.jpg)

ref: [The dynamics of chromatin architecture in brain development and function. 2021](https://doi.org/10.1016/j.gde.2020.12.008)

### 参考基因组

- [NCBI-Genome](https://www.ncbi.nlm.nih.gov/data-hub/genome/)
- [Ensembl](https://asia.ensembl.org/index.html)
- [GENECODE](https://www.gencodegenes.org/)
- [UCSC](https://genome.ucsc.edu/)
- [REPBASE](https://www.girinst.org/)  #重复序列

基因组注释统计信息：  
人类（Homo sapiens）：[Human assembly and gene annotation](https://asia.ensembl.org/Homo_sapiens/Info/Annotation)  
小鼠（Mus musculus）：[Mouse assembly and gene annotation](https://asia.ensembl.org/Mus_musculus/Info/Annotation)

### 数据表示方法

[List of file formats](https://en.wikipedia.org/wiki/List_of_file_formats)

#### 字符形式

- fasta

[FASTA format](https://en.wikipedia.org/wiki/FASTA_format)  
[K-mer](https://en.wikipedia.org/wiki/K-mer)

- fastq

[FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format)  
[SAM format](https://en.wikipedia.org/wiki/SAM_(file_format))  

- sequence features

[DiProDB: a database for dinucleotide properties](https://doi.org/10.1093/nar/gkn597)
[Seq2Feature: a comprehensive web-based feature extraction tool](https://doi.org/10.1093/bioinformatics/btz432)

|Properties I|Properties II (details)|Description|
|--|--|--|
|[Physicochemical properties](https://www.iitm.ac.in/bioinfo/SBFE/physico_dna.html)|Stacking energy, Enthalpy, Entropy, Flexibility shift, Flexibility_slide, Free energy, Melting Temperature, Mobility to bend towards major groove, Mobility to bend towards minor groove, Probability contacting nucleosome core, Rise stiffness, Roll stiffness, Shift stiffness, Slide stiffness, Tilt stiffness, Twist stiffness|理化性质|
|[Conformational properties](https://www.iitm.ac.in/bioinfo/SBFE/conform_dna.html)|Bend, Rise, Roll, Inclination, Major Groove Depth, Major Groove Distance, Major Groove Size, Major Groove Width, Minor Groove Depth, Minor Groove Distance, Minor Groove Size, Minor Groove Width, Shift, Propeller Twist, Slide, Tilt, Tip, Twist|构象性质|
|Nucleotide content|Adenine content, Cytosine content, GC content, Guanine content, Keto (GT) content, Purine (AG) content, Thymine content, Pyrimidine (CT)|碱基含量|

- gene features

![Eukaryotic protein-coding gene](https://upload.wikimedia.org/wikipedia/commons/5/54/Gene_structure_eukaryote_2_annotated.svg)
![Prokaryotic operon of protein-coding genes](https://upload.wikimedia.org/wikipedia/commons/f/fd/Gene_structure_prokaryote_2_annotated.svg)

- toolkits

[bio.tools](https://bio.tools/)  
[Sequence Manipulation Suite](https://www.bioinformatics.org/sms2/)  
[Seq2Feature](https://www.iitm.ac.in/bioinfo/SBFE/index.html)

#### 数字形式

- One-Hot

One-hot方法，将DNA序列直观的表示为0-1的矩阵，如图所示：
![]({{ '/assets/img/mds/dna_onehot.png' | relative_url }})

[DeepSEA](https://www.nature.com/articles/nmeth.3547)  
[DeepBind](https://www.nature.com/articles/nbt.3300)  
[Basset](https://genome.cshlp.org/content/26/7/990.short)

- Embedding

词向量表示方法，借鉴自然语言处理的思路，先将DNA序列切分为一定长度的k-mer，然后把k-mer的序列片段视为词语，一段序列视为句子，训练词向量模型，可参考[博客](https://jalammar.github.io/illustrated-word2vec/)理解。
![]({{ '/assets/img/mds/dna_embedding.png' | relative_url }})

[word2vec](https://arxiv.org/pdf/1411.2738.pdf)  
[Transformer](https://arxiv.org/pdf/1706.03762.pdf)  

[dna2vec](https://arxiv.org/pdf/1701.06279.pdf)  
[kmer2vec](https://doi.org/10.1089/cmb.2021.0536)  
[DNABERT](https://doi.org/10.1093/bioinformatics/btab083)  

---

`The real voyage of discovery consists not in seeking new lands but seeing with new eyes.` --Marcel Proust, 1923, La Prisonierre
