---
layout: post
title: Deep learning 在基因组学中的应用
tags: [deep learning, genomics, papers]
---

### 综述

Google学术搜索关键词 deep learning; genomics; review

- [A primer on deep learning in genomics. 2018](https://www.nature.com/articles/s41588-018-0295-5)

>云平台和软件资源 [DeepLearning Resources](https://www.nature.com/articles/s41588-018-0295-5/tables/1)

- [Deep learning: new computational modelling techniques for genomics. 2019](https://www.nature.com/articles/s41576-019-0122-6)

>Applications of genomics:  
>
>1. finding associations between genotype and phenotype
>2. discovering biomarkers for patient stratification
>3. predicting the function of genes
>4. charting biochemically active genomic regions such as transcriptional enhancers
>
>Machine learning algorithms are designed to automatically detect patterns in data. A central issue is that classification performance depends heavily on the quality and the relevance of handcrafted features. Deep learning, a subdiscipline of machine learning, addresses this issue by embedding the computation of features into the machine learning model itself to yield end-to-end models
>![fig1](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41576-019-0122-6/MediaObjects/41576_2019_122_Fig2_HTML.png?as=webp)
>CNNs appplications:
>
>1. classifying transcription factor binding sites
>2. predicting molecular phenotypes such as chromatin features
>3. DNA contact maps
>4. DNA methylation
>5. gene expression
>6. translation effiency
>7. RBP binding
>8. microRNA (miRNA) targets
>9. predict the specificity of guide RNA
>10. denoise ChIP–seq
>11. enhance Hi-C data resolution
>12. predict the laboratory of origin from DNA sequences
>13. call genetic variants
>14. ...
>
>??? Modelling molecular phenotypes from the linear DNA sequence, albeit a crude approximation of the chromatin, can be improved by allowing for long-range dependencies and allowing the model to implicitly learn aspects of the 3D organization, such as promoter–enhancer looping.
>
>RNNs applications:
>
>1. predicting single-cell DNA methylation states
>2. RBP binding
>3. transcription factor binding
>4. DNA accessibility
>5. miRNA biology ...
>
>GCNs applications:
>
>1. derive new features of proteins from protein–protein interaction networks
>2. modelling polypharmacy side effects
>3. predict various molecular properties including solubility, drug efficacy and photovoltaic efficiency
>4. predicting binarized gene expression given the expression of other genes
>5. classification of cancer subtypes
>6. ...
>
>??? In simple models such as linear models, the parameters of the model often measure the contribution of an input feature to prediction. Therefore, they can be directly interpreted in cases where the input features are relatively independent. By contrast, the parameters of a deep neural network are difficult to interpret because of their redundancy and nonlinear relationship with the output.
>
>**Feature importance scores** can be divided into two main categories on the basis of whether they are computed using input perturbations or using backpropagation.
>
>1. Perturbation-based approaches systematically perturb the input features and observe the change in the output.
>2. Backpropagation-based approaches, to the contrary, are much more computationally efficient.
>
>Advantages:
>
>1. end-to-end learning
>2. deal with multimodal data effectively
>3. abstraction of the mathematical details

- [Deep Learning for Genomics: A Concise Overview. 2018](https://arxiv.org/abs/1802.00810)

>基因组学历程：DNA双螺旋结构-1953，人类基因组计划-2001，基因组研究相关项目（FANTOM-2001，ENCODE-2013， Roadmap Epigenomics-2015等）
>
>基因组学和传统遗传学：Genomic research aims to understand the genomes of different species. It studies the roles assumed by multiple genetic factors and the way they interact with the surrounding environment under different conditions. In contrast to genetics that deals with limited number of specific genes, genomics takes a global view that involves the entirety of genes possessed by an organism.
>
>|Models or strategies|Researches|Description|
>|--|--|--|
>|CNN|DeepBind, DeepSEA, Basset|When applying convolutional neural networks in genomic, since deep learning models are always over-parameterized, simply changing the network depth would not account for much improvement of model performance. Researchers should pay more attention to particular techniques that can be used in CNNs, such as the kernel size, the number of feature map, the design of pooling or convolution kernels, the choice of window size of input DNA sequences, etc., or include prior genomic information if possible.|
>|RNN|ProLanGO, DeepNano, DanQ, |-|
>|Autoencoder|-|Now they have proved successful for feature extraction because of being able to learn a compact representation of input through the encode-decode procedure. When applying autoencoders, one should be aware that the better reconstruction accuracy does not necessarily lead to model improvement.|
>|CNN-RNN|Deep GDashboard, |-|
>|Transfer Learning & multitask learning|PEDLA, TFImpute|Transfer learning is such a framework that allows deep learning to adapt the previously-trained model to exploit a new but relevant problem more effectively|
>|Multi-view learning|gRNM, |Multi-view learning can be achieved by, for example, concatenating features, ensemble methods, or multi-modal learning.|
>
>??? Applications of deep learning in genomic problems have fully proven its power. Although the pragmatism of deep learning is surprisingly successful, this method suffers from lacking the physical transparency to be well interpreted so as to better assist the understanding of genomic problems.
>
>**Genomics applications:**
>
>1. Gene Expression  
>1.1 Characterization  
>1.2 Prediction
>2. Regulatory Genomics  
>2.1 Promoters and Enhancers  
>2.2 Splicing  
>2.3 Transcription Factors and RNA-binding Proteins
>3. Functional Genomics  
>3.1 Mutations and Functional Activities  
>3.2 Subcellular Localization
>4. Structural Genomics  
>4.1 Structural Classification of Proteins  
>4.2 Protein Secondary Structure  
>4.3 Protein Tertiary Structure and Quality Assessment  
>4.4 Contact Map
>
>|应用领域|研究对象|模型|相关文章|
>|--|--|--|--|
>|[Gene expression](https://en.wikipedia.org/wiki/Gene_expression) |Characterization||[D.Urda et al. 2017](https://link.springer.com/chapter/10.1007/978-3-319-59147-6_5), [Jie Tan et al. 2017](https://www.sciencedirect.com/science/article/pii/S2405471217302314?via%3Dihub), [Padideh Danaee et al. 2017](https://www.worldscientific.com/doi/abs/10.1142/9789813207813_0022), [Aman Gupta et al. 2015](https://ieeexplore.ieee.org/document/7359871), [Lujia Chen et al. 2016](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0852-1), [Gregory P. Way et al. 2017](https://arxiv.org/abs/1711.04828), [Ayse Dincer at al. 2018](https://www.biorxiv.org/content/10.1101/278739v1), [Hossein Sharifi-Noghabi et al. 2018](https://www.biorxiv.org/content/10.1101/276055v2), [Haohan Wang et al. 2017](https://www.inderscienceonline.com/doi/10.1504/IJDMB.2017.085711), |
>||Prediction||[Yifei Chen et al. 2016](https://academic.oup.com/bioinformatics/article/32/12/1832/1743989), [Rui Xie et al. 2017](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-017-4226-0), [Ritambhara Singh et al.](https://academic.oup.com/bioinformatics/article/32/17/i639/2450757), |
>|Regulatory Genomics|[Promoter](https://en.wikipedia.org/wiki/Promoter_(genetics))||[Ramzan Kh.Umarov et al. 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0171410), [Shashank Singh et al. 2019](https://link.springer.com/article/10.1007/s40484-019-0154-0), [Sean Whalen et al. 2016](https://www.nature.com/articles/ng.3539), |
>||[Enhancer](https://en.wikipedia.org/wiki/Enhancer_(genetics))||[Dikla Cohn et al. 2018](https://www.biorxiv.org/content/10.1101/264200v2), [Feng Liu et al. 2016](https://www.nature.com/articles/srep28517#Sec2), [Xu Min et al. 2017](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1878-3), [Bite Yang et al. 2017](https://academic.oup.com/bioinformatics/article/33/13/1930/3002766), |
>||[non-coding DNA](https://en.wikipedia.org/wiki/Non-coding_DNA)|||
>||[microRNA, miRNA](https://en.wikipedia.org/wiki/MicroRNA)|||
>||[Transcription factor, TF](https://en.wikipedia.org/wiki/Transcription_factor)||[Babak Alipanahi et al. 2015](https://www.nature.com/articles/nbt.3300), [Dexiong Chen et al. 2017](https://www.biorxiv.org/content/10.1101/217257v1), [Haoyang Zeng et al. 2016](https://academic.oup.com/bioinformatics/article/32/12/i121/2240609), [Jack Lanchantin et al. 2016](https://arxiv.org/abs/1608.03644), [Avanti Shrikumar et al. 2017](https://www.biorxiv.org/content/10.1101/103663v1), |
>||[Alternative splicing, AS](https://en.wikipedia.org/wiki/Alternative_splicing)||[Anupama Jha et al.](https://academic.oup.com/bioinformatics/article/33/14/i274/3953982)|
>||[RNA-binding protein, RBP](https://en.wikipedia.org/wiki/RNA-binding_protein)|||
>||[DNA methylation](https://en.wikipedia.org/wiki/DNA_methylation)|||
>||DNA accessbility|||
>|[Functional Genomics](https://en.wikipedia.org/wiki/Functional_genomics)|[Mutations](https://en.wikipedia.org/wiki/Mutation)||[David R. Kelley et al. 2016](https://genome.cshlp.org/content/26/7/990.long), [Jian Zhou et al. 2015](https://www.nature.com/articles/nmeth.3547), [Adam J. Riesselman et al. 2018](https://www.nature.com/articles/s41592-018-0138-4), [Gabriel E Hoffman et al. 2019](https://academic.oup.com/nar/article/47/20/10597/5572574?login=true)|
>||Subcellular Localization|||
>|[Structural Genomics](https://en.wikipedia.org/wiki/Structural_genomics)|[Protein structure](https://en.wikipedia.org/wiki/Protein_structure)||[John Jumper et al. 2021](https://www.nature.com/articles/s41586-021-03819-2)|
>||Contact map|||
>

- [Deep learning models in genomics; are we there yet? 2020](https://www.sciencedirect.com/science/article/pii/S2001037020303068)

>Major limitations of the DL models in the genomics area:
>
>1. Model interpretation (the black box)
>2. The curse of dimensionality
>3. Imbalanced classes. Transfer learning can provide a solution to tackle the class imbalanced problem.
>4. Heterogeneity of data
>5. Parameters and hyper-parameters tuning
>
>![fig1](https://ars.els-cdn.com/content/image/1-s2.0-S2001037020303068-gr2_lrg.jpg)

- [Deep learning for computational biology. 2016](https://www.embopress.org/doi/full/10.15252/msb.20156651)

>![fig1](https://www.embopress.org/cms/asset/4562d390-ceb3-483a-83b6-8e265fac3483/msb156651-fig-0001-m.jpg)

- [Interpretation of deep learning in genomics and epigenomics. 2020](https://academic.oup.com/bib/article/22/3/bbaa177/5894987)

>**Interpretability is defined as**:  
>The ability to explain or to present in understandable terms to a human. -- Doshi-Velez  
>The science of comprehending what a model did. -- Gilpin  
>which is the first step toward **explainability**.
>
>A classification of common DNN interpretation approaches
>![fig1](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bib/22/3/10.1093_bib_bbaa177/1/m_bbaa177f1.jpeg?Expires=1668848289&Signature=DSWAsKZN82AIJLzh7~AlUDd03UNQzWG~FkJB79uNQQuWp-7hrfNh26A0imSbpi6yjWOkUR0SvPMYRzMBaJDt780ZnfRRgyLqEnodK2KmvYkTy8JFS5T7pWusnuBoxwwdo3Sg5F4DGvqNh~VhcSdzBM7hrB6i4FmnVga643PzbVfObqTpLI7L8VCnsZnaP4DP36IrQvgpu4dkgPDWnep1BVxMQEDkj2AtlJwjNFVeWFaxGXUz8aSBVHd4hsPdwSnkZC7nEZSDWt3I-YGFDQhnZC9SC~56W1HPni4hKhtQ4eG3-xGU8nqD7r-NoeitoCSMJnRDMsUbl13F6oW0eWgbXQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)
>![fig2](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bib/22/3/10.1093_bib_bbaa177/1/m_bbaa177f2.jpeg?Expires=1668848289&Signature=CPisbj2ctLRUWOa1CjcFav3aGBWhZQlpGZSkHjtpnuyzO7z0wNyVq4nU5XySuk6f7ZOCNGO2BP8aTFJ-Kq5G1JI9~Hd~LYLpNBRkZRV9z~E43zi0CsYmRiDf1mVCSVM8eNqFg2lTIvguXhuY0pZTX8YD6OSPsj23zwDLRDzkxUa~0-y7yHgk-69Cv9vkt7mxaRO4h86TgAvvg7OKzX4dVYoMSC2BZFxeVuS2OzyY7whcFhxpDQ6bkfU5Ij~vm5s-kqK5po~9k3ZiuBRn7EwTkwc9QhQVQbodfkU0RK~Hv9br1zFTl046nYFMF3wb1BBt~yxEEbhFOZrlxofhfUmW4A__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

- [Application of deep learning in genomics. 2020](https://link.springer.com/article/10.1007/s11427-020-1804-5)

- [Deep learning. 2015](https://www.nature.com/articles/nature14539)

>Multilayer neural networks and backpropagation
>![fig1](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fnature14539/MediaObjects/41586_2015_BFnature14539_Fig1_HTML.jpg?as=webp)

---

### 文章阅读

[DeepSEA]
[Basset]

[DeepFIGV]

[DeepAccess]

- 基因表达

[Histone modification levels are predictive for gene expression. 2010. /Rosa Karlić, Martin Vingron/](https://doi.org/10.1073/pnas.0909344107)

[A statistical framework for modeling gene expression using chromatin features and application to modENCODE datasets. 2011. /Chao Cheng, Chong Shou & Mark Gerstein/](https://doi.org/10.1186/gb-2011-12-2-r15)

[Modeling gene expression using chromatin features in various cellular contexts. 2012. /Xianjun Dong, Ewan Birney & Zhiping Weng/](https://doi.org/10.1186/gb-2012-13-9-r53)

[DeepChrome: deep-learning for predicting gene expression from histone modifications. 2016. /Ritambhara Singh, Yanjun Qi/](https://doi.org/10.1093/bioinformatics/btw427)

>Drawbacks in the previous studies:
>
>1. they rely on multiple models to separate prediction and combinatorial analysis
>2. For input features, some of them take the average value of histone modification signal from the gene region and fail to capture the subtle differences among signal distributions of histone modifications
>
>'binning' approach:
>that is, a large region surrounding the gene transcription start site (TSS) is converted into consecutive smaller bins.
>![fig1](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/32/17/10.1093_bioinformatics_btw427/3/btw427f1p.jpeg?Expires=1669103982&Signature=nnIwBA4khPD5QrkHwi~tRap~EEJKaZF-~R6mmTL8p3Oc2lvYoHmOa74GZhOjKUPaHRdF7efb2sneFgNqIhHOSCD8ud6o6vmMDZOP6wu7Vas7kuY23AxfNglHlD38RBpRQvhWdE3zMhpC6gYKK38hWhV7i2M~0sJzgVqB0TSZ-~SW~tY1fxEnAmTSB56OUT0J4PhCdQGgUf~oXF5xYBjfb~dxt74fcdNXVEqMpaysxzeold2eca0bhz4AIJQcSHPtbLkyIF6eVE5pq~j5EWHY0s8MHocyzIsxiUfpURlHOC1d-9LOTOm2CltiA1HW8G8GGIYCu0W7OoRB6AgYBn2ILg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)
>
>[Roadmap Epigenome Project, REMC](http://www.roadmapepigenomics.org/)
>
>Model:
>![fig2](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/32/17/10.1093_bioinformatics_btw427/3/btw427f2p.jpeg?Expires=1669104557&Signature=Z1dOIXKExKMkU-VshNEFLQCjC2DFW66J24Kc1vZpqbUs1NHilurCKn-X2A~LTAKvQpESXhJTCkq7TcmtbAQ5mQ8gEWmW-ea~JcFKWz2jgF8NLA35KAna3tTUA5Mtm~virgZd6e1KyquKFV0SR9DsVg5NKcQNAJhZ6ojhavJJnNgQzh1sg6HIJYCpmhR95ONbZ82TQEmrtdrbq1Ii5dNQWJfPWFq3KEnCxdMXm6rx89mHv4Q4ZJq~hJlXmvzhY~9KDUKamhh1MzJPbFZSj7TlEbqZb1lLU9VI6B-J5RLAyEEutP36Xg~~GD5SfnTMOzBsDVxgGcl7fZDSd9lVjnZ7FA__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

[DeepDiff: DEEP-learning for predicting DIFFerential gene expression from histone modifications](https://doi.org/10.1093/bioinformatics/bty612)

>Challenges:
>
>1. Genome-wide HM signals are spatially structured and may have long-range dependency.
>2. The core aim is to understand what the relevant HM factors are and how they work together to control differential expression.
>3. Since the fundamental goal of such analysis is to understand how HMs affect gene regulation, it requires the modeling techniques to provide a degree of interpretability and allowing for automatically discovering what features are essential for predictions.
>4. There exist a small number of genes exhibiting a significant change of gene expression (differential patterns) across two human cell types like A and B. This makes the prediction task using differential gene expression as outputs much harder than predicting gene expression directly in a single condition like A alone or B alone.
>
>Input feature:
>![fig1](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/34/17/10.1093_bioinformatics_bty612/2/bty612f1.jpeg?Expires=1669105514&Signature=J9hWb1uFq0en8HGBAdYtgoXN-OtVDPg9Z0de5d90E9~uimiW8ZuMyLrEijVF~D~X0Kfj40wb9Qdc6p5eMCwStN8lGTv7JH2wPJpN6nKldS3pVJD9Hjqn5y-x3OyecWfpO2HCd0ZLbRaCGYAJJj7~U2L2abj0lq5~2qKpYH18LD-D9xtVN13rMzLDC-V5PfAdAyTkEBLzCUjIJ7~LFoja2bixUjpdgU~SUnzRGmOoRME76bOaJmMKNjPxGRl5ASb2-96XiQS~ToKawL~aoaFxARjYXbfwNuqbcR-zN84xR4T8URKeve1Elc~~KhN06RvGoxaSgObufpaIw8jV9SBrSg__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)
>Model:
>![fig2](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/34/17/10.1093_bioinformatics_bty612/2/bty612f2.jpeg?Expires=1669105750&Signature=zDq7St55FjFMmKUaMHFlL7lgb4ts67pC5tz9mIs2hgIqNR0ZJuz21t4eY3xNbxQz4DxNff55QyRVyIfM3PQwMG6p6xjQaBEOtYfV3hxTee8N82gCZR~CXQ4ytQ~QQaoBEPmp9U1TCijm4JUpa63OeW3br0-1eTRzyooUPgUrctAoUTGtvPCeTCXeN45SPsq-Eucwg0akPni5-72XWZRJsxcEmFbc8U8lTOsBHOhQKBdMuDeJZNnBsdAf6ir-6xcvXGRFaQeMWbiniXfpdv3uQqoAvo~dknEOiDSGltpuNF87jTgigohm~W0mbWUX4UR-rIDpWvFnObVXlh1jShE0nw__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

[Deep learning decodes the principles of differential gene expression](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7363043/)

[Gene Expression Classification Based on Deep Learning](https://ieeexplore.ieee.org/abstract/document/9019357)

[Gene expression inference with deep learning. 2016](https://academic.oup.com/bioinformatics/article/32/12/1832/1743989)

- 转录因子

[Improving representations of genomic sequence motifs in convolutional networks with exponential activations](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8315445/)

[Lisa: inferring transcriptional regulators through integrative modeling of public chromatin accessibility and ChIP-seq data. 2020. /Qian Qin, X. Shirley Liu/](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1934-6)

>We developed [Lisa](http://lisa.cistrome.org/) to predict the transcriptional regulators (TRs) of differentially expressed or co-expressed gene sets. Based on the input gene sets, Lisa first uses histone mark ChIP-seq and chromatin accessibility profiles to construct a chromatin model related to the regulation of these genes. Using TR ChIP-seq peaks or imputed TR binding sites, Lisa probes the chromatin models using in silico deletion to find the most relevant TRs. Applied to gene sets derived from targeted TF perturbation experiments, Lisa boosted the performance of imputed TR cistromes and outperformed alternative methods in identifying the perturbed TRs.
>
>Transcriptional regulators (TRs), which include transcription factors (TFs) and chromatin regulators (CRs), play essential roles in controlling normal biological processes and are frequently implicated in disease.
>
>??? To infer the TRs that regulate a query gene set derived from differential or correlated gene expression analyses in humans or mice.
>
>算法迭代：  
>**MARGE** builds a classifier based on H3K27ac ChIP-seq RPs from the Cistrome DB to discriminate the genes in a query differentially expressed gene set from a set of background genes.  
>**BART** extends MARGE, to predict the TRs that regulate the query gene set through an analysis of the predicted cis-regulatory elements.  
>**Lisa** (epigenetic Landscape In Silico deletion Analysis and the second descendent of MARGE), a more accurate method of integrating H3K27ac ChIP-seq and DNase-seq with TR ChIP-seq or imputed TR binding sites to predict the TRs that regulate a query gene set.
>
>![fig1](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13059-020-1934-6/MediaObjects/13059_2020_1934_Fig1_HTML.png?as=webp)
>
>Changes in H3K27ac ChIP-seq and DNase-seq associated with cell state perturbations are often a matter of degree rather than switch-like. 细胞状态转变，涉及染色质状态的改变（组蛋白修饰 or 开放程度）只是程度改变而非有无的改变。

[Modeling cis-regulation with a compendium of genome-wide histone H3K27ac profiles. 2016. /Su Wang, X. Shirley Liu/](https://genome.cshlp.org/content/26/10/1417)

>reading

[Accurate prediction of cell type-specific transcription factor binding. 2019. /Jens Keilwagen, Jan Grau/](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1614-y)

可借鉴其研究思路和写作方式。

>Prediction of cell type-specific, in vivo transcription factor binding sites is one of the central challenges in regulatory genomics. Here, we present our approach that earned a shared first rank in the “ENCODE-DREAM in vivo Transcription Factor Binding Site Prediction Challenge” in 2017. In post-challenge analyses, we benchmark the influence of different feature sets and find that chromatin accessibility and binding motifs are sufficient to yield state-of-the-art performance. Finally, we provide 682 lists of predicted peaks for a total of 31 transcription factors in 22 primary cell types and tissues and a user-friendly version of our approach, Catchitt, for download.
>
>[ENCODE-DREAM challenge](https://www.synapse.org/#!Synapse:syn6131484/wiki/402026)
>
>

[FactorNet: A deep learning framework for predicting cell type specific
transcription factor binding from nucleotide-resolution sequential data. 2018. /DanielQuang, XiaohuiXie/](https://www.sciencedirect.com/science/article/pii/S1046202318303293)

[DeepBind]