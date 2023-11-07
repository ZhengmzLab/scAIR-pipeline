# scAIR-ScSmOP

scAIR pipeline will integrated into ScSmOP (Single cell Single Molecule Multiple Omics Pipeline).

Kai Jing, Yewen Xu, Yang Yang, Pengfei Yin, Duo Ning, Guangyu Huang, Yuqing Deng, Gengzhan Chen, Guoliang Li, Simon Zhongyuan Tian, Meizhen Zheng, **ScSmOP: a universal computational pipeline for single-cell single-molecule multiomics data analysis**, Briefings in Bioinformatics, Volume 24, Issue 6, November 2023, bbad343, https://doi.org/10.1093/bib/bbad343

-----------
scAIR read structure
RNA part:
  Read 1: 10× Barcode(16bp) - UMI(12bp) - Additional base(122bp)
  Read 2: cDNA fragment(150bp)
ATAC part:
  Read 1: DNA fragment(150bp)
  Read 2: Space(8bp) - 10× Barcode(16bp)
  Read 3: DNA fragment(150bp)
RNA 10× Barcode whitelist: 10× Genomics Single Cell Multiome ATAC + Gene Expression RNA barcode (737K-arc-v1-scrna.txt)
ATAC 10× Genomics Single Cell Multiome ATAC + Gene Expression ATAC barcode reverse complementary (reverse complementary of 737K-arc-v1-scatac.txt).

