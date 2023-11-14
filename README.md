# scAIR-pipeline is specifically coded for scAIR (single cell ATAC/Interactome/RNA) methods, based on ScSmOP.

Code: Kai Jing

---

Kai Jing, Yewen Xu, Yang Yang, Pengfei Yin, Duo Ning, Guangyu Huang, Yuqing Deng, Gengzhan Chen, Guoliang Li, Simon Zhongyuan Tian, Meizhen Zheng, **ScSmOP: a universal computational pipeline for single-cell single-molecule multiomics data analysis**, Briefings in Bioinformatics, Volume 24, Issue 6, November 2023, bbad343, https://doi.org/10.1093/bib/bbad343

-----------
## scAIR-pipeline description

### scAIR read structure
**RNA part:**

  Read 1: 10× Barcode(16bp) - UMI(12bp) - Additional base(122bp)

  Read 2: cDNA fragment(150bp)

**ATAC part:**

  Read 1: DNA fragment(150bp)

  Read 2: Space(8bp) - 10× Barcode(16bp)

  Read 3: DNA fragment(150bp)

*RNA 10× Barcode whitelist: 10× Genomics Single Cell Multiome ATAC + Gene Expression RNA barcode (737K-arc-v1-scrna.txt)
ATAC 10× Barcode whitelist: 10× Genomics Single Cell Multiome ATAC + Gene Expression ATAC barcode reverse complementary (reverse complementary of 737K-arc-v1-scatac.txt).*

### Run

Processing scAIR-seq, specify FASTQ through -1 DNA R1 -2 DNA R2 -3 DNA R3 -4 RNA R1 -5 RNA R2
```
  ~/ScSmOP/scsmop.sh -t scair -1 R1.fq.gz -2 R2.fq.gz -3 R3.fq.gz -4 R4.fq.gz -5 R5.fq.gz -l 3000 -r ~/RefGenome/refdata-gex-GRCh38-2020-A-STAR/ -b ~/RefGenome/bwa_dm3_index/dm3.fa -s ~/ScSmOP/ChromSize/hg38.size.txt
```

### Process

scAIR data were processed for RNA part and ATAC-Hi-C part separately:

**RNA-seq part**

Reads were processed using ScSmOP with parameter ‘-t scrna_10x_v3 -c 10x_scair-rna_config.json’.

**ATAC-Hi-C part**

Barcodes were identified using ScSmOP ***BARP idb*** with config file ‘10x_scair-atac-hic_config.json’. Reads with valid barcode were mapped to reference genome with ***bwa -SP5M***, contacts were identified using ***pairtools parse --add-columns algn_ref_span --min-mapq 30***. Contacts were further separated into self-ligation fragments if the genome span of contact is less than certain bp and inter-ligation fragments (long range interaction) otherwise. Fragments carrying same barcode with same start position were considered as duplications and only one was retained for downstream analysis using custom scripts ***Parsepairs.py***. ***juicertools*** was used to generate interaction heatmap for visualization. Barcodes with same order within the barcode whitelist for RNA-seq part and ATAC-Hi-C part were considered from same droplet, RNA BARP ID and ATAC-Hi-C BARP ID were correlated based on the concept, ***PickUpAtacBarcode.py***. 

### Output files

- RNAResult
  - 01.BarcodeIden
  - 02.ReadAlign
  - 04.QualityAssess
- DNAResult
  - 01.BarcodeIden
  - 02.ReadAlign
  - 03.GroupAndRefine
    - De-duplicated interaction pairs ([pairs format](https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md))
    
    SCAIR.Dedup.Pairs.gz
      ```
      ## pairs format v1.0.0
      #shape: upper triangle
      #genome_assembly: unknown
      #chromsize: 21 46709983
      #samheader: @SQ	SN:21	LN:46709983
      #samheader: @PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:/home/kasen/ScSmOP_test-main/Tools/bwa mem -SP5M -t 1 -o SCAIR_DNA.bam /home/kasen/testdatascsmop-master/hg38_chr21/genome.fa ../01.BarcodeIden/SCAIR_NDNR_1.fq ../01.BarcodeIden/SCAIR_NDNR_3.fq
      #samheader: @PG	ID:pairtools_parse	PN:pairtools_parse	CL:/home/kasen/ScSmOP_test-main/Tools/pairtools parse -c /home/kasen/testdatascsmop-master/hg38_chr21/genome.fa -o SCAIR.All.pairs.gz --cmd-out /home/kasen/ScSmOP_test-main/Tools/pigz --add-columns algn_ref_span --min-mapq 30 ../02.ReadAlign/SCAIR_DNA.bam	PP:bwa	VN:1.0.2
      #columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type ligatetype fppos1 fppos2 algn_ref_span1 algn_ref_span2 barpID dupMark
      #samheader: @PG	ID:parsepairs	PN:parsepairs	VN:0.1.1	CL:/home/kasen/ScSmOP_test-main/PythonScript/ParsePairs.py -i SCAIR.All.pairs.gz -s 3000 -3 500 -5 0 -o SCAIR
      #command: /home/kasen/ScSmOP_test-main/PythonScript/ParsePairs.py -i SCAIR.All.pairs.gz -s 3000 -3 500 -5 0 -o SCAIR
      A00174:200:HJMHMDSX7:3:1101:29957:1830|||BC:GATTCACCTTTGTTGA|||CELL_1239	21	22819050	21	23807472	-	-	UU	pair	22819050	23807472	150	145CELL_1239	nondup
      A00174:200:HJMHMDSX7:3:1101:1443:3082|||BC:CTGGCTTGATTAGCGA|||CELL_222	21	21348295	21	31768176	-	-	UU	pair	21348301	31768181	48	126	CELL_222	nondup
      A00174:200:HJMHMDSX7:3:1101:9914:4570|||BC:GGTCATTCTTATTGAC|||CELL_4389	21	17301970	21	29357629	-	-	UU	pair	17301970	29357633	150	146	CELL_4389	nondup
      A00174:200:HJMHMDSX7:3:1101:9471:25848|||BC:CCATTGTACTCATAAG|||CELL_2532	21	32623165	21	42949956	+	+	RU	pair	32623150	42949955	136	96CELL_2532	nondup
      A00174:200:HJMHMDSX7:3:1101:20455:31814|||BC:TCTTGCAACCCTCCTT|||CELL_1506	21	17124086	21	28164364	+	-	UU	pair	17124085	28164364	152	150CELL_1506	nondup
      A00174:200:HJMHMDSX7:3:1102:10285:13604|||BC:AACCCTCACGATAGAA|||CELL_7539	21	32618474	21	38573684	-	+	UU	pair	32618474	38573683	138	150CELL_7539	nondup
      A00174:200:HJMHMDSX7:3:1102:10782:15687|||BC:GTGATCCTGGCTTTCT|||CELL_5415	21	13575629	21	13579268	+	-	RU	pair	13575628	13579268	150	68CELL_5415	nondup
      A00174:200:HJMHMDSX7:3:1102:16306:31610|||BC:GATTGCTCTTAAAGCG|||CELL_2668	21	39927369	21	42740798	-	-	UU	pair	39927371	42740798	91	96CELL_2668	nondup
      A00174:200:HJMHMDSX7:3:1103:2058:4429|||BC:GGTGAGGGACCGTTGT|||CELL_32536	21	37590314	21	40957043	+	+	UU	pair	37590313	40957042	150	148CELL_32536	nondup
      A00174:200:HJMHMDSX7:3:1103:14606:7686|||BC:GTATTGAGATCCGTCA|||CELL_16175	21	18151914	21	22652368	-	+	RU	pair	18151914	22652367	150	52CELL_16175	nondup
      A00174:200:HJMHMDSX7:3:1103:29740:11819|||BC:ATGGCTAACTTCGTCA|||CELL_20564	21	9959758	21	40214454	+	+	UU	pair	9959752	40214442	58	73	CELL_20564nondup
      A00174:200:HJMHMDSX7:3:1103:28637:16830|||BC:AGACTAGTGACCGCTA|||CELL_6051	21	24122381	21	24125559	-	-	UU	pair	24122381	24125559	150	150CELL_6051	nondup
      A00174:200:HJMHMDSX7:3:1103:20175:17237|||BC:ATGACTTTGCATTTGG|||CELL_5129	21	26731638	21	28776724	+	-	UR	pair	26731626	28776724	65	88CELL_5129	nondup
      A00174:200:HJMHMDSX7:3:1103:12581:36307|||BC:GCACCTTACATGCAAC|||CELL_646	21	42932929	21	42948908	+	+	UU	pair	42932928	42948897	149	140CELL_646	nondup
      A00174:200:HJMHMDSX7:3:1104:31385:14544|||BC:CTGTTACACAAGCTCA|||CELL_1259	21	8592642	21	10624712	-	+	UU	pair	8592642	10624711	150	150	CELL_1259	nondup
      A00174:200:HJMHMDSX7:3:1104:21187:15139|||BC:AGGTTCCTGACCTTGG|||CELL_1165	21	16761781	21	18205382	-	-	UR	pair	16761781	18205382	77	103CELL_1165	nondup
      ```
    - De-duplicated fragments including pairs (marked ligation type column as pair)
    
    SCAIR.Dedup.Frags.gz
      ```
      ## pairs format v1.0.0
      #shape: upper triangle
      #genome_assembly: unknown
      #chromsize: 21 46709983
      #samheader: @SQ	SN:21	LN:46709983
      #samheader: @PG	ID:bwa	PN:bwa	VN:0.7.17-r1188	CL:/home/kasen/ScSmOP_test-main/Tools/bwa mem -SP5M -t 1 -o SCAIR_DNA.bam /home/kasen/testdatascsmop-master/hg38_chr21/genome.fa ../01.BarcodeIden/SCAIR_NDNR_1.fq ../01.BarcodeIden/SCAIR_NDNR_3.fq
      #samheader: @PG	ID:pairtools_parse	PN:pairtools_parse	CL:/home/kasen/ScSmOP_test-main/Tools/pairtools parse -c /home/kasen/testdatascsmop-master/hg38_chr21/genome.fa -o SCAIR.All.pairs.gz --cmd-out /home/kasen/ScSmOP_test-main/Tools/pigz --add-columns algn_ref_span --min-mapq 30 ../02.ReadAlign/SCAIR_DNA.bam	PP:bwa	VN:1.0.2
      #columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type ligatetype fppos1 fppos2 algn_ref_span1 algn_ref_span2 barpID dupMark
      #samheader: @PG	ID:parsepairs	PN:parsepairs	VN:0.1.1	CL:/home/kasen/ScSmOP_test-main/PythonScript/ParsePairs.py -i SCAIR.All.pairs.gz -s 3000 -3 500 -5 0 -o SCAIR
      #command: /home/kasen/ScSmOP_test-main/PythonScript/ParsePairs.py -i SCAIR.All.pairs.gz -s 3000 -3 500 -5 0 -o SCAIR
      A00174:200:HJMHMDSX7:3:1101:8874:1016|||BC:CCACTTAACTAGAAGG|||CELL_22	!	0	21	17122086	-	-	MU	self	0	17122086	130	132	CELL_22	nondup
      A00174:200:HJMHMDSX7:3:1101:31367:1078|||BC:GTAACGCGATTAAGCG|||CELL_136	21	8568825	21	8568967	+	-	UU	self	8568824	8568967	138	143	CELL_136	nondup
      A00174:200:HJMHMDSX7:3:1101:2908:1141|||BC:CTTTAGGTGCTCACTC|||CELL_201	!	0	21	32663384	-	-	MU	self	0	32663384	141	150	CELL_201	nondup
      A00174:200:HJMHMDSX7:3:1101:6994:1141|||BC:GGATATTGAGCTTATC|||CELL_205	21	8819540	21	8819794	+	-	UU	self	8819539	8819794	110	134	CELL_205	nondup
      A00174:200:HJMHMDSX7:3:1101:2772:1157|||BC:CTTGAAGGACAGCCTT|||CELL_223	!	0	21	33617487	-	+	MU	self	0	33617486	140	151	CELL_223	nondup
      A00174:200:HJMHMDSX7:3:1101:9805:1188|||BC:AGCTTACCTCAAGGTA|||CELL_289	!	0	21	12281710	-	+	NR	self	0	12281709	0	99	CELL_289	nondup
      A00174:200:HJMHMDSX7:3:1101:15121:1188|||BC:GCTACTCACACGCGGT|||CELL_295	21	10586155	21	10586517	+	-	UU	self	10586154	10586517	150	150	CELL_295	nondup
      A00174:200:HJMHMDSX7:3:1101:2392:1219|||BC:CTTGAAGGACAGCCTT|||CELL_223	!	0	21	33617487	-	+	MU	self	0	33617486	140	151	CELL_223	nondup
      A00174:200:HJMHMDSX7:3:1101:23981:1219|||BC:TAGCCACCTACCTGCG|||CELL_348	21	16333556	21	16333682	+	-	UU	self	16333555	16333682	127	127	CELL_348	nondup
      A00174:200:HJMHMDSX7:3:1101:3414:1235|||BC:ATGCTTATGCAATCGT|||CELL_356	21	9292826	21	9293164	+	-	UU	self	9292825	9293164	150	150	CELL_356	nondup
      A00174:200:HJMHMDSX7:3:1101:14714:1235|||BC:TGTCAATCTCCACAAC|||CELL_362	!	0	21	10637704	-	-	NR	self	0	10637713	0	140	CELL_362	nondup
      A00174:200:HJMHMDSX7:3:1101:17833:1251|||BC:GCGACTCCTAGTAAGT|||CELL_382	!	0	21	44610247	-	+	MU	self	0	44610246	53	47	CELL_382	nondup
      A00174:200:HJMHMDSX7:3:1101:28031:1282|||BC:GCCTGGTCTCGCTTGA|||CELL_428	!	0	21	23359301	-	-	MU	self	0	23359301	73	41	CELL_428	nondup
      A00174:200:HJMHMDSX7:3:1101:13304:1297|||BC:CGGATTAGATTGTTTG|||CELL_443	21	13850645	21	13850907	+	-	UU	self	13850644	13850915	150	142	CELL_443	nondup
      A00174:200:HJMHMDSX7:3:1101:13801:1313|||BC:TTACTCCCTATGGCTT|||CELL_465	21	34964907	21	34965028	+	-	UU	self	34964906	34965028	122	122	CELL_465	nondup
      A00174:200:HJMHMDSX7:3:1101:7536:1329|||BC:TCTCTATTGCGACTTA|||CELL_482	!	0	21	44908545	-	-	NR	self	0	44908559	0	127	CELL_482	nondup
      A00174:200:HJMHMDSX7:3:1101:7365:1376|||BC:CAAGGGTTGGATTCAT|||CELL_554	!	0	21	32622866	-	+	NU	self	0	32622854	0	50	CELL_554	nondup
      A00174:200:HJMHMDSX7:3:1101:3685:1391|||BC:ATTCCCAGATATTGGG|||CELL_579	21	46572430	21	46572806	+	-	UU	self	46572429	46572806	150	150	CELL_579	nondup
      A00174:200:HJMHMDSX7:3:1101:12961:1391|||BC:TGGTGCGGAAATTACG|||CELL_585	!	0	21	13200611	-	+	MU	self	0	13200610	150	150	CELL_585	nondup
      A00174:200:HJMHMDSX7:3:1101:7238:1407|||BC:CCTTTCTGATTGAAGC|||CELL_610	!	0	21	9105822	-	-	NR	self	0	9105829	0	44	CELL_610	nondup
      A00174:200:HJMHMDSX7:3:1101:10312:1438|||BC:CCACTTAACCAGTTGG|||CELL_652	!	0	21	45315093	-	-	MU	self	0	45315113	153	116	CELL_652	nondup
      A00174:200:HJMHMDSX7:3:1101:8115:1485|||BC:CATGCAAGATGATGAC|||CELL_740	21	19971609	21	19971925	+	-	UU	self	19971608	19971925	150	150	CELL_740	nondup
      A00174:200:HJMHMDSX7:3:1101:22381:1517|||BC:GTTGCTGTGTTAGCTG|||CELL_797	!	0	21	12281720	-	+	NU	self	0	12281719	0	58	CELL_797	nondup
      A00174:200:HJMHMDSX7:3:1101:7247:1548|||BC:TACTCACACCTTGTTG|||CELL_836	!	0	21	43000082	-	+	MU	self	0	43000081	82	84	CELL_836	nondup
      A00174:200:HJMHMDSX7:3:1101:29197:1548|||BC:CACCATATGCAAGGAC|||CELL_850	!	0	21	32020128	-	-	MR	self	0	32020128	70	71	CELL_850	nondup
      A00174:200:HJMHMDSX7:3:1101:31096:1548|||BC:AGCTTGTCTCTAACTA|||CELL_854	21	8563344	21	8563411	+	-	UU	self	8563343	8563420	68	64	CELL_854	nondup
      A00174:200:HJMHMDSX7:3:1101:13313:1564|||BC:AGGTGCACTTGAGGGC|||CELL_866	21	38569588	21	38569679	+	-	UU	self	38569587	38569679	93	92	CELL_866	nondup
      A00174:200:HJMHMDSX7:3:1101:14922:1564|||BC:TAGGAGTCTCACCCTC|||CELL_869	!	0	21	27057499	-	-	MU	self	0	27057518	138	116	CELL_869	nondup
      A00174:200:HJMHMDSX7:3:1101:12897:1595|||BC:TGGTGCGGAAATTACG|||CELL_585	!	0	21	13200611	-	+	MU	self	0	13200610	150	150	CELL_585	nondup
      A00174:200:HJMHMDSX7:3:1101:6325:1642|||BC:ATGACTTTGAAACTGG|||CELL_982	!	0	21	13156201	-	-	MU	self	0	13156201	150	150	CELL_982	nondup
      ```
    - De-duplicated fragment BAM file (including pair marked with tag LT [ligate type])

    SCAIR.PS.Dedup.bam
      ```
      A00174:200:HJMHMDSX7:3:1101:2989:1000|||BC:TAATCCGTGTGACCCG|||CELL_1    77      *       0       0       *       *       0       0       NCTGAAGGATTTTCAGCATTTGGGTTATTGAACTGATCATTTGCCCATCATTAAGAGGGTTGTCTATTAATCTTTCCGTTGAACAGTCTTTGCTTCTGGATCTTAAACTTTTAATTAAGAAAGTATGTAAGTTATTAGCTATATTGGTCA #FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFF:FFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFF  AS:i:0  XS:i:0  Yt:Z:NN LT:A:N
      A00174:200:HJMHMDSX7:3:1101:2989:1000|||BC:TAATCCGTGTGACCCG|||CELL_1    141     *       0       0       *       *       0       0       ATGTAAGATAGGATACAAGTTTGGTCTTTTCTTCTCGATTGAATATCTTTTGTATTCGAAACCGGTAAGATTGCAGATGATTTGCACATAATGGTGACTAAACCAGTTGGTAAATACCGACCTTTAACCACGAAAATGATTCTGATTCAC  FFFFF:FFFFFFFFFFFFFF:FFFFFFFFFF:FFFFFF,FFFF:FFFFF:FFFFFFFFFFFF,FF:FFFF:FFFFFFFFFF:FFFFFFF::FFFFFFFFFFFFFFF::FFFFFFFFFFFFFFFF:FFF,FF:,FFFFF:FFF:FFFFFFF  AS:i:0  XS:i:0  Yt:Z:NN LT:A:N
      ```
    - DNA Cell ID and RNA Cell ID correlate table
    
    SCAIR.RNA.DNA.correlate.tsv
      ```
      #/RNAResult/01.BarcodeIden/SCAIR.CELL	/01.BarcodeIden/SCAIR.CELL
      CELL_45056	CELL_14311
      CELL_7635	CELL_17417
      CELL_26850	-
      CELL_3038	CELL_15741
      CELL_26618	-
      -	CELL_42717
      CELL_14645	CELL_30437
      CELL_20584	-
      -	CELL_10827
      ```
      ***Note that Cell ID for RNA part and ATAC part are different, if '-' means no such Cell ID found in the type of genomic material.***
  - 04.QualityAssess

Files' detailed description not mentioned here can found at [ScSmOP wiki](https://github.com/ZhengmzLab/ScSmOP/wiki/ScSmOP-Standard-Output).

