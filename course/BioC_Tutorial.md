Bioconductor Tutorial
========================================================
css: Rpress.css
author: MRC LMS Bioinformatics Core
date:https://lmsbioinformatics.github.io/MRCLMSBioinfo/LMStraining.html
width: 1440
height: 1100
autosize: true
font-import: <link href='http://fonts.googleapis.com/css?family=Slabo+27px' rel='stylesheet' type='text/css'>
font-family: 'Slabo 27px', serif;
<!-- css:style.css -->


Overview
========================================================

* [Introduction to Bioconductor](#/Intro)
* [GenomicRanges](#/GRanges)
* [Reading Sequence alignemnts](#/BAM)
* [Annotations](#/Annotation)


Bioconductor
========================================================
id: Intro

Bioconductor (BioC) is an open source, open development software project to provide tools for the analysis and comprehension of high-throughput genomics data 

- Started in 2001
- Gained popularity through microarray analysis packages
- Colloborative effort by developers across the world
- Distributed as R packages
- Two releases every year


Objectives
========================================================
- Provide access to powerful statistical and graphical methods for analysing genomics data.
- Facilitate retrieval and integration of annotation data (GenBank, GO, Entrez Gene, PubMed).
- Allow the rapid development of extensible, interoperable, and scalable software.
- Promote high-quality documentation and reproducible research.
- Provide training in computational and statistical methods.

www.bioconductor.org
========================================================
![BioC webpage](./BioC.png)


Bioconductor Release 3.9 - it works with R version 3.6.0
========================================================

- Software (1741)

    + Provides implementation of analysis methods
    
- AnnotationData (948)

    + mapping between microarray probe, gene, pathway, gene ontology, homology and other annotations
    + Representations of GO, KEGG and other annotations, and can easily access NCBI, Biomart, UCSC and other sources
    
- ExperimentData (371)

    + code, data and documentation for specific experiments or projects
    
- Workflows (27)

    + Common analysis work flows for genomics data

Installating BioC Packages
========================================================

Installing & Using a package for R version < 3.5.0 , i.e. R.3.4.0:


```r
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
library("GenomicRanges")
```

Installing & Using a package for R version after 3.5.0, i.e. R.3.6.0:


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
library("GenomicRanges")
```

BioC packages for Sequencing Data Analysis
========================================================
class: small-code
- <b>Data structures</b>: IRanges, GenomicRanges, Biostrings, BSgenome
- <b>Input/Ouput</b>: ShortRead, Rsamtools, GenomicAlignments and rtracklayer (GTF,GFF,BED)
- <b>Annotation</b>: GenomicFeatures, BSgenome, biomaRt, TxDb.\*, org.\*
- <b>Alignment</b>: Rsubread, Biostrings
- <b>Accessing Database</b>: SRAdb & GEOquery 
- <b>ChIP-seq peak identification, motif discovery and annotation</b>: ChIPQC, chipseq, ChIPseqR, ChIPpeakAnno, DiffBind, rGADEM, BayesPeak, MotifDb, SeqLogo.
- <b>RNA-seq and Differential expression analysis</b>: Rsubread, GenomicAlignments, edgeR, DESeq2, DEXseq, goseq
- <b>SNP</b>: snpStats, SeqVarTools, GGtools
- <b>Work-flows</b>: ReportingTools, easyRNASeq, ArrayExpressHTS, oneChannelGUI


Finding Help
========================================================

- ?function 

```r
library("limma")
?lmFit
```
- Browse Vignettes browseVignettes(package="limma")
- Bioconductor Course Materials: http://www.bioconductor.org/help/course-materials/
- Bioconductor Support Forum -  https://support.bioconductor.org/

Prints version information about R and all loaded packages

```r
sessionInfo()
```


GenomicRanges
=========================================================
type:section
id: GRanges

GenomicRanges
========================================================

- Genomic Ranges provides data structure for efficiently storing genomic coordinates
    + Collection of genes coordinates
    + Transcription factor binding sites (ChIP-Seq peaks)
    + Collection of aligned sequencing reads
- Builds on top of Interval Ranges (IRanges) package and lays foundation for sequencing analysis. 
- IRanges are collection of integer interval and GenomicRanges extends IRanges by including chromosome and strand.
- Provides collection of functions for accessing and manipulating Genomic coordinates
- Use cases: Identifying TF binding overlap, counting sequencing reads overlap with a gene
- Main classes: GRanges and GRangesList



Run Length Encoding (Rle)
========================================================
- Run length encoding is a data compression technique
- Efficiently encoding the redundant information



```r
# Orginial vector
# chr1, chr2, chr2, chr2, chr1, chr1, chr3, chr3
library(GenomicRanges)
chr <- Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 2))
chr
```

```
character-Rle of length 8 with 4 runs
  Lengths:      1      3      2      2
  Values : "chr1" "chr2" "chr1" "chr3"
```

The above Rle can be interpreted as a run of length 1 of chr1, followed by run length of 3 of chr2, followed by run length of 2 of chr1 and followed by run length of 2 of chr3



Constructing GRanges object
========================================================

GRanges class represents a collection of genomic features with single start and end location on the genome.  GRanges object can be cretated using <b>GRanges</b> function.


```r
library("GenomicRanges")
gr1 <- GRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
               ranges = IRanges(start=11:20, end = 50:59, names = head(letters,10)),
               strand = Rle(c("-", "+", "-", "+", "-"), c(1,2, 2, 3, 2)),
               score = 1:10, GC = runif(10,0,1))
gr1
```

```
GRanges object with 10 ranges and 2 metadata columns:
    seqnames    ranges strand |     score                 GC
       <Rle> <IRanges>  <Rle> | <integer>          <numeric>
  a     chr1     11-50      - |         1 0.0831755038816482
  b     chr2     12-51      + |         2 0.0912680625915527
  c     chr2     13-52      + |         3  0.803204068914056
  d     chr2     14-53      - |         4  0.522232784423977
  e     chr1     15-54      - |         5  0.730534711154178
  f     chr1     16-55      + |         6  0.230501104611903
  g     chr3     17-56      + |         7  0.982684309827164
  h     chr3     18-57      + |         8  0.244172841543332
  i     chr3     19-58      - |         9  0.666269501904026
  j     chr3     20-59      - |        10  0.848706887569278
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

Constructing GRanges object
========================================================

- The coordinates in GRanges are 1-based and left-most (start of a read will always be left-most coordinate of the read regardless of which strand the read aligned to).
- Additional data stored beyond genomic coordinates (separated by “|”) are called metadata. In our case metadata column contains score and GC content. 
- Metadata columns are optional and can be extracted from GRanges object using <b>mcols</b> function.


```r
mcols(gr1)
```

```
DataFrame with 10 rows and 2 columns
      score                 GC
  <integer>          <numeric>
a         1 0.0831755038816482
b         2 0.0912680625915527
c         3  0.803204068914056
d         4  0.522232784423977
e         5  0.730534711154178
f         6  0.230501104611903
g         7  0.982684309827164
h         8  0.244172841543332
i         9  0.666269501904026
j        10  0.848706887569278
```

Constructing GRanges object from data frame
========================================================


```r
mm9genes <- read.table("./data/mm9Genes.txt",sep="\t",header=T)
head(mm9genes)
```

```
    chr     start       end strand                ens        Symbol
1  chr9 105729415 105731415      - ENSMUSG00000043719        Col6a6
2 chr18  43636450  43638450      - ENSMUSG00000043424        Gm9781
3  chr7  70356455  70358455      - ENSMUSG00000030525        Chrna7
4  chrX 100818093 100820093      - ENSMUSG00000086370 B230206F22Rik
5 chr12  82881157  82883157      - ENSMUSG00000042724        Map3k9
6  chr7 152124774 152126774      - ENSMUSG00000070348         Ccnd1
```

```r
mm9genes.GR <- GRanges(seqnames=mm9genes$chr,
                       ranges=IRanges(start=mm9genes$start,end=mm9genes$end),
                       strand=mm9genes$strand,
                       ENSID=mm9genes$ens,
                       Symbol=mm9genes$Symbol)
```

Another option: <b>makeGRangesFromDataFrame()</b> 

Converting GRanges object to data frame: <b>as.data.frame(GRanges)</b>

Constructing GRanges object from data frame
========================================================

```r
head(mm9genes.GR)
```

```
GRanges object with 6 ranges and 2 metadata columns:
      seqnames              ranges strand |              ENSID
         <Rle>           <IRanges>  <Rle> |           <factor>
  [1]     chr9 105729415-105731415      - | ENSMUSG00000043719
  [2]    chr18   43636450-43638450      - | ENSMUSG00000043424
  [3]     chr7   70356455-70358455      - | ENSMUSG00000030525
  [4]     chrX 100818093-100820093      - | ENSMUSG00000086370
  [5]    chr12   82881157-82883157      - | ENSMUSG00000042724
  [6]     chr7 152124774-152126774      - | ENSMUSG00000070348
             Symbol
           <factor>
  [1]        Col6a6
  [2]        Gm9781
  [3]        Chrna7
  [4] B230206F22Rik
  [5]        Map3k9
  [6]         Ccnd1
  -------
  seqinfo: 21 sequences from an unspecified genome; no seqlengths
```


Please try to use makeGRangesFromDataFrame to create a GRanges object from mm9genes (1/2)
========================================================

```r
?makeGRangesFromDataFrame
```

### Usage
    makeGRangesFromDataFrame(df,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)

Please try to use makeGRangesFromDataFrame to create a GRanges object from mm9genes (2/2)
========================================================

```r
mm9genes.GR2 <- makeGRangesFromDataFrame(mm9genes)
head(mm9genes.GR2)
```

```
GRanges object with 6 ranges and 0 metadata columns:
      seqnames              ranges strand
         <Rle>           <IRanges>  <Rle>
  [1]     chr9 105729415-105731415      -
  [2]    chr18   43636450-43638450      -
  [3]     chr7   70356455-70358455      -
  [4]     chrX 100818093-100820093      -
  [5]    chr12   82881157-82883157      -
  [6]     chr7 152124774-152126774      -
  -------
  seqinfo: 21 sequences from an unspecified genome; no seqlengths
```

```r
mm9genes.GR3 <- makeGRangesFromDataFrame(mm9genes,
                                         keep.extra.columns=T)
head(mm9genes.GR3)
```

```
GRanges object with 6 ranges and 2 metadata columns:
      seqnames              ranges strand |                ens
         <Rle>           <IRanges>  <Rle> |           <factor>
  [1]     chr9 105729415-105731415      - | ENSMUSG00000043719
  [2]    chr18   43636450-43638450      - | ENSMUSG00000043424
  [3]     chr7   70356455-70358455      - | ENSMUSG00000030525
  [4]     chrX 100818093-100820093      - | ENSMUSG00000086370
  [5]    chr12   82881157-82883157      - | ENSMUSG00000042724
  [6]     chr7 152124774-152126774      - | ENSMUSG00000070348
             Symbol
           <factor>
  [1]        Col6a6
  [2]        Gm9781
  [3]        Chrna7
  [4] B230206F22Rik
  [5]        Map3k9
  [6]         Ccnd1
  -------
  seqinfo: 21 sequences from an unspecified genome; no seqlengths
```


GenomicRangesList
========================================================

- To represent hierarchical structured data, ex: Exons in a transcript
- List-like data structure
- Each element of the list is GRanges instance


```r
gr2 <- GRanges(seqnames = Rle(c("chr1", "chr3","chr2", "chr1", "chr3"), c(1, 2,1, 2, 4)),
                ranges = IRanges(start=55:64, end = 94:103, names = letters[11:20]),
                strand = Rle(c("+", "-", "+", "-"), c(1, 4, 3, 2)),
                score = 1:10, GC = runif(10,0,1))

GRL <- GRangesList("Peak1" = gr1, "Peak2" = gr2)
```



Operations on GenomicRanges
========================================================

<b>Accessors:</b> seqnames, start, end, ranges, strand, width, names, mcols, length

<b>Extraction:</b> GR[i], GRL[[i]], head, tail

<b>Set operations:</b> reduce, disjoin

<b>Overlaps:</b> findOverlaps, subsetByOverlaps, countOverlaps, nearest, precede, follow

<b>Arithmetic:</b> shift, resize, distance, distanceToNearest



Operations on GenomicRanges
========================================================
<div align="center">
<img src="./GenomicRanges.jpg" alt="GenomicRanges" height="700" width="1000">
</div>


Operations on GenomicRanges
========================================================

```r
head(ranges(gr1))
```

```
IRanges object with 6 ranges and 0 metadata columns:
        start       end     width
    <integer> <integer> <integer>
  a        11        50        40
  b        12        51        40
  c        13        52        40
  d        14        53        40
  e        15        54        40
  f        16        55        40
```

```r
start(gr1)
```

```
 [1] 11 12 13 14 15 16 17 18 19 20
```

```r
width(gr1)
```

```
 [1] 40 40 40 40 40 40 40 40 40 40
```

```r
length(gr1)
```

```
[1] 10
```

Operations on GenomicRanges
========================================================
Subsetting GRanges

```r
gr1[seqnames(gr1)=="chr1"]
```

```
GRanges object with 3 ranges and 2 metadata columns:
    seqnames    ranges strand |     score                 GC
       <Rle> <IRanges>  <Rle> | <integer>          <numeric>
  a     chr1     11-50      - |         1 0.0831755038816482
  e     chr1     15-54      - |         5  0.730534711154178
  f     chr1     16-55      + |         6  0.230501104611903
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

Merge overlapping genomic ranges within the same GRanges object

```r
reduce(gr1)
```

```
GRanges object with 6 ranges and 0 metadata columns:
      seqnames    ranges strand
         <Rle> <IRanges>  <Rle>
  [1]     chr1     16-55      +
  [2]     chr1     11-54      -
  [3]     chr2     12-52      +
  [4]     chr2     14-53      -
  [5]     chr3     17-57      +
  [6]     chr3     19-59      -
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```


Finding overlapping regions
========================================================
- One of the common tasks in sequencing data analysis
- Ex: Identifying transcription factor binding sites overlap with promoters
- <b>findOverlaps</b> function finds intervals overlap between two GRanges object.
- Usage: <b>findOverlaps(query,subject)</b>

```r
gr1_overlaps <- findOverlaps(gr1,gr2,ignore.strand=F)
gr1_overlaps
```

```
Hits object with 5 hits and 0 metadata columns:
      queryHits subjectHits
      <integer>   <integer>
  [1]         6           1
  [2]         9           2
  [3]         9           3
  [4]        10           2
  [5]        10           3
  -------
  queryLength: 10 / subjectLength: 10
```

Output of findOverlaps is a 'Hits' object indicating which of the query and subject intervals overlap.

Finding overlapping regions - 1
========================================================
Convert the 'Hits' object to 2 column matrix using <b>as.matrix()</b>. Values in the first column are indices of the query and values in second column are indices of the subject.


```r
gr1_overlaps.m <- as.matrix(gr1_overlaps)
gr1[gr1_overlaps.m[,"queryHits"], ]
```

```
GRanges object with 5 ranges and 2 metadata columns:
    seqnames    ranges strand |     score                GC
       <Rle> <IRanges>  <Rle> | <integer>         <numeric>
  f     chr1     16-55      + |         6 0.230501104611903
  i     chr3     19-58      - |         9 0.666269501904026
  i     chr3     19-58      - |         9 0.666269501904026
  j     chr3     20-59      - |        10 0.848706887569278
  j     chr3     20-59      - |        10 0.848706887569278
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

Finding overlapping regions - 2.1
========================================================
Or use *queryHits()*


```r
gr1[queryHits(gr1_overlaps)]
```

```
GRanges object with 5 ranges and 2 metadata columns:
    seqnames    ranges strand |     score                GC
       <Rle> <IRanges>  <Rle> | <integer>         <numeric>
  f     chr1     16-55      + |         6 0.230501104611903
  i     chr3     19-58      - |         9 0.666269501904026
  i     chr3     19-58      - |         9 0.666269501904026
  j     chr3     20-59      - |        10 0.848706887569278
  j     chr3     20-59      - |        10 0.848706887569278
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

Finding overlapping regions - 2.2
========================================================
Or use *queryHits()*


```r
gr1[unique(queryHits(gr1_overlaps))]
```

```
GRanges object with 3 ranges and 2 metadata columns:
    seqnames    ranges strand |     score                GC
       <Rle> <IRanges>  <Rle> | <integer>         <numeric>
  f     chr1     16-55      + |         6 0.230501104611903
  i     chr3     19-58      - |         9 0.666269501904026
  j     chr3     20-59      - |        10 0.848706887569278
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

Finding overlapping regions
========================================================

<b>subsetByOverlaps</b> extracts the query intervals overlap with subject intervals

```r
subsetByOverlaps(gr1,gr2,ignore.strand=F)
```

```
GRanges object with 3 ranges and 2 metadata columns:
    seqnames    ranges strand |     score                GC
       <Rle> <IRanges>  <Rle> | <integer>         <numeric>
  f     chr1     16-55      + |         6 0.230501104611903
  i     chr3     19-58      - |         9 0.666269501904026
  j     chr3     20-59      - |        10 0.848706887569278
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

Finding overlapping regions - %in% and %over%
========================================================

Alternate ways of find overlapping regions

```r
gr1[gr1 %in% gr2]
```

```
GRanges object with 0 ranges and 2 metadata columns:
   seqnames    ranges strand |     score        GC
      <Rle> <IRanges>  <Rle> | <integer> <numeric>
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```

```r
gr1[gr1 %over% gr2]
```

```
GRanges object with 3 ranges and 2 metadata columns:
    seqnames    ranges strand |     score                GC
       <Rle> <IRanges>  <Rle> | <integer>         <numeric>
  f     chr1     16-55      + |         6 0.230501104611903
  i     chr3     19-58      - |         9 0.666269501904026
  j     chr3     20-59      - |        10 0.848706887569278
  -------
  seqinfo: 3 sequences from an unspecified genome; no seqlengths
```
Note the difference between %in% and %over%

Finding overlapping regions
========================================================
Other interesting functions: <b>nearest()</b> and <b>distanceToNearest()</b>

```r
distanceToNearest(gr1,gr2,ignore.strand=F)
```

```
Hits object with 8 hits and 1 metadata column:
      queryHits subjectHits |  distance
      <integer>   <integer> | <integer>
  [1]         1           5 |         8
  [2]         4           4 |         4
  [3]         5           5 |         4
  [4]         6           1 |         0
  [5]         7           7 |         4
  [6]         8           7 |         3
  [7]         9           3 |         0
  [8]        10           3 |         0
  -------
  queryLength: 10 / subjectLength: 10
```

precede, follow
========================================================
Find nearest range in gr2 that precede or follow each range in gr1

```r
precede(gr1,gr2)
```

```
 [1] NA NA NA NA NA  6  7  7 NA NA
```

```r
follow(gr1,gr2)
```

```
 [1]  5 NA NA  4  5 NA NA NA  9  9
```
precede, follow returns the index of previous/next ranges. Overlapping ranges are excluded

Set operations
========================================================
- Operations between individual ranges within two GRanges object
- punion, pintersect, psetdiff

Example: Finding number of overlapping bases between TFBS and promoters (using pintersect).


Counting overlapping regions
========================================================
<b>countOverlaps()</b> tabulates number of subject intervals overlap with each interval in query, ex: counting number of sequencing reads overlap with genes in RNA-Seq


```r
countOverlaps(gr1,gr2,ignore.strand=T) # note the strand!
```

```
a b c d e f g h i j 
0 0 0 0 0 1 1 2 2 2 
```

Computing Coverage
========================================================
<b>coverage</b> calculates how many ranges overlap with individual positions in the genome. <b>coverage</b> function returns the coverage as Rle instance.

```r
coverage(gr1)
```

```
RleList of length 3
$chr1
integer-Rle of length 55 with 6 runs
  Lengths: 10  4  1 35  4  1
  Values :  0  1  2  3  2  1

$chr2
integer-Rle of length 53 with 6 runs
  Lengths: 11  1  1 38  1  1
  Values :  0  1  2  3  2  1

$chr3
integer-Rle of length 59 with 8 runs
  Lengths: 16  1  1  1 37  1  1  1
  Values :  0  1  2  3  4  3  2  1
```

Coverage can be exported as BigWig, Bedgraph, Wiggle and other formats to visualise in genome browsers.

Time for Exercises!
========================================================
* [Exercises Part1](./Bioconductor_Exercises_Part1.html)
<br><br>

* [Exercises Part1 Solutions](./Bioconductor_Exercises_Part1_solutions.html)



Reading Sequence Alignments
=========================================================
type:section
id: BAM


File formats in NGS
========================================================
- FASTQ (Fasta with quality)
- SAM/BAM - Sequence Alignment/Map format
- BED
- Wiggle/bedgraph

* [Genomic File Formats](https://lmsbioinformatics.github.io/LMS_genomic_formats/)

<b>FASTQ</b>

```r
@ERR590398.1 HWI-ST1146:148:C2FVTACXX:8:1101:1172:2059
NAAAATGCATATTCCTAGCATACTTCCCAAACATACTGAATTATAATCTC
+ERR590398.1 HWI-ST1146:148:C2FVTACXX:8:1101:1172:2059
A1BDFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJIJJ
```


BAM/SAM
========================================================
BAM/SAM consist two sections: Header and Alignment

Header section: 

Meta data (reference genome, aligner), starts with “@”

BAM/SAM
========================================================
Alignment section:

<b>QNAME:</b> ID of the read (“query”)<br>
<b>FLAG:</b>  alignment flags<br>
<b>RNAME:</b> ID of the reference (typically: chromosome name)<br>
<b>POS:</b>   Position in reference (1-based, left side)<br>
<b>MAPQ:</b>  Mapping quality (as Phred score)<br>
<b>CIGAR:</b> Alignment description (mismatch, gaps etc.)<br>
<b>RNEXT:</b> Mate/next read reference sequence name<br>
<b>MPOS:</b>  Mate/next read position<br>
<b>TLEN:</b>  observed Template Length<br>
<b>SEQ:</b>   sequence of the read<br>
<b>QUAL:</b>  quality string of the read<br>




Reading Sequence alignments (BAM/SAM)
========================================================
Methods for reading BAM/SAM

- <b>readAligned</b> from ShortRead package 
    – Accept multiple formats – BAM, export
    - Reads all files in a directory
    - Reads base call qualities, chromosome, position, and strand
- <b>scanBam</b> from Rsamtools package
    - scanBam reads BAM files into list structure
    - Options to select what fields and which records to import using <b>ScanBamParam</b>
- <b>readGAlignments</b> from GenomicAlignments package


Reading Sequence alignments (BAM/SAM)
========================================================

```r
library("Rsamtools")
BamFile <- system.file("extdata", "ex1.bam", package="Rsamtools")
which <- IRangesList(seq1=IRanges(1000, 2000),seq2=IRanges(c(100, 1000), c(1000, 2000)))
what <- c("rname", "strand", "pos", "qwidth", "seq")
param <- ScanBamParam(which=which, what=what)
bamReads <-  scanBam(BamFile, param=param)

length(bamReads) # Each element of the list corresponds to a range speciﬁed by the which argument in ScanBamParam
```

```
[1] 3
```

```r
names(bamReads[[1]]) # elements speciﬁed by the what and tag arguments to ScanBamParam
```

```
[1] "rname"  "strand" "pos"    "qwidth" "seq"   
```



Reading Sequence alignments (BAM/SAM)
========================================================

<b>ScanBamParam:</b>


```r
# Constructor
ScanBamParam(flag = scanBamFlag(), what = scanBamWhat(), which=which)

# Constructor helpers
scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA, 
    hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
    isFirstMateRead = NA, isSecondMateRead = NA, isSecondaryAlignment = NA,
    isNotPassingQualityControls = NA, isDuplicate = NA)
```

Reading Sequence alignments (BAM/SAM)
========================================================

- GenomicRanges package defines the GAlignments class – a specialised class for storing set of genomic alignments (ex: sequencing data) 
- Only BAM support now – future version may include other formats
- The readGAlignments function takes an additional argument, <b>param</b> allowing the user to customise which genomic regions and which fields to read from BAM
-<b>param</b> can be constructed using </b>ScanBamParam</b> function


```r
library(GenomicAlignments)
SampleAlign <- readGAlignments(BamFile)
```



Reading Sequence alignments (BAM/SAM)
========================================================
We can also customise which regions to read

```r
region <- IRangesList(seq1=IRanges(1000, 2000),seq2=IRanges(1000, 2000))
param1 <- ScanBamParam(what=c("rname", "pos", "cigar","qwidth"),which=region)
SampleAlign1 <- readGAlignments(BamFile,param=param1)
```



Time for Exercises!
========================================================
* [Exercises Part2](./Bioconductor_Exercises_Part2.html)
<br><br>

* [Exercises Part2 Solutions](./Bioconductor_Exercises_Part2_solutions.html)



Annotation Packages
=========================================================
type:section
id: Annotation

BioC Annotation Packages
========================================================
Annotation packages can be broadly classified in to gene-centric and genome-centric.

<b>Gene-centric annotation packages (AnnotationDbi):</b>

-  Organism level packages: contains gene annotation for entire organism. Follows “org.XX.YY.db” pattern (Ex: org.Hs.eg.db)
- General System biology data: KEGG.db (association between pathways and genes), GO.db (Gene ontology term and genes) and ReactomeDb
- Platform level packages: Annotation for a specific platform (ex: hgu133a.db for Affymetrix HGU133A microarray). 




BioC Annotation Packages
========================================================

<b>Genomic-centric GenomicFeatures packages:</b>

- TranscriptDB (TxDB) packages contains genomic coordiantes for transcripts specific to a genome build, ex: TxDb.Hsapiens.UCSC.hg19.knownGene. These packages allow access to various features on transcriptome, including exons, genes and transcripts coordinates

<b>Web-based annotation services:</b>

- biomaRt provides interface to query web-based `biomart' resource for genes, sequence, SNPs, and etc. 


AnnotationDbi Accessor Functions
========================================================
- <b>columns </b>What kind of annotation available in AnnotationDb object.
- <b>keytypes </b>Displays which type of identifiers can be passed in to <b>select</b> function. 
- <b>keys </b> returns keys (index) for the database contained in the AnnotationDb object. Used along with <b>keytypes</b> in <b>select</b> function to retrieve interested annotation
- <b>select</b> will retrieve the annotation data as a data.frame based on the supplied keys, keytypes and columns. 

<b>Note:</b> Package name = Annotation object 

We will explore how to retrieve annotation from gene-centric organism level annotation package (org.Hs.eg.db)



Accessing annotation from org.Hs.eg.db
========================================================

Load the package and list the contents

```r
library("org.Hs.eg.db")
```

```r
columns(org.Hs.eg.db)
 [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
 [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
 [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
[13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
[17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
[21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
[25] "UNIGENE"      "UNIPROT"     
```

To know more about the above identifier types

```r
help(SYMBOL)
```

Accessing annotation from org.Hs.eg.db
========================================================

Which keytypes can be used to query this database? <b>keytypes</b> (What is the difference between <b>columns</b> and <b>keytypes</b>?)

```r
keytypes(org.Hs.eg.db)
 [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT" 
 [5] "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"    
 [9] "EVIDENCEALL"  "GENENAME"     "GO"           "GOALL"       
[13] "IPI"          "MAP"          "OMIM"         "ONTOLOGY"    
[17] "ONTOLOGYALL"  "PATH"         "PFAM"         "PMID"        
[21] "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
[25] "UNIGENE"      "UNIPROT"     
```



Accessing annotation from org.Hs.eg.db
========================================================
If we want to extract few identifiers of a particular keytype, we can use <b>keys</b> function

```r
head(keys(org.Hs.eg.db, keytype="SYMBOL"))
[1] "A1BG"  "A2M"   "A2MP1" "NAT1"  "NAT2"  "NATP" 
```

We can extract other annotations for a particular identifier using <b>select</b> function

```r
select(org.Hs.eg.db, keys = "A1BG", keytype = "SYMBOL", columns = c("SYMBOL", "GENENAME", "ENSEMBL") )
  SYMBOL               GENENAME         ENSEMBL
1   A1BG alpha-1-B glycoprotein ENSG00000121410
```


Annotating results - Example
========================================================
How can we annotate our results (ex: RNA-Seq differential expression analysis results)? 

First, we will load an example results:

```r
load(system.file("extdata", "resultTable.Rda", package="AnnotationDbi"))
head(resultTable)
             logConc     logFC LR.statistic       PValue          FDR
100418920  -9.639471 -4.679498     378.0732 3.269307e-84 2.613484e-80
100419779 -10.638865 -4.264830     291.1028 2.859424e-65 1.142912e-61
100271867 -11.448981 -4.009603     222.3653 2.757135e-50 7.346846e-47
100287169 -11.026699 -3.486593     206.7771 6.934967e-47 1.385953e-43
100287735 -11.036862  3.064980     204.1235 2.630432e-46 4.205535e-43
100421986 -12.276297 -4.695736     190.5368 2.427556e-43 3.234314e-40
```
Rownames of the above dataframe are "Entrez gene identifiers" (human). We will extract gene symbol for these Entrez gene identifiers from org.Hs.eg.db package using select fucntion.




Annotating results - Example
========================================================

```r
SYM <- select(org.Hs.eg.db, keys = rownames(resultTable), keytype = "ENTREZID", columns = "SYMBOL")
head(SYM)
```

```
   ENTREZID       SYMBOL
1 100418920     ACTG1P13
2 100419779 LOC100419779
3 100271867      MPVQTL1
4 100287169         <NA>
5 100287735      TTTY13B
6 100421986         <NA>
```

```r
NewResult <- merge(resultTable,SYM,by.x=0,by.y=1)
head(NewResult)
```

```
  Row.names   logConc     logFC LR.statistic       PValue          FDR
1 100127888 -10.57050  2.758937     182.8937 1.131473e-41 1.130624e-38
2 100131223 -12.37808 -4.654318     179.2331 7.126423e-41 6.329847e-38
3 100271381 -12.06340  3.511937     188.4824 6.817155e-43 7.785191e-40
4 100271867 -11.44898 -4.009603     222.3653 2.757135e-50 7.346846e-47
5 100287169 -11.02670 -3.486593     206.7771 6.934967e-47 1.385953e-43
6 100287735 -11.03686  3.064980     204.1235 2.630432e-46 4.205535e-43
        SYMBOL
1  SLCO4A1-AS1
2 LOC100131223
3      RPS28P8
4      MPVQTL1
5         <NA>
6      TTTY13B
```


TranscriptDB (TxDb) packages
========================================================
- TxDb packages provide access genomic coordinates to various transcript-related features from UCSC and Biomart data sources.
- TxDb objects contains relationship between mRNA transcripts, exons, CDS and their associated identifiers
- TxDb packages follows specific naming scheme, ex: TxDb.Mmusculus.UCSC.mm9.knownGene

We will explore TxDb package for Mouse mm9 genome from UCSC. We will first install the package and load in to our current working space.

To install this package, before R < (version "3.5") and enter:


```r
source("http://bioconductor.org/biocLite.R")
biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene")
```

To install this package, start R (version "3.6") and enter:


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("TxDb.Mmusculus.UCSC.mm9.knownGene")
```

TranscriptDb - TxDb.Mmusculus.UCSC.mm9.knownGene
========================================================
By default, the annotation object will have same name as package name. Create an alias for convenience.

```r
library("TxDb.Mmusculus.UCSC.mm9.knownGene")
txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
```

GenomicFeatures also provides functions to create TxDb objects directly from UCSC and Biomart databases: <b> makeTxDbFromBiomart</b> and <b>makeTxDbFromUCSC</b>


```r
library(GenomicFeatures)
txdb <- makeTxDbFromUCSC(genome = "mm9", tablename = "knownGene")
```



Since TxDb are inherited from AnnotationDb object, we can use <b>columns, keys, select,</b> and <b>keytypes</b> functions.


```r
keys <- c("100009600", "100009609", "100009614")
select(txdb,keys=keys,columns=c("GENEID","TXNAME"),keytype="GENEID")
     GENEID     TXNAME
1 100009600 uc009veu.1
2 100009609 uc012fog.1
3 100009614 uc011xhj.1
```

TranscriptDb - TxDb.Mmusculus.UCSC.mm9.knownGene
========================================================
Most common operations performed on TxDb objects are retrieving exons, transcripts and CDS genomic coordinates.  The functions <b>genes, transcripts, exons</b>, and <b>cds</b> return the coordinates for the group as GRanges objects.


```r
TranscriptRanges <- transcripts(txdb)
TranscriptRanges[1:3]
GRanges object with 3 ranges and 2 metadata columns:
      seqnames          ranges strand |     tx_id     tx_name
         <Rle>       <IRanges>  <Rle> | <integer> <character>
  [1]     chr1 4797974-4832908      + |         1  uc007afg.1
  [2]     chr1 4797974-4836816      + |         2  uc007afh.1
  [3]     chr1 4847775-4887990      + |         3  uc007afi.2
  -------
  seqinfo: 35 sequences (1 circular) from mm9 genome
ExonRanges <- exons(txdb)
ExonRanges[1:2]
GRanges object with 2 ranges and 1 metadata column:
      seqnames          ranges strand |   exon_id
         <Rle>       <IRanges>  <Rle> | <integer>
  [1]     chr1 4797974-4798063      + |         1
  [2]     chr1 4798536-4798567      + |         2
  -------
  seqinfo: 35 sequences (1 circular) from mm9 genome
```



TxDb.Mmusculus.UCSC.mm9.knownGene
========================================================
TxDb package also provides interface to discover how genomic features are related to each other. Ex:  Access all transcripts or exons associated to a gene. Such grouping can be achieved by <b>transcriptsBy</b>, <b>exonsBy</b>, and <b>cdsBy</b> functions. The results are returned as GRangesList objects.


```r
Transcripts <- transcriptsBy(txdb, by = "gene")
Transcripts[1:2]
GRangesList object of length 2:
$100009600 
GRanges object with 1 range and 2 metadata columns:
      seqnames            ranges strand |     tx_id     tx_name
         <Rle>         <IRanges>  <Rle> | <integer> <character>
  [1]     chr9 20866837-20872369      - |     28943  uc009veu.1

$100009609 
GRanges object with 1 range and 2 metadata columns:
      seqnames            ranges strand | tx_id    tx_name
  [1]     chr7 92088679-92112519      - | 23717 uc012fog.1

-------
seqinfo: 35 sequences (1 circular) from mm9 genome
```

TranscriptDb - TxDb.Mmusculus.UCSC.mm9.knownGene
========================================================

Other interesting functions: <b>intronsByTranscript, fiveUTRsByTranscript</b> and <b>threeUTRsByTranscript</b>

GenomicFeatures also provides functions to create TxDb objects directly from UCSC and Biomart databases: <b> makeTxDbFromBiomart</b> and <b>makeTxDbFromUCSC</b>

```r
USCmm9KnownGene <- makeTxDbFromUCSC(genome = "mm9", tablename = "knownGene")
```

Save the annotation object and label them appropriately to facilitate reproducible research:

```r
saveDb(txdb,file="Mouse_ucsc_mm9_20190919.sqlite")
txdb <- loadDb("Mouse_ucsc_mm9_20190919.sqlite")
```

Another way of creating TxDb: </b>makeTxDbFromGFF()</b>

Annotations from the web – biomaRt
========================================================
biomaRt package oﬀers access to biomart based online annotation resources (marts). Each mart has several datasets. <b>getBM</b> function can be used to retrieve annotation from the biomarts. Use the following functions to find values for the arguments in <b>getBM</b> 

<b>listMarts():</b>  list the available biomart resources<br>
<b>useMart():</b>    select the mart<br>
<b>listDatasets():</b>  available dataset in the select biomart resource<br>
<b>useDataset():</b> select a dataset in the select mart<br>
<b>listAttributes():</b> available annotation attributes for the selected dataset<br>
<b>listFiltersList():</b> available filters for the selected dataset<br>

biomaRt - Example
========================================================


```r
library("biomaRt")
marts <- listMarts(host="www.ensembl.org") # List available marts
marts[1,]
```

```
               biomart          version
1 ENSEMBL_MART_ENSEMBL Ensembl Genes 97
```

```r
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", host = "sep2015.archive.ensembl.org") # select ensembl
ens_datasets <- listDatasets(ensembl) # list datasets 
ens_human <- useDataset("hsapiens_gene_ensembl",mart=ensembl) # select human dataset
ens_human_Attr <- listAttributes(ens_human) # list available annotation
ens_human_filters <- listFilters(ens_human) # list availabel filters
```

biomaRt - Example
========================================================

Extract genomic coordinates, ensembl gene identifier and gene symbol for genes in chromsome X

```r
chrXGenes <- getBM(attributes = c("chromosome_name","start_position","end_position","ensembl_gene_id","strand","external_gene_name"), filter="chromosome_name",values="X", mart=ens_human)
head(chrXGenes)
```

```
  chromosome_name start_position end_position ensembl_gene_id strand
1               X        3717149      3717600 ENSG00000182888      1
2               X      154620096    154620171 ENSG00000273773      1
3               X      154622428    154623500 ENSG00000213338     -1
4               X       68070520     68071767 ENSG00000228328     -1
5               X      154632470    154633182 ENSG00000272681     -1
6               X       98251679     98252543 ENSG00000236243      1
  external_gene_name
1      RP11-558O12.1
2       RP5-865E18.2
3             ATF4P1
4      RP11-516A11.1
5            FAM223A
6            RPL6P29
```


Session Info
========================================================

```r
sessionInfo()
```

```
R version 3.6.1 (2019-07-05)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS High Sierra 10.13.6

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] biomaRt_2.40.4                         
 [2] TxDb.Mmusculus.UCSC.mm9.knownGene_3.2.2
 [3] GenomicFeatures_1.36.4                 
 [4] org.Hs.eg.db_3.8.2                     
 [5] AnnotationDbi_1.46.1                   
 [6] GenomicAlignments_1.20.1               
 [7] SummarizedExperiment_1.14.1            
 [8] DelayedArray_0.10.0                    
 [9] BiocParallel_1.18.1                    
[10] matrixStats_0.54.0                     
[11] Biobase_2.44.0                         
[12] Rsamtools_2.0.0                        
[13] Biostrings_2.52.0                      
[14] XVector_0.24.0                         
[15] GenomicRanges_1.36.0                   
[16] GenomeInfoDb_1.20.0                    
[17] IRanges_2.18.1                         
[18] S4Vectors_0.22.0                       
[19] BiocGenerics_0.30.0                    
[20] knitr_1.24                             

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.2             compiler_3.6.1         pillar_1.4.2          
 [4] prettyunits_1.0.2      progress_1.2.2         bitops_1.0-6          
 [7] tools_3.6.1            zlibbioc_1.30.0        digest_0.6.20         
[10] zeallot_0.1.0          bit_1.1-14             memoise_1.1.0         
[13] evaluate_0.14          RSQLite_2.1.2          tibble_2.1.3          
[16] lattice_0.20-38        pkgconfig_2.0.2        rlang_0.4.0           
[19] Matrix_1.2-17          DBI_1.0.0              curl_4.0              
[22] xfun_0.9               GenomeInfoDbData_1.2.1 rtracklayer_1.44.2    
[25] httr_1.4.1             stringr_1.4.0          hms_0.5.0             
[28] vctrs_0.2.0            bit64_0.9-7            grid_3.6.1            
[31] R6_2.4.0               XML_3.98-1.20          blob_1.2.0            
[34] magrittr_1.5           backports_1.1.4        assertthat_0.2.1      
[37] stringi_1.4.3          RCurl_1.95-4.12        crayon_1.3.4          
```

Time for Exercises!
========================================================
* [Exercises Part3](./Bioconductor_Exercises_Part3.html)
<br><br>

* [Exercises Part3 Solutions](./Bioconductor_Exercises_Part3_solutions.html)


