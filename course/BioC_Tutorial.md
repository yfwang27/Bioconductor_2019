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
    seqnames    ranges strand |     score                GC
       <Rle> <IRanges>  <Rle> | <integer>         <numeric>
  a     chr1     11-50      - |         1 0.827230074210092
  b     chr2     12-51      + |         2 0.242001093691215
  c     chr2     13-52      + |         3 0.502610171213746
  d     chr2     14-53      - |         4 0.496603355975822
  e     chr1     15-54      - |         5 0.318836323916912
  f     chr1     16-55      + |         6 0.606191824423149
  g     chr3     17-56      + |         7  0.80108821648173
  h     chr3     18-57      + |         8 0.191388507606462
  i     chr3     19-58      - |         9 0.732713744742796
  j     chr3     20-59      - |        10 0.552968603093177
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
      score                GC
  <integer>         <numeric>
a         1 0.827230074210092
b         2 0.242001093691215
c         3 0.502610171213746
d         4 0.496603355975822
e         5 0.318836323916912
f         6 0.606191824423149
g         7  0.80108821648173
h         8 0.191388507606462
i         9 0.732713744742796
j        10 0.552968603093177
```

Constructing GRanges object from data frame
========================================================


```r
mm9genes <- read.table("mm9Genes.txt",sep="\t",header=T)
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
    seqnames    ranges strand |     score                GC
       <Rle> <IRanges>  <Rle> | <integer>         <numeric>
  a     chr1     11-50      - |         1 0.827230074210092
  e     chr1     15-54      - |         5 0.318836323916912
  f     chr1     16-55      + |         6 0.606191824423149
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
  f     chr1     16-55      + |         6 0.606191824423149
  i     chr3     19-58      - |         9 0.732713744742796
  i     chr3     19-58      - |         9 0.732713744742796
  j     chr3     20-59      - |        10 0.552968603093177
  j     chr3     20-59      - |        10 0.552968603093177
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
  f     chr1     16-55      + |         6 0.606191824423149
  i     chr3     19-58      - |         9 0.732713744742796
  i     chr3     19-58      - |         9 0.732713744742796
  j     chr3     20-59      - |        10 0.552968603093177
  j     chr3     20-59      - |        10 0.552968603093177
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
  f     chr1     16-55      + |         6 0.606191824423149
  i     chr3     19-58      - |         9 0.732713744742796
  j     chr3     20-59      - |        10 0.552968603093177
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
  f     chr1     16-55      + |         6 0.606191824423149
  i     chr3     19-58      - |         9 0.732713744742796
  j     chr3     20-59      - |        10 0.552968603093177
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
  f     chr1     16-55      + |         6 0.606191824423149
  i     chr3     19-58      - |         9 0.732713744742796
  j     chr3     20-59      - |        10 0.552968603093177
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














































```
processing file: BioC_Tutorial.Rpres
Loading required package: stats4
Loading required package: BiocGenerics
Loading required package: parallel

Attaching package: 'BiocGenerics'

The following objects are masked from 'package:parallel':

    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB

The following objects are masked from 'package:stats':

    IQR, mad, sd, var, xtabs

The following objects are masked from 'package:base':

    anyDuplicated, append, as.data.frame, basename, cbind,
    colnames, dirname, do.call, duplicated, eval, evalq, Filter,
    Find, get, grep, grepl, intersect, is.unsorted, lapply, Map,
    mapply, match, mget, order, paste, pmax, pmax.int, pmin,
    pmin.int, Position, rank, rbind, Reduce, rownames, sapply,
    setdiff, sort, table, tapply, union, unique, unsplit, which,
    which.max, which.min

Loading required package: S4Vectors

Attaching package: 'S4Vectors'

The following object is masked from 'package:base':

    expand.grid

Loading required package: IRanges
Loading required package: GenomeInfoDb
Loading required package: Biostrings
Loading required package: XVector

Attaching package: 'Biostrings'

The following object is masked from 'package:base':

    strsplit

Quitting from lines 445-456 (BioC_Tutorial.Rpres) 
Error: The RangesList() constructor is defunct. Please coerce to
  IRangesList instead e.g. do 'as(list(x1, x2), "IRangesList")'
  instead of 'RangesList(x1, x2)'. Alternatively, you can use the
  IRangesList() constructor e.g. 'IRangesList(x1, x2,
  compress=FALSE)'. See '?IRangesList' for more information.
Execution halted
```
