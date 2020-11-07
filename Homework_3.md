# Homework 3
## Pavan Nayak
## EE282


## Summarize Genome Assembly
### File Integrity of gzipped fasta file

Use the following command to download the chromosome fasta file from Flybase:

```
wget ftp://ftp.flybase.net/genomes/dmel/current/fasta/dmel-all-chromosome-r6.36.fasta.gz
```

Next, we look at the md5sum.txt file stored on the flybase download page showing:

fbd2855a20c3610050ff249dd975821 dmel-all-chromosome-r6.36.fasta.gz

upon downloading this file, we use the command:

```
md5sum dmel-all-chromosome-r6.36.fasta.gz
```
and obtain the same md5sum number:


` fbd2855a20c3610050ff249dd975821 dmel-all-chromosome-r6.36.fasta.gz `

This validates the file integrity of the gzipped fasta file after download


### Calculate Summaries of the Genome

We use the commands:

```
 gunzip dmel-all-chromosome-r6.36.fasta.gz
 faSize dmel-all-chromosome-r6.36.fasta
 ```

 and return:

```
143726002 bases (1152978 N's 142573024 real 142573024 upper 0 lower) in 1870 sequences in 1 files
Total size: mean 76858.8 sd 1382100.2 min 544 (211000022279089) max 32079331 (3R) median 1577
N count: mean 616.6 sd 6960.7
U count: mean 76242.3 sd 1379508.4
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real

```
Hence, the total number of nucleotides (bases) in the chromosome fasta file is 143726002, the number of N's is 1152970, and the number of sequences is 1870

Alternatively, we can use the h3_fasta_summary.sh script to perform the entire pipeline in one go, where the md5sum check is performed by the diff -a function. This function returns no output if files are identical.

```
#! usr/bin/env bash

wget ftp://ftp.flybase.net/genomes/dmel/current/fasta/dmel-all-chromosome-r6.36.fasta.gz
wget ftp://ftp.flybase.net/genomes/dmel/current/fasta/md5sum.txt

grep -E chromosome md5sum.txt > check1.txt
md5sum dmel-all-chromosome-r6.36.fasta.gz > check2.txt

diff -a check1.txt check2.txt

gunzip dmel-all-chromosome-r6.36.fasta.gz
faSize dmel-all-chromosome-r6.36.fasta

```

## Summarize an Annotation
### File Integrity
We can download the gtf file using:

`wget ftp://ftp.flybase.net/genomes/dmel/current/gtf/dmel-all-r6.36.gtf.gz `

and using md5sum compare it against the provided md5sum on the flybase website.

using md5sum on dmel-all-r6.36.gtf.gz returns:

`9085d2f3d2449fc1f6159015511240b8  dmel-all-r6.36.gtf.gz`

and when comparing to the md5sum.txt file provided, it is the same:

`9085d2f3d2449fc1f6159015511240b8`


### Compile a Report Summarizing the Annotation

#### Total number of features of each type, sorted from most common to the least common

Total number of features can be determined using a pipeline of commands as follows. to unzip the gtf file after downloading, we use:

`gunzip dmel-all-r6.36.gtf.gz`

and on the resultant file we can use the following pipeline:

```
gawk -F '\t' '{print $3}' dmel-all-r6.36.gtf | sort | uniq -c | sort -k1,1nr > unique_col3_sorted.txt
```
This pipeline prints the 3rd column of the gtf file (the column of features) and sorts them, then counts the number of unique values, and finally performs a sort in reverse numerical order (i.e. descending order) on the new first column (the number of counts). Finally the output is placed in a file called unique_col3_sorted.txt

The resulting summary is as follows:

```
189268 exon
162578 CDS
46664 5UTR
33629 3UTR
30812 start_codon
30754 stop_codon
30728 mRNA
17875 gene
3047 ncRNA
485 miRNA
366 pseudogene
312 tRNA
300 snoRNA
262 pre_miRNA
115 rRNA
32 snRNA

```
#### Total number of genes per chromosome arm (X, Y, 2L, 2R, 3L, 3R, 4)

This can be achieved using a similar pipeline on the gtf file, with some modifications:

```
$gawk -F '\t' '$3=="gene"{print $1}'  dmel-all-r6.36.gtf | grep -P '^[4XY23][LR]?$' | sort | uniq -c > chromosome_gene_count.txt
```

This pipeline prints only the rows in column 1 (column of Chromosomes) where the respective values in column 3 equal exactly the string "gene". This filters out all other values, including pseudogenes. Next, grep is used to find only values with 4, X, Y, 2, or 3 in the first position, and optionally L or R in the second position with an end of line anchor tag. Lastly the resultant output is sorted and unique values are found to give the count of genes per chromosome arm and is output to a file called chromosome_gene_count.txt.

The resulting output is:

```

   3516 2L
   3653 2R
   3486 3L
   4225 3R
    114 4
   2691 X
    113 Y

```

Alternatively, you can use the script provided below and in the 'code' directory titled h3_gtf_summary.sh to perform all of this in one go:

```
#! usr/bin/env bash

wget ftp://ftp.flybase.net/genomes/dmel/current/gtf/dmel-all-r6.36.gtf.gz
wget ftp://ftp.flybase.net/genomes/dmel/current/gtf/md5sum.txt

md5sum dmel-all-r6.36.gtf.gz > check3.txt

diff -a check3.txt md5sum.txt

gunzip dmel-all-r6.36.gtf.gz

gawk -F '\t' '{print $3}' dmel-all-r6.36.gtf | sort | uniq -c | sort -k1,1nr > unique_col3_sorted.txt

gawk -F '\t' '$3=="gene"{print $1}'  dmel-all-r6.36.gtf | grep -P '^[4XY23][LR]?$' | sort | uniq -c > chromosome_gene_count.txt


```

