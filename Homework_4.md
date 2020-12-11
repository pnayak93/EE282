# Homework 4
## Pavan Nayak
## EE282

## Summarize Partitions of a Genome Assembly

To begin our summary, we download the gzipped chromosome fasta file from flybase.org and perform an md5sum using the script:

```
#! usr/bin/env bash

wget ftp://ftp.flybase.net/genomes/dmel/current/fasta/dmel-all-chromosome-r6.36.fasta.gz
wget ftp://ftp.flybase.net/genomes/dmel/current/fasta/md5sum.txt

grep -E chromosome md5sum.txt > check1.txt
md5sum dmel-all-chromosome-r6.36.fasta.gz > check2.txt

diff -a check1.txt check2.txt

gunzip dmel-all-chromosome-r6.36.fasta.gz

```

Once the checksum is complete, and the fasta file is unzipped, we can use bioawk to partition the file into sequences greater than 100kb, and sequences less than or equal to 100kb using the following block of code:

```

bioawk -c fastx '{ if(length($seq) > 100000) { print ">"$name; print $seq }}' dmel-all-chromosome-r6.36.fasta > greater_than_100kb.fasta

$bioawk -c fastx '{ if(length($seq) <= 100000) { print ">"$name; print $seq }}' dmel-all-chromosome-r6.36.fasta > less_than_100kb.fasta

```

Then we can summarize the partitions for number of Nucleotides, N's, and Sequences using:

```
faSize greater_than_100kb.fasta
faSize less_than_100kb.fasta

```

The resulting summary for sequences greater than 100kb look like so:

```
137547960 bases (490385 N's 137057575 real 137057575 upper 0 lower) in 7 sequences in 1 files
Total size: mean 19649708.6 sd 12099037.5 min 1348131 (4) max 32079331 (3R) median 23542271
N count: mean 70055.0 sd 92459.2
U count: mean 19579653.6 sd 12138278.9
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real

```
The resulting summary for sequences less than or equal to 100kb looks like so:

```
6178042 bases (662593 N's 5515449 real 5515449 upper 0 lower) in 1863 sequences in 1 files
Total size: mean 3316.2 sd 7116.2 min 544 (211000022279089) max 88768 (Unmapped_Scaffold_8_D1580_D1567) median 1567
N count: mean 355.7 sd 1700.6
U count: mean 2960.5 sd 6351.5
L count: mean 0.0 sd 0.0
%0.00 masked total, %0.00 masked real

```
## Plots of the following for for all sequences ≤ 100kb and all sequences > 100kb

After having partitioned the sequences into greater than and less than/equal to 100kb files, we can use the following code to print sequence length and GC content, as well as sort the length from largest to smallest to be used in the Cumulative sequence size plot

```

bioawk -c fastx ' { print length($seq) "\t" gc($seq) } ' less_than_100kb.fasta | sort -k1,1rn > len_gc_less.txt

bioawk -c fastx ' { print length($seq) "\t" gc($seq) } ' greater_than_100kb.fasta | sort -k1,1rn > len_gc_greater.txt

```

1. Sequence length distribution:

    After we have partitioned the assembly into greater than and less than/equal to 100kb and returned sequence length and gc content using bioawk (file names len_gc_greater.txt and len_gc_less.txt) we can use WinSCP or Webdrive to transfer the files locally to plot in Rstudio.

    We find the file and convert to dataframe using:
    
    ```
    library(ggplot2)
    len_gc_less <- read.delim(file.choose(), header=FALSE)
    colnames(len_gc_less) <- c("length", "gc")
    len_gc_greater <- read.delim(file.choose(), header=FALSE)
    colnames(len_gc_greater) <- c("length", "gc")
    ```
    then plot the histogram for the sequences using:

    ```
    ggplot(len_gc_less, aes(x=length))+geom_histogram(color="darkblue", fill="lightblue", bins=70)
    ```
    ```
    ggplot(len_gc_greater, aes(x=length))+geom_histogram(color="darkblue", fill="lightblue", bins=8)
    ```
    We can then use Webdrive or WinSCP to transfer the image back onto the HPC. They have been placed under /output/figures under the names:

    len_gc_greater_length_seq.png

    ![histogram of seq lengths > 100kb](/output/figures/len_gc_greater_length_seq.png)

    len_gc_less_length_seq.png

    ![histogram of seq lengths <= 100kb](/output/figures/len_gc_less_length_seq.png)


2. Sequence GC% distribution

    Since we have already created the dataframes in the previous step in R with the length and GC content, we can simply plot the histograms in the same fashion using ggplot
    ```
    ggplot(len_gc_less, aes(x=gc))+geom_histogram(color="darkblue", fill="lightblue", bins=100)
    ```
    ```
    ggplot(len_gc_greater, aes(x=gc))+geom_histogram(color="darkblue", fill="lightblue", bins=7)
    ```
    the files are saved in the /output/figures folder under the names
    len_gc_greater_gc.png

    ![histogram of seq GC content > 100kb](/output/figures/len_gc_greater_gc.png)

    len_gc_less_gc.png
    ![histogram of seq GC content <= 100kb](/output/figures/len_gc_less_gc.png)

3. Cumulative sequence size sorted from largest to smallest sequences

    Here we use the plotCDF command to get the distribution functions for each, and the commands are as follows:

    ```
    plotCDF <(cut -f 1 len_gc_less.txt) disp_less.png
    plotCDF <(cut -f 1 len_gc_greater.txt) disp_greater.png
    ```
    And the CDFs can be displayed using:

    ```
    display disp_less.png
    ```
    ![CDF of <= 100kb](/output/figures/disp_less.png)

    and
    
    ```
    display disp_greater.png
    ```
    ![CDF of > 100kb](/output/figures/disp_greater.png)

    The CDF png files are contained within the output/figures folder

# Genome Assembly

## Assemble a genome from MinION reads

to download reads, overlap reads, and construct assembly we can use the following code:

download reads:

```
wget https://hpc.oit.uci.edu/~solarese/ee282/iso1_onp_a2_1kb.fastq.gz
gunzip iso1_onp_a2_1kb.fastq.gz 
```

overlap reads after using srun to request 32 cores:

```
minimap -t 32 -Sw5 -L100 -m0  iso1_onp_a2_1kb.fastq iso1_onp_a2_1kb.fastq\
| gzip -1 > iso1.paf.gz
```

construct assembly:

```
miniasm -f iso1_onp_a2_1kb.fastq iso1.paf.gz > iso1reads.gfa
```

## Assembly Assessment

1. To calculate the N50 of the assembly, we can use the following code:

```
n50 () {
  bioawk -c fastx ' { print length($seq); n=n+length($seq); } END { print n; } ' $1 \
  | sort -rn \
  | gawk ' NR == 1 { n = $1 }; NR > 1 { ni = $1 + ni; } ni/n > 0.5 { print $1; exit; } '
}

awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' iso1reads.gfa \
| tee >(n50 /dev/stdin > n50.txt) \
| fold -w 60 \
> unitigs.fa
```
and then open the n50.txt file using
```
less n50.txt
```
to get an N50 of:

4494246

2. We want to find the length of sequences in the unitig.fa file to compare against the scaffold and contig assemblies, so we first do:

```
bioawk -c fastx ' { print length($seq) } ' unitigs.fa | sort -rn > unitig_seq_length.sizes.txt
```

To compare against the scaffold and contig assemblies from flybase we can run the following script (assuming that the plotCDF, createProject, and faSplitbyN scripts are within the current directory)


```
createProject pipeline .
cd pipeline

r6url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/215/GCA_000001215.4_Release_6_plus_ISO1_MT/GCA_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz"

trusequrl="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/705/575/GCA_000705575.1_D._melanogaster_TruSeq_synthetic_long-read_assembly/GCA_000705575.1_D._melanogaster_TruSeq_synthetic_long-read_assembly_genomic.fna.gz"

wget -O - -q $trusequrl \
| tee data/raw/ISO1.truseq.ctg.fa.gz \
| gunzip -c \
| faSize -detailed /dev/stdin \
| sort -rnk 2,2 \
| tee data/processed/ISO1.truseq.ctg.sorted.namesizes.txt \
| cut -f 2 \
> data/processed/ISO1.truseq.ctg.sorted.sizes.txt

wget -O - -q $r6url \
| tee data/raw/ISO1.r6.scaff.fa.gz \
| gunzip -c \
| tee >(faSize -detailed /dev/stdin \
        | sort -rnk 2,2 \
        | tee data/processed/ISO1.r6.scaff.sorted.namesizes.txt \
        | cut -f 2 \
        > data/processed/ISO1.r6.scaff.sorted.sizes.txt) \
| faSplitByN /dev/stdin /dev/stdout 10 \
| tee >(gzip -c > data/raw/ISO1.r6.ctg.fa.gz) \
| faSize -detailed /dev/stdin \
| sort -rnk 2,2 \
| tee data/processed/ISO1.r6.ctg.sorted.namesizes.txt \
| cut -f 2 \
> data/processed/ISO1.r6.ctg.sorted.sizes.txt
```
we then move our unitig seq length file from earlier into the pipeline\data\processed folder using

```
cd ..
mv unitig_seq_length.sizes.txt pipeline\data\processed
```
then plot all of the seq lengths for all contigs and assemblies in one CDF using:

```
plotCDF pipeline/data/processed/*.sizes.txt /dev/stdout \
| tee output/figures/CDF.png \
| display
```
In this CDF.png plot, the solid line represents the sequence lengths of the r6 contig assembly, the dashed line represents the sequence lengths of the r6 scaffold assembly, the dot-and-dashed line represents the sequence lengths of our minION assembly, and the dotted line represents the truseq synthetic longreads assembly.

The final plot is named CDF.png, in the directory :

```
/data/raw/pipeline/output/figures
```
![CDF of <= 100kb](/data/raw/pipeline/output/figures)



3. Running BUSCO on the minion assembly, we can use srun to request 32 cores, then run the following command on our unitigs.fa file, located in /data/raw

```
busco -c 32 -i unitigs.fa -l diptera_odb10 -o unitig_busco -m genome
```
this returns the output as short_summary.specific.diptera_odb10.unitig_busco.txt in the data/raw/unitig_busco folder

opening this txt file using less, we get:

```
# BUSCO version is: 4.1.4
# The lineage dataset is: diptera_odb10 (Creation date: 2020-08-05, number of species: 56, number of BUSCOs: 3285)
# Summarized benchmarking in BUSCO notation for file unitigs.fa
# BUSCO was run in mode: genome

        ***** Results: *****

        C:0.2%[S:0.2%,D:0.0%],F:2.0%,M:97.8%,n:3285
        7       Complete BUSCOs (C)
        7       Complete and single-copy BUSCOs (S)
        0       Complete and duplicated BUSCOs (D)
        66      Fragmented BUSCOs (F)
        3212    Missing BUSCOs (M)
        3285    Total BUSCO groups searched
```
If we want to compare this against the flybase contig assembly, scaffold assembly, and truseq assemblies, we can run an batch job on hpc3 using the following code, (assuming you have busco installed and your virtual environment activated):

```
#!/bin/bash

#SBATCH --job-name=fly-busco-ee282      ## Name of the job.
#SBATCH -A ecoevo282     ## account to charge
#SBATCH -p standard          ## partition/queue name
#SBATCH --nodes=1            ## (-N) number of nodes to use
#SBATCH --ntasks=1           ## (-n) number of tasks to launch
#SBATCH --cpus-per-task=32    ## number of cores the job needs
#SBATCH --error=fly-busco.err ## error log file


busco -c 31 -i /data/homezvol1/pnayak/myrepos/ee282/data/raw/pipeline/data/raw/ISO1.r6.ctg.fa -l diptera_odb10 -o flybase_ctg_busco -m genome
busco -c 31 -i /data/homezvol1/pnayak/myrepos/ee282/data/raw/pipeline/data/raw/ISO1.r6.scaff.fa -l diptera_odb10 -o flybase_scaff_busco -m genome
busco -c 31 -i /data/homezvol1/pnayak/myrepos/ee282/data/raw/pipeline/data/raw/ISO1.truseq.ctg.fa -l diptera_odb10 -o truseq_ctg_busco -m genome

```
The resultant output is a folder for each assembly.

The BUSCO score for the r6 contig assembly is as follows:

```

# BUSCO version is: 4.1.4
# The lineage dataset is: diptera_odb10 (Creation date: 2020-08-05, number of species: 56, number of BUSCOs: 3285)
# Summarized benchmarking in BUSCO notation for file /data/homezvol1/pnayak/myrepos/ee282/data/raw/pipeline/data/raw/ISO1.r6.ctg.fa
# BUSCO was run in mode: genome

        ***** Results: *****

        C:99.5%[S:99.1%,D:0.4%],F:0.2%,M:0.3%,n:3285
        3269    Complete BUSCOs (C)
        3255    Complete and single-copy BUSCOs (S)
        14      Complete and duplicated BUSCOs (D)
        5       Fragmented BUSCOs (F)
        11      Missing BUSCOs (M)
        3285    Total BUSCO groups searched

```

The BUSCO score for the r6 scaffold assembly is as follows:

```
# BUSCO version is: 4.1.4
# The lineage dataset is: diptera_odb10 (Creation date: 2020-08-05, number of species: 56, number of BUSCOs: 3285)
# Summarized benchmarking in BUSCO notation for file /data/homezvol1/pnayak/myrepos/ee282/data/raw/pipeline/data/raw/ISO1.r6.scaff.fa
# BUSCO was run in mode: genome

        ***** Results: *****

        C:99.5%[S:99.1%,D:0.4%],F:0.2%,M:0.3%,n:3285
        3268    Complete BUSCOs (C)
        3254    Complete and single-copy BUSCOs (S)
        14      Complete and duplicated BUSCOs (D)
        5       Fragmented BUSCOs (F)
        12      Missing BUSCOs (M)
        3285    Total BUSCO groups searched

```
And the BUSCO score for the truseq assembly is:

```
# BUSCO version is: 4.1.4
# The lineage dataset is: diptera_odb10 (Creation date: 2020-08-05, number of species: 56, number of BUSCOs: 3285)
# Summarized benchmarking in BUSCO notation for file /data/homezvol1/pnayak/myrepos/ee282/data/raw/pipeline/data/raw/ISO1.truseq.ctg.fa
# BUSCO was run in mode: genome

        ***** Results: *****

        C:99.1%[S:97.2%,D:1.9%],F:0.3%,M:0.6%,n:3285
        3255    Complete BUSCOs (C)
        3192    Complete and single-copy BUSCOs (S)
        63      Complete and duplicated BUSCOs (D)
        11      Fragmented BUSCOs (F)
        19      Missing BUSCOs (M)
        3285    Total BUSCO groups searched
```

This comparison shows that our minION assembly is fairly poor compared to the flybase assemblies and truseq assembly

