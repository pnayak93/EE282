# Final Project Writeup
## Pavan Nayak
## EE282

# Introduction

Upon completion of almost any Next Generation Sequencing (NGS) experiment, scientists are usually presented with huge amounts of data. These data often present themselves in the form of vast gene lists, for example differential expression analyses in RNA-seq experiments, from which it is difficult to quickly and effectively derive actionable, testable hypotheses in the lab. Many tools exist to try and solve this problem, such as KEGG pathways and GO-term analysis, but often these analyses can lead researchers to vague results and unclear directions as to which genes to select for practical drill-down and further functional study, especially when such a large group of seemingly unrelated or understudied genes in the given context can arise from the experiment itself. The simplest first step to begin to understand whether a single given gene in a gene list from an NGS experiment output is important for the context that a researcher is working on, would be to simply search up the gene and the context (e.g. ‘klf2a heart valve development’) on a research database, usually something like NCBI Pubmed. This is a simple task, albeit slightly time-consuming, when presented with one, ten, or even fifty genes, but when an NGS experiment outputs genes in the range of hundreds to thousands, the problem begs for a simple solution. To help address this problem, I have created a simple group of scripts that:

1) (Written in Bash) Performs a brief analysis of read quality of fastq files by utilizing a bash script which identifies and returns the total number of reads, total number of poor quality reads as defined by reads containing 10 consecutive N’s, and the percentage of poor quality reads in each fastq file.

2) (Written in Python) Takes lists of gene names (either conventional names or ENSEMBL-ID) in excel file format that are output from a typical RNA-seq or other NGS experiment, and searches these genes, along with a user input search term, against the pubmed research database, and returns the number of article hits for each search in an excel file.

3) (Written in R) Takes the output from script 2 as input to create a bar plot of the genes to get a better visual representation of which gene search hits on pubmed.


# Methods

The script pipeline begins with a brief read-quality analysis, written in Bash. For each fastq file in the working directory, the total number of reads is calculated using the cat command piped into the “wc -l “ and divided by 4, to account for fastq file format showing one read sequence per four lines. This calculation is done using a pipe into the bc calculator tool in linux. Next, the number of reads with 10 consecutive N’s is calculated using the grep command and piped into another wc -l. Lastly the percentage of bad reads is calculated by assigning each of the two previously calculated values to variables, then dividing the number of bad reads by the total reads and multiplied by 100, again piped into the bc calculator. These values are appended to the end of an output .txt file. The script repeats for each fastq file in the folder by using a for-loop and wildcard expression (*.fastq). This script is entitled run_all_reads.sh and is placed within the code/scripts directory. The script needs to be placed in the directory with all fastq files and can be run using the command:

```
sh run_all_reads.sh
```
The script itself looks as follows:

```
#!/bin/bash
for filename in *.fastq
do
echo  ${filename} >>number_of_Ns.txt
echo "Total number of reads:">>number_of_Ns.txt
T=`echo $(cat ${filename} |wc -l)/4|bc`
echo $T>>number_of_Ns.txt
echo "Number of reads with 10 consecutive N's">>number_of_Ns.txt
N=`echo $(grep NNNNNNNNNN ${filename} | wc -l)`
echo $N>>number_of_Ns.txt
echo "percentage of bad reads in ${filename}:">>number_of_Ns.txt
P=`echo $(echo "scale=7;($N/$T)*100"|bc)`
echo $P>>number_of_Ns.txt
done
```

The next script in the pipeline is a python command-line program, entitled ncbi-scraper-with-ensembl-id-input.py (in the code/scripts folder) which takes an excel file containing a list of gene names, either conventional name or ensembl id, and searches each gene along with a user-input search term against the pubmed database. This is done using multiple python packages and some basic scripting. In brief, the program first reads in user specified path to the input file using sys.argv values inbuilt in python, and converts the excel file into a dataframe, then into a list for easier manipulation. The program uses conditionals to check if the genes in the dataframe are in ENSEMBL-ID form or not, then continues onwards. If the list is entirely ENSEMBL-ID, then the program first loops through the list, and does a GET HTTP request, using the python requests package, against the ENSEMBL REST API (http://rest.ensembl.org/) to return the conventional gene name in another column next to the gene list at the correct respective position. The program then uses the Entrez eutils API (https://www.ncbi.nlm.nih.gov/books/NBK25501/) to perform GET requests for each gene along with the search term. If the list of genes are already in the conventional gene name format, the program skips to the second step automatically. This dataframe is finally converted back into an excel file and, after some formatting changes, is output. For aesthetic/user interface purposes, the tqdm package was used to wrap all for-loops, providing a progress bar for the program’s run. The entire program is converted into a simple command line interface with “help” details for each argument by using the argparse python package. Assuming all dependencies are installed in your virtual environment (bs4, pandas, requests, tqdm, xlrd, openpyxl), the script can be run (with my sample data sets) with the commands:

```
python ncbi-scraper-with-ensembl-id-input.py gene-list-top-499.xlsx tendon "output-file-name.xlsx"
```
for a sample gene list input already in conventional gene format and:

```
python ncbi-scraper-with-ensembl-id-input.py test-ens-gene-list-short2.xlsx tendon "output-file-name2.xlsx"
```
for a sample gene list in ENSEMBL ID format.

After the excel file is output from the python program, an R script utilizes the readxl library with file.choose() command to locally to allow the user to select the excel file of their output from the python program and convert the excel file into an R dataframe. Next, the ggplot2 library is used to create a barplot to visualize the number of research papers written per search, along with some formatting adjustments for readability improvement. This script needs to be run on Rstudio locally due to issues I faced with package management for R libraries on the hpc3.


# Results

For pipeline testing purposes, I used a sample dataset of 100bp paired-end reads from an RNA-seq experiment I previously conducted looking at FACS isolated tendon fibroblasts (tenocytes) from zebrafish embryos comparing between 36hpf and 48 hpf. These timepoints correspond to embryonic developmental stages prior to, and after, the onset of active muscle contraction respectively. The first bash program outputs a text file called number_of_Ns.txt, the output format of which is shown in Figure 1. Each fastq file has a bad read percentage less than 1% indicating that the reads from this particular sequencing experiment may be of good quality.

Figure 1:

![fastq_read_quality_analysis](/code/scripts/bash_output.png)

Figure 2 displays the input and output of the python program where the inputs are conventional gene names, and Figure 3 shows the input and output of the python program where the inputs are ENSEMBL IDs. It appears the ncbi API sometimes gives large numbers of articles written when the gene is not fully known, or given obscure gene names by the ENSEMBL API, hence the top genes having hits of 122. In future, a filtering function  could be written which removes unknown/unnamed genes to streamline the process and prevent unwanted data. Both the input and output of this data are excel files, which can be passed next into the R script.

Figure 2: (input (top) and output (bottom) of python program with gene list with conventional gene names)

![gene-list-input-nonens](/code/scripts/python_nonens_input.png)

![gene-list-output-nonens](/code/scripts/python_nonens_output.png)

Figure 3: (input (top) and output (bottom) of python program with gene list with ENSEMBL ID names)

![gene-list-input-ens](/code/scripts/python_ens_input.png)

![gene-list-output-ens](/code/scripts/python_ens_output.png)

The R script outputs a plot, shown in Figure 4 and Figure 5, displaying visually the genes with the most research articles written for the given search term. Ignoring the high hits from unknown/unnamed genes, we can see the top hit genes for the current search term of "tendon" being tnmd, with 114 hits, pcna with 60 hits, tgfb3 with 53 hits, lum with 16 hits, etc... These are genes with known roles in tendon development and repair, showing that the scripts indeed provide a "narrowing down" of the gene list from 500 genes (in the non-ensembl gene list example) to ~40 interesting genes in this context. We can also validate that the number of hits is indeed true by actually searching a few of these terms, and checking the number of hits on our browser. An image of this is shown in Figure 6.

Figure 4: output of R barplot accepting as input the output from previous python program for 500 genes:

![R-barplot-500 genes](/code/scripts/r-plot-500.png)

Figure 5: output of R barplot accepting as input the output from previous python program for ensembl IDs with ~10 genes:

![R-barplot-500 genes](/code/scripts/r-plot-short.png)

Figure 6: validation of gene list hits against pubmed articles web search

![POC](/code/scripts/POC.PNG)

# Discussion

In this analysis, I have demonstrated a group of programs for determination of read quality and drill-down of relevant genes from the output of a NGS experiment. The read quality analysis was completed using bash scripting and piping tools, whereas the "narrow-down" script pipeline was written as API calling functions in python and a plotting function in R. In our sample data, the read-quality analysis bash script returned all fastq files as having bad-read percentages less than 1% with "bad reads" being defined as reads with 10 consecutive N's within the fastq files. The python script and R script combination reduced the scope of the output of our RNA-seq gene list, which had resulted from a differential expression pipeline, to a more user-friendly list of genes relevant to the research terms we were interested in. To further improve this pipeline in the future a few changes could be made. For the read-quality script, I could rewrite the loop to check in both fastq files of a single paired-end sample run to check for poor quality reads, and make sure not to double count reads between a sample, such that the program returns poor quality read counts PER SAMPLE instead of per fastq file. I could add in a feature in the python program searching the ENSEMBL API for orthologous genes for each gene in multiple species, and then search all of these orthologuess with a search term against the pubmed API, to account for orthologues which may been studied in a particular context in one species, but not in the species of the original dataset. This would allow for easier detection of potentially understudied genes within the NGS output gene list for a particular research context.A bash script could also be written to pipe the input gene list through the python file and utilize stdout/stdin to pipe the output xlsx file as input into the R-plot to seamlessly generate the data and plot in a single script. Additionally, a further direction could be to integrate this pipeline into the back end of a web application GUI to make the tool easier and more accessible to use for scientists without background in programming/command line usage.
