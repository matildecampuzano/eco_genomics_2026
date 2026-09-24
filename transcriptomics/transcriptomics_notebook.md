# Transcriptomics Notebook

**Course:** Intro to Ecological Genomics - Fall 2026

**Name:** Matilde Campuzano

------------------------------------------------------------------------

## 9.15.2026 - Setting up lab notebook and learning markdown

-   Setting up transcriptomics notebook

-   Learn how to take notes in markdown

-   Push notes to github

**Working Directory:**

`/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/transcriptomics`

**Input Files:**

`none`

**Output Files:**

`/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/transcriptomics/transcriptomics_notebook.md`

**Programs and dependencies:**

-   `R version 4.5.1`

-   `R-Studio`

**Scripts**

`none`

**Code:**

``` r
print("Hello World")
```

**Table**:

| Col1 | Col2 | Col3 |
|------|------|------|
|      |      |      |
|      |      |      |
|      |      |      |

**Image:**

![](images/markdown-syntax-cheatsheet.webp)

**Notes/Observations:**

oh 67

**Next Steps:**

------------------------------------------------------------------------

## 9.17.2026 - Diving into code

-   Discussion

-   Research background info and questions

-   Learn how to unzip fast.qc files, etc.

**Working Directory:**

`/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/transcriptomics`

**Commands Used:**

`zcat filename | head -n 4`

-   opened a gzipped file, head opens the first part, -n 4 puts in number first 4 lines

`cd /filepath/`

-   change directory to filepath that you put in

`ls ll`

-   List, list long

`zcat .fq.gz(file) | wc -l`

-   open gzipped file, wc means count, -l means lines

`history`

-   shows all my previous commands of the day

**Notes/Observations:**

remember tab to finish commands yayy

------------------------------------------------------------------------

## 9.22.2026 - Gene Expression Analysis

-   Set up Rstudio working environment

-   Set up directories in Transcriptomics directory

-   Made mydata directory secret from github

-   Copied data to import into DeSeq2

**Working Directory:**

`/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/transcriptomics`

**Script Used:**

`~/Projects/eco_genomics_2026/transcriptomics/myscripts/ahud_DESeq2_inclas.R`

**Commands Learned:**

`git status`

-   See whether git is up to date

`git branch`

-   ??

`git pull/push`

-   duh

`mkdir filename`

-   make a directory within whatever directory I'm in

`cp * /filepath`

-   copy stuff into filepath

`touch ./filepath/filename.txt`

-   make file(filename.txt) in filepath

`rm filename.txt`

-   remove filename.txt

**Final Plot:**

![](myresults/PCA_allGens.png){width="1000"}

------------------------------------------------------------------------

## 9.22.2026 - Gene Expression Analysis

-   Set up Rstudio working environment

-   Set up directories in Transcriptomics directory

-   Made mydata directory secret from github

-   Copied data to import into DeSeq2

**Working Directory:**

`/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/transcriptomics`

-   mine :)

`/gpfs1/cl/biol3990/Transcriptomics/CountsMatrix`

-   Where we got our data in our class filing cabinet!

**Script Used:**

`~/Projects/eco_genomics_2026/transcriptomics/myscripts/ahud_DESeq2_inclas.R`

**Commands Learned:**
