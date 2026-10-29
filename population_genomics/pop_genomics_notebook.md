# Population Genomics Notebook

**Course:** Intro to Ecological Genomics - Fall 2026

**Name:** Matilde Campuzano

------------------------------------------------------------------------

## 10.15.2026 - Setting up lab notebook and learning markdown

-   Set up population genomics notebook
-   Put data in popgen again
-   Learned how to look in a VCF

**Working Directory:**

`/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/population_genomics`

-   mine

`/gpfs1/cl/biol3990/PopulationGenomics/VCF/.vcf`

-   classes data

**Code:**

``` r
print("Hello World")
```

**Image:**

![](images/markdown-syntax-cheatsheet.webp)

**Notes/Observations:**

-   SRA in NIH has a bajill sequences free, lots of data
-   0 is reference allele, 1 is alternate one, . is missing
-   VCF column meanings:
    -   #CHROM Chromosome (for us, scaffold) where the variants are found
    -   POS Position of the variant on that chromosome/scaffold ID Variant identifier; usually . (empty) unless you work on human genetics
    -   REF Reference allele
    -   ALT Alternative allele(s)
    -   QUAL Phred-scaled quality of the variant call
    -   FILTER PASS if the site passed all filters, otherwise the name of the filter it failed
    -   INFO Site-level information, e.g. total depth and allele counts; each tag is defined in a ##INFO header line
    -   FORMAT The order of the individual-level fields that follow, e.g. GT:PL:DP:AD:GQ
-   When reading variant record given output:
    -   Each individual has one entry, with its fields in the order given by FORMAT. For the first individual, 0/1:0,9,107:3:3,0:127 means:
    -   GT: 0/0: Genotype: 0/ = REF allele, /0 = ALT allele, so 0/0 is homozygous for the reference (A/A)
    -   PL: 0,9,107: Phred-scaled likelihoods of the genotypes 0/0, 0/1 and 1/1; the most likely genotype is scaled to 0 and higher values are less likely
    -   DP: 3: Number of reads covering this site in this individual
    -   AD: 3,0: Allelic depth: 3 reads support the REF allele and 0 support the ALT allele
    -   GQ: 127: Genotype quality: Phred-scaled confidence in the GT call

**Next Steps:**

------------------------------------------------------------------------

## 10.20.2026 - Population Structure

-   Ran batch job

**Working Directory:**

`/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/population_genomics`

-   mine

**Commands:**

`DATA_DIR=/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/population_genomics/mydata`

-   Do that every time, labels data so we can zoom over there

`VCF=Lottia_gigantea_flt3_sub0.01.vcf.gz`

-   name ts

**Code:**

``` r
print("Hello World")
```

**Image:**

**Notes/Observations:**

-   Population structure things we look at:

    -   Genetic diversity

        -   PCA - also helps to exclude extreme outliers, look at data/quality assurance)

        -   Admixture

        -   \*filter out linkage disequilibrium (SNP's inherited together) which is based on r which is correlation coefficient by removing one of them! and the other represents both in the data

    -   Population Structure

        -   MAF- Minor Allele Frequence, filter out samples w/ weird allele not seen(rare or mistake)

------------------------------------------------------------------------

## 10.22.2026 - Genetic Diversity

-   

**Working Directory:**

`/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/population_genomics`

-   mine

**Commands:**

**Code:**

``` r
print("Hello World")
```

**Image:**

**Notes/Observations:**

-   Looking at Genetic Diversity:

    -   Within Individual:

        -   Diversity between structure/chromosomes, what are more heterozygous etc.

    -   Within Population:

        -   Fitness, projecting dynamics?

    -   Within Species:

        -   Diversification of populations within species

-   Metrics for measuring:

    -   Fst -\> fixation index(population idfferentiation)

    -   Heterozygosity (using Hardy-Weinberg Equilibrium)

    -   Ne -\> Genetic, Nc -\> Consensus

    -   SNP's (Nucleotide diversity aka pi)

    -   Fis -\> Inbreeding coefficient

        -   if F\>0=Inbred

        -   H0 is het. indiv/total indiv.

    -   Allele richness (\# unique alleles per locus)

    -   Tajima'sD= How many segregating sites do we have

        -   (pi-Vaterson's theta)/(Var\*sqrt(pi-theta))

            -   pi is pairwise nucleotide diversity -\>pairwise makes it more sensitive to rare alleles

            -   theta is (segregating sites)/(sample size coefficients)

        -   Should equal 0 if expected

        -   \<0 is excess of rare mutations

            -   Could come from sensitive sweep or rapid population expansion

        -   \>0 is excess of common mutations

            -   Balancing selection, more heterozygosity

------------------------------------------------------------------------

## 10.27.2026 - Genetic Diversity & Selection Sweeps

-   Heterozygosity

-   Fst

-   etc.

**Working Directory:**

`/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/population_genomics`

-   mine

**Commands:**

**Code:**

``` r
print("Hello World")
```

**Image:**

**Notes/Observations:**

-   biooracle has all the environmental data, but its averages so can't show seasonal nuance, etc.

-   looked at windows of Fst, see where it is highest in the genome between two groups

------------------------------------------------------------------------

## 10.29.2026 - Selection Sweeps

-   BayPass

**Working Directory:**

`/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/population_genomics`

-   mine

**Commands:**

**Code:**

``` r
print("Hello World")
```

**Image:**

**Notes/Observations:**

-   looked at windows of Fst, see where it is highest in the genome between two groups

-   BayPass- give it prior information that it can build off of

    -   Can solve problems of Fst( only two populations at a time, ignore shared history of population, only show differentiation in locus, not why)

    -   Omega: pairwise differences between sites, which represents neutral population structure drift etc.

        -   Corrects for demographic effects on genomic structure, like genetic drift, non-random mating to pick out NATURAL SELECTION

    -   Uses alelle counts instead of allele frequency

    -   uses Xtx based on omega

    -   Run Baypass 5+ times and then average it to make sure it isn't pulling out super random stuff

    -   Higher Bays factor is significant, anything above 20
