
##THIS IS ALL LINUX

#first look
zcat Lottia_gigantea_biallelic.vcf.gz | head -n 20

#Read header etc####
module load gcc
module load bcftools

bcftools --version   # record the versions you use in your methods
bcftools -h          # -h or --help usually gets you the manual

bcftools head Lottia_gigantea_pl_flt3.renamed.vcf.gz

#put output into grep to remove big string of contig lines
bcftools head Lottia_gigantea_pl_flt3.renamed.vcf.gz | grep -v "##contig="

#Defines tags
bcftools view -h Lottia_gigantea_pl_flt3.renamed.vcf.gz | grep "^##INFO"
bcftools view -h Lottia_gigantea_pl_flt3.renamed.vcf.gz | grep "^##FORMAT"

#recosnturct how file was made and structured
bcftools view -h Lottia_gigantea_pl_flt3.renamed.vcf.gz | grep "Command="

#list samples####
bcftools query -l Lottia_gigantea_pl_flt3.renamed.vcf.gz

# count the samples by piping into wc -l (counts lines)
bcftools query -l Lottia_gigantea_pl_flt3.renamed.vcf.gz | wc -l

#Reading Variant record####
# just the first record
bcftools view -H Lottia_gigantea_pl_flt3.renamed.vcf.gz | head -n 1

# the first 5 records; less -S stops long lines wrapping (press q to exit)
bcftools view -H Lottia_gigantea_pl_flt3.renamed.vcf.gz | head -n 5 | less -S

# -v selects variant types: here SNPs only
bcftools view -H -v snps Lottia_gigantea_pl_flt3.renamed.vcf.gz | head -n 5 | less -S

# what do you notice when we ask for indels?
bcftools view -H -v indels Lottia_gigantea_pl_flt3.renamed.vcf.gz | head -n 5 | less -S

#Viewing a region####
ls Lottia_gigantea_pl_flt3.renamed.vcf.gz*   # is there a .csi or .tbi index?
  
 # all SNPs between positions 8,000 and 9,000 on scaffold KB201661.1
  bcftools view -H -r KB201661.1:8000-9000 Lottia_gigantea_pl_flt3.renamed.vcf.gz

# for a large region, pipe into head
bcftools view -H -r KB201661.1:8000-9000 Lottia_gigantea_pl_flt3.renamed.vcf.gz \
| head -n 5

#Pulling out specific fields####
# scaffold, position, alleles and every individual's genotype for the first 5 SNPs
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
Lottia_gigantea_pl_flt3.renamed.vcf.gz | head -n 5

#Summary statistics with vcftools####
#calculate the proportion of missing data and the mean depth for each individual:
module load vcftools

# proportion of missing genotypes per individual -> missing_data.imiss
vcftools --gzvcf Lottia_gigantea_flt3_sub0.01.vcf.gz --missing-indv --out missing_data

# mean depth per individual -> depth_data.idepth
vcftools --gzvcf Lottia_gigantea_flt3_sub0.01.vcf.gz --depth --out depth_data

#OK NOW THIS IS R#
#PLot ts####
library(ggplot2)

##REPLACE WITH YOUR OWN PATH##
setwd("/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/population_genomics")

##Notice here the mydata sub directory!
imiss  <- read.delim("mydata/missing_data.imiss")
idepth <- read.delim("mydata/depth_data.idepth")

p1<-ggplot(imiss, aes(x = F_MISS)) +
  geom_histogram(bins = 30) +
  labs(x = "Proportion of missing genotypes", y = "Number of individuals")

p2<-ggplot(idepth, aes(x = MEAN_DEPTH)) +
  geom_histogram(bins = 30) +
  labs(x = "Mean depth", y = "Number of individuals")
p1
p2
###REPLACE WITH YOUR PATH TO MY RESULTS
OUTDIR="/gpfs1/home/m/c/mcampuza/Projects/eco_genomics_2026/population_genomics/myresults"

ggsave(p1, file = paste0(OUTDIR,"/missing_data.png"))
ggsave(p2, file = paste0(OUTDIR,"/depth_data.png"))
