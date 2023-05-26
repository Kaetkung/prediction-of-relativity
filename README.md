# prediction-of-relativity

This file contain 2 steps
1: run bash script for generate "IBD segment" and "PI_HAT"(percent DNA sharing) -> complete_PHASE_PIHAT_IBD_SEGMENT_4.sh
2: run python for relativeness prediction -> predictive.py or predictive.ipynb

########################################################################
describe of bash script
This script use for performing genetic data processing and analysis using a variety of tools like PLINK, Beagle, bcftools, and Refined IBD. Here's a step-by-step overview:
- Convert VCF to PLINK binary format: The script starts by converting a VCF (Variant Call Format) file, which contains genetic data, into the PLINK binary format. This is a format used by PLINK, a free, open-source whole genome association analysis toolset. The converted output is saved as df1.
- List Duplicate Variants: The script then identifies duplicate genetic variants in the data and saves this information.
- Perform Quality Control and Generate Compressed VCF: This step removes genetic variants based on missing genotype rates, minor allele frequencies, and individual missingness. It also excludes any duplicate variants identified in the previous step. The cleaned data is then converted back to the VCF format, compressed and saved as afterQC.vcf.gz.
- Index the Compressed VCF File: The cleaned and compressed VCF file is then indexed using bcftools, a set of tools for interacting with variant call (VCF) and binary variant call (BCF) formats.
- Split VCF File by Chromosome: The script splits the cleaned VCF file by chromosome, resulting in individual VCF files for each of the 22 autosomes.
- Phase Each Chromosome using Beagle without Imputation: Phasing is the process of predicting which alleles originate from the mother and which come from the father. The Beagle software is used to phase each individual chromosome file. Note that imputation is turned off; imputation would estimate missing genotypes, but here we're just phasing.
- Index Phased VCF Files: Each of the phased chromosome files is then indexed using bcftools.
- Concatenate Phased Chromosome Files: All of the phased and indexed chromosome files are then concatenated (joined) into one single file.
- Convert VCF to PLINK Binary Format Again: The concatenated VCF file is then converted back into the PLINK binary format, resulting in df2.
- Estimate Pairwise Identity by Descent (IBD): This step uses PLINK to estimate pairwise IBD, which is a measure of the genetic relationship between individuals.
- Concatenate Genetic Map: This step concatenates (joins together) all the chromosome-specific genetic map files into one combined genetic map.
- Generate IBD Segments: The Refined IBD tool is then used to generate IBD segments. An IBD segment is a stretch of DNA that is identical in two or more individuals because they have inherited it from a common ancestor.
- Merge IBD Segments: The last step merges these IBD segments using merge-ibd-segments tool.
########################################################################

To execute this script, you need to have several specific files in your current directory. Here's the list:
- start.vcf.gz: The input VCF file that contains the original genetic data.
- plink: The executable binary for PLINK software. You might want to ensure that the file has the necessary permissions to be run as an executable.
- beagle.22Jul22.46e.jar: The Beagle software's Java JAR file for phasing genetic data.
- refined-ibd.17Jan20.102.jar: The Java JAR file for the Refined IBD software.
- merge-ibd-segments.17Jan20.102.jar: The Java JAR file for the IBD segment merging tool.
- chr${i}.1kg.phase3.v5a.b37.bref3: The reference panel for Beagle, where ${i} is each number from 1 to 22, so you should have 22 files in total. These files contain phased haplotype information from the 1000 Genomes Project.
- plink.chr${i}.GRCh37.map: The genetic map files used in phasing and IBD detection, where ${i} is each number from 1 to 22. This implies that you should have 22 files in total.
########################################################################

describe of predictive script(python)
This Python script performs an analysis on genetic data. It takes in genotype data, call rates, and segments, and then filters, merges, and clusters the data based on different criteria to predict the degrees of relationships between different samples. The script can be divided into several sections.

1. Data loading and initial cleaning:
The script starts by importing necessary libraries and reading in data from the files 'plink.genome' and 'callrate_real.csv'. The data is then sorted, filtered, and merged based on certain criteria. 'FID1' and 'FID2' columns are dropped, and certain character cleaning operations on 'IID1' and 'IID2' columns are commented out. Call rates are filtered and outliers are removed.

2. Segment file loading and merging:
Next, the script reads a segment file ('merged_ibd_segments.ibd') and merges it with the processed data from the previous steps. The columns 'IID1' and 'IID2' are sorted in both dataframes to ensure that they align correctly when merging.

3. Data transformation for analysis:
The script calculates the average 'cM' (a measure of genetic distance) for each pair of samples and the number of shared segments for each pair. Then it applies certain rules to define the relationship between the samples based on the PI_HAT value, which is a measure of genetic relatedness.

4. Clustering:
After defining relationships for identical twins, parent-child pairs, full siblings, etc., the script applies k-means clustering on the remaining samples. It uses a combination of PI_HAT, shared_segment, and avg_cM for clustering. The result of clustering is visualized using a 3D scatter plot.

5. Relationship definition and output:
Finally, the script defines relationships for the remaining samples based on the clusters from the k-means algorithm, combines all of the processed data into one dataframe, and outputs it as a .csv file named 'dfx.csv'. It also filters the result by a specific sample ID and saves it as 'testdf.csv'.
########################################################################

To execute this python script, you need to have this files
- plink.genome(% DNA sharing or PI_HAT) from step 1
- merged_ibd_segments.ibd(IBD segment) from step 1
