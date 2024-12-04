# Polarization of ancestral and derived polymorphic alleles for inferences of extended haplotype homozygosity in human populations

---

**Overview**

The Polaris package consists of two distinct programs written in the C-Python language: (i) Panderas and (ii) Panderas_Plots that run in a command-line Terminal window on Mac OS X and Linux machines. This document will guide users through the process of applying our software to recode alleles and generate an associated genetic map (*.hap and *map, respectively) for direct use as input files by selscan and HaploSweep. In addition, we provide detailed instructions on how to apply the Panderas_Plots software to create Manhattan plots of standardized statistics and EHH line graphs using an array of command-line options.

---

**Table of Contents**

- [Installation of Panderas](#installation-of-panderas)
- [Requirements for Panderas](#requirements-for-panderas)
- [Ancestral allele file](#ancestral-allele-file)
  - [Basic Format](#basic-format-1)
  - [Example Format](#example-format-1)
- [Reference Genetic Map File](#reference-genetic-map-file)
  - [Basic Format](#basic-format-2)
  - [Example Format](#example-format-2)
- [Run the Panderas software](#run-the-panderas-software)
  - [Command-Line Usage](#command-line-usage)
  - [Example of Usage](#example-of-usage)
  - [Help Command](#help-command)
- [Output Files](#output-files)

---

**Panderas**

---

**Overview**

The Panderas software classifies ancestral and derived alleles on chromosomes in phased vcf files for analysis by two haplotype-based programs, selscan and HaploSweep.

---

## Installation of Panderas

The Polaris package can be downloaded to users’ local computers by clicking on Polaris under "Releases" in our Github repository `(https://github.com/alisi1989/Polaris)`. Alternatively, users can download the Polaris package from DropBox `(https://www.dropbox.com/scl/fo/2s6hm1owv26wivrn29saw/AGkE6YM91jp_DzJlgNScaSY?rlkey=pxe1i6ytbcsup7jv6h1q5fq3m&st=b1dyusit&dl=0)`. To unzip the file from GitHub (the same file from DropBox will be uncompressed), type `"unzip Polaris_main.zip"` at the command line prompt `(typically indicated by a ">" symbol)`, and the uncompressed `"Polaris_Main"` folder will appear. At the command-line prompt in the Terminal window, users will change the current working directory to the "Panderas" directory where the software is located.

Users will then need to ensure the software is executable with the following command:

`chmod +x Panderas`

---

## Requirements for Panderas

1.    Variant call format file
The input for the Panderas software consists of phased bi-allelic genotypes in a variant call format (vcf) file. The vcf file should have the following format: 

a)    It should be phased. This task can be accomplished with a number of publicly available programs, such as Beagle, Shapeit4, Eagle, or similar phasing software.

b)    It should be gzipped (e.g., *.vcf.gz) and properly indexed (e.g., *.vcf.gz.tbi).

c)    All variants must have an rs identifier (no missing or empty identifiers are allowed). Variants with missing rs identifiers (e.g., '.') could cause issues in downstream analyses.

d)    For vcf files created using build hg38, chromosomes should be annotated with "chr" followed by the chromosome number (e.g., "chr1", "chr2",…"chr22").

In cases where the vcf file contains variants with missing rs identifiers (e.g., "."), users can run our custom script, provided in the Fill_rsidentifier folder in our GitHub repository, to replace missing identifiers with unique identifiers (e.g., rs1, rs2, rs3, etc.).

`Command-line Usage: ` 

```bash
./Fill_._rsunknown.py –-input [argument] --output [argument]
```

where users will provide: 1) the path to the gzipped input vcf file after the "--input" flag; and 2) the path to where the gzipped output vcf file will be saved after the "--output" flag. 

`Example of Usage:`

```bash
./Fill_._rsunknown.py –-input path_to_folder/Dataset-Test/ --output path_to_folder/Output/
```

## Ancestral allele file

Given the space limitations on GitHub, the Homo_sapiens_hg38_reference file from Ensembl for each chromosome (https://ftp.ensembl.org/pub/release-112/fasta/ancestral_alleles/) is stored separately in DropBox (https://www.dropbox.com/scl/fo/0du8z7xoeqs5qqr73qk2s/ADuGddGTZXwUIeZkmQZ-geM?rlkey=chre7oc2fcwpl4bi9p2gab4jb&st=5igiuhrx&dl=0) and are available for downloaded. If users wish to incorporate a different ancestral allele reference file, they will need to make sure the file is in a tab-delimited format (please see below). 


### `Basic Format`

```bash
Position  Allele
```

### `Example of Format:`

```
16050075  G
16050103  A
16055683  T
```

## Reference Genetic Map File

The PLINK genetic map files (i.e., reference maps from https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/) are provided in DropBox `(https://www.dropbox.com/scl/fo/b3a9z16sqjvksprqudlpg/ACv4xl3Nk9HHZP0m4Em0CVI?rlkey=hn23fb1j5vie5d84ss3hpf8nv&st=rnm2ve4n&dl=0)`. Genetic map files should contain the genetic position and physical position in white space delimited columns. 

`Basic Format:`

```
GeneticDistance PhysicalPosition
```

  
`Example of Format:`

```
0.123456 16050000  
0.125678 16060000  
0.128901 16070000  
```


## Run the Panderas software

The Panderas software is provided as a standalone executable file. Users can run this program as follows:

`Command-Line Usage`

```
./Panderas -v [argument] -o [argument] -c [argument] –-ancestor [argument] –-genetic-map [argument]
```

where users must provide: 1) the name of the gzipped input vcf file with the extension (e.g., input_file.vcf.gz) after the `"-v"` flag; 2) the prefix of an output filename with no file extension (e.g., "output_file") after the `"-o"` flag; 3) the chromosome number to extract using the `-c` flag; 4) the name of the ancestral allele file with the *.txt extension after the `"-–ancestor"` flag;  5) the name of the reference genetic map file (e.g., mod_genetic_map_GRCh38_chr2.txt) after the `"-–genetic-map"` flag

The runtime for generating the *.hap and *.map files for alleles on Chromosome 2 in the Finnish population was 15 minutes. Please note that this runtime could vary across populations and chromosomes. 

Panderas also provides error messages to help identify issues with input files and/or command-line arguments.

`Example of Usage`

```
./Panderas -v Example/Finnish_1KG.vcf.gz -c chr2 --ancestor Ancestor-Alleles/ancestor_alleles_chr2.txt --genetic-map Genetic-Map/genetic_map_chr2.txt -o Output/Finnish_chr2
```

Additional information about these arguments also can be accessed using the following command: 

./Panderas.py --help

Finally, regarding the genetic map, the Panderas software will perform a linear interpolation to determine the genetic position (centimorgans, cM) of loci in a given dataset based on the physical position (i.e., genomic coordinate) and the corresponding genetic position in the generic PLINK map file (i.e., the reference). Specifically, Panderas will search for the two physical and genetic positions in the PLINK reference file that flank a query locus in a given dataset and calculate the genetic position for that query locus using information in the reference.

---

## Output Files

After running Panderas, the following output files will be generated:

```
`*.hap:` The haplotype file contains the polarized alleles (0/1) for each haplotype.

`*.map:` The genetic map file contains the interpolated genetic positions.
```

The output files are now ready for direct use as input file by selscan and/or HaploSweep. More specific information on the usage of selscan and/or HaploSweep can be found at https://github.com/szpiech/selscan and https://github.com/ChenHuaLab/HaploSweep, respectively. 

The normalized output files from these programs will serve as input files for Panderas_Plots. If users wish to familiarize themselves with Panderas using the example *.hap and *.map files on our GitHub repository, these files must be uncompressed prior to use.



**Panderas_Plots**

Table of Contents
•    Overview
•    Installation
•    Requirements
•    Manhattan Plotting 
o    Basic-Command-Line Options
•    Example of Usage
o    Manhattan plot with threshold line and SNP annotations
•    EHH Plotting
o    Basic-Command-Line Options
o    Example of Usage
    Plotting EHH from Selscan Output
    Plotting EHH from HaploSweep Output
•    License
•    Contact


**Overview** 
The Panderas_Plot software offers an array of options to visualize and annotate the numerical output of selscan and HaploSweep. Specifically, it creates Manhattan plots of standardized statistics and EHH line graphs showing the decay of haplotype homozygosity on chromosomes around ancestral and derived alleles.

**Installation**
The Polaris package can be downloaded to users’ local computers by 
clicking on Polaris under "Releases" in our Github repository (https://github.com/alisi1989/Polaris). Alternatively, users can download the Polaris package from DropBox (https://www.dropbox.com/scl/fo/2s6hm1owv26wivrn29saw/AGkE6YM91jp_DzJlgNScaSY?rlkey=pxe1i6ytbcsup7jv6h1q5fq3m&st=b1dyusit&dl=0). To unzip the file from GitHub (the same file from DropBox will be uncompressed), type "unzip Polaris_main.zip" at the command line prompt (typically indicated by a ">" symbol), and the uncompressed "Polaris_Main" folder will appear. Using the command-line interface in the Terminal window, users will change the current working directory to the "Panderas_Plots" directory where the software is located.

Users will need to ensure that the software is executable by running the following command:

chmod +x Panderas_Plots

**Requirements**

The normalized output files from selscan and/or HaploSweep (e.g., Finnish_1KG_chr2.ihs.out.100bins.norm from selscan and Finnish_norm.txt from HaploSweep) for iHS calculation. Please note that these input files must be uncompressed prior to use. 


**Manhattan Plotting**

To visualize haplotype-based statistics along a given chromosome, users will execute Panderas_Plots as follows:

`Basic Command-Line Usage: `

```
./Panderas_Plots.py manhattan –-alg [argument] --input [argument] --chr [argument] --output [argument] 
```

where users will provide: 1) the "manhattan" function before the "--alg" flag; 2) the algorithm used to generate the haplotype-based statistic (either selscan or haplosweep in lowercase letters) after the "--alg" flag;  3) the input file from selscan or HaploSweep with the file extension (e.g., Finnish_1KG_chr2.ihs.out.100bins.norm) after the "--input" flag; 4) the chromosome analyzed (e.g., "2" or "chr 2") after the "--chr" flag; and 5) the output filename (with either a pdf, svg, eps, or png extension at the end of the filename) after the "--output" flag. 

Panderas_Plots also provides error messages to help identify issues with input files or command-line arguments.

`Example of Usage: `

<pre>
./Panderas_Plots manhattan --alg selscan --input selscan_ihs.out --chr 2 --output selscan_manhattan.pdf
</pre>

Notably, Panderas_Plots also offers an array of options to visualize and annotate the numerical output of selscan and/or HaploSweep. For example, 





Manhattan plot with threshold line and SNP annotations

<pre>
./Panderas_Plots manhattan --alg haplosweep --input haplosweep_ihs.out --chr 2 --output haplosweep_manhattan.pdf --threshold-line 2.0 --color-line red --rs-annot 2.0 --label-annot y --statistic iHS
<pre>

where users can select options to create a threshold line and annotate SNPs in addition to the required arguments described above. In particular, users can specify: 1) an iHS Y-value with the "--threshold-line" flag for outlier values (e.g., 2); 2) a specific color of the threshold line using the "--color-line" flag (e.g., red); 3) threshold for annotating statistics with rs identifiers in absolute numbers (e.g., 2, 3, 4) with "--rs-annot" flag; and 4) specify (yes or no) if users want to add labels (rs identifiers) to dots with the "--label-annot" flag. Because labels can overlap, this flag gives users the option to remove the labels for a cleaner plot, if they prefer.

Other noteworthy features include: i) changing the size and/or color of each dot in Manhattan plots with the --size–dots flag; and ii) creating a list of genes harboring outlier iHS values in a *.txt file using "--gene–annot". 

For this latter feature ("--gene–annot"), if users wish to include gene annotations in the Manhattan plots, they will need to prepare a separate gene annotation text file in tab-delimited format:

chr    start    end    gene

a)    chr: Chromosome number (e.g., chr2).
b)    start: Start position of the gene (e.g., 135800000).
c)    end: End position of the gene (e.g., 135850000).
d)    gene: Gene name (e.g., MCM6).

Although we provide a reference file to annotate by gene name (HG38_UCSC_refGene_filtered .txt), users may want to use another reference file. If so, this reference file should be formatted as above. However, it is our recommendation that users utilize the gene reference file that we provide.

Users can also access a more complete list of plotting features with the following command:

./Panderas_Plots.py --help

 

EHH Plotting

To visualize the decay of haplotype homozygosity along chromosomes, users will execute Panderas_Plots as follows:

`Basic Command-Line Usage: `

```
./Panderas_Plots ehh --alg {selscan,haplosweep} --input INPUT --output OUTPUT 
```

where users will provide: 1) the "ehh" function before the "--alg" flag; 2) the algorithm used to generate the haplotype-based statistic (either selscan or haplosweep in lowercase letters) after the "--alg" flag 3) the input created by selscan or HaploSweep (with the extension) after the "--input" flag; 4) the user-specified output filename (with either a pdf, svg, eps, or png extension at the end of the filename) after the "--output" flag; and. 

Panderas_Plots will also provide error messages to help identify issues with input files or command-line arguments.

`Examples of Usage: `

Plotting EHH from Selscan Output

```
./Panderas_Plots ehh --alg selscan --input selscan_ehh.out --output selscan_ehh_plot.pdf
```

Plotting EHH from HaploSweep Output

```
./Panderas_Plots ehh --alg haplosweep --input haplosweep_ehh.out --output haplosweep_ehh_plot.pdf 
```

INCLUDE IMAGE HERE

Recommendations and Notes
•    Consistent Chromosome Naming: Ensure that chromosome identifiers are consistent across your input files (e.g., "chr2" or "2").
•    Highlighting SNPs: To highlight specific SNPs, users are required to create a text file (snps_to_highlight.txt) with one SNP identifier per line. The filename will be entered after the "--snps-to-highlight" flag.

`License`

This project is licensed under the MIT License.

`Contact`

For questions or comments, please contact:
1.    Alessandro Lisi, alisi@usc.edu or on GitHub (alisi1989)
2.    Michael C. Campbell, mc44680@usc.edu or on GitHub (mc44680)
