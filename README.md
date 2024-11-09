**VCF-to-Hap-Map_Polarized**

In this README file, we present VCF-to-Hap-Map_Polarized_Fast, a tool designed to process phased VCF files and generate haplotype and map files, with polarization based on ancestral alleles. This document will guide users through the process of preparing the input files, running the tool, and understanding the output.
1. Prepare Input Datasets

**Phased VCF File**

'Requirement:' The VCF file must be phased using tools like Shapeit4, Eagle, or similar phasing software.
'Format:' The VCF should be gzipped (.vcf.gz) and properly indexed if necessary.
'Important:' Ensure that all variants have valid rsIDs (no missing or empty IDs). Variants with missing IDs (e.g., .) may cause issues in downstream analyses.
'Chromosome Annotation:' For the hg38 build, chromosomes should be annotated with chr followed by the chromosome number (e.g., chr1, chr2, ..., chr22).

**Ancestral Alleles File**

Create an ancestral alleles file in tab-delimited format. The file should contain positions and their corresponding ancestral alleles.
The Homo sapiens ancestral alleles are already provided in the repository. Soure: Ensemble Homo sapiens hg38. 


'Position  Allele'


'Example:'

<pre>
16050075  G
16050103  A
16055683  T
</pre>

'Genetic Map File'

Prepare a genetic map file containing physical positions and their corresponding genetic distances.

'Format:'

<pre>
PhysicalPosition  GeneticDistance
  </pre>
  
'Example:'

<pre>
16050000  0.123456
16060000  0.125678
16070000  0.128901
</pre>

2. Run VCF-to-Hap-Map_Polarized_Fast

The tool is provided as a standalone executable. You can run it directly with the appropriate command-line arguments.
Command-Line Usage
usage: vcf_to_hap_map_polarized_fast [-h] -v VCF -o OUTPUT -c CHROMOSOME --ancestor ANCESTOR --genetic-map GENETIC_MAP
-v, --vcf: Input phased VCF file (gzipped).
-o, --output: Base name for the output files.
-c, --chromosome: Chromosome to extract (e.g., chr22).
--ancestor: Ancestral alleles file.
--genetic-map: Genetic map file for interpolation.
-h, --help: Show the help message and exit.
Example Command
./vcf_to_hap_map_polarized_fast -v your_phased_data.vcf.gz -c chr22 --ancestor ancestor_alleles_chr22.txt --genetic-map genetic_map_chr22.txt -o Output/Baka_chr22
3. Output Files

After running the tool, the following output files will be generated:
<output>.hap: The haplotype file containing the polarized genotypes (0/1) for each haplotype.
<output>.map: The genetic map file with interpolated genetic distances.
4. Recommendations and Notes

Filtering Missing rsIDs: Ensure your VCF file does not contain variants with missing IDs. Use bcftools to filter out such variants.
Example:
bcftools view -e 'ID=="."' -Oz -o filtered.vcf.gz input.vcf.gz
Chromosome Consistency: Make sure that the chromosome naming in your VCF file matches that in your genetic map file (e.g., chr22).
Phasing: If your VCF file is not phased, use a phasing tool like Shapeit4 or Eagle before running this software.
5. License

This project is licensed under the MIT License.
6. Contact

For questions or comments, please contact:
Your Name
Email: your.email@example.com
GitHub: yourusername
