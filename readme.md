# Artemis

Artemis is a simple bioinformatics pipeline which scans a reference genome for PAM sites, to then annotate intersecting variants in a VCF file with sequence information for gRNA design.

Steps:
- Load and process genome files (supports gzipped files)
- Scan for PAM sites in the sequence
- Intersect PAM sites with an optional exclusion VCF file (e.g. 1000 genomes)
- Sort and merge PAM-associated seed sites
- Intersect seed regions with the input VCF file (e.g. clinvar)
- Annotate variants with CRISPR target sequences in vcf format

For more information and citations: https://doi.org/10.1016/j.crmeth.2024.100912 

# Dependencies

- bedtools
- pybedtools
- pysam

for the examples below:

- samtools
- bcftools

## Installation

To install Artemis, clone the repository and use `setup.py`:

```sh
git clone https://github.com/jasperlinthorst/artemis.git
cd artemis
python [setup.py] install
```
or use pip

```sh

python3 -m pip install git+https://github.com/jasperlinthorst/ARTEMIS --user

```

## Download
For published use-case download GRCh38, clinvar and 1000 genomes data:

```sh
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20181203_biallelic_SNV/ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz
```

## Usage

```sh
./artemis --help
```

To run the Artemis pipeline, use the following command:  

### Arguments
reference_genome: Path to the human reference genome file (required)  
vcf_file: Path to the VCF file, e.g., a ClinVar VCF file (required)  

### Options
--exclude <kg>: Exclude PAMs that intersect with variants in this VCF (e.g., 1000 genomes data)  
-o <outputfile>: Where to write the annotated output VCF file (default: stdout)  
--seedsize <size>: Specify the size of the seed region (default: 5)  
--pamregex <regex>: Specify regular expression that matches PAM sites, first group for positive strand and second group for negative strand matches, (default Cas12a: "(?=(TTT[ACG]))|(?=([CGT]AAA))")  

NOTE, make sure vcf and reference genomes use the same version and contig annotations.  

## Example

```sh
artemis Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz clinvar.vcf.gz > clinvar.cas12a.annotated.vcf
```

NOTE, the above will intersect all ~60 million PAM sites in the human genome with all variants in clinvar and thus wil take long to run and uses large amounts of memory. Better to subset the data, which can be done from the commandline directly. If the reference is bgzipped and indexed using samtools, and vcf files are index using bcftools, we can intersect PAMs on chromosome 22 with pathogenic SNVs on chromosome 22 in clinvar like this:

```sh
artemis <(samtools faidx Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz 22) <(bcftools view clinvar.vcf.gz -i 'INFO/CLNSIG = "Pathogenic" && CLNVC="single_nucleotide_variant"' 22) -o clinvar.cas12a.chr22.pathogenic.snv.vcf
```

To exclude PAM sites that are potentially interrupted by common variants:

```sh
artemis <(samtools faidx Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz 22) <(bcftools view clinvar.vcf.gz -i 'INFO/CLNSIG = "Pathogenic" && CLNVC="single_nucleotide_variant"' 22) --excl <(bcftools view -i 'INFO/AF>=0.01 & INFO/AF<=0.99' ALL.wgs.shapeit2_integrated_v1a.GRCh38.20181129.sites.vcf.gz chr22 | sed -s 's/chr//1') -o clinvar.cas12a.chr22.pathogenic.snv.excl1kg.maf1p.vcf
```

To use different seed size definitions and match PAMs for other Cas proteins use --seedsize and --pamregex arguments. E.g to scan for variants in the human genome that are within 8bp of a SpCas9 PAM:

```sh
artemis Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz <(bcftools view clinvar.vcf.gz -i 'INFO/CLNSIG = "Pathogenic" && CLNVC="single_nucleotide_variant"' 22) --seedsize 8 --pamregex "(?=([ATCG]GG))|(?=(CC[ATCG]))" -o clinvar.cas9.chr22.pathogenic.snv.vcf

```

or for Cas12b

```sh
artemis Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz <(bcftools view clinvar.vcf.gz -i 'INFO/CLNSIG = "Pathogenic" && CLNVC="single_nucleotide_variant"' 22) --seedsize 6 --pamregex "(?=(TT[CT]))|(?=([GA]AA))" -o clinvar.cas12b.chr22.pathogenic.snv.vcf

```

etc.

The resulting vcf file can be imported into genome browsers like IGV, contains all targetable variants and the following additional annotations with respect to CRISPR gRNA spacer sequence and the nearest PAM site:

- SNP_position: Position of SNP relative to PAM site
- CRISPR_ref_target_sequence: CRISPR target sequence to address the reference (WT) allele
- CRISPR_alt_target_sequence: CRISPR target sequence to address the alternative (pathogenic) allele.'
- STRAND: Whether PAM site is found on the positive (+) or negative (-) strand.'

## Contributing
Contributions are welcome! Please open an issue or submit a pull request on GitHub.