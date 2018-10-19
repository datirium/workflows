cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement


inputs:

  genome:
    type: string
    label: "Genome"
    doc: "Used by BioWardrobe to set genome"

  genome_label:
    type: string?
    label: "Genome label"
    doc: "Genome label is used by web-ui to show label"
    'sd:preview':
      position: 1

  genome_description:
    type: string?
    label: "Genome description"
    doc: "Genome description is used by web-ui to show description"
    'sd:preview':
      position: 2

  genome_details:
    type: string?
    label: "Genome details"
    doc: "Genome details"
    'sd:preview':
      position: 3

  fasta:
    type: File
    label: "Genome FASTA file"
    format: "http://edamontology.org/format_1929"
    doc: "Reference genome FASTA file"

  fasta_ribosomal:
    type: File?
    label: "Ribosomal DNA sequence FASTA file"
    format: "http://edamontology.org/format_1929"
    doc: "Ribosomal DNA sequence FASTA file"

  fasta_mitochondrial:
    type: File?
    label: "Mitochondrial chromosome sequence FASTA file"
    format: "http://edamontology.org/format_1929"
    doc: "Mitochondrial chromosome sequence FASTA file"

  effective_genome_size:
    type: string
    label: "Effective genome size"
    doc: "MACS2 effective genome size: hs, mm, ce, dm or number, for example 2.7e9"

  input_annotation_gtf:
    type: File
    label: "GTF input file"
    format: "http://edamontology.org/format_2306"
    doc: "Annotation input file"

  annotation_tab:
    type: File
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file"

  genome_sa_index_n_bases:
    type: int?
    label: "length of the SA pre-indexing string"
    doc: |
      For small genomes, the parameter --genomeSAindexNbases must to be scaled down, with a typical value of
      min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal to 9,
      for 100 kiloBase genome, this is equal to 7.

      default: 14

      int: length (bases) of the SA pre-indexing string.
      Typically between 10 and 15. Longer strings will use much more memory,
      but allow faster searches. For small genomes, the parameter –genomeSAindexNbases
      must be scaled down to min(14, log2(GenomeLength)/2 - 1).
    'sd:layout':
      advanced: true

  genome_sa_index_n_bases_mitochondrial:
    type: int?
    default: 6
    label: "length (mitochondrial) of the SA pre-indexing string"
    doc: |
      For small genomes, the parameter --genomeSAindexNbases must to be scaled down, with a typical value of
      min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal to 9,
      for 100 kiloBase genome, this is equal to 7.

      default: 14

      int: length (bases) of the SA pre-indexing string.
      Typically between 10 and 15. Longer strings will use much more memory,
      but allow faster searches. For small genomes, the parameter –genomeSAindexNbases
      must be scaled down to min(14, log2(GenomeLength)/2 - 1).
    'sd:layout':
      advanced: true

  genome_chr_bin_n_bits:
    type: int?
    label: "Genome Chr Bin NBits"
    doc: |
      If you are using a genome with a large (>5,000) number of references (chrosomes/scaﬀolds),
      you may need to reduce the --genomeChrBinNbits to reduce RAM consumption.
      The following scaling is recommended: --genomeChrBinNbits = min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)]).
      For example, for 3 gigaBase genome with 100,000 chromosomes/scaﬀolds, this is equal to 15.

      default: 18

      int: =log2(chrBin), where chrBin is the size of the bins for genome storage:
      each chromosome will occupy an integer number of bins.
      For a genome with large number of contigs, it is recommended to scale this parameter
      as min(18, log2[max(GenomeLength/NumberOfReferences,ReadLength)]).
    'sd:layout':
      advanced: true

  genome_sa_sparse_d:
    type: int?
    label: "Genome SA sparse (Use 2 to decrease RAM usage)"
    doc: |
      default: 1

      int>0: suffux array sparsity, i.e. distance between indices: use bigger
      numbers to decrease needed RAM at the cost of mapping speed reduction
    'sd:layout':
      advanced: true

  limit_genome_generate_ram:
    type: long?
    label: "Genome Generate RAM (31G default)"
    inputBinding:
      position: 1
      prefix: --limitGenomeGenerateRAM
    doc: |
      31000000000
      int>0: maximum available RAM (bytes) for genome generation
    'sd:layout':
      advanced: true

  threads:
    type: int?
    label: "Number of threads to run tools"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  star_indices:
    type: Directory
    label: "STAR indices folder"
    doc: "Folder which includes all STAR generated indices folder"
    outputSource: star_generate_indices/indices

  bowtie_indices:
    type: Directory
    label: "Bowtie indices folder"
    doc: "Folder which includes all Bowtie generated indices folder"
    outputSource: bowtie_generate_indices/indices

  ribosomal_indices:
    type: Directory
    label: "Ribosomal DNA indices folder"
    doc: "Ribosomal DNA Bowtie generated indices folder"
    outputSource: ribosomal_generate_indices/indices

  mitochondrial_indices:
    type: Directory
    label: "Mitochondrial chromosome index folder"
    doc: "Mitochondrial chromosome index folder"
    outputSource: mitochondrial_generate_indices/indices

  annotation_gtf:
    type: File
    label: "GTF input file"
    format: "http://edamontology.org/format_2306"
    doc: "Annotation input file"
    outputSource: input_annotation_gtf

  annotation:
    type: File
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file"
    outputSource: annotation_tab

  genome_size:
    type: string
    label: "Effective genome size"
    doc: "MACS2 effective genome size: hs, mm, ce, dm or number, for example 2.7e9"
    outputSource: effective_genome_size

  chrom_length:
    type: File
    label: "Chromosome length file"
    format: "http://edamontology.org/format_2330"
    outputSource: star_generate_indices/chr_name_length
    doc: "Chromosome length file"


steps:

  star_generate_indices:
    run: ../tools/star-genomegenerate.cwl
    in:
      genome_fasta_files: fasta
      sjdb_gtf_file: input_annotation_gtf
      genome_sa_sparse_d: genome_sa_sparse_d
      limit_genome_generate_ram: limit_genome_generate_ram
      genome_sa_index_n_bases: genome_sa_index_n_bases
      genome_chr_bin_n_bits: genome_chr_bin_n_bits
      threads: threads
    out: [indices, chr_name_length]

  mitochondrial_generate_indices:
    run: ../tools/star-genomegenerate.cwl
    in:
      genome_fasta_files: fasta_mitochondrial
      sjdb_gtf_file: input_annotation_gtf
      genome_sa_sparse_d: genome_sa_sparse_d
      genome_sa_index_n_bases: genome_sa_index_n_bases_mitochondrial
      threads: threads
    out: [indices]

  bowtie_generate_indices:
    run: ../tools/bowtie-build.cwl
    in:
      fasta_file: fasta
      index_base_name: genome
    out: [indices]

  ribosomal_generate_indices:
    run: ../tools/bowtie-build.cwl
    in:
      fasta_file: fasta_ribosomal
      index_base_name: genome
    out: [indices]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "Generate genome indices for STAR & bowtie"
s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/genome-indices.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium, LLC"
  s:email: mailto:support@datirium.com
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45226"
    s:streetAddress: "3559 Kroger Ave"

doc: >
  Generates genome indices for STAR v2.5.3a (03/17/2017) & bowtie v1.2.0 (12/30/2016).

s:about: |
  Creates indices for:
   * [STAR](https://github.com/alexdobin/STAR) v2.5.3a (03/17/2017) PMID: [23104886](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
   * [bowtie](http://bowtie-bio.sourceforge.net/tutorial.shtml) v1.2.0 (12/30/2016)

  It performs the following steps:

  1. `STAR --runMode genomeGenerate` to generate indices, based on [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) and [GTF](http://mblab.wustl.edu/GTF2.html) input files, returns results as an array of files
  2. Outputs indices as [Direcotry](http://www.commonwl.org/v1.0/CommandLineTool.html#Directory) data type
  3. Separates *chrNameLength.txt* file from Directory output
  4. `bowtie-build` to generate indices requires genome [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) file as input, returns results as a group of main and secondary files

