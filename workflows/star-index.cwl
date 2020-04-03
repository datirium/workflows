cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement


'sd:metadata':
  - "../metadata/indices-header.cwl"


inputs:

  genome:
    type: string
    label: "Genome type"
    doc: "Genome type, such as mm10, hg19, hg38, etc"

  fasta_file:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "Reference genome FASTA file"
    doc: "Reference genome FASTA file. Includes all chromosomes"

  annotation_gtf_file:
    type: File?
    format: "http://edamontology.org/format_2306"
    label: "GTF annotation file"
    doc: "GTF annotation file"

  genome_sa_sparse_d:
    type: int?
    default: 2
    label: "Suffix array sparsity (use 2 to decrease needed RAM)"
    doc: |
      Suffix array sparsity, i.e. distance between indices: use bigger
      numbers to decrease needed RAMat the cost of mapping speed reduction"
    'sd:layout':
      advanced: true

  genome_sa_index_n_bases:
    type: int?
    label: "Length of SA pre-indexing string"
    doc: |
      Length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more memory,
      but allow faster searches. For small genomes, the parameter –genomeSAindexNbases must be scaled down to
      min(14, log2(GenomeLength)/2 - 1). For example, for 1 megaBase genome, this is equal to 9, for 100 kiloBase genome,
      this is equal to 7.
      default: 14
    'sd:layout':
      advanced: true

  genome_chr_bin_n_bits:
    type: int?
    label: "Number of bins allocated for each chromosome"
    doc: |
      If you are using a genome with a large (>5,000) number of references (chrosomes/scaﬀolds), you may need to reduce the
      --genomeChrBinNbits to reduce RAM consumption. For a genome with large number of contigs, it is recommended to scale
      this parameter as min(18, log2[max(GenomeLength/NumberOfReferences,ReadLength)]).
      default: 18
    'sd:layout':
      advanced: true

  limit_genome_generate_ram:
    type: long?
    label: "Limit maximum available RAM (bytes) for genome generation"
    doc: "Maximum available RAM (bytes) for genome generation. Default 31000000000"
    'sd:layout':
      advanced: true

  threads:
    type: int?
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  indices_folder:
    type: Directory
    outputSource: star_generate_indices/indices_folder
    label: "STAR indices"
    doc: "STAR generated indices folder"

  chrom_length_file:
    type: File
    format: "http://edamontology.org/format_2330"
    outputSource: star_generate_indices/chrom_length
    label: "Chromosome length file"
    doc: "Chromosome length file"
    
  stdout_log:
    type: File
    label: "STAR stdout log"
    doc: "STAR generated stdout log"
    outputSource: star_generate_indices/stdout_log

  stderr_log:
    type: File
    label: "STAR stderr log"
    doc: "STAR generated stderr log"
    outputSource: star_generate_indices/stderr_log


steps:

  star_generate_indices:
    run: ../tools/star-genomegenerate.cwl
    in:
      genome_dir: genome
      genome_fasta_files: fasta_file
      sjdb_gtf_file: annotation_gtf_file
      genome_sa_sparse_d: genome_sa_sparse_d
      genome_sa_index_n_bases: genome_sa_index_n_bases
      genome_chr_bin_n_bits: genome_chr_bin_n_bits
      limit_genome_generate_ram: limit_genome_generate_ram
      threads: threads
    out:
    - indices_folder
    - chrom_length
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "Build STAR indices"
label: "Build STAR indices"
s:alternateName: "Generates indices for STAR v2.5.3a (03/17/2017)"

s:downloadUrl: https://raw.githubusercontent.com/barski-lab/workflows/master/workflows/scidap/star-index.cwl
s:codeRepository: https://github.com/barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898
      - class: s:Person
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0001-9102-5681


# doc:
#   $include: ../descriptions/star-index.md


doc: |
  Workflow runs [STAR](https://github.com/alexdobin/STAR) v2.5.3a (03/17/2017) PMID: [23104886](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
  to build indices for reference genome provided in a single FASTA file as fasta_file input and GTF annotation file from annotation_gtf_file input.
  Generated indices are saved in a folder with the name that corresponds to the input genome.
