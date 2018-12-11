cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

inputs:

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
    label: "FASTA input file"
    format: "http://edamontology.org/format_1929"
    doc: "Reference genome input FASTA file"

  annotation_gtf:
    type: File
    label: "GTF input file"
    format: "http://edamontology.org/format_2306"
    doc: "Annotation input file"

  annotation_tab:
    type: File
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file"

  genome_sa_sparse_d:
    type: int?
    default: 2
    label: "Use 2 to decrease needed RAM for STAR"
    doc: |
      int>0: suffux array sparsity, i.e. distance between indices: use bigger
      numbers to decrease needed RAM at the cost of mapping speed reduction
    'sd:layout':
      advanced: true

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

  limit_genome_generate_ram:
    type: long?
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
    doc: "Folder which includes all STAR generated indices files"
    outputSource: star_generate_indices/indices

  annotation:
    type: File
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file"
    outputSource: annotation_tab

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
      genome_dir:
        default: "."
      genome_fasta_files: fasta
      sjdb_gtf_file: annotation_gtf
      genome_sa_sparse_d: genome_sa_sparse_d
      limit_genome_generate_ram: limit_genome_generate_ram
      genome_sa_index_n_bases: genome_sa_index_n_bases
      genome_chr_bin_n_bits: genome_chr_bin_n_bits
      threads: threads
    out: [indices, chr_name_length]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "Generate genome index STAR RNA"
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

doc: >
  Generates indices for STAR v2.5.3a (03/17/2017).

s:about: |
  Workflow makes indices for [STAR](https://github.com/alexdobin/STAR) v2.5.3a (03/17/2017) PMID: [23104886](https://www.ncbi.nlm.nih.gov/pubmed/23104886).

  It performs the following steps:
  1. Runs `STAR --runMode genomeGenerate` to generate indices, based on [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) and [GTF](http://mblab.wustl.edu/GTF2.html) input files, returns results as an array of files
  2. Transforms array of files into [Direcotry](http://www.commonwl.org/v1.0/CommandLineTool.html#Directory) data type
  3. Separates *chrNameLength.txt* file as an output

