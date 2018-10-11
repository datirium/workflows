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

  fasta:
    type: File
    label: "FASTA input file"
    format: "http://edamontology.org/format_1929"
    doc: "Reference genome input FASTA file"

  annotation_gtf:
    type: File?
    label: "GTF input file"
    format: "http://edamontology.org/format_2306"
    doc: "Annotation input file"

  annotation_tab:
    type: File?
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file"

  genome_sa_sparse_d:
    type: int?
    label: "Use 2 to decrease needed RAM for STAR"
    doc: |
      int>0: suffux array sparsity, i.e. distance between indices: use bigger
      numbers to decrease needed RAM at the cost of mapping speed reduction
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

  bowtie_indices:
    type: Directory
    label: "Bowtie indices folder"
    doc: "Folder which includes all Bowtie generated indices files"
    outputSource: bowtie_generate_indices/indices

  annotation:
    type: File?
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file"
    outputSource: annotation_tab

  genome_size:
    type: string?
    label: "Effective genome size"
    doc: "MACS2 effective genome size: hs, mm, ce, dm or number, for example 2.7e9"

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
#      genome_dir:
#        default: "./star"
      genome_fasta_files: fasta_input_file
      sjdb_gtf_file: annotation_gtf_file
      genome_sa_sparse_d: genome_sa_sparse_d
      limit_genome_generate_ram:
        default: 15000000000
      threads: threads
    out: [indices, chr_name_length]

  bowtie_generate_indices:
    run: ../tools/bowtie-build.cwl
    in:
      fasta_file: fasta_input_file
      index_base_name: genome
    out: [indices]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "genome-indices"
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
  Creates indices for STAR v2.5.3a (03/17/2017) & bowtie v1.2.0 (12/30/2016).

s:about: |
  Creates indices for:
   * [STAR](https://github.com/alexdobin/STAR) v2.5.3a (03/17/2017) PMID: [23104886](https://www.ncbi.nlm.nih.gov/pubmed/23104886)
   * [bowtie](http://bowtie-bio.sourceforge.net/tutorial.shtml) v1.2.0 (12/30/2016)

  It performs the following steps:

  1. Runs `STAR --runMode genomeGenerate` to generate indices, based on [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) and [GTF](http://mblab.wustl.edu/GTF2.html) input files, returns results as an array of files
  2. Transforms array of files into [Direcotry](http://www.commonwl.org/v1.0/CommandLineTool.html#Directory) data type
  3. Separates *chrNameLength.txt* file as an output

  1. Executes `bowtie-build` to generate indices requires genome [FASTA](http://zhanglab.ccmb.med.umich.edu/FASTA/) file as input, returns results as a group of main and secondary files
  2. Transforms results from the previous step into Direcotry data type

