cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement


inputs:

  alias:
    type: string
    label: "Sample short name/Alias:"
    'sd:localLabel': true
    doc: |
      Short name for the analysis.
    sd:preview:
      position: 1

  reference_fasta:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "Reference genome FASTA file to index:"
    'sd:localLabel': true
    doc: "FASTA file of the reference genome that will be indexed. May be compressed with gzip. Must end in '.gz', '.fasta', '.fa', or '.fna'.\n\n
      NOTE: Sequence names (string after the '>') in the transcriptome FASTA must match column 3 (Chrom/gene id/name) of the annotation TSV."
    sd:preview:
      position: 2

  input_annotation_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Annotation file (gff, gtf, tsv):"
    'sd:localLabel': true
    doc: "TSV file containing gene annotations for the reference genome (from kallisto index upstream).\n\n
      Required columns (include headers as row 1 of TSV):\n
      \t1. RefseqId\n
      \t2. GeneId\n
      \t3. Chrom (gene/transcript id/name)\n
      \t4. TxStart\n
      \t5. TxEnd\n
      \t6. Strand\n\n
      NOTE: Sequence names (string after the >) in the transcriptome FASTA must match column 3 (Chrom) of the annotation TSV."
    sd:preview:
      position: 3

  threads:
    type: int?
    default: 10
    label: "Threads:"
    'sd:localLabel': true
    doc: |
      Number of threads to use for steps that support multithreading.


outputs:

  index_file:
    type: File
    label: "Kallisto index file."
    outputSource: kallisto_index/kallisto_index

  log_file_stdout:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout logfile"
    outputSource: kallisto_index/log_file_stdout
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  log_file_stderr:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr logfile"
    outputSource: kallisto_index/log_file_stderr     


steps:

  kallisto_index:
    run: ../tools/kallisto-index.cwl
    in:
      ref_genome_fasta: reference_fasta
      annotation_tsv: input_annotation_file
      threads: threads
    out: [kallisto_index, log_file_stdout, log_file_stderr]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Kallisto index pipeline"
label: "Kallisto index pipeline"
s:alternateName: "Kallisto index pipeline"

s:downloadUrl: https://github.com/datirium/workflows/tree/master/workflows/workflows/kallisto-index.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X


doc: |
  This workflow indexes the input reference FASTA with kallisto, and generates a kallisto index file (.kdx).
  This index sample can then be used as input into the kallisto transcript-level quantification workflow
  (kallisto-quant-pe.cwl), or others that may include this workflow as an upstream source.

  ### __Inputs__
   - FASTA file of the reference genome that will be indexed
   - number of threads to use for multithreading processes
  
  ### __Outputs__
   - kallisto index file (.kdx).
   - stdout log file (output in Overview tab as well)
   - stderr log file

  ### __Data Analysis Steps__
  1. cwl calls dockercontainer robertplayer/scidap-kallisto to index reference FASTA with `kallisto index`, generating a kallisto index file.

  ### __References__
    -   Bray, N. L., Pimentel, H., Melsted, P. & Pachter, L. Near-optimal probabilistic RNA-seq quantification, Nature Biotechnology 34, 525-527(2016), doi:10.1038/nbt.3519
