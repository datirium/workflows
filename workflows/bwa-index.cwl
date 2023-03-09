cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement


'sd:upstream':
  genome_indices: "genome-indices.cwl"


inputs:

  alias:
    type: string
    label: "Sample short name/Alias"
    sd:preview:
      position: 1

  reference_fasta:
    type: File
    'sd:upstreamSource': "genome_indices/fasta_output"
    label: "Reference genome FASTA file to index:"
    'sd:localLabel': true
    doc: |
      FASTA file of the reference genome that will be indexed.
    sd:preview:
      position: 2


outputs:

  index_directory:
    type: Directory
    label: "Directory containing the original FASTA, faidx, dict, and bwa index files."
    outputSource: index_reference/bwa_index

  log_file_stdout:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout logfile"
    outputSource: index_reference/log_file_stdout
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  log_file_stderr:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr logfile"
    outputSource: index_reference/log_file_stderr     


steps:

  index_reference:
    run: ../tools/bwa-index.cwl
    in:
      ref_genome_fasta: reference_fasta
    out: [bwa_index, log_file_stdout, log_file_stderr]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "BWA index pipeline"
label: "BWA index pipeline"
s:alternateName: "BWA index pipeline"

s:downloadUrl: https://github.com/datirium/workflows/tree/master/workflows/workflows/bwa-index.cwl
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
  This workflow indexes the input reference FASTA with bwa, and generates faidx and dict file using samtools.
  This index sample can then be used as input into the germline variant calling workflow, or others that may
  include this workflow as an upstream source.

  ### __Inputs__
   - FASTA file of the reference genome that will be indexed.
  
  ### __Outputs__
   - Directory containing the original FASTA, faidx, dict, and bwa index files.
   - stdout log file (output in Overview tab as well)
   - stderr log file

  ### __Data Analysis Steps__
  1. cwl calls dockercontainer robertplayer/scidap-gatk4 to index reference FASTA with bwa, and generates faidx and dict files using samtools

  ### __References__
    - Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics, 25(14), 1754–1760.
