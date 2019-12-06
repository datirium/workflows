cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  first_biological_condition:
    - "chipseq-se.cwl"
    - "chipseq-pe.cwl"
    - "trim-chipseq-se.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "trim-atacseq-pe.cwl"
  second_biological_condition:
    - "chipseq-se.cwl"
    - "chipseq-pe.cwl"
    - "trim-chipseq-se.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "trim-atacseq-pe.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  read_files_cond_1:
    type: File[]
    format: "http://edamontology.org/format_2572"
    label: "Biological condition 1 samples"
    doc: "Read files for condition 1. Minimim 2 files in BAM format"
    'sd:upstreamSource': "first_biological_condition/bambai_pair"
    'sd:localLabel': true

  read_files_cond_2:
    type: File[]
    format: "http://edamontology.org/format_2572"
    label: "Biological condition 2 samples"
    doc: "Read files for condition 2. Minimim 2 files in BAM format"
    'sd:upstreamSource': "second_biological_condition/bambai_pair"
    'sd:localLabel': true

  peak_files_cond_1:
    type: File[]
    format: "http://edamontology.org/format_3468"
    label: "Biological condition 1 samples"
    doc: "XLS peak files for condition 1 from MACS2. Minimim 2 files. Order corresponds to read_files_cond_1"
    'sd:upstreamSource': "first_biological_condition/macs2_called_peaks"

  peak_files_cond_2:
    type: File[]
    format: "http://edamontology.org/format_3468"
    label: "Biological condition 2 samples"
    doc: "XLS peak files for condition 2 from MACS2. Minimim 2 files. Order corresponds to read_files_cond_2"
    'sd:upstreamSource': "second_biological_condition/macs2_called_peaks"

  name_cond_1:
    type: string?
    default: "condition_1"
    label: "Condition 1 name, single word with letters and numbers only"
    doc: "Condition 1 name, single word with letters and numbers only"
    'sd:layout':
      advanced: true

  name_cond_2:
    type: string?
    default: "condition_2"
    label: "Condition 2 name, single word with letters and numbers only"
    doc: "Condition 2 name, single word with letters and numbers only"
    'sd:layout':
      advanced: true


outputs:

  diffbind_report_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Differential binding analysis results"
    doc: "Differential binding analysis results exported as TSV"
    outputSource: diffbind/diffbind_report_file
    'sd:visualPlugins':
      - syncfusiongrid:
          tab: 'Differential Peak Calling'
          Title: 'Differential Binding Analysis Results'

  diffbind_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "diffbind stdout log"
    doc: "diffbind stdout log"
    outputSource: diffbind/stdout_log

  diffbind_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "diffbind stderr log"
    doc: "diffbind stderr log"
    outputSource: diffbind/stderr_log


steps:

  diffbind:
    run: ../tools/diffbind.cwl
    in:
      read_files_cond_1: read_files_cond_1
      read_files_cond_2: read_files_cond_2
      peak_files_cond_1: peak_files_cond_1
      peak_files_cond_2: peak_files_cond_2
      name_cond_1: name_cond_1
      name_cond_2: name_cond_2
      peakformat:
        default: "macs"
    out:
      - diffbind_report_file
      - stdout_log
      - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

label: "DiffBind - Differential Binding Analysis of ChIP-Seq Peak Data"
s:name: "DiffBind - Differential Binding Analysis of ChIP-Seq Peak Data"
s:alternateName: "Compute differentially bound sites from multiple ChIP-seq experiments using affinity (quantitative) data."

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/diffbind.cwl
s:codeRepository: https://github.com/datirium/workflows
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

doc: |
  Differential Binding Analysis of ChIP-Seq Peak Data
  ---------------------------------------------------
  
  DiffBind processes ChIP-Seq data enriched for genomic loci where specific protein/DNA binding occurs, including peak sets identified by ChIP-Seq peak callers and
  aligned sequence read datasets. It is designed to work with multiple peak sets simultaneously, representing different ChIP experiments (antibodies, transcription
  factor and/or histone marks, experimental conditions, replicates) as well as managing the results of multiple peak callers.

  For more information please refer to:
  -------------------------------------
  Ross-Innes CS, Stark R, Teschendorff AE, Holmes KA, Ali HR, Dunning MJ, Brown GD, Gojis O, Ellis IO, Green AR, Ali S, Chin S, Palmieri C, Caldas C, Carroll JS (2012).
  “Differential oestrogen receptor binding is associated with clinical outcome in breast cancer.” Nature, 481, -4.
  
