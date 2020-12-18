cwlVersion: v1.0
class: Workflow


requirements:
- class: InlineJavascriptRequirement


'sd:upstream':
  rnaseq_sample:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-pe-smarter-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"
    - "trim-quantseq-mrnaseq-se.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  rpkm_genes:
    type: File[]
    format: "http://edamontology.org/format_3475"
    label: "RNA-Seq experiments"
    doc: |
      CSV/TSV file with RPKM grouped by gene name
    'sd:upstreamSource': "rnaseq_sample/rpkm_genes"
    'sd:localLabel': true

  sample_names:
    type: string[]
    label: "RNA-Seq experiments"
    doc: |
      Aliases for RNA-Seq experiments to be used as the names for the reported
      columns in the merged file. Order corresponds to the rpkm_genes input
    'sd:upstreamSource': "rnaseq_sample/alias"
    'sd:localLabel': true


outputs:

  merged_file:
    type: File
    label: "Merged gene expression file"
    format: "http://edamontology.org/format_3475"
    doc: |
      Merged by GeneId, Chrom, TxStart, TxEnd and Strand gene expression files
      with reported and renamed Rpkm columns.
    outputSource: feature_merge/merged_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'RPKM Gene Expression'
        Title: 'Merged RPKM Gene Expression Results'

  feature_merge_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Gene expression merge stdout log"
    doc: "Gene expression merge stdout log"
    outputSource: feature_merge/stdout_log

  feature_merge_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Gene expression merge stderr log"
    doc: "Gene expression merge stderr log"
    outputSource: feature_merge/stderr_log


steps:

  feature_merge:
    run: ../tools/feature-merge.cwl
    in:
      feature_files: rpkm_genes
      feature_aliases: sample_names
      mergeby:
        default: ["GeneId", "Chrom", "TxStart", "TxEnd", "Strand"]
      report:
        default: "Rpkm"
    out:
    - merged_file
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/version/9.0/schemaorg-current-http.rdf

s:name: "Gene expression merge - combines RPKM gene expression from several experiments"
label: "Gene expression merge - combines RPKM gene expression from several experiments"
s:alternateName: "Gene expression merge - combines RPKM gene expression from several experiments"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/fastqc.cwl
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
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Gene expression merge - combines RPKM gene expression from several experiments
  ===================================================================================

  Workflows merges RPKM gene expression from several experiments based on the values
  from GeneId, Chrom, TxStart, TxEnd and Strand columns. Reported RPKM columns are
  renamed based on the experiments names.
