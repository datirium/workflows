cwlVersion: v1.0
class: Workflow


'sd:upstream':
  first_chipseq_sample:
    - "chipseq-se.cwl"
    - "chipseq-pe.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-chipseq-se.cwl"
  second_chipseq_sample:
    - "chipseq-se.cwl"
    - "chipseq-pe.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-chipseq-se.cwl"

inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  peak_file_first:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "First TSV peak file"
    doc: "TSV peak file, formatted as iaintersect output"
    'sd:upstreamSource': "first_chipseq_sample/iaintersect_result"
    'sd:localLabel': true

  peak_file_second:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Second TSV peak file"
    doc: "TSV peak file, formatted as iaintersect output"
    'sd:upstreamSource': "second_chipseq_sample/iaintersect_result"
    'sd:localLabel': true

  bam_file_first:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "First BAM file"
    doc: "BAM alignment file"
    'sd:upstreamSource': "first_chipseq_sample/bambai_pair"
    'sd:localLabel': true

  bam_file_second:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Second BAM file"
    doc: "BAM alignment file"
    'sd:upstreamSource': "second_chipseq_sample/bambai_pair"
    'sd:localLabel': true

  fragment_size_first:
    type: int?
    label: "First fragment size"
    doc: "Fragment size, int"
    default: 150

  fragment_size_second:
    type: int?
    label: "Second fragment size"
    doc: "Fragment size, int"
    default: 150

outputs:

  common_peak_merged_file:
    type: File
    label: "MAnorm resutls, TSV"
    format: "http://edamontology.org/format_3475"
    doc: "MAnorm generated list of common peaks"
    outputSource: manorm/common_peak_merged_file
    'sd:visualPlugins':
    - syncfusiongrid:
        Title: 'MAnorm results'

steps:

  manorm:
    in:
      peak_file_first: peak_file_first
      peak_file_second: peak_file_second
      bam_file_first: bam_file_first
      bam_file_second: bam_file_second
      fragment_size_first: fragment_size_first
      fragment_size_second: fragment_size_second
    out: [common_peak_merged_file]
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: DockerRequirement
        dockerPull: biowardrobe2/manorm:v0.0.1
      inputs:
        peak_file_first:
          type: File
          inputBinding:
            position: 5
        peak_file_second:
          type: File
          inputBinding:
            position: 6
        bam_file_first:
          type: File
          inputBinding:
            position: 7
        bam_file_second:
          type: File
          inputBinding:
            position: 9
        fragment_size_first:
          type: int
          inputBinding:
            position: 10
        fragment_size_second:
          type: int
          inputBinding:
            position: 11
      outputs:
        common_peak_merged_file:
          type: File
          outputBinding:
            glob: "MAnorm_result_commonPeak_merged.xls"
      baseCommand: ["run_manorm.sh"]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "MAnorm: a robust model for quantitative comparison of ChIP-Seq data sets"
label: "MAnorm: a robust model for quantitative comparison of ChIP-Seq data sets"
s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/manorm.cwl
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
  MAnorm: a robust model for quantitative comparison of ChIP-Seq data sets

s:about: |
  MAnorm: a robust model for quantitative comparison of ChIP-Seq data sets
  =====================================

  ChIP-Seq is widely used to characterize genome-wide binding patterns of transcription factors and other chromatin-associated proteins. Although comparison of ChIP-Seq data sets is critical for understanding cell type-dependent and cell state-specific binding, and thus the study of cell-specific gene regulation, few quantitative approaches have been developed. Here, we present a simple and effective method, MAnorm, for quantitative comparison of ChIP-Seq data sets describing transcription factor binding sites and epigenetic modifications. The quantitative binding differences inferred by MAnorm showed strong correlation with both the changes in expression of target genes and the binding of cell type-specific regulators.
