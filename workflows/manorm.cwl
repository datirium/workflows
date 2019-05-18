cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


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
  genome_indices:
    - "genome-indices.cwl"


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

  annotation_file:
    type: File
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file"
    'sd:upstreamSource': "genome_indices/annotation"
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

  output_filename:
    type: string?
    label: "MAnorm output TSV filename"
    doc: "MAnorm output TSV filename"
    default: "manorm_common_peak.tsv"
    'sd:layout':
      advanced: true

outputs:

  common_peak_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "MAnorm common peak resutls, TSV"
    doc: "MAnorm generated list of common peaks with assigned genes"
    outputSource: restore_columns/restored_peak_file
    'sd:visualPlugins':
    - syncfusiongrid:
        Title: 'MAnorm Common Peak Results'


steps:

  manorm:
    in:
      peak_file_first: peak_file_first
      peak_file_second: peak_file_second
      bam_file_first: bam_file_first
      bam_file_second: bam_file_second
      fragment_size_first: fragment_size_first
      fragment_size_second: fragment_size_second
    out: [common_peak_file]
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
        common_peak_file:
          type: File
          outputBinding:
            glob: "MAnorm_result_commonPeak_merged.xls"
      baseCommand: ["run_manorm.sh"]

  filter_columns:
    in:
      peak_file: manorm/common_peak_file
      script:
        default: >
          cat $0 | grep -v "start" | awk
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t"NR"\t0\t0\t0\t0"}' > `basename $0`
    out: [filtered_peak_file]
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap:v0.0.3
      inputs:
        script:
          type: string
          inputBinding:
            position: 1
        peak_file:
          type: File
          inputBinding:
            position: 2
      outputs:
        filtered_peak_file:
          type: File
          outputBinding:
            glob: "*"
      baseCommand: [bash, '-c']

  assign_genes:
      in:
        peak_file: filter_columns/filtered_peak_file
        annotation_file: annotation_file
        output_filename: output_filename
        promoter_bp:
          default: 1000
      out: [peaks_and_genes_file]
      run:
        cwlVersion: v1.0
        class: CommandLineTool
        requirements:
        - class: DockerRequirement
          dockerPull: biowardrobe2/iaintersect:v0.0.2
        inputs:
          peak_file:
            type: File
            inputBinding:
              position: 1
              prefix: --in=
              separate: false
          annotation_file:
            type: File
            inputBinding:
              position: 2
              prefix: --a=
              separate: false
          output_filename:
            type: string
            inputBinding:
              position: 3
              prefix: --out=
              separate: false
          promoter_bp:
            type: int
            inputBinding:
              position: 4
              prefix: --promoter=
              separate: false
        outputs:
          peaks_and_genes_file:
            type: File
            outputBinding:
              glob: $(inputs.output_filename)
        baseCommand: [iaintersect]

  restore_columns:
    in:
      peak_files: [assign_genes/peaks_and_genes_file, manorm/common_peak_file]
      script:
        default: |
          cat $0 | grep -v "start" | sort -k 11n > sorted_iaintersect_result.tsv
          cat $1 | grep -v "start" > manorm_result.tsv
          echo -e "refseq_id\tgene_id\ttxStart\ttxEnd\tstrand\tchrom\tstart\tend\tlength\tregion\tdescription\t#raw_read_1\t#raw_read_2\tM_value_rescaled\tA_value_rescaled\t-log10(p-value)" > `basename $0`;
          cat sorted_iaintersect_result.tsv | paste - manorm_result.tsv | cut -f 1-9,15,19-24 >> `basename $0`
          rm sorted_iaintersect_result.tsv manorm_result.tsv
    out: [restored_peak_file]
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      requirements:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap:v0.0.3
      inputs:
        script:
          type: string
          inputBinding:
            position: 1
        peak_files:
          type: File[]
          inputBinding:
            position: 2
      outputs:
        restored_peak_file:
          type: File
          outputBinding:
            glob: "*"
      baseCommand: [bash, '-c']


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
