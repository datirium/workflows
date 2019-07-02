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
  genome_indices:
    - "genome-indices.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  bambai_pair_cond_1:
    type: File[]
    format: "http://edamontology.org/format_2572"
    label: "Biological condition 1"
    doc: "Coordinate sorted BAM alignment and index BAI files for the first biological condition"
    'sd:upstreamSource': "first_biological_condition/bambai_pair"
    'sd:localLabel': true

  bambai_pair_cond_2:
    type: File[]
    format: "http://edamontology.org/format_2572"
    label: "Biological condition 2"
    doc: "Coordinate sorted BAM alignment and index BAI files for the second biological condition"
    'sd:upstreamSource': "second_biological_condition/bambai_pair"
    'sd:localLabel': true

  chrom_length_file:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Chromosome length file"
    doc: "Chromosome length file"
    'sd:upstreamSource': "genome_indices/chrom_length"

  merge_peaks:
    type: boolean?
    label: "Merge peaks closer than fragment size"
    doc: "Merge peaks which have a distance less than the estimated mean fragment size (recommended for histone data)"
    'sd:layout':
      advanced: true

  housekeeping_genes_bed_file:
    type: File?
    format: "http://edamontology.org/format_3003"
    label: "Housekeeping genes file"
    doc: "Define housekeeping genes (BED format) used for normalizing"
    'sd:layout':
      advanced: true

  deadzones_bed_file:
    type: File?
    format: "http://edamontology.org/format_3003"
    label: "Dead zones file"
    doc: "Define blacklisted genomic regions avoided for analysis"
    'sd:layout':
      advanced: true

  pvalue_cutoff:
    type: float?
    label: "P-value cutoff for peak detection"
    doc: "P-value cutoff for peak detection. Call only peaks with p-value lower than cutoff. [default: 0.1]"
    'sd:layout':
      advanced: true

  extension_size:
    type:
      - "null"
      - int[]
    label: "Read's extension size"
    doc: |
      Read's extension size for BAM files (comma separated list for each BAM file in config file).
      If option is not chosen, estimate extension sizes
    'sd:layout':
      advanced: true


outputs:

  diffpeaks_bed_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Estimated differential peaks"
    doc: "Estimated differential peaks"
    outputSource: thor/diffpeaks_bed_file
    'sd:visualPlugins':
    - igvbrowser:
        id: 'igvbrowser'
        type: 'bed'
        name: "Differential peaks"
        height: 120

  cond_1_bigwig_file:
    type: File[]
    format: "http://edamontology.org/format_3006"
    label: "First biological condition ChIP-seq signals"
    doc: "Postprocessed ChIP-seq signals from the first biological condition samples"
    outputSource: thor/cond_1_bigwig_file

  cond_2_bigwig_file:
    type: File[]
    format: "http://edamontology.org/format_3006"
    label: "Second biological condition ChIP-seq signals"
    doc: "Postprocessed ChIP-seq signals from the second biological condition samples"
    outputSource: thor/cond_2_bigwig_file

  thor_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "rgt-THOR log"
    doc: "rgt-THOR log"
    outputSource: thor/stdout_log


steps:

  thor:
    run: ../tools/rgt-thor.cwl
    in:
      bambai_pair_cond_1: bambai_pair_cond_1
      bambai_pair_cond_2: bambai_pair_cond_2
      chrom_length_file: chrom_length_file
      merge_peaks: merge_peaks
      housekeeping_genes_bed_file: housekeeping_genes_bed_file
      deadzones_bed_file: deadzones_bed_file
      pvalue_cutoff: pvalue_cutoff
      extension_size: extension_size
    out:
      - diffpeaks_bed_file
      - cond_1_bigwig_file
      - cond_2_bigwig_file
      - stdout_log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "THOR - differential peak calling of ChIP-seq signals with replicates"
label: "THOR - differential peak calling of ChIP-seq signals with replicates"
s:alternateName: "THOR is an HMM-based approach to detect and analyze differential peaks in two sets of ChIP-seq data from distinct biological conditions with replicates"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/rgt-thor.cwl
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
  What is THOR?
  --------------

  THOR is an HMM-based approach to detect and analyze differential peaks in two sets of ChIP-seq data
  from distinct biological conditions with replicates. THOR performs genomic signal processing, peak
  calling and p-value calculation in an integrated framework.

  For more information please refer to:
  -------------------------------------

  Allhoff, M., Sere K., Freitas, J., Zenke, M., Costa, I.G. (2016), Differential Peak Calling of ChIP-seq
  Signals with Replicates with THOR, Nucleic Acids Research, epub gkw680.
