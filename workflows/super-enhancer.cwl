#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement


inputs:

  islands_file:
    type: File
    label: "Input TSV file"
    format: "http://edamontology.org/format_3475"
    doc: "Input tab-delimited file (output from iaintersect.cwl or macs2-callpeak-biowardrobe-only.cwl tool)"

  islands_control_file:
    type: File?
    label: "Input TSV file (control)"
    format: "http://edamontology.org/format_3475"
    doc: "Control tab-delimited file (output from iaintersect.cwl or macs2-callpeak-biowardrobe-only.cwl tool)"

  bam_file:
    type: File
    label: "Input TSV file (control)"
    format: "http://edamontology.org/format_3475"
    doc: "Indexed bamfile to rank super-enhancer"

  genome_type:
    type: string
    label: "Genome type"
    doc: "Genome type"

  stitch_distance:
    type: int
    label: "Stitching distance"
    doc: "Linking distance for stitching"

  tss_distance:
    type: int
    label: "TSS distance"
    doc: "Distance from TSS to exclude. 0 = no TSS exclusion"


outputs:

  gff_directory:
    type: Directory
    outputSource: rose/gff_directory
    label: "gff_directory"
    doc: "gff_directory"

  mapped_gff_directory:
    type: Directory
    outputSource: rose/mapped_gff_directory
    label: "mapped_gff_directory"
    doc: "mapped_gff_directory"

  stitched_enhancer_region_map:
    type: File?
    outputSource: rose/stitched_enhancer_region_map
    label: "stitched_enhancer_region_map"
    doc: "stitched_enhancer_region_map"

  all_enhancers_table:
    type: File?
    outputSource: rose/all_enhancers_table
    label: "all_enhancers_table"
    doc: "all_enhancers_table"

  super_enhancers_table:
    type: File?
    outputSource: rose/super_enhancers_table
    label: "super_enhancers_table"
    doc: "super_enhancers_table"

  enhancers_with_super_bed:
    type: File?
    outputSource: rose/enhancers_with_super_bed
    label: "enhancers_with_super_bed"
    doc: "enhancers_with_super_bed"

  plot_points_pic:
    type: File?
    outputSource: rose/plot_points_pic
    label: "plot_points_pic"
    doc: "plot_points_pic"

  gateway_enhancers_bed:
    type: File?
    outputSource: rose/gateway_enhancers_bed
    label: "gateway_enhancers_bed"
    doc: "gateway_enhancers_bed"

  gateway_super_enhancers_bed:
    type: File?
    outputSource: rose/gateway_super_enhancers_bed
    label: "gateway_super_enhancers_bed"
    doc: "gateway_super_enhancers_bed"


steps:

  makegff:
    run: ../tools/makegff.cwl
    in:
      islands_file: islands_file
      islands_control_file: islands_control_file
    out: [gff_file]

  rose:
    run: ../tools/rose.cwl
    in:
      binding_sites_file: makegff/gff_file
      bam_file: bam_file
      genome_type: genome_type
      stitch_distance: stitch_distance
      tss_distance: tss_distance
    out:
    - gff_directory
    - mapped_gff_directory
    - stitched_enhancer_region_map
    - all_enhancers_table
    - super_enhancers_table
    - enhancers_with_super_bed
    - plot_points_pic
    - gateway_enhancers_bed
    - gateway_super_enhancers_bed


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "super-enhancer"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/workflows/super-enhancer.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
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
  Workflow
  Both islands_file and islands_control_file should be produced by the same tool

s:about: |
  Workflow