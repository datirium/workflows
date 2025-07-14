cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/rose:v0.0.2
inputs:
  binding_sites_file:
    type: File
    inputBinding:
      position: 5
      prefix: -i
    doc: GFF or BED file of binding sites used to make enhancers
  bam_file:
    type: File
    inputBinding:
      position: 6
      prefix: -r
    secondaryFiles:
    - .bai
    doc: Indexed bamfile to rank enhancer by
  annotation_file:
    type: File
    inputBinding:
      position: 7
      prefix: -g
    doc: TSV genome annotation file
  stitch_distance:
    type: int
    inputBinding:
      position: 8
      prefix: -s
    doc: Linking distance for stitching
  tss_distance:
    type: int
    inputBinding:
      position: 9
      prefix: -t
    doc: Distance from TSS to exclude. 0 = no TSS exclusion
outputs:
  gff_directory:
    type: Directory
    outputBinding:
      glob: gff
  mapped_gff_directory:
    type: Directory
    outputBinding:
      glob: mappedGFF
  stitched_enhancer_region_map:
    type: File?
    outputBinding:
      glob: '*STITCHED_TSS_DISTAL_ENHANCER_REGION_MAP.txt'
  all_enhancers_table:
    type: File?
    outputBinding:
      glob: '*AllEnhancers.table.txt'
  super_enhancers_table:
    type: File?
    outputBinding:
      glob: '*SuperEnhancers.table.txt'
  enhancers_with_super_bed:
    type: File?
    outputBinding:
      glob: '*Enhancers_withSuper.bed'
  plot_points_pic:
    type: File?
    outputBinding:
      glob: '*Plot_points.png'
  gateway_enhancers_bed:
    type: File?
    outputBinding:
      glob: '*Gateway_Enhancers.bed'
  gateway_super_enhancers_bed:
    type: File?
    outputBinding:
      glob: '*Gateway_SuperEnhancers.bed'
baseCommand:
- ROSE_main
- -o
- ./
doc: |
  Tool runs ROSE to get Super Enhancers regions
  -b and -c arguments are not supported
label: rose
