cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement


inputs:

  sam_file:
    type: File
    label: "SAM file"
    doc: "SAM file with reads from both strains mapped to concatenated genome"

  chrom_length_file:
    type: File
    label: "Chromosome length file for reference genome"
    doc: "Chromosome length file for reference genome"

  hal_file:
    type: File
    label: "HAL file"
    doc: "HAL file that includes current and reference strain information"

  current_strain_name:
    type: string
    label: "Current strain name"
    doc: "Current strain name"

  reference_strain_name:
    type: string
    label: "Reference strain name"
    doc: "Reference strain name to be projected to"

  mapped_reads_number:
    type: int
    label: "Mapped to concatenated genome reads number"
    doc: "Mapped to concatenated genome reads number to calculate scaling factor"

  output_file_prefix:
    type: string
    label: "Prefix for all generated output files"
    doc: "Corresponds to UID"

  threads:
    type: int?
    default: 2
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"


outputs:

  bambai_pair:
    type: File
    outputSource: samtools_sort_index/bam_bai_pair
    label: "BAM mapped to strain genome, not projected to reference genome"
    doc: "Coordinate sorted BAM file mapped to the current strain genome, not projected to reference genome"

  bigwig_file:
    type: File
    label: "Strain specific bigWig file"
    doc: "Generated bigWig file for the current strain, projected to reference genome"
    outputSource: bedgraph_to_bigwig/bigwig_file


steps:

  filter_sam:
    run: ../tools/custom-bash.cwl
    in:
      input_file: sam_file
      script:
        default: 'cat "$0" | grep "$1" | sed "s/$1/chr/g"  > `basename $0`'
      param:
        source: current_strain_name
        valueFrom: $(self+"_")
    out: [output_file]

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: filter_sam/output_file
      sort_output_filename:
        source: [output_file_prefix, current_strain_name]
        valueFrom: $(self.join("_")+".bam")
      threads: threads
    out: [bam_bai_pair]

  bam_to_bedgraph:
    run: ../tools/bedtools-genomecov.cwl
    in:
      input_file: samtools_sort_index/bam_bai_pair
      depth:
        default: "-bg"
      split:
        default: true
      mapped_reads_number: mapped_reads_number
    out: [genome_coverage_file]

  sort_bedgraph:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: bam_to_bedgraph/genome_coverage_file
      key:
        default: ["1,1","2,2n"]
    out: [sorted_file]

  project_bedgraph:
    run: ../tools/halliftover.cwl
    in:
      input_bed_file: sort_bedgraph/sorted_file
      hal_file: hal_file
      source_genome_name: current_strain_name
      target_genome_name: reference_strain_name
    out: [projected_bed_file]

  filter_projected_bedgraph:
    run: ../tools/custom-bash.cwl
    in:
      input_file: project_bedgraph/projected_bed_file
      script:
        default: 'cat "$0" | grep "$1" > `basename $0`'
      param:
        default: "chr"
    out: [output_file]

  sort_filtered_bedgraph:
    run: ../tools/custom-bedops.cwl
    in:
      input_file: filter_projected_bedgraph/output_file
      script:
        default: 'sort-bed "$0" > `basename $0`'
    out: [output_file]

  remove_overlaps:
    run: ../tools/custom-bedops.cwl
    in:
      input_file: sort_filtered_bedgraph/output_file
      script:
        default: bedops --partition "$0" | bedmap --echo --echo-map-id --delim "\t" - "$0" | awk '{n=split($4,reads,";"); sum=0; for(i=1;i<=n;i++) sum+=reads[i]; $4=sum; print $0 }' > `basename $0`
    out: [output_file]

  bedgraph_to_bigwig:
    run: ../tools/ucsc-bedgraphtobigwig.cwl
    in:
      bedgraph_file: remove_overlaps/output_file
      chrom_length_file: chrom_length_file
    out: [bigwig_file]