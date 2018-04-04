cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement


'sd:metadata':
  - "https://raw.githubusercontent.com/Barski-lab/workflows/master/metadata/rnaseq-header.cwl"


inputs:

  fastq_file:
    type: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format"

  star_indices_folder:
    type: Directory
    label: "STAR indices folder for reference genome"
    doc: "Path to STAR generated indices for reference genome"

  insilico_star_indices_folder:
    type: Directory
    label: "STAR indices folder for insilico genome"
    doc: "Path to STAR generated indices for insilico genome"

  chrom_length_file:
    type: File
    label: "Chromosome length file for reference genome"
    format: "http://edamontology.org/format_2330"
    doc: "Chromosome length file for reference genome"

  strain1:
    type: string
    label: "I strain name"
    doc: "First strain name"

  strain2:
    type: string
    label: "II strain name"
    doc: "Second strain name"

  strain1_refmap_file:
    type: File
    label: "I strain refmap file"
    format: "http://edamontology.org/format_2330"
    doc: "Refmap file generated while making strain specific indices for the first strain"

  strain2_refmap_file:
    type: File
    label: "II strain refmap file"
    format: "http://edamontology.org/format_2330"
    doc: "Refmap file generated while making strain specific indices for the second strain"

  threads:
    type: int?
    default: 2
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"


outputs:

  strain1_bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "I strain bigWig file"
    doc: "Generated bigWig file for the first strain"
    outputSource: strain1_mea_createtracks/bigwig_file

  strain2_bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "II strain bigWig file"
    doc: "Generated bigWig file for the second strain"
    outputSource: strain2_mea_createtracks/bigwig_file

  reference_bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "Reference bigWig file"
    doc: "Generated BigWig file for the reference genome"
    outputSource: reference_bam_to_bigwig/bigwig_file

  insilico_star_final_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "STAR final log for insilico genome"
    doc: "STAR Log.final.out for insilico genome"
    outputSource: mea_alignreads_star_se_pe/insilico_star_final_log

  insilico_star_out_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR log out for insilico genome"
    doc: "STAR Log.out for insilico genome"
    outputSource: mea_alignreads_star_se_pe/insilico_star_out_log

  insilico_star_progress_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR progress log for insilico genome"
    doc: "STAR Log.progress.out for insilico genome"
    outputSource: mea_alignreads_star_se_pe/insilico_star_progress_log

  insilico_star_stdout_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR stdout log for insilico genome"
    doc: "STAR Log.std.out for insilico genome"
    outputSource: mea_alignreads_star_se_pe/insilico_star_stdout_log

  reference_star_final_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "STAR final log for reference genome"
    doc: "STAR Log.final.out for reference genome"
    outputSource: mea_alignreads_star_se_pe/reference_star_final_log

  reference_star_out_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR log out for reference genome"
    doc: "STAR Log.out for reference genome"
    outputSource: mea_alignreads_star_se_pe/reference_star_out_log

  reference_star_progress_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR progress log for reference genome"
    doc: "STAR Log.progress.out for reference genome"
    outputSource: mea_alignreads_star_se_pe/reference_star_progress_log

  reference_star_stdout_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR stdout log for reference genome"
    doc: "STAR Log.std.out for reference genome"
    outputSource: mea_alignreads_star_se_pe/reference_star_stdout_log


steps:

  extract_fastq:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file
    out: [fastq_file]

  mea_alignreads_star_se_pe:
    run: ./mea-alignreads-star-se-pe.cwl
    in:
      fastq_files: extract_fastq/fastq_file
      insilico_star_indices_folder: insilico_star_indices_folder
      reference_star_indices_folder: star_indices_folder
      strain1: strain1
      strain2: strain2
      threads: threads
    out:
    - strain1_sorted_bam
    - strain2_sorted_bam
    - reference_sorted_bam
    - reference_uniquely_mapped_reads_number
    - insilico_star_final_log
    - insilico_star_out_log
    - insilico_star_progress_log
    - insilico_star_stdout_log
    - reference_star_final_log
    - reference_star_out_log
    - reference_star_progress_log
    - reference_star_stdout_log

  strain1_mea_createtracks:
    run: ./mea-createtracks.cwl
    in:
      bam_file: mea_alignreads_star_se_pe/strain1_sorted_bam
      reference_uniquely_mapped_reads_number: mea_alignreads_star_se_pe/reference_uniquely_mapped_reads_number
      reference_chrom_length_file: chrom_length_file
      refmap_file: strain1_refmap_file
    out: [bigwig_file]

  strain2_mea_createtracks:
    run: ./mea-createtracks.cwl
    in:
      bam_file: mea_alignreads_star_se_pe/strain2_sorted_bam
      reference_uniquely_mapped_reads_number: mea_alignreads_star_se_pe/reference_uniquely_mapped_reads_number
      reference_chrom_length_file: chrom_length_file
      refmap_file: strain2_refmap_file
    out: [bigwig_file]

  reference_bam_to_bigwig:
    run: bam-bedgraph-bigwig.cwl
    in:
      bam_file: mea_alignreads_star_se_pe/reference_sorted_bam
      chrom_length_file: chrom_length_file
      mapped_reads_number: mea_alignreads_star_se_pe/reference_uniquely_mapped_reads_number
    out: [bigwig_file]