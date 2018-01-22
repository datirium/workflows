cwlVersion: v1.0
class: Workflow


inputs:

  sra_input_file:
    type: File
  illumina_adapters_file:
    type: File
  threads:
    type: int?


outputs:

  fastq:
    type: File
    outputSource: trimmomatic/output_read1_trimmed_file


steps:

  sra_to_fastq:
    run: ../tools/fastq-dump.cwl
    in:
      sra_file: sra_input_file
    out: [fastq_file_1]

  fastqc:
    run: ../tools/fastqc.cwl
    in:
      fastq_file: sra_to_fastq/fastq_file_1
    out: [summary_file]

  fastqc_results_trigger:
    run: ../expressiontools/fastqc-results-trigger.cwl
    in:
      summary_file: fastqc/summary_file
    out: [trigger]

  trimmomatic:
    run: ../tools/trimmomatic.cwl
    in:
      input_read1_fastq_file: sra_to_fastq/fastq_file_1
      input_adapters_file: illumina_adapters_file
      trigger: fastqc_results_trigger/trigger
      end_mode:
        default: "SE"
      illuminaclip:
        default: '2:30:15'
      threads: threads
    out: [output_read1_trimmed_file]

