cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


inputs:

  sra_input_file:
    type: File
  illumina_adapters_file:
    type: File
  threads:
    type: int?


outputs:

  upstream_fastq:
    type: File
    outputSource: trimmomatic/upstream_trimmed_file
  downstream_fastq:
    type: File
    outputSource: trimmomatic/downstream_trimmed_file


steps:

  sra_to_fastq:
    run: ../tools/fastq-dump.cwl
    in:
      sra_file: sra_input_file
      split_files:
        default: true
    out: [fastq_file_1, fastq_file_2]

  fastqc_1:
    run: ../tools/fastqc.cwl
    in:
      fastq_file: sra_to_fastq/fastq_file_1
    out: [summary_file]

  fastqc_2:
    run: ../tools/fastqc.cwl
    in:
      fastq_file: sra_to_fastq/fastq_file_2
    out: [summary_file]

  fastqc_results_trigger_1:
    run: ../expressiontools/fastqc-results-trigger.cwl
    in:
      summary_file: fastqc_1/summary_file
    out: [trigger]

  fastqc_results_trigger_2:
    run: ../expressiontools/fastqc-results-trigger.cwl
    in:
      summary_file: fastqc_2/summary_file
    out: [trigger]

  trimmomatic:
    run: ../tools/trimmomatic.cwl
    in:
      fastq_file_upstream: sra_to_fastq/fastq_file_1
      fastq_file_downstream: sra_to_fastq/fastq_file_2
      adapters_file: illumina_adapters_file
      trigger:
        source: [fastqc_results_trigger_1/trigger, fastqc_results_trigger_2/trigger]
        valueFrom: |
          ${
              return self[0] && self[1];
          }
      lib_type:
        default: "PE"
      illuminaclip_step_param:
        default: '2:30:15'
      threads: threads
    out: [upstream_trimmed_file, downstream_trimmed_file]