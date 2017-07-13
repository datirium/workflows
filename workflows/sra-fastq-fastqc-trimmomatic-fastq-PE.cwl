cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
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
    outputSource: trimmomatic/output_read1_trimmed_file
  downstream_fastq:
    type: File
    outputSource: trimmomatic/output_read2_trimmed_paired_file

steps:

  sra_to_fastq:
    run: ../tools/fastq-dump.cwl
    in:
      input_file: sra_input_file
      split_files:
        default: true
    out: [output_file_1, output_file_2]

  fastqc_1:
    run: ../tools/fastqc.cwl
    in:
      fastq_file: sra_to_fastq/output_file_1
    out: [summary_file]

  fastqc_2:
    run: ../tools/fastqc.cwl
    in:
      fastq_file: sra_to_fastq/output_file_2
    out: [summary_file]


  pick_file_from_array_1:
    run: ../expressiontools/get-file-by-name.cwl
    in:
      input_files: fastqc_1/summary_file
      basename_regex:
        source: sra_to_fastq/output_file_1
        valueFrom: |
          ${
            return self.basename.split(".")[0] + "_1.fastq"
          }
    out: [selected_file]

  pick_file_from_array_2:
    run: ../expressiontools/get-file-by-name.cwl
    in:
      input_files: fastqc_2/summary_file
      basename_regex:
        source: sra_to_fastq/output_file_2
        valueFrom: |
          ${
            return self.basename.split(".")[0] + "_2.fastq"
          }
    out: [selected_file]

  fastqc_results_trigger_1:
    run: ../expressiontools/fastqc-results-trigger.cwl
    in:
      summary: pick_file_from_array_1/selected_file
    out: [trigger]

  fastqc_results_trigger_2:
    run: ../expressiontools/fastqc-results-trigger.cwl
    in:
      summary: pick_file_from_array_2/selected_file
    out: [trigger]

  trimmomatic:
    run: ../tools/trimmomatic.cwl
    in:
      input_read1_fastq_file: sra_to_fastq/output_file_1
      input_read2_fastq_file: sra_to_fastq/output_file_2
      input_adapters_file: illumina_adapters_file
      trigger:
        source: [fastqc_results_trigger_1/trigger, fastqc_results_trigger_2/trigger]
        valueFrom: |
          ${
              return self[0] && self[1];
          }
      end_mode:
        default: "PE"
      illuminaclip:
        default: '2:30:15'
      threads: threads
    out: [output_read1_trimmed_file, output_read2_trimmed_paired_file]