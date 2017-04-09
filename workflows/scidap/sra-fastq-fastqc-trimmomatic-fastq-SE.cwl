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
  fastq:
    type: File
    outputSource: trimmomatic/output_read1_trimmed_file

steps:
  sra_to_fastq:
    run: ../../tools/fastq-dump.cwl
    in:
      input_file: sra_input_file
    out: [output_file_1]

  fastqc:
    run: ../../tools/fastqc.cwl
    in:
      fastq_file: sra_to_fastq/output_file_1
    out: [summary_file]

  pick_file_from_array:
    run: ../../expressiontools/get-file-by-name.cwl
    in:
     input_files: fastqc/summary_file
     basename_regex:
        source: sra_to_fastq/output_file_1
        valueFrom: |
          ${
            return self.basename.split(".")[0] + ".fastq"
          }
    out: [selected_file]

  fastqc_results_trigger:
    run: ../../expressiontools/fastqc-results-trigger.cwl
    in:
      summary: pick_file_from_array/selected_file
    out: [trigger]

  trimmomatic:
    run: ../../tools/trimmomatic.cwl
    in:
      input_read1_fastq_file: sra_to_fastq/output_file_1
      input_adapters_file: illumina_adapters_file
      trigger: fastqc_results_trigger/trigger
      end_mode:
        default: "SE"
      illuminaclip:
        default: '2:30:15'
      threads: threads
    out: [output_read1_trimmed_file]

