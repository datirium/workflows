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

outputs:
  fastq:
    type: File
    outputSource: trimmomatic/output_read1_trimmed_file

steps:
  sra_to_fastq:
    run: ../../tools/fastq-dump.cwl
    in:
      inputFiles: sra_input_file
    out: [outputFile]

  fastqc:
    run: ../../tools/fastqc.cwl
    in:
      fastqFile: sra_to_fastq/outputFile
    out: [summary_file]

  fastqc_results_trigger:
    run: ../../expressiontools/fastqc-results-trigger.cwl
    in:
      summary: fastqc/summary_file
    out: [trigger]

  trimmomatic:
    run: ../../tools/trimmomatic.cwl
    in:
      input_read1_fastq_file: sra_to_fastq/outputFile
      input_adapters_file: illumina_adapters_file
      trigger: fastqc_results_trigger/trigger
      end_mode:
        default: "SE"
      illuminaclip:
        default: '2:30:15'
    out: [output_read1_trimmed_file]

