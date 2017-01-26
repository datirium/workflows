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
  rsem_reference_name_dir:
    type: Directory
  rsem_reference_name:
    type: string
  rsem_star:
    type: boolean?
  rsem_bowtie2:
    type: boolean?
  chrLengthFile:
    type: File


outputs:
  rsem_isoform_results:
    type: File
    outputSource: rsem_calculate_expression/isoform_results
  rsem_gene_results:
    type: File
    outputSource: rsem_calculate_expression/gene_results
  rsem_genome_sorted_bam_bai_pair:
    type: File
    outputSource: rsem_calculate_expression/genome_sorted_bam_bai_pair
  bigwig_outfile:
    type: File
    outputSource: bamToBigwig/outfile

steps:

  sra_to_fastq:
    run: ../../tools/fastq-dump.cwl
    in:
      inputFiles: sra_input_file
      split_files:
        default: true
    out: [outputFile, outputFile_2]

  fastqc_1:
    run: ../../tools/fastqc.cwl
    in:
      fastqFile: sra_to_fastq/outputFile
    out: [summary_file]

  fastqc_2:
    run: ../../tools/fastqc.cwl
    in:
      fastqFile: sra_to_fastq/outputFile_2
    out: [summary_file]

  parse_fastqc_results_1:
    run: ../../tools/parse-fastqc-results.cwl
    in:
      fastqc_report: fastqc_1/summary_file
    out: [output]

  parse_fastqc_results_2:
    run: ../../tools/parse-fastqc-results.cwl
    in:
      fastqc_report: fastqc_2/summary_file
    out: [output]

  trimmomatic:
    run: ../../tools/trimmomatic.cwl
    in:
      input_read1_fastq_file: sra_to_fastq/outputFile
      input_read2_fastq_file: sra_to_fastq/outputFile_2
      input_adapters_file: illumina_adapters_file
      end_mode:
        default: "PE"
      illuminaclip:
        default: '2:30:15'
    out: [output_read1_trimmed_file, output_read2_trimmed_paired_file]

  rsem_calculate_expression:
    run: ../../tools/rsem-calculate-expression.cwl
    in:
      upstream_read_file:
        source: [sra_to_fastq/outputFile, trimmomatic/output_read1_trimmed_file, parse_fastqc_results_1/output, parse_fastqc_results_2/output]
        valueFrom: |
          ${
            if (self.slice(2)[0].basename == "fail.txt" || self.slice(3)[0].basename == "fail.txt"){
              return self.slice(1,2);
            } else {
              return self.slice(0,1);
            }
          }
      downstream_read_file:
        source: [sra_to_fastq/outputFile_2, trimmomatic/output_read2_trimmed_paired_file, parse_fastqc_results_1/output, parse_fastqc_results_2/output]
        valueFrom: |
          ${
            if (self.slice(2)[0].basename == "fail.txt" || self.slice(3)[0].basename == "fail.txt"){
              return self.slice(1,2);
            } else {
              return self.slice(0,1);
            }
          }
      reference_name_dir: rsem_reference_name_dir
      reference_name: rsem_reference_name
      star: rsem_star
      bowtie2: rsem_bowtie2
      sort_bam_by_coordinate:
        default: true
      output_genome_bam:
        default: true
      sample_name:
        source: sra_input_file
        valueFrom: |
          ${
            return self.basename.split('.')[0];
          }
    out: [isoform_results, gene_results, genome_sorted_bam_bai_pair]

  bamtoolsStats:
    run: ../../tools/bamtools-stats.cwl
    in:
      inputFiles: rsem_calculate_expression/genome_sorted_bam_bai_pair
    out: [mappedreads]

  bamToBigwig:
    run: bam-genomecov-bigwig.cwl
    in:
      input: rsem_calculate_expression/genome_sorted_bam_bai_pair
      genomeFile: chrLengthFile
      mappedreads: bamtoolsStats/mappedreads
      bigWig:
        source: sra_input_file
        valueFrom: |
          ${
            return self.basename.split('.')[0]+".bigwig";
          }
    out: [outfile]
