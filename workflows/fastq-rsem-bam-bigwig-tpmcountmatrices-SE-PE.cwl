cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement

inputs:
  upstream_fastq:
    type: File
  downstream_fastq:
    type: File?
  rsem_reference_name_dir:
    type: Directory
  chrLengthFile:
    type: File
  aligner_type:
    type:
      name: "aligner_type"
      type: enum
      symbols: ["bowtie","star","bowtie2"]
  threads:
    type: int?

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
  bam_quality_log:
    type: File
    outputSource: bamtoolsStats/stats_log


steps:

  rsem_calculate_expression:
    run: ../tools/rsem-calculate-expression.cwl
    in:
      upstream_read_file: upstream_fastq
      downstream_read_file: downstream_fastq
      reference_name_dir: rsem_reference_name_dir
      star:
        source: aligner_type
        valueFrom: |
          ${
           if (self == "star"){
             return true;
           } else {
             return false;
           }
          }
      bowtie2:
        source: aligner_type
        valueFrom: |
          ${
           if (self == "bowtie2"){
             return true;
           } else {
             return false;
           }
          }
      sort_bam_by_coordinate:
        default: true
      output_genome_bam:
        default: true
      sample_name:
        source: upstream_fastq
        valueFrom: |
          ${
            return self.basename.split('.')[0];
          }
      num_threads: threads
    out: [isoform_results, gene_results, genome_sorted_bam_bai_pair]

  bamtoolsStats:
    run: ../tools/bamtools-stats.cwl
    in:
      input_files: rsem_calculate_expression/genome_sorted_bam_bai_pair
    out: [mappedreads, stats_log]

  bamToBigwig:
    run: bam-genomecov-bigwig.cwl
    in:
      input: rsem_calculate_expression/genome_sorted_bam_bai_pair
      genomeFile: chrLengthFile
      mappedreads: bamtoolsStats/mappedreads
      bigWig:
        source: upstream_fastq
        valueFrom: |
          ${
            return self.basename.split('.')[0]+".bigwig";
          }
    out: [outfile]