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
#  rsem_reference_name:
#    type: string?
  rsem_aligner_type:
    type:
      name: "aligner_type"
      type: enum
      symbols: ["bowtie","star","bowtie2"]
  chrLengthFile:
    type: File

outputs:
  rsem_isoform_results:
    type: File
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/rsem_isoform_results
  rsem_gene_results:
    type: File
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/rsem_gene_results
  rsem_genome_sorted_bam_bai_pair:
    type: File
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/rsem_genome_sorted_bam_bai_pair
  bigwig_outfile:
    type: File
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/bigwig_outfile
  isoforms_tpm_matrix:
    type: File
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/isoforms_tpm_matrix
  isoforms_counts_matrix:
    type: File
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/isoforms_counts_matrix
  genes_tpm_matrix:
    type: File
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/genes_tpm_matrix
  genes_counts_matrix:
    type: File
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/genes_counts_matrix
  fastq:
    type: File
    outputSource: sra_fastqc_trimmomatic_fastq_SE/fastq
  bam_quality_log:
    type: File
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/bam_quality_log

steps:

  sra_fastqc_trimmomatic_fastq_SE:
    run: sra-fastq-fastqc-trimmomatic-fastq-SE.cwl
    in:
      sra_input_file: sra_input_file
      illumina_adapters_file: illumina_adapters_file
    out: [fastq]

  fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE:
    run: fastq-rsem-bam-bigwig-tpmcountmatrices-SE-PE.cwl
    in:
      upstream_fastq: sra_fastqc_trimmomatic_fastq_SE/fastq
      rsem_reference_name_dir: rsem_reference_name_dir
#      rsem_reference_name: rsem_reference_name
      aligner_type:
        sour: rsem_aligner_type
        valueFrom: $(self)
      chrLengthFile: chrLengthFile
    out: [rsem_isoform_results, rsem_gene_results, rsem_genome_sorted_bam_bai_pair, bigwig_outfile, isoforms_tpm_matrix, isoforms_counts_matrix, genes_tpm_matrix, genes_counts_matrix, bam_quality_log]
