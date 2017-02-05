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
    label: "SRA archive file"
    format: "http://edamontology.org/format_3698"
    doc: "SDLoadFile"
  illumina_adapters_file:
    type: File
    label: "Illumina adapters file"
    format: "http://edamontology.org/format_1929"
    doc: "SDIlluminaAdapters"
  rsem_reference_name_dir:
    type: Directory
    label: "RSEM references folder"
    doc: "SDRsemRef"
#  rsem_reference_name:
#    type: string?
  rsem_star:
    type: boolean?
  rsem_bowtie2:
    type: boolean?
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
  upstream_fastq:
    type: File
    outputSource: sra_fastqc_trimmomatic_fastq_PE/upstream_fastq
  downstream_fastq:
    type: File
    outputSource: sra_fastqc_trimmomatic_fastq_PE/downstream_fastq

steps:

  sra_fastqc_trimmomatic_fastq_PE:
    run: sra-fastq-fastqc-trimmomatic-fastq-PE.cwl
    in:
      sra_input_file: sra_input_file
      illumina_adapters_file: illumina_adapters_file
    out: [upstream_fastq, downstream_fastq]

  fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE:
    run: fastq-rsem-bam-bigwig-tpmcountmatrices-SE-PE.cwl
    in:
      upstream_fastq: sra_fastqc_trimmomatic_fastq_PE/upstream_fastq
      downstream_fastq: sra_fastqc_trimmomatic_fastq_PE/downstream_fastq
      rsem_reference_name_dir: rsem_reference_name_dir
#      rsem_reference_name: rsem_reference_name
      rsem_star: rsem_star
      rsem_bowtie2: rsem_bowtie2
      chrLengthFile: chrLengthFile
    out: [rsem_isoform_results, rsem_gene_results, rsem_genome_sorted_bam_bai_pair, bigwig_outfile, isoforms_tpm_matrix, isoforms_counts_matrix, genes_tpm_matrix, genes_counts_matrix]
