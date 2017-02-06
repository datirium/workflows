class: Workflow
cwlVersion: v1.0
requirements:
  - class: ScatterFeatureRequirement
  - class: SubworkflowFeatureRequirement

inputs:
  sra_input_files:
    type: File[]
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
    type: File[]
    outputSource: rsem_rna_se/rsem_isoform_results
  rsem_gene_results:
    type: File[]
    outputSource: rsem_rna_se/rsem_gene_results
  rsem_genome_sorted_bam_bai_pair:
    type: File[]
    outputSource: rsem_rna_se/rsem_genome_sorted_bam_bai_pair
  bigwig_outfile:
    type: File[]
    outputSource: rsem_rna_se/bigwig_outfile
  isoforms_tpm_matrix:
    type: File[]
    outputSource: rsem_rna_se/isoforms_tpm_matrix
  isoforms_counts_matrix:
    type: File[]
    outputSource: rsem_rna_se/isoforms_counts_matrix
  genes_tpm_matrix:
    type: File[]
    outputSource: rsem_rna_se/genes_tpm_matrix
  genes_counts_matrix:
    type: File[]
    outputSource: rsem_rna_se/genes_counts_matrix
  fastq:
    type: File[]
    outputSource: rsem_rna_se/fastq

steps:
  rsem_rna_se:
    run: sra-fastq-fastqc-trimmomatic-fastq-rsem-bam-bigwig-tpmcountmatrices-SE.cwl
    in:
      sra_input_file: sra_input_files
      illumina_adapters_file: illumina_adapters_file
      rsem_reference_name_dir: rsem_reference_name_dir
#      rsem_reference_name: rsem_reference_name
      rsem_aligner_type:
        source: rsem_aligner_type
        valueFrom: $(self)
      chrLengthFile: chrLengthFile
    scatter: sra_input_file
    out: [rsem_isoform_results,rsem_gene_results, rsem_genome_sorted_bam_bai_pair, bigwig_outfile, isoforms_tpm_matrix, isoforms_counts_matrix, genes_tpm_matrix, genes_counts_matrix, fastq]
