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
    type: File[]
    outputSource: rsem_rna_pe/rsem_isoform_results
  rsem_gene_results:
    type: File[]
    outputSource: rsem_rna_pe/rsem_gene_results
  rsem_genome_sorted_bam_bai_pair:
    type: File[]
    outputSource: rsem_rna_pe/rsem_genome_sorted_bam_bai_pair
  bigwig_outfile:
    type: File[]
    outputSource: rsem_rna_pe/bigwig_outfile
  isoforms_tpm_matrix:
    type: File[]
    outputSource: rsem_rna_pe/isoforms_tpm_matrix
  isoforms_counts_matrix:
    type: File[]
    outputSource: rsem_rna_pe/isoforms_counts_matrix
  genes_tpm_matrix:
    type: File[]
    outputSource: rsem_rna_pe/genes_tpm_matrix
  genes_counts_matrix:
    type: File[]
    outputSource: rsem_rna_pe/genes_counts_matrix

steps:
  rsem_rna_pe:
    run: rsem-rna-pe.cwl
    in:
      sra_input_file: sra_input_files
      illumina_adapters_file: illumina_adapters_file
      rsem_reference_name_dir: rsem_reference_name_dir
      rsem_reference_name: rsem_reference_name
      rsem_star: rsem_star
      rsem_bowtie2: rsem_bowtie2
      chrLengthFile: chrLengthFile
    scatter: sra_input_file
    out: [rsem_isoform_results, rsem_gene_results, rsem_genome_sorted_bam_bai_pair, bigwig_outfile, isoforms_tpm_matrix, isoforms_counts_matrix, genes_tpm_matrix, genes_counts_matrix]
