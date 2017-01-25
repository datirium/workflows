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

steps:
  rsem_rna_se:
    run: rsem-rna-se.cwl
    in:
      sra_input_file: sra_input_files
      illumina_adapters_file: illumina_adapters_file
      rsem_reference_name_dir: rsem_reference_name_dir
      rsem_reference_name: rsem_reference_name
      rsem_star: rsem_star
      rsem_bowtie2: rsem_bowtie2
      chrLengthFile: chrLengthFile
    scatter: sra_input_file
    out: [rsem_isoform_results,rsem_gene_results, rsem_genome_sorted_bam_bai_pair, bigwig_outfile]
