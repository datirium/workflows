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
  rsem_aligner_type:
    type:
      name: "aligner_type"
      type: enum
      symbols: ["bowtie","star","bowtie2"]
    label: "Aligner STAR|Bowtie|Bowtie2"
    doc: "SDSelectAligner"
  chrLengthFile:
    type: File
    label: "Chromosome length file"
    doc: "SDSelectChrLengthFile"

outputs:
  rsem_isoform_results:
    type: File
    label: "Isoform level expression estimates"
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/rsem_isoform_results
  rsem_gene_results:
    type: File
    label: "Gene level expression estimates"
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/rsem_gene_results
  rsem_genome_sorted_bam_bai_pair:
    type: File
    label: "Coordinate sorted BAM(+BAI) file"
    format: http://edamontology.org/format_2572
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/rsem_genome_sorted_bam_bai_pair
  bigwig_outfile:
    type: File
    label: "BigWig file"
    format: http://edamontology.org/format_3006
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/bigwig_outfile
  isoforms_tpm_matrix:
    type: File
    label: "Isoform/TPM matrix file"
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/isoforms_tpm_matrix
  isoforms_counts_matrix:
    type: File
    label: "Isoform/ReadsCount matrix file"
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/isoforms_counts_matrix
  genes_tpm_matrix:
    type: File
    label: "Gene/TPM matrix file"
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/genes_tpm_matrix
  genes_counts_matrix:
    type: File
    label: "Gene/ReadsCount matrix file"
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/genes_counts_matrix
  upstream_fastq:
    type: File
    label: "Upstream FASTQ file"
    format: http://edamontology.org/format_1930
    outputSource: sra_fastqc_trimmomatic_fastq_PE/upstream_fastq
  downstream_fastq:
    type: File
    label: "Downstream FASTQ file"
    format: http://edamontology.org/format_1930
    outputSource: sra_fastqc_trimmomatic_fastq_PE/downstream_fastq
  bam_quality_log:
    type: File
    label: "BAM statistics file"
    outputSource: fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE/bam_quality_log


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
      aligner_type:
        source: rsem_aligner_type
        valueFrom: $(self)
      chrLengthFile: chrLengthFile
    out: [rsem_isoform_results, rsem_gene_results, rsem_genome_sorted_bam_bai_pair, bigwig_outfile, isoforms_tpm_matrix, isoforms_counts_matrix, genes_tpm_matrix, genes_counts_matrix, bam_quality_log]



$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "sra-fastq-fastqc-trimmomatic-fastq-rsem-bam-bigwig-tpmcountmatrices-PE"
s:downloadUrl: https://github.com/SciDAP/workflows/blob/master/workflows/scidap/sra-fastq-fastqc-trimmomatic-fastq-rsem-bam-bigwig-tpmcountmatrices-PE.cwl
s:codeRepository: https://github.com/SciDAP/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

s:about: >
  Current workflow should be used only with the paired-end RNA-Seq data obtained from Illumina sequencing machine. It performs the following steps:
  1. Convert input SRA file into the pair of FASTQ files (run fastq-dump)
  3. Check the quality of generated upstream and downstream FASTQ files (run fastqc)
  4. If fastqc returned FAIL for "Per base sequence quality" or "Per sequence quality scores" or "Overrepresented sequences" for upstream or downstream
     fastq file, run Trimmomatic for both of the files. If files passed all of the fastqc test parameters, continue to work with the original FASTQ files.
  5. Run rsem-calculate-expression utility. Return isoform and gene expression files, a pair of coordinate sorted BAM and index BAI files.
  6. Calculate basic statistics for BAM file (bamtools stats). Return file with calculated statistics.
  7. Generate BigWig file on the base of BAM file
  8. Calculate isoform TPM, isoform raw reads count, gene TPM, gene raw reads count matrices
