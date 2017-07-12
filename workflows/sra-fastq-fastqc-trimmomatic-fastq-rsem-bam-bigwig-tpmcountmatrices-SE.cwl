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
  rsem_aligner_type:
    type:
      name: "aligner_type"
      type: enum
      symbols: ["bowtie","star","bowtie2"]
  chrLengthFile:
    type: File
  threads:
    type: int?
    label: "Number of threads to run tools"
    doc: "SDInput"

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
      threads: threads
    out: [fastq]

  fastq_rsem_bam_bigwig_tpmcountmatrices_SE_PE:
    run: fastq-rsem-bam-bigwig-tpmcountmatrices-SE-PE.cwl
    in:
      upstream_fastq: sra_fastqc_trimmomatic_fastq_SE/fastq
      rsem_reference_name_dir: rsem_reference_name_dir
      aligner_type:
        source: rsem_aligner_type
        valueFrom: $(self)
      chrLengthFile: chrLengthFile
      threads: threads
    out: [rsem_isoform_results, rsem_gene_results, rsem_genome_sorted_bam_bai_pair, bigwig_outfile, bam_quality_log]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "sra-fastq-fastqc-trimmomatic-fastq-rsem-bam-bigwig-tpmcountmatrices-PE"
s:downloadUrl: https://github.com/SciDAP/workflows/blob/master/workflows/scidap/sra-fastq-fastqc-trimmomatic-fastq-rsem-bam-bigwig-tpmcountmatrices-SE.cwl
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
  Current workflow should be used only with the single-end RNA-Seq data obtained from Illumina sequencing machine. It performs the following steps:
  1. Convert input SRA file into FASTQ file (run fastq-dump)
  3. Check the quality of generated FASTQ file (run fastqc)
  4. If fastqc returned FAIL for "Per base sequence quality" or "Per sequence quality scores" or "Overrepresented sequences" or "Adapter Content" for fastq file,
     run Trimmomatic. If file passed all of the fastqc test parameters, continue to work with the original FASTQ file.
  5. Run rsem-calculate-expression utility. Return isoform and gene expression files, a pair of coordinate sorted BAM and index BAI files.
  6. Calculate basic statistics for BAM file (bamtools stats). Return file with calculated statistics.
  7. Generate BigWig file on the base of BAM file