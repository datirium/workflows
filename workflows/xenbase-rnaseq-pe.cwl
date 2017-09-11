cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement


inputs:
  sra_input_file:
    type: File
  illumina_adapters_file:
    type: File
  rsem_indices_folder:
    type: Directory
  chr_length_file:
    type: File
  threads:
    type: int?


outputs:

  rsem_isoforms:
    type: File
    outputSource: fastq_to_bigwig/rsem_isoforms
  rsem_genes:
    type: File
    outputSource: fastq_to_bigwig/rsem_genes
  bam_file:
    type: File
    outputSource: fastq_to_bigwig/bam_file
  bamtools_log:
    type: File
    outputSource: fastq_to_bigwig/bamtools_log
  bed:
    type: File
    outputSource: fastq_to_bigwig/bed
  bigwig:
    type: File
    outputSource: fastq_to_bigwig/bigwig


steps:

  sra_to_fastq:
    run: xenbase-sra-to-fastq-pe.cwl
    in:
      sra_input_file: sra_input_file
      illumina_adapters_file: illumina_adapters_file
      threads: threads
    out:
      - upstream_fastq
      - downstream_fastq

  fastq_to_bigwig:
    run: xenbase-fastq-rsem-bigwig-se-pe.cwl
    in:
      upstream_fastq: sra_to_fastq/upstream_fastq
      downstream_fastq: sra_to_fastq/downstream_fastq
      rsem_indices_folder: rsem_indices_folder
      chr_length_file: chr_length_file
      threads: threads
    out:
      - rsem_isoforms
      - rsem_genes
      - bam_file
      - bamtools_log
      - bed
      - bigwig


$namespaces:
  s: http://schema.org/
  dct: http://purl.org/dc/terms/
  foaf: http://xmlns.com/foaf/0.1/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "xenbase-rnaseq-pe"
s:downloadUrl: https://raw.githubusercontent.com/SciDAP/workflows/master/workflows/xenbase-rnaseq-pe.cwl
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

dct:creator:
  "@id": "http://orcid.org/0000-0002-6486-3898"
  foaf:name: "Michael Kotliar"
  foaf:mbox: "mailto:michael.kotliar@cchmc.org"

doc: >
  XenBase workflow for analysing RNA-Seq paired-end data

s:about: |
  1. Convert input SRA file into pair of upsrtream and downstream FASTQ files (run fastq-dump)
  2. Analyze quality of FASTQ files (run fastqc with each of the FASTQ files)
  3. If any of the following fields in fastqc generated report is marked as failed for at least one of input FASTQ files:
        "Per base sequence quality",
        "Per sequence quality scores",
        "Overrepresented sequences",
        "Adapter Content",
    - trim adapters (run trimmomatic)
  4. Align original or trimmed FASTQ files to reference genome, calculate genes and isoforms expression (run RSEM)
  5. Count mapped reads number in sorted BAM file (run bamtools stats)
  6. Generate genome coverage BED file (run bedtools genomecov)
  7. Sort genearted BED file (run sort)
  8. Generate genome coverage bigWig file from BED file (run bedGraphToBigWig)