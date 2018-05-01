cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement


inputs:

  fastq_file:
    type: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format"

  star_indices_folder:
    type: Directory
    label: "STAR indices folder for reference genome"
    doc: "Path to STAR generated indices for reference genome"

  insilico_star_indices_folder:
    type: Directory
    label: "STAR indices folder for insilico genome"
    doc: "Path to STAR generated indices for insilico genome"

  chrom_length_file:
    type: File
    label: "Chromosome length file for reference genome"
    format: "http://edamontology.org/format_2330"
    doc: "Chromosome length file for reference genome"

  strain1:
    type: string
    label: "I strain name"
    doc: "First strain name"

  strain2:
    type: string
    label: "II strain name"
    doc: "Second strain name"

  strain1_chain_file:
    type: File
    label: "I strain chain file"
    format: "http://edamontology.org/format_2330"
    doc: "Chain file to project strain I to reference genome"

  strain2_chain_file:
    type: File
    label: "II strain chain file"
    format: "http://edamontology.org/format_2330"
    doc: "Chain file to project strain II to reference genome"

  threads:
    type: int?
    default: 2
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"


outputs:

  strain1_bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "I strain bigWig file"
    doc: "Generated bigWig file for the first strain"
    outputSource: allele_vcf_alignreads_se_pe/strain1_bigwig

  strain2_bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "II strain bigWig file"
    doc: "Generated bigWig file for the second strain"
    outputSource: allele_vcf_alignreads_se_pe/strain2_bigwig

  reference_bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "Reference bigWig file"
    doc: "Generated BigWig file for the reference genome"
    outputSource: allele_vcf_alignreads_se_pe/reference_bigwig

  strain1_bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Strain I coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file for strain I"
    outputSource: allele_vcf_alignreads_se_pe/strain1_bambai_pair

  strain2_bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Strain II coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file for strain II"
    outputSource: allele_vcf_alignreads_se_pe/strain2_bambai_pair

  reference_bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Reference coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file for reference genome"
    outputSource: allele_vcf_alignreads_se_pe/reference_bambai_pair

  insilico_star_final_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "STAR final log for insilico genome"
    doc: "STAR Log.final.out for insilico genome"
    outputSource: allele_vcf_alignreads_se_pe/insilico_star_final_log

  insilico_star_out_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR log out for insilico genome"
    doc: "STAR Log.out for insilico genome"
    outputSource: allele_vcf_alignreads_se_pe/insilico_star_out_log

  insilico_star_progress_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR progress log for insilico genome"
    doc: "STAR Log.progress.out for insilico genome"
    outputSource: allele_vcf_alignreads_se_pe/insilico_star_progress_log

  insilico_star_stdout_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR stdout log for insilico genome"
    doc: "STAR Log.std.out for insilico genome"
    outputSource: allele_vcf_alignreads_se_pe/insilico_star_stdout_log

  reference_star_final_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "STAR final log for reference genome"
    doc: "STAR Log.final.out for reference genome"
    outputSource: allele_vcf_alignreads_se_pe/reference_star_final_log

  reference_star_out_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR log out for reference genome"
    doc: "STAR Log.out for reference genome"
    outputSource: allele_vcf_alignreads_se_pe/reference_star_out_log

  reference_star_progress_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR progress log for reference genome"
    doc: "STAR Log.progress.out for reference genome"
    outputSource: allele_vcf_alignreads_se_pe/reference_star_progress_log

  reference_star_stdout_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR stdout log for reference genome"
    doc: "STAR Log.std.out for reference genome"
    outputSource: allele_vcf_alignreads_se_pe/reference_star_stdout_log

steps:

  extract_fastq:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file
    out: [fastq_file]

  allele_vcf_alignreads_se_pe:
    run: ./allele-vcf-alignreads-se-pe.cwl
    in:
      fastq_files: extract_fastq/fastq_file
      insilico_star_indices_folder: insilico_star_indices_folder
      reference_star_indices_folder: star_indices_folder
      reference_chrom_length_file: chrom_length_file
      strain1: strain1
      strain2: strain2
      strain1_chain_file: strain1_chain_file
      strain2_chain_file: strain2_chain_file
      threads: threads
    out:
    - strain1_bambai_pair
    - strain2_bambai_pair
    - reference_bambai_pair
    - strain1_bigwig
    - strain2_bigwig
    - reference_bigwig
    - insilico_star_final_log
    - insilico_star_out_log
    - insilico_star_progress_log
    - insilico_star_stdout_log
    - reference_star_final_log
    - reference_star_out_log
    - reference_star_progress_log
    - reference_star_stdout_log