cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement

'sd:metadata':
  - "https://raw.githubusercontent.com/Barski-lab/workflows/master/metadata/rnaseq-header.cwl"

inputs:

# General inputs

  fastq_file:
    type: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format"

  star_indices_folder:
    type: Directory
    label: "STAR indices folder"
    'sd:parent': "https://raw.githubusercontent.com/Barski-lab/workflows/master/workflows/star-index.cwl"
    doc: "Path to STAR generated indices"

  bowtie_indices_folder:
    type: Directory
    label: "BowTie Ribosomal Indices"
    'sd:parent': "https://raw.githubusercontent.com/Barski-lab/workflows/master/workflows/bowtie-index.cwl"
    doc: "Path to Bowtie generated indices"

  chrom_length_file:
    type: File
    label: "Chromosome length file"
    format: "http://edamontology.org/format_2330"
    'sd:parent': "https://raw.githubusercontent.com/Barski-lab/workflows/master/workflows/star-index.cwl"
    doc: "Chromosome length file"

  annotation_file:
    type: File
    label: "Annotation file"
    format:
      - "http://edamontology.org/format_2306"
      - "http://edamontology.org/format_3475"
    'sd:parent': "https://raw.githubusercontent.com/Barski-lab/workflows/master/workflows/star-index.cwl"
    doc: "GTF or TAB-separated annotation file"

# Advanced inputs

  exclude_chr:
    type: string?
    'sd:layout':
      advanced: true
    label: "Chromosome to be excluded in rpkm calculation"
    doc: "Chromosome to be excluded in rpkm calculation"

  clip_3p_end:
    type: int?
    default: 0
    'sd:layout':
      advanced: true
    label: "Clip from 3p end"
    doc: "Number of bases to clip from the 3p end"

  clip_5p_end:
    type: int?
    default: 0
    'sd:layout':
      advanced: true
    label: "Clip from 5p end"
    doc: "Number of bases to clip from the 5p end"

# System dependent

  threads:
    type: int?
    default: 2
    'sd:layout':
      advanced: true
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"

outputs:

  bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated BigWig file"
    outputSource: bam_to_bigwig/bigwig_file

  star_final_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "STAR final log"
    doc: "STAR Log.final.out"
    outputSource: star_aligner/log_final

  star_out_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR log out"
    doc: "STAR Log.out"
    outputSource: star_aligner/log_out

  star_progress_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR progress log"
    doc: "STAR Log.progress.out"
    outputSource: star_aligner/log_progress

  star_stdout_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR stdout log"
    doc: "STAR Log.std.out"
    outputSource: star_aligner/log_std

  star_sj_log:
    type: File?
    format: "http://edamontology.org/format_2330"
    label: "STAR sj log"
    doc: "STAR SJ.out.tab"
    outputSource: star_aligner/log_sj

  fastx_statistics:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "FASTQ statistics"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats/statistics_file

  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index/bam_bai_pair

  bowtie_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Bowtie alignment log"
    doc: "Bowtie alignment log file"
    outputSource: bowtie_aligner/log_file

  rpkm_isoforms:
    type: File
    format: "http://edamontology.org/format_3752"
    label: "RPKM, grouped by isoforms"
    doc: "Calculated rpkm values, grouped by isoforms"
    outputSource: rpkm_calculation/isoforms_file

  get_stat_log:
    type: File?
    label: "Bowtie, STAR and GEEP combined log"
    format: "http://edamontology.org/format_2330"
    doc: "Processed and combined Bowtie & STAR aligner and GEEP logs"
    outputSource: get_stat/output_file

  trim_report:
    type: File
    label: "TrimGalore report"
    doc: "TrimGalore generated log"
    outputSource: trim_fastq/report_file

steps:

  extract_fastq:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file
    out: [fastq_file]

  trim_fastq:
    run: ../tools/trimgalore.cwl
    in:
      input_file: extract_fastq/fastq_file
      dont_gzip:
        default: true
      length:
        default: 30
    out:
      - trimmed_file
      - report_file

  rename:
    run: ../tools/rename.cwl
    in:
      source_file: trim_fastq/trimmed_file
      target_filename:
        source: extract_fastq/fastq_file
        valueFrom: $(self.basename)
    out:
      - target_file

  star_aligner:
    run: ../tools/star-alignreads.cwl
    in:
      readFilesIn: rename/target_file
      genomeDir: star_indices_folder
      outFilterMultimapNmax:
        default: 1
      outFilterMismatchNmax:
        default: 5
      alignSJDBoverhangMin:
        default: 1
      seedSearchStartLmax:
        default: 15
      clip3pNbases: clip_3p_end
      clip5pNbases: clip_5p_end
      threads: threads
    out:
      - aligned_file
      - log_final
      - uniquely_mapped_reads_number
      - log_out
      - log_progress
      - log_std
      - log_sj

  fastx_quality_stats:
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename/target_file
    out: [statistics_file]

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: star_aligner/aligned_file
      sort_output_filename:
        source: rename/target_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'.bam')
      threads: threads
    out: [bam_bai_pair]

  bam_to_bigwig:
    run: bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
      chrom_length_file: chrom_length_file
      mapped_reads_number: star_aligner/uniquely_mapped_reads_number
#     fragmentsize is not set (STAR gives only read length). It will be calculated automatically by bedtools genomecov.
    out: [bigwig_file]

  bowtie_aligner:
    run: ../tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: rename/target_file
      indices_folder: bowtie_indices_folder
      clip_3p_end: clip_3p_end
      clip_5p_end: clip_5p_end
      v:
        default: 3
      m:
        default: 1
      best:
        default: true
      strata:
        default: true
      sam:
        default: true
      threads: threads
    out: [log_file]

  rpkm_calculation:
    run: ../tools/geep.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
      annotation_file: annotation_file
      rpkm_threshold:
        default: 0.001
      exclude_chr: exclude_chr
      threads: threads
    out: [isoforms_file]

  get_stat:
      run: ../tools/python-get-stat-rnaseq.cwl
      in:
        star_log: star_aligner/log_final
        bowtie_log: bowtie_aligner/log_file
        rpkm_isoforms: rpkm_calculation/isoforms_file
      out: [output_file]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "trim-rnaseq-se"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/workflows/trim-rnaseq-se.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
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
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898
      - class: s:Person
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0001-9102-5681

doc: |
  Runs RNA-Seq BioWardrobe basic analysis with single-end data file.

s:about: |
  The original [BioWardrobe's](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465)
  **RNA-Seq** basic analysis for a **single-end** experiment.
  A corresponded input [FASTQ](http://maq.sourceforge.net/fastq.shtml) file has to be provided.

  Current workflow should be used only with the single-end RNA-Seq data. It performs the following steps:
  1. Trim adapters from input FASTQ file
  2. Use STAR to align reads from input FASTQ file according to the predefined reference indices; generate unsorted BAM file and alignment statistics file
  3. Use fastx_quality_stats to analyze input FASTQ file and generate quality statistics file
  4. Use samtools sort to generate coordinate sorted BAM(+BAI) file pair from the unsorted BAM file obtained on the step 1 (after running STAR)
  5. Generate BigWig file on the base of sorted BAM file
  6. Map input FASTQ file to predefined rRNA reference indices using Bowtie to define the level of rRNA contamination; export resulted statistics to file
  7. Calculate isoform expression level for the sorted BAM file and GTF/TAB annotation file using GEEP reads-counting utility; export results to file










