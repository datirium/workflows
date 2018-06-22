cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement

'sd:metadata':
- "https://raw.githubusercontent.com/datirium/workflows/master/metadata/rnaseq-header.cwl"

'sd:upstream':
  star_index: "https://raw.githubusercontent.com/datirium/workflows/master/workflows/star-index.cwl"
  star_index_mitochondrial: "https://raw.githubusercontent.com/datirium/workflows/master/workflows/star-index.cwl"
  bowtie_index: "https://raw.githubusercontent.com/datirium/workflows/master/workflows/bowtie-index.cwl"

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
    'sd:upstreamSource': "star_index/indices_folder"
    doc: "Path to STAR generated indices"

  star_indices_folder_mitochondrial:
    type: Directory
    label: "STAR indices mitochondrial folder"
    'sd:upstreamSource': "star_index_mitochondrial/indices_folder"
    doc: "Path to STAR generated indices for mitochondrial dna"

  bowtie_indices_folder:
    type: Directory
    label: "BowTie Ribosomal Indices"
    'sd:upstreamSource': "bowtie_index/indices_folder"
    doc: "Path to Bowtie generated indices"

  chrom_length_file:
    type: File
    label: "Chromosome length file"
    format: "http://edamontology.org/format_2330"
    'sd:upstreamSource': "star_index/chrom_length"
    doc: "Chromosome length file"

  annotation_file:
    type: File
    label: "Annotation file"
    format:
      - "http://edamontology.org/format_2306"
      - "http://edamontology.org/format_3475"
    'sd:upstreamSource': "star_index/annotation_file"
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

  bigwig_upstream:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated upstream BigWig file"
    outputSource: bam_to_bigwig_upstream/bigwig_file

  bigwig_downstream:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated downstream BigWig file"
    outputSource: bam_to_bigwig_downstream/bigwig_file

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

  bam_merged_index:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: merge_original_and_mitochondrial_index/bam_bai_pair

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

steps:

  extract_fastq:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file
    out: [fastq_file]

  star_aligner:
    run: ../tools/star-alignreads.cwl
    in:
      readFilesIn: extract_fastq/fastq_file
      genomeDir: star_indices_folder
      outFilterMultimapNmax:
        default: 1
      outFilterMismatchNmax:
        default: 5
      alignSJDBoverhangMin:
        default: 1
      seedSearchStartLmax:
        default: 15
      outReadsUnmapped:
        default: "Fastx"
      clip3pNbases: clip_3p_end
      clip5pNbases: clip_5p_end
      threads: threads
    out:
      - aligned_file
      - unmapped_mate_1_file
      - log_final
      - uniquely_mapped_reads_number
      - log_out
      - log_progress
      - log_std
      - log_sj

  star_aligner_mitochondrial:
    run: ../tools/star-alignreads.cwl
    in:
      readFilesIn: star_aligner/unmapped_mate_1_file
      genomeDir: star_indices_folder_mitochondrial
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
      input_file: extract_fastq/fastq_file
    out: [statistics_file]

  samtools_sort_index_mitochondrial:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: star_aligner_mitochondrial/aligned_file
      sort_output_filename:
        source: extract_fastq/fastq_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'_mitochondrial.bam')
      threads: threads
    out: [bam_bai_pair]

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: star_aligner/aligned_file
      sort_output_filename:
        source: extract_fastq/fastq_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'_sorted.bam')
      threads: threads
    out: [bam_bai_pair]

  merge_original_and_mitochondrial:
    run: ../tools/samtools-merge.cwl
    in:
      output_name:
        source: extract_fastq/fastq_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'_merged.bam')
      input: [ samtools_sort_index/bam_bai_pair, samtools_sort_index_mitochondrial/bam_bai_pair ]
    out: [output]

  merge_original_and_mitochondrial_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: merge_original_and_mitochondrial/output
      sort_output_filename:
        source: extract_fastq/fastq_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'.bam')
      threads: threads
    out: [bam_bai_pair]

  bam_to_bigwig_upstream:
    run: bam-bedgraph-bigwig.cwl
    in:
      bam_file: merge_original_and_mitochondrial/output
      chrom_length_file: chrom_length_file
      mapped_reads_number: star_aligner/uniquely_mapped_reads_number
      bigwig_filename:
        source: extract_fastq/fastq_file
        valueFrom: |
          ${
            let root = self.basename.split('.').slice(0,-1).join('.');
            let ext = "_upstream.bigWig";
            return (root == "")?self.basename+ext:root+ext;
          }
      strand:
        default: '+'
    out: [bigwig_file]

  bam_to_bigwig_downstream:
    run: bam-bedgraph-bigwig.cwl
    in:
      bam_file: merge_original_and_mitochondrial/output
      chrom_length_file: chrom_length_file
      mapped_reads_number:
        source: star_aligner/uniquely_mapped_reads_number
        valueFrom: $(-self)
      bigwig_filename:
        source: extract_fastq/fastq_file
        valueFrom: |
          ${
            let root = self.basename.split('.').slice(0,-1).join('.');
            let ext = "_downstream.bigWig";
            return (root == "")?self.basename+ext:root+ext;
          }
      strand:
        default: '-'
    out: [bigwig_file]

  bowtie_aligner:
    run: ../tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: extract_fastq/fastq_file
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
      bam_file: merge_original_and_mitochondrial/output
      annotation_file: annotation_file
      dutp:
        default: true
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

s:name: "rnaseq-se-dutp-mitochondrial"
s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/rnaseq-se-dutp-mitochondrial.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium, LLC"
  s:member:
  - class: s:Person
    s:name: Artem BArski
    s:email: mailto:Artem.Barski@datirum.com
  - class: s:Person
    s:name: Andrey Kartashov
    s:email: mailto:Andrey.Kartashov@datirium.com
    s:sameAs:
    - id: http://orcid.org/0000-0001-9102-5681

doc: |
  RNA-Seq strand specific mitochondrial workflow for single-read experiment based on BioWardrobe's basic analysis.

s:about: |
  Slightly changed original [BioWardrobe's](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465)
  **RNA-Seq** basic analysis for **strand specific single-read** experiment.
  An additional steps were added to map data to mitochondrial chromosome only and then merge the output.

  Experiment files in [FASTQ](http://maq.sourceforge.net/fastq.shtml) format either compressed or not can be used.

  Current workflow should be used only with single-read strand specific RNA-Seq data. It performs the following steps:
  1. `STAR` to align reads from input FASTQ file according to the predefined reference indices; generate unsorted BAM file and alignment statistics file
  2. `fastx_quality_stats` to analyze input FASTQ file and generate quality statistics file
  3. `samtools sort` to generate coordinate sorted BAM(+BAI) file pair from the unsorted BAM file obtained on the step 1 (after running STAR)
  5. Generate BigWig file on the base of sorted BAM file
  6. Map input FASTQ file to predefined rRNA reference indices using Bowtie to define the level of rRNA contamination; export resulted statistics to file
  7. Calculate isoform expression level for the sorted BAM file and GTF/TAB annotation file using `GEEP` reads-counting utility; export results to file









