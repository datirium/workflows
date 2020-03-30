cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_root = function(basename) {
          return basename.split('.').slice(0,1).join('.');
      };

'sd:metadata':
- "../metadata/rnaseq-header.cwl"


'sd:upstream':
  genome_indices: "genome-indices.cwl"

inputs:

# General inputs

  star_indices_folder:
    type: Directory
    label: "STAR indices folder"
    'sd:upstreamSource': "genome_indices/star_indices"
    doc: "Path to STAR generated indices"

  star_indices_folder_mitochondrial:
    type: Directory
    label: "STAR indices mitochondrial folder"
    'sd:upstreamSource': "genome_indices/mitochondrial_indices"
    doc: "Path to STAR generated indices for mitochondrial dna"

  bowtie_indices_folder:
    type: Directory
    label: "BowTie Ribosomal Indices"
    'sd:upstreamSource': "genome_indices/ribosomal_indices"
    doc: "Path to Bowtie generated indices"

  chrom_length_file:
    type: File
    label: "Chromosome length file"
    format: "http://edamontology.org/format_2330"
    'sd:upstreamSource': "genome_indices/chrom_length"
    doc: "Chromosome length file"

  annotation_file:
    type: File
    label: "Annotation file"
    format:
      - "http://edamontology.org/format_2306"
      - "http://edamontology.org/format_3475"
    'sd:upstreamSource': "genome_indices/annotation"
    doc: "GTF or TAB-separated annotation file"

  fastq_file_upstream:
    type: File
    label: "FASTQ 1 input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, received after paired end sequencing"

  fastq_file_downstream:
    type: File
    label: "FASTQ 2 input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, received after paired end sequencing"

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
    doc: "Generated BigWig file"
    outputSource: bam_to_bigwig_upstream/bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "BigWig 1"
        height: 120

  bigwig_downstream:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated BigWig file"
    outputSource: bam_to_bigwig_downstream/bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "BigWig 2"
        height: 120

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

  fastx_statistics_upstream:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "FASTQ 1 statistics"
    doc: "fastx_quality_stats generated FASTQ 1 quality statistics file"
    outputSource: fastx_quality_stats_upstream/statistics_file
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'FASTQ 1 Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'FASTQ 1 Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  fastx_statistics_downstream:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "FASTQ 2 statistics"
    doc: "fastx_quality_stats generated FASTQ 2 quality statistics file"
    outputSource: fastx_quality_stats_downstream/statistics_file
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'FASTQ 2 Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'FASTQ 2 Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  bam_merged_index:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: merge_original_and_mitochondrial_index/bam_bai_pair
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        optional: true
        type: 'alignment'
        format: 'bam'
        name: "BAM Track"
        displayMode: "SQUISHED"

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

  rpkm_genes:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "RPKM, grouped by gene name"
    doc: "Calculated rpkm values, grouped by gene name"
    outputSource: group_isoforms/genes_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Gene Expression'
        Title: 'RPKM, grouped by gene name'

  rpkm_common_tss:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "RPKM, grouped by common TSS"
    doc: "Calculated rpkm values, grouped by common TSS"
    outputSource: group_isoforms/common_tss_file

  get_stat_log:
    type: File?
    label: "YAML formatted combined log"
    format: "http://edamontology.org/format_3750"
    doc: "YAML formatted combined log"
    outputSource: get_stat/collected_statistics_yaml

  get_stat_markdown:
    type: File?
    label: "Markdown formatted combined log"
    format: "http://edamontology.org/format_3835"
    doc: "Markdown formatted combined log"
    outputSource: get_stat/collected_statistics_md
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  get_formatted_stats:
    type: File?
    label: "Bowtie, STAR and GEEP mapping stats"
    format: "http://edamontology.org/format_2330"
    doc: "Processed and combined Bowtie & STAR aligner and GEEP logs"
    outputSource: get_stat/collected_statistics_tsv
    'sd:visualPlugins':
    - tableView:
        vertical: true
        tab: 'Overview'
    'sd:preview':
      'sd:visualPlugins':
      - pie:
          colors: ['#b3de69', '#99c0db', '#fdc381', '#fb8072']
          data: [$2, $3, $4, $5]

  bam_statistics_report:
    type: File
    label: "BAM statistics report"
    format: "http://edamontology.org/format_2330"
    doc: "BAM statistics report (right after alignment and sorting)"
    outputSource: get_bam_statistics/log_file

  insert_size_report:
    type: File
    label: "Insert size distribution report"
    format: "http://edamontology.org/format_3475"
    doc: "Insert size distribution report (right after alignment and sorting)"
    outputSource: get_bam_statistics/ext_is_section
    'sd:visualPlugins':
    - scatter:
        tab: 'QC Plots'
        Title: 'Insert Size Distribution'
        xAxisTitle: 'Insert size'
        yAxisTitle: 'Pairs total'
        colors: ["#4b78a3"]
        height: 500
        data: [$1, $2]


steps:

  extract_fastq_upstream:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_upstream
    out: [fastq_file]

  extract_fastq_downstream:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_downstream
    out: [fastq_file]

  fastx_quality_stats_upstream:
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: extract_fastq_upstream/fastq_file
    out: [statistics_file]

  fastx_quality_stats_downstream:
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: extract_fastq_downstream/fastq_file
    out: [statistics_file]

  star_aligner:
    run: ../tools/star-alignreads.cwl
    in:
      readFilesIn: [extract_fastq_upstream/fastq_file, extract_fastq_downstream/fastq_file]
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
      - log_final
      - unmapped_mate_1_file
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

  samtools_sort_index_mitochondrial:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: star_aligner_mitochondrial/aligned_file
      sort_output_filename:
        source: extract_fastq_upstream/fastq_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'_mitochondrial.bam')
      threads: threads
    out: [bam_bai_pair]

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: star_aligner/aligned_file
      sort_output_filename:
        source: extract_fastq_upstream/fastq_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'_sorted.bam')
      threads: threads
    out: [bam_bai_pair]

  merge_original_and_mitochondrial:
    run: ../tools/samtools-merge.cwl
    in:
      output_filename:
        source: extract_fastq_upstream/fastq_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'_merged.bam')
      alignment_files: [ samtools_sort_index/bam_bai_pair, samtools_sort_index_mitochondrial/bam_bai_pair ]
    out: [merged_alignment_file]

  merge_original_and_mitochondrial_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: merge_original_and_mitochondrial/merged_alignment_file
      sort_output_filename:
        source: extract_fastq_upstream/fastq_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'.bam')
      threads: threads
    out: [bam_bai_pair]

  bam_to_bigwig_upstream:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: merge_original_and_mitochondrial_index/bam_bai_pair
      chrom_length_file: chrom_length_file
      mapped_reads_number: star_aligner/uniquely_mapped_reads_number
      bigwig_filename:
        source: extract_fastq_upstream/fastq_file
        valueFrom: |
          ${
            var root = self.basename.split('.').slice(0,-1).join('.');
            var ext = "_upstream.bigWig";
            return (root == "")?self.basename+ext:root+ext;
          }
      strand:
        default: '+'
    out: [bigwig_file]

  bam_to_bigwig_downstream:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: merge_original_and_mitochondrial_index/bam_bai_pair
      chrom_length_file: chrom_length_file
      mapped_reads_number:
        source: star_aligner/uniquely_mapped_reads_number
        valueFrom: $(-self)
      bigwig_filename:
        source: extract_fastq_upstream/fastq_file
        valueFrom: |
          ${
            var root = self.basename.split('.').slice(0,-1).join('.');
            var ext = "_downstream.bigWig";
            return (root == "")?self.basename+ext:root+ext;
          }
      strand:
        default: '-'
    out: [bigwig_file]

  bowtie_aligner:
    run: ../tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: extract_fastq_upstream/fastq_file
      downstream_filelist: extract_fastq_downstream/fastq_file
      indices_folder: bowtie_indices_folder
      clip_3p_end: clip_3p_end
      clip_5p_end: clip_5p_end
      v:
        default: 3
      m:
        default: 1
      sam:
        default: true
      threads: threads
    out: [log_file]

  rpkm_calculation:
    run: ../tools/geep.cwl
    in:
      bam_file: merge_original_and_mitochondrial_index/bam_bai_pair
      annotation_file: annotation_file
      dutp:
        default: true
      rpkm_threshold:
        default: 0.001
      exclude_chr: exclude_chr
      threads: threads
    out: [isoforms_file]

  group_isoforms:
    run: ../tools/group-isoforms.cwl
    in:
      isoforms_file: rpkm_calculation/isoforms_file
    out:
      - genes_file
      - common_tss_file

  get_bam_statistics:
    run: ../tools/samtools-stats.cwl
    in:
      bambai_pair: samtools_sort_index/bam_bai_pair
      output_filename:
        source: samtools_sort_index/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_bam_statistics_report.txt")
    out: [log_file, ext_is_section]

  get_stat:
      run: ../tools/collect-statistics-rna-seq.cwl
      in:
        star_alignment_report: star_aligner/log_final
        bowtie_alignment_report: bowtie_aligner/log_file
        bam_statistics_report: get_bam_statistics/log_file
        isoforms_file: rpkm_calculation/isoforms_file
        paired_end:
          default: true
      out: [collected_statistics_yaml, collected_statistics_tsv, collected_statistics_md]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "RNA-Seq pipeline paired-end stranded mitochondrial"
label: "RNA-Seq pipeline paired-end stranded mitochondrial"
s:alternateName: "RNA-Seq strand specific mitochondrial workflow for pair-end experiment based on BioWardrobe's basic analysis"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/rnaseq-pe-dutp-mitochondrial.cwl
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


doc:
  $include: ../descriptions/rnaseq-pe-dutp-mitochondrial.md
