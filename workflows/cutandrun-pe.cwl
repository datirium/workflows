cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_root = function(basename) {
          return basename.split('.').slice(0,1).join('.');
      };
  

'sd:metadata':
  - "../metadata/chipseq-header.cwl"

'sd:upstream':
  genome_indices: "genome-indices.cwl"
  genome_indices_spikein: "genome-indices.cwl"


inputs:

  indices_folder:
    type: Directory
    'sd:upstreamSource': "genome_indices/bowtie_indices"
    label: "Sample genome index for alignment and peak calling:"
    'sd:localLabel': true
    doc: "Sample genome index for alignment and peak calling"

  indices_folder_for_spikein:
    type: Directory
    'sd:upstreamSource': "genome_indices_spikein/bowtie_indices"
    label: "Spike-in genome index for normalization:"
    'sd:localLabel': true
    doc: "Spike-in genome index for normalization"

  annotation_file:
    type: File
    'sd:upstreamSource': "genome_indices/annotation"
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: "Tab-separated annotation file"

  genome_size:
    type: string
    'sd:upstreamSource': "genome_indices/genome_size"
    label: "Effective genome size"
    doc: "effective genome size: hs, mm, ce, dm or number, for example 2.7e9"

  chrom_length:
    type: File
    'sd:upstreamSource': "genome_indices/chrom_length"
    label: "Chromosomes length file"
    format: "http://edamontology.org/format_2330"
    doc: "Chromosomes length file"

  fastq_file_upstream:
    type:
      - File
      - type: array
        items: File
    label: "FASTQ file for R1:"
    'sd:localLabel': true
    format: "http://edamontology.org/format_1930"
    doc: "Read1 data in a FASTQ format, received after paired end sequencing"

  fastq_file_downstream:
    type:
      - File
      - type: array
        items: File
    label: "FASTQ file for R2:"
    'sd:localLabel': true
    format: "http://edamontology.org/format_1930"
    doc: "Read2 data in a FASTQ format, received after paired end sequencing"

  clip_3p_end:
    type: int?
    default: 0
    'sd:layout':
      advanced: true
    label: "Number of bases to clip from the 3p end:"
    'sd:localLabel': true
    doc: "Number of bases to clip from the 3p end"

  clip_5p_end:
    type: int?
    default: 0
    'sd:layout':
      advanced: true
    label: "Number of bases to clip from the 3p end:"
    'sd:localLabel': true
    doc: "Number of bases to clip from the 5p end"

  remove_duplicates:
    type: boolean?
    default: false
    'sd:layout':
      advanced: true
    label: "Call samtools rmdup to remove duplicates from sorted BAM file?"
    'sd:localLabel': true
    doc: "Calls samtools rmdup to remove duplicates from sorted BAM file"

  promoter_dist:
    type: int?
    default: 1000
    'sd:layout':
      advanced: true
    label: "Max distance (bp) from gene TSS (in both directions) overlapping which the peak will be assigned to the promoter region:"
    'sd:localLabel': true
    doc: "Max distance (bp) from gene TSS (in both directions) overlapping which the peak will be assigned to the promoter region:"

  upstream_dist:
    type: int?
    default: 20000
    'sd:layout':
      advanced: true
    label: "Max distance (bp) from the promoter (only in upstream directions) overlapping which the peak will be assigned to the upstream region:"
    'sd:localLabel': true
    doc: "Max distance (bp) from the promoter (only in upstream directions) overlapping which the peak will be assigned to the upstream region:"    

  threads:
    type: int?
    default: 2
    'sd:layout':
      advanced: true
    label: "Number of threads for steps that support multithreading:"
    'sd:localLabel': true
    doc: "Number of threads for steps that support multithreading"


outputs:

  unaligned_fastq:
    type:
      - "null"
      - File[]
    format: "http://edamontology.org/format_1930"
    label: "Unaligned FASTQ file(s)"
    doc: "Unaligned FASTQ file(s)"
    outputSource: bowtie_aligner/unaligned_fastq

  multimapped_fastq:
    type:
      - "null"
      - File[]
    format: "http://edamontology.org/format_1930"
    label: "Multimapped FASTQ file(s)"
    doc: "Multimapped FASTQ file(s)"
    outputSource: bowtie_aligner/multimapped_fastq

  fastx_statistics_upstream:
    type: File
    label: "FASTQ 1 statistics"
    format: "http://edamontology.org/format_2330"
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
    label: "FASTQ 2 statistics"
    format: "http://edamontology.org/format_2330"
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

  bowtie_log:
    type: File
    label: "BOWTIE alignment log"
    format: "http://edamontology.org/format_2330"
    doc: "BOWTIE generated alignment log"
    outputSource: bowtie_aligner/log_file

  samtools_rmdup_log:
    type: File
    label: "Remove duplicates log"
    format: "http://edamontology.org/format_2330"
    doc: "Samtools rmdup generated log"
    outputSource: samtools_rmdup/rmdup_log

  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index_after_rmdup/bam_bai_pair
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        optional: true
        type: 'alignment'
        format: 'bam'
        name: "BAM Track"
        displayMode: "SQUISHED"

  bam_statistics_report:
    type: File
    label: "BAM statistics report (original)"
    format: "http://edamontology.org/format_2330"
    doc: "BAM statistics report (right after alignment and sorting)"
    outputSource: get_bam_statistics/log_file

  bam_statistics_spikein_report:
    type: File
    label: "BAM statistics report for spike in alignment (original)"
    format: "http://edamontology.org/format_2330"
    doc: "BAM statistics report for spike in alignment (right after alignment and sorting)"
    outputSource: get_spikein_bam_statistics/log_file

  bam_statistics_report_after_filtering:
    type: File
    label: "BAM statistics report (after filtering)"
    format: "http://edamontology.org/format_2330"
    doc: "BAM statistics report (after all filters applied)"
    outputSource: get_bam_statistics_after_filtering/log_file

  insert_size_report_after_filtering:
    type: File
    label: "Insert size distribution report (after filtering)"
    format: "http://edamontology.org/format_3475"
    doc: "Insert size distribution report (after all filters applied)"
    outputSource: get_bam_statistics_after_filtering/ext_is_section
    'sd:visualPlugins':
    - scatter:
        tab: 'QC Plots'
        Title: 'Insert Size Distribution (after filtering)'
        xAxisTitle: 'Insert size'
        yAxisTitle: 'Pairs total'
        colors: ["#4b78a3"]
        height: 500
        data: [$1, $2]
        comparable: "isdp"

  trim_report_upstream:
    type: File
    label: "TrimGalore report Upstream"
    doc: "TrimGalore generated log for FASTQ 1"
    outputSource: trim_fastq/report_file

  trim_report_downstream:
    type: File
    label: "TrimGalore report Downstream"
    doc: "TrimGalore generated log for FASTQ 2"
    outputSource: trim_fastq/report_file_pair

  preseq_estimates:
    type: File?
    label: "Preseq estimates"
    format: "http://edamontology.org/format_3475"
    doc: "Preseq estimated results"
    outputSource: preseq/estimates_file
    'sd:visualPlugins':
    - scatter:
        tab: 'QC Plots'
        Title: 'Preseq Estimates'
        xAxisTitle: 'Total reads count'
        yAxisTitle: 'Expected distinct reads count'
        colors: ["#4b78a3"]
        height: 500
        data: [$1, $2]
        comparable: "preseq"

  bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated BigWig file"
    outputSource: bam_to_bigwig/bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "BigWig Track"
        height: 120

  bigwig_scaled:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "scaled BigWig file"
    doc: "Generated scaled BigWig file"
    outputSource: bam_to_bigwig_scaled/bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "BigWig Track Scaled"
        height: 120

  fc_sorted_bed:
    type: File
    format: "http://edamontology.org/format_3583"
    label: "bedgraph file"
    doc: "Generated bedgraph file"
    outputSource: fragment_counts/sorted_bed 

  fc_sorted_bed_scaled:
    type: File
    format: "http://edamontology.org/format_3583"
    label: "bedgraph file"
    doc: "Generated bedgraph file"
    outputSource: fragment_counts/sorted_bed_scaled

  fc_log_file_stderr:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr logfile"
    outputSource: fragment_counts/log_file_stderr

  fc_log_file_stdout:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout logfile"
    outputSource: fragment_counts/log_file_stdout     

  norm_peak_log_stderr:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "seacr logfile"
    doc: "stderr from seacr command"
    outputSource: seacr_callpeak/log_file_stderr

  norm_peak_log_stdout:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "seacr logfile stdout"
    doc: "stdout from seacr command"
    outputSource: seacr_callpeak/log_file_stdout

  stats_for_vis_md:
    type: File?
    label: "Markdown formatted combined log"
    format: "http://edamontology.org/format_3835"
    doc: "Markdown formatted combined log"
    outputSource: stats_for_vis/modified_file_md
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  stats_for_vis_tsv:
    type: File?
    label: "Bowtie & Samtools Rmdup combined formatted log"
    format: "http://edamontology.org/format_3475"
    doc: "Processed and combined Bowtie aligner and Samtools rmdup formatted log"
    outputSource: stats_for_vis/modified_file_tsv
    'sd:visualPlugins':
    - tableView:
        vertical: true
        tab: 'Overview'
    'sd:preview':
      'sd:visualPlugins':
      - pie:
          colors: ['#b3de69', '#99c0db', '#fb8072', '#fdc381']
          data: [$2, $3, $4, $5]

  stats_for_vis_yaml:
    type: File?
    label: "YAML formatted combined log"
    format: "http://edamontology.org/format_3750"
    doc: "YAML formatted combined log"
    outputSource: stats_for_vis/modified_file_yaml

  stats_for_vis_log_stdout:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stats logfile stdout"
    doc: "stdout from stats_for_vis step"
    outputSource: stats_for_vis/log_file_stdout

  stats_for_vis_log_stderr:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stats logfile stderr"
    doc: "stderr from stats_for_vis step"
    outputSource: stats_for_vis/log_file_stderr

  macs2_called_peaks:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "bedgraph file of peaks from seacr"
    doc: "Bed file of enriched regions called by seacr (from normalized bigwig) in macs2's bed format. *peaks removed where col2>col3"
    outputSource: convert_bed_to_xls/output_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'bed'
        name: "Called Peaks"
        height: 120

  annotated_peaks:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "gene annotated peaks file"
    doc: "nearest gene annotation per peak"
    outputSource: island_intersect/result_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Annotated Peaks'
        Title: 'sparse enriched peak list with nearest gene annotation'


steps:

  extract_fastq_upstream:
    label: "Loading unmapped sequence data for read 1"
    doc: |
      Most DNA cores and commercial NGS companies return unmapped sequence data in FASTQ format.
      The data can be uploaded from users computer, downloaded directly from an ftp server of
      the core facility by providing a URL or from GEO by providing SRA accession number.
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_upstream
      output_prefix:
        default: "merged_R1"
    out: [fastq_file]

  extract_fastq_downstream:
    label: "Loading unmapped sequence data for read 2"
    doc: |
      Most DNA cores and commercial NGS companies return unmapped sequence data in FASTQ format.
      The data can be uploaded from users computer, downloaded directly from an ftp server of
      the core facility by providing a URL or from GEO by providing SRA accession number.
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file_downstream
      output_prefix:
        default: "merged_R2"
    out: [fastq_file]

  trim_fastq:
    label: "Adapter trimming"
    doc: |
      For libraries sequenced on the Illumina platform it’s recommended to remove adapter sequences
      from the reads. If adapters are not trimmed there is a high risk of reads being unmapped to a
      reference genome. This becomes particularly important when the reads are long and the fragments
      are short - resulting in sequencing adapters at the end of read. If adapter trimming will cause
      all the reads become too short (<30bp), this step will be skipped.
    run: ../tools/trimgalore.cwl
    in:
      input_file: extract_fastq_upstream/fastq_file
      input_file_pair: extract_fastq_downstream/fastq_file
      dont_gzip:
        default: true
      length:
        default: 30
      trim1:
        default: true
      paired:
        default: true
    out:
      - trimmed_file
      - trimmed_file_pair
      - report_file
      - report_file_pair

  bypass_trim:
    run: ../tools/bypass-trimgalore-pe.cwl
    in:
      original_fastq_file_1: extract_fastq_upstream/fastq_file
      trimmed_fastq_file_1: trim_fastq/trimmed_file
      trimming_report_file_1: trim_fastq/report_file
      original_fastq_file_2: extract_fastq_downstream/fastq_file
      trimmed_fastq_file_2: trim_fastq/trimmed_file_pair
      trimming_report_file_2: trim_fastq/report_file_pair
      min_reads_count:
        default: 100  # any small number should be good, as we are catching the case when trimgalore discarded all reads
    out:
      - selected_fastq_file_1
      - selected_report_file_1
      - selected_fastq_file_2
      - selected_report_file_2

  rename_upstream:
    run: ../tools/rename.cwl
    in:
      source_file: bypass_trim/selected_fastq_file_1
      target_filename:
        source: extract_fastq_upstream/fastq_file
        valueFrom: $(self.basename)
    out:
      - target_file

  rename_downstream:
    run: ../tools/rename.cwl
    in:
      source_file: bypass_trim/selected_fastq_file_2
      target_filename:
        source: extract_fastq_downstream/fastq_file
        valueFrom: $(self.basename)
    out:
      - target_file

  fastx_quality_stats_upstream:
    label: "Quality control of unmapped sequence data for read 1"
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename_upstream/target_file
    out: [statistics_file]

  fastx_quality_stats_downstream:
    label: "Quality control of unmapped sequence data for read 2"
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename_downstream/target_file
    out: [statistics_file]

  bowtie_aligner:
    label: "Alignment to reference genome"
    doc: |
      Aligns reads to the reference genome keeping only uniquely mapped reads with
      less than 3 mismatches.
    run: ../tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: rename_upstream/target_file
      downstream_filelist: rename_downstream/target_file
      indices_folder: indices_folder
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
      unaligned_prefix:
        default: "unaligned_reads"
      multimapped_prefix:
        default: "multimapped_reads"
      threads: threads
      q:
        default: true
      X:
        default: 500
    out: [sam_file, log_file, unaligned_fastq, multimapped_fastq]

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: bowtie_aligner/sam_file
      threads: threads
    out: [bam_bai_pair]

  clean_sam_headers_for_preseq:
    run: ../tools/samtools-clean-headers.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
    out: [preseq_bam]

  preseq:
    label: "Sequencing depth estimation"
    doc: |
      Estimates the complexity of the sequencing library, evaluates how many reads can
      be expected from the additional sequencing of the same experiment.
    run: ../tools/preseq-lc-extrap.cwl
    in:
      bam_file: clean_sam_headers_for_preseq/preseq_bam
      pe_mode:
        default: true
      extrapolation:
        default: 1000000000
    out: [estimates_file]

  samtools_rmdup:
    label: "PCR duplicates removal"
    doc: |
      Removes potential PCR duplicates. This step is used to remove reads overamplified
      in PCR. Unfortunately, it may also remove "good" reads. We do not recommend to
      remove duplicates unless the library is heavily duplicated.
    run: ../tools/samtools-rmdup.cwl
    in:
      trigger: remove_duplicates
      bam_file: samtools_sort_index/bam_bai_pair
    out: [rmdup_output, rmdup_log]

  samtools_sort_index_after_rmdup:
    run: ../tools/samtools-sort-index.cwl
    in:
      trigger: remove_duplicates
      sort_input: samtools_rmdup/rmdup_output
      threads: threads
    out: [bam_bai_pair]

  bowtie_aligner_spikein:
    label: "Alignment to spike-in reference genome"
    doc: |
      Aligns reads to the spike-in reference genome keeping only uniquely mapped reads with
      less than 3 mismatches.
    run: ../tools/bowtie-alignreads.cwl
    in:
      upstream_filelist:
          source: bowtie_aligner/unaligned_fastq
          valueFrom: $(self[0])
      downstream_filelist:
          source: bowtie_aligner/unaligned_fastq
          valueFrom: $(self[1])
      indices_folder: indices_folder_for_spikein
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
      unaligned_prefix:
        default: "unaligned_reads_spikein"
      multimapped_prefix:
        default: "multimapped_reads_spikein"
      threads: threads
      q:
        default: true
      X:
        default: 500
    out: [sam_file, log_file, unaligned_fastq, multimapped_fastq]

  samtools_sort_index_spikein:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: bowtie_aligner_spikein/sam_file
      threads: threads
    out: [bam_bai_pair]

  get_spikein_bam_statistics:
    label: "Quality control of aligned sequence data to spike-in reference (default is E. coli)"
    doc: |
      Calculates alignment statistics, such as reads mapped/unmapped, average
      read length and quality score, etc.
    run: ../tools/samtools-stats.cwl
    in:
      bambai_pair: samtools_sort_index_spikein/bam_bai_pair
      output_filename:
        source: samtools_sort_index_spikein/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_bam_statistics_spikein_report.txt")
    out: [log_file, reads_mapped]

  get_scale_from_spikein:
    label: "return scaling factor from spike-in aligned read count (x)"
    doc: |
      scaling factor = (10000/x)
    run: ../tools/get_scale.cwl
    in:
      reads_mapped: get_spikein_bam_statistics/reads_mapped
    out: [scaling_factor]

  get_bam_statistics:
    label: "Quality control of aligned sequence data"
    doc: |
      Calculates alignment statistics, such as reads mapped/unmapped, average
      read length and quality score, etc.
    run: ../tools/samtools-stats.cwl
    in:
      bambai_pair: samtools_sort_index/bam_bai_pair
      output_filename:
        source: samtools_sort_index/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_bam_statistics_report.txt")
    out: [log_file]

  get_bam_statistics_after_filtering:
    run: ../tools/samtools-stats.cwl
    in:
      bambai_pair: samtools_sort_index_after_rmdup/bam_bai_pair
      output_filename:
        source: samtools_sort_index_after_rmdup/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_bam_statistics_report_after_filtering.txt")
    out: [log_file, ext_is_section]

  bam_to_bigwig:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index_after_rmdup/bam_bai_pair
      chrom_length_file: chrom_length
      pairchip:
        default: true
    out: [bigwig_file]

  bam_to_bigwig_scaled:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index_after_rmdup/bam_bai_pair
      chrom_length_file: chrom_length
      scale: get_scale_from_spikein/scaling_factor
      pairchip:
        default: true
      bigwig_filename:
        source: samtools_sort_index_after_rmdup/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_scaled.bigWig")
    out: [bigwig_file]

  fragment_counts:
    run: ../tools/bedtools-fragmentcounts.cwl
    in:
      bam_file: samtools_sort_index_after_rmdup/bam_bai_pair
      scale: get_scale_from_spikein/scaling_factor
      output_prefix:
        source: samtools_sort_index_after_rmdup/bam_bai_pair
        valueFrom: $(get_root(self.basename))
    out: [sorted_bed, sorted_bed_scaled, log_file_stderr, log_file_stdout]

  seacr_callpeak:
    run: ../tools/seacr.cwl
    in:
      treatment_bedgraph: fragment_counts/sorted_bed_scaled
      numeric_threshold:
        default: 0.01
      norm_control_to_treatment:
        default: "non"
      peakcalling_mode:
        default: "stringent"
      output_prefix:
        source: fragment_counts/sorted_bed_scaled
        valueFrom: $(get_root(self.basename)+"_scaled")
    out: [peak_tsv_file, log_file_stderr, log_file_stdout]
    doc: |
      Input is normalized depth data using spike-in mapped reads (E. coli by default).
      Scaling factor (sf) for seq library normalization:
            sf=(C/[mapped reads]) where C is a constant (10000 used here)
      Henikoff protocol, Section 16: https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=16#step-4A3D8C70DC3011EABA5FF3676F0827C5)

  get_stat:
    run: ../tools/collect-statistics-cutandrun.cwl
    in:
      trimgalore_report_fastq_1: bypass_trim/selected_report_file_1
      trimgalore_report_fastq_2: bypass_trim/selected_report_file_2
      bowtie_alignment_report: bowtie_aligner/log_file
      bam_statistics_report: get_bam_statistics/log_file
      bam_statistics_after_filtering_report: get_bam_statistics_after_filtering/log_file
      seacr_called_peaks: seacr_callpeak/peak_tsv_file
      paired_end:
        default: True
      output_prefix:
        source: seacr_callpeak/peak_tsv_file
        valueFrom: $(get_root(self.basename))
    out: [collected_statistics_yaml, collected_statistics_tsv, collected_statistics_md, mapped_reads]
    doc: |
      Statistics pulled from alignment reports and other logs, as well as spike-in normalized peak calling.

  stats_for_vis:
    run: ../tools/collect-statistics-frip.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
      seacr_called_peaks_norm: seacr_callpeak/peak_tsv_file
      collected_statistics_md: get_stat/collected_statistics_md
      collected_statistics_tsv: get_stat/collected_statistics_tsv
      collected_statistics_yaml: get_stat/collected_statistics_yaml
      spikein_reads_mapped: get_spikein_bam_statistics/reads_mapped
    out: [modified_file_md, modified_file_tsv, modified_file_yaml, log_file_stdout, log_file_stderr]

  convert_bed_to_xls:
    run: ../tools/custom-bash.cwl
    in:
      input_file: seacr_callpeak/peak_tsv_file
      script:
        default: >
          cat $0 | awk -F'\t'
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {if($3>$2){print $1"\t"$2"\t"$3"\t"$3-$2+1"\t"$5"\t"$4"\t0\t0\t0\tpeak_"NR}}' > `basename $0`
    out:
    - output_file
    doc: |
      This step also removes rows where start (col2) > end (col3). Still investigating why SEACR produces peak coordinates in this way.

  island_intersect:
    label: "Peak annotation"
    doc: |
      Assigns nearest genes to peaks to explore the biological implication of the open
      chromatin binding sites.
    run: ../tools/iaintersect.cwl
    in:
      input_filename: convert_bed_to_xls/output_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
    out: [result_file, log_file]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "cutandrun paired-end workflow"
label: "cutandrun paired-end workflow"
s:alternateName: "Cut & Run basic analysis workflow for a paired-end experiment with Trim Galore"

s:downloadUrl: https://github.com/datirium/workflows/tree/master/workflows/workflows/cutandrun-pe.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:robert.player@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X


doc: |
  A basic analysis workflow for a **paired-read** Cut & Run/TAG experiment.
    
  This workflow utilizes the tool SEACR (Sparse Enrichment Analysis of Cut & Run data) which calls enriched
  regions in target data by selecting the top 1% of regions by AUC. This workflow is loosely based on the
  [CUT-RUNTools-2.0 pipeline](https://github.com/fl-yu/CUT-RUNTools-2.0) pipeline, and the ChIP-Seq pipeline
  from [BioWardrobe](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465)
  was used as a CWL template.

  User-provided inputs required to run this workflow include the [FASTQ](http://maq.sourceforge.net/fastq.shtml)
  input files for R1 (upstream) and R2 (downstream) reads, selection from the dropdown of a reference genome the
  sample library was prepared from, and selection from the dropdown of a spike-in genome (default should be E.
  coli K-12).

  Intermediate and final downloadable outputs include coordinate sorted BAM file with associated BAI file,
  quality statistics for both R1/R2 input FASTQ files, read pileup/coverage in BigWig format, and SEACR-called
  peak BED files (raw and spike-in normalized).


  ### Data Analysis  
  This pipeline includes the following steps:
  1. Trimming the adapters with TrimGalore. This step is particularly important
          when the reads are long and the fragments are short-resulting in
          sequencing adapters at the end of read. If adapter is not removed the
          read will not map. TrimGalore can recognize standard adapters, such as
          Illumina or Nextera/Tn5 adapters.
  2. Quality control of trimmed, unmapped sequence data
  3. (Optional) trimming adapters on 5' or 3' end by the specified number of
          bases.
  4. Mapping reads with BowTie. Only uniquely mapped reads with less than 3
          mismatches are used in the downstream analysis. Results are saved as a
          .sam file, which is then sorted and indexed. Final outputs are in
          bam/bai format, which are also used to extrapolate effects of
          additional sequencing based on library complexity.
  5.  (Optional) Removal of duplicates (reads/pairs of reads mapping to exactly
          same location). This step is used to remove reads overamplified in PCR.
          Unfortunately, it may also remove "good" reads. We usually do not
          remove duplicates unless the library is heavily duplicated.
  6.  Call enriched regions using SEACR. Note, there is no native way of
          specifying a custom calculated read extension length for SEACR as there
          is for the ChIP-Seq peak calling tool MACS2. The only output of SEACR
          are the called peaks in a single BED format file.
  7.  Generation of read, alignment, and peak statistics, as well as BigWig
          coverage files for display on the browser. The coverage shows the number
          of fragments at each base in the genome normalized to the number of
          millions of mapped reads.


  ### Additional Processing Details & Notes
  _Trim Galore_ is a wrapper around [Cutadapt](https://github.com/marcelm/cutadapt) and
  [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to consistently apply adapter and quality
  trimming to FASTQ files, with extra functionality for RRBS data. This workflow starts with running
  fastx_quality_stats (step: fastx_quality_stats_upstream, fastx_quality_stats_downstream) from FASTX-Toolkit
  to calculate quality statistics for both R1 and R2 input FASTQ files. Bowtie is run in parallel to align reads
  from input FASTQ files to the selected reference genome (step: bowtie_aligner). The output of this step is an
  unsorted SAM file that is then sorted and indexed by samtools sort and index (step: samtools_sort_index).
  Depending on workflow input parameters, indexed and sorted BAM files could be processed by samtools rmdup
  (step: samtools_rmdup) to remove all detected instances of read duplication. In the case when removing
  duplicates is not necessary, this step returns the original input BAM and BAI files without any processing. If
  the duplicates were removed, the following step (step: samtools_sort_index_after_rmdup) reruns samtools sort
  and index with the BAM and BAI files, if not, the step returns the original unchanged input files.
