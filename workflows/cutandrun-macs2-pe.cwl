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
  

'sd:upstream':
  genome_indices: "genome-indices.cwl"
  genome_indices_spikein: "genome-indices.cwl"
  control_file: "cutandrun-macs2-pe.cwl"


inputs:

  alias:
    type: string
    label: "Sample short name"
    sd:preview:
      position: 1

  cells:
    type: string
    label: "Cell type/name/id"
    sd:preview:
      position: 2

  conditions:
    type: string
    label: "Experimental condition(s)"
    sd:preview:
      position: 3

  catalog:
    type: string?
    label: "Catalog #"

  indices_folder:
    type: Directory
    'sd:upstreamSource': "genome_indices/bowtie_indices"
    label: "Primary genome index for peak calling:"
    'sd:localLabel': true
    doc: "Preprocessed genome index of sample organism for primary alignment and peak calling."

  indices_folder_for_spikein:
    type: Directory
    'sd:upstreamSource': "genome_indices_spikein/bowtie_indices"
    label: "Secondary genome index for spike-in normalization:"
    'sd:localLabel': true
    doc: "preprocessed genome index of spike-in organism for secondary alignment (of unaligned
      reads from primary alignment) and spike-in normalization, default should be E. coli K-12"

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

  control_file:
    type: File?
    default: null
    'sd:upstreamSource': "control_file/bambai_pair"
    'sd:localLabel': true
    label: "Use sample as a control (optional)"
    format: "http://edamontology.org/format_2572"
    doc: "If skipping upstream (i.e. not using an upstream sample as input here), please upload a BAM format file."

  broad_peak:
    type: boolean?
    default: False
    label: "Select to call broad peaks"
    doc: "Sets broad peak calling for MACS2 (toggle on for histone binding libraries)"
    
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

  exp_fragment_size:
    type: int?
    default: 150
    'sd:layout':
      advanced: true
    label: "Expected fragment size"
    doc: "Expected fragment size for MACS2"

  force_fragment_size:
    type: boolean?
    default: false
    'sd:layout':
      advanced: true
    label: "Force fragment size"
    doc: "Force MACS2 to use exp_fragment_size"

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
    label: "Number of bases to clip from the 5p end:"
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

  fragment_length_filter:
    type:
    - "null"
    - type: enum
      name: "Fragment Length Filter"
      symbols:
      - default_below_1000
      - histones_130_to_300
      - TF_below_130
    default: "default_below_1000"
    'sd:layout':
      advanced: true
    label: "Fragment Length Filter will retain fragments between set ranges for peak analysis."
    'sd:localLabel': true
    doc: "Fragment Length Filter options: 1) default_below_1000 retains fragments <1000 bp, 2) histones_130_to_300 retains fragments
      between 130-300 bp, and 3) TF_below_130 retains fragments <130 bp."

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

  multimapped_fastq_files:
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
        Title: 'FASTQ 1 Base frequency plot (post trim)'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'FASTQ 1 Quality Control (post trim)'
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
        Title: 'FASTQ 2 Base frequency plot (post trim)'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'FASTQ 2 Quality Control (post trim)'
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

  bigwig_raw:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated BigWig file from filtered bam"
    outputSource: bam_to_bigwig/bigwig_file

  bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "scaled BigWig file"
    doc: "Generated SCALED BigWig file from filtered bam, used as input for diffbind"
    outputSource: bam_to_bigwig_scaled/bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "Scaled BigWig"
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
    type: File?
    label: "Called peaks"
    format: "http://edamontology.org/format_3468"
    doc: "XLS file to include information about called peaks"
    outputSource: macs2_callpeak/peak_xls_file

  macs2_narrow_peaks:
    type: File?
    label: "Narrow peaks"
    format: "http://edamontology.org/format_3613"
    doc: "Contains the peak locations together with peak summit, pvalue and qvalue"
    outputSource: macs2_callpeak/narrow_peak_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Narrow peaks"
        displayMode: "COLLAPSE"
        height: 40

  macs2_broad_peaks:
    type: File?
    label: "Broad peaks"
    format: "http://edamontology.org/format_3614"
    doc: "Contains the peak locations together with peak summit, pvalue and qvalue"
    outputSource: macs2_callpeak/broad_peak_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Broad peaks"
        displayMode: "COLLAPSE"
        height: 40

  macs2_peak_summits:
    type: File?
    label: "Peak summits"
    format: "http://edamontology.org/format_3003"
    doc: "Contains the peak summits locations for every peaks"
    outputSource: macs2_callpeak/peak_summits_file

  macs2_moder_r:
    type: File?
    label: "MACS2 generated R script"
    format: "http://edamontology.org/format_2330"
    doc: "R script to produce a PDF image about the model based on your data"
    outputSource: macs2_callpeak/moder_r_file

  macs2_gapped_peak:
    type: File?
    label: "Gapped peaks"
    format: "http://edamontology.org/format_3586"
    doc: "Contains both the broad region and narrow peaks"
    outputSource: macs2_callpeak/gapped_peak_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Gapped peaks"
        displayMode: "COLLAPSE"
        height: 40

  macs2_log:
    type: File?
    label: "MACS2 log"
    format: "http://edamontology.org/format_2330"
    doc: "MACS2 output log"
    outputSource: macs2_callpeak/macs_log

  macs2_fragment_stat:
    type: File?
    label: "FRAGMENT, FRAGMENTE, ISLANDS"
    format: "http://edamontology.org/format_2330"
    doc: "fragment, calculated fragment, islands count from MACS2 results"
    outputSource: macs2_callpeak/macs2_stat_file

  estimated_fragment_size:
    type: int
    label: "Estimated fragment size"
    doc: "Estimated fragment size for downstream analyses"
    outputSource: macs2_callpeak/macs2_fragments_calculated

  iaintersect_result:
    type: File?
    format: "http://edamontology.org/format_3475"
    label: "gene annotated peaks file"
    doc: "nearest gene annotation per peak"
    outputSource: island_intersect/result_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Annotated Peaks'
        Title: 'Peak list with nearest gene annotation'

  atdp_log:
    type: File
    label: "ATDP log"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density generated log"
    outputSource: average_tag_density/log_file

  atdp_result:
    type: File
    label: "ATDP results"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density generated results"
    outputSource: average_tag_density/result_file
    'sd:visualPlugins':
    - scatter:
        tab: 'QC Plots'
        Title: 'Average Tag Density'
        xAxisTitle: 'Distance From TSS (bases)'
        yAxisTitle: 'Average Tag Density (per bp)'
        colors: ["#b3de69"]
        height: 500
        data: [$1, $2]
        comparable: "atdp"


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
      Returns a float scaling factor equal to (10000/x), where x is the
      number of spike-in aligned reads. This scale is then used for normalizing
      the bam-to-bedgraph that is used as input to the `bam-bedgraph-bigwig.cwl`
      tool (in step `bam_to_bigwig_scaled`).
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

  fragment_counts:
    run: ../tools/bedtools-fragmentcounts.cwl
    in:
      bam_file: samtools_sort_index_after_rmdup/bam_bai_pair
      scale: get_scale_from_spikein/scaling_factor
      output_prefix:
        source: samtools_sort_index_after_rmdup/bam_bai_pair
        valueFrom: $(get_root(self.basename))
      chrom_length_file: chrom_length
      fragment_length_filter:
        source: fragment_length_filter
        valueFrom: $(self)
    out: [sorted_bed, sorted_bed_scaled, log_file_stderr, log_file_stdout]
    doc: |
      Formatting alignment file to account for fragments based on PE bam.
      Output is a filtered (see fragment_length_filter input) and scaled
      (normalized) bed file to be used as input for peak calling.

  macs2_callpeak:
    label: "Peak detection"
    doc: |
      Identifies enriched with aligned reads genome areas. Those areas correspond to the
      transcription factor binding sites.
    run: ../tools/macs2-callpeak-biowardrobe-only.cwl
    in:
      treatment_file: fragment_counts/sorted_bed_scaled
      control_file: control_file
      nolambda:
        source: control_file
        valueFrom: $(!self)
      genome_size: genome_size
      mfold:
        default: "4 40"
      verbose:
        default: 3
      nomodel: force_fragment_size
      extsize: exp_fragment_size
      bw: exp_fragment_size
      broad: broad_peak
      call_summits:
        source: broad_peak
        valueFrom: $(!self)
      keep_dup:
        default: auto
      q_value:
        default: 0.05
      format_mode:
        default: AUTO
      buffer_size:
        default: 10000
    out:
    - peak_xls_file
    - narrow_peak_file
    - peak_summits_file
    - broad_peak_file
    - moder_r_file
    - gapped_peak_file
    - treat_pileup_bdg_file
    - control_lambda_bdg_file
    - macs_log
    - macs2_stat_file
    - macs2_fragments_calculated

  filter_fragment_lengths:
    run: ../tools/samtools-filter-fragmentlengths.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
      fragment_length_filter:
        source: fragment_length_filter
        valueFrom: $(self)
    out: [log_file_stdout, log_file_stderr, filtered_bam]
    doc: |
      Filtering bam file before scaling for IGV.

  samtools_sort_index_filtered:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: filter_fragment_lengths/filtered_bam
      threads: threads
    out: [bam_bai_pair]

  bam_to_bigwig:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index_filtered/bam_bai_pair
      chrom_length_file: chrom_length
      pairchip:
        default: true
    out: [bigwig_file]

  bam_to_bigwig_scaled:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index_filtered/bam_bai_pair
      chrom_length_file: chrom_length
      scale: get_scale_from_spikein/scaling_factor
      pairchip:
        default: true
      bigwig_filename:
        source: samtools_sort_index_filtered/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_scaled.bigWig")
    out: [bigwig_file]
    
  get_bam_statistics_after_filtering:
    run: ../tools/samtools-stats.cwl
    in:
      bambai_pair: samtools_sort_index_filtered/bam_bai_pair
      output_filename:
        source: samtools_sort_index_filtered/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_bam_statistics_report_after_filtering.txt")
    out: [log_file, ext_is_section]

  get_stat:
    run: ../tools/collect-statistics-cutandrun.cwl
    in:
      trimgalore_report_fastq_1: bypass_trim/selected_report_file_1
      trimgalore_report_fastq_2: bypass_trim/selected_report_file_2
      bowtie_alignment_report: bowtie_aligner/log_file
      bam_statistics_report: get_bam_statistics/log_file
      bam_statistics_after_filtering_report: get_bam_statistics_after_filtering/log_file
      cutandrun_called_peaks: macs2_callpeak/peak_xls_file
      paired_end:
        default: True
      output_prefix:
        source: macs2_callpeak/peak_xls_file
        valueFrom: $(get_root(self.basename))
    out: [collected_statistics_yaml, collected_statistics_tsv, collected_statistics_md, mapped_reads]
    doc: |
      Statistics pulled from alignment reports and other logs, as well as spike-in normalized peak calling.

  stats_for_vis:
    run: ../tools/collect-statistics-frip.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
      called_peaks_norm: macs2_callpeak/peak_xls_file
      collected_statistics_md: get_stat/collected_statistics_md
      collected_statistics_tsv: get_stat/collected_statistics_tsv
      collected_statistics_yaml: get_stat/collected_statistics_yaml
      spikein_reads_mapped: get_spikein_bam_statistics/reads_mapped
      peak_caller:
        default: "MACS2"
    out: [modified_file_md, modified_file_tsv, modified_file_yaml, log_file_stdout, log_file_stderr]

  island_intersect:
    label: "Peak annotation"
    doc: |
      Assigns nearest genes to peaks to explore the biological implication of the open
      chromatin binding sites.
    run: ../tools/iaintersect.cwl
    in:
      input_filename: macs2_callpeak/peak_xls_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
    out: [result_file, log_file]

  average_tag_density:
    label: "Read enrichment around genes TSS"
    doc: |
      Generates average tag density plot around genes TSS as a lot of cis-regulatory
      elements are close to the TSS of their targets.
    run: ../tools/atdp.cwl
    in:
      input_file: samtools_sort_index_after_rmdup/bam_bai_pair
      annotation_filename: annotation_file
      fragmentsize_bp: macs2_callpeak/macs2_fragments_calculated
      avd_window_bp:
        default: 5000
      avd_smooth_bp:
        default: 50
      ignore_chr:
        default: chrM
      double_chr:
        default: "chrX chrY"
      avd_heat_window_bp:
        default: 200
      mapped_reads: get_stat/mapped_reads
    out: [result_file, log_file]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "CUT&RUN/TAG MACS2 pipeline paired-end"
label: "CUT&RUN/TAG MACS2 pipeline paired-end"
s:alternateName: "CUR&RUN and CUT&TAG basic sparse enrichment analysis workflow for a paired-end experiment with Trim Galoreand MACS2"

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
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X


doc: |
  A basic analysis workflow for paired-read CUT&RUN and CUT&TAG sequencing experiments. These sequencing library prep methods are ultra-sensitive chromatin mapping technologies compared to the ChIP-Seq methodology. Its primary benefits include 1) length filtering, 2) a higher signal-to-noise ratio, and 3) built-in normalization for between sample comparisons.

  This workflow utilizes the tool MACS2 which calls enriched regions in the target sequence data by identifying the top regions by area under a poisson distribution (of the alignment pileup). This workflow is loosely based on the [CUT-RUNTools-2.0 pipeline](https://github.com/fl-yu/CUT-RUNTools-2.0) pipeline, and the ChIP-Seq pipeline from [BioWardrobe](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465) was used as a CWL template.

  ### __Inputs__
  *General Info (required\*):*
  - Experiment short name/Alias* - a unique name for the sample (e.g. what was used on tubes while processing it)
  - Cells* - sample cell type or organism name
  - Conditions* - experimental condition name
  - Catalog # - catalog number for cells from vender/supplier
  - Primary [genome index](https://scidap.com/tutorials/basic/genome-indices) for peak calling* - preprocessed genome index of sample organism for primary alignment and peak calling
  - Secondary [genome index](https://scidap.com/tutorials/basic/genome-indices) for spike-in normalization* - preprocessed genome index of spike-in organism for secondary alignment (of unaligned reads from primary alignment) and spike-in normalization, default should be E. coli K-12
  - FASTQ file for R1* - read 1 file of a pair-end library
  - FASTQ file for R2* - read 2 file of a pair-end library

  *Advanced:*
  - 
  - Number of bases to clip from the 3p end - used by bowtie aligner to trim <int> bases from 3' (right) end of reads
  - Number of bases to clip from the 5p end - used by bowtie aligner to trim <int> bases from 5' (left) end of reads
  - Call samtools rmdup to remove duplicates from sorted BAM file? - toggle on/off to remove duplicate reads from analysis
  - Fragment Length Filter will retain fragments between set base pair (bp) ranges for peak analysis - drop down menu
      - `default_below_1000` retains fragments <1000 bp
      - `histones_130_to_300` retains fragments between 130-300 bp
      - `TF_below_130` retains fragments <130 bp
  - Max distance (bp) from gene TSS (in both directions) overlapping which the peak will be assigned to the promoter region - default set to `1000`
  - Max distance (bp) from the promoter (only in upstream directions) overlapping which the peak will be assigned to the upstream region - default set to `20000`
  - Number of threads for steps that support multithreading - default set to `2`

  ### __Outputs__
  Intermediate and final downloadable outputs include:
  - IGV with gene, BigWig (raw and normalized), and stringent peak tracks
  - quality statistics and visualizations for both R1/R2 input FASTQ files
  - coordinate sorted BAM file with associated BAI file for primary alignment
  - read pileup/coverage in BigWig format (raw and normalized)
  - cleaned bed files (containing fragment coordinates), and spike-in normalized peak-called BED files (also includes "narrow" and "broad" peaks).
  - stringent peak call bed file with nearest gene annotations per peak

  ### __Data Analysis Steps__
  1. Trimming the adapters with TrimGalore.
      - This step is particularly important when the reads are long and the fragments are short - resulting in sequencing adapters at the ends of reads. If adapter is not removed the read will not map. TrimGalore can recognize standard adapters, such as Illumina or Nextera/Tn5 adapters.
  2. Generate quality control statistics of trimmed, unmapped sequence data
  3. (Optional) Clipping of 5' and/or 3' end by the specified number of bases.
  4. Mapping reads to primary genome index with Bowtie.
      - Only uniquely mapped reads with less than 3 mismatches are used in the downstream analysis. Results are then sorted and indexed. Final outputs are in bam/bai format, which are also used to extrapolate effects of additional sequencing based on library complexity.
  5. (Optional) Removal of duplicates (reads/pairs of reads mapping to exactly the same location).
      - This step is used to remove reads overamplified during amplification of the library. Unfortunately, it may also remove "good" reads. We usually do not remove duplicates unless the library is heavily duplicated.
  6. Mapping unaligned reads from primary alignment to secondary genome index with Bowtie.
      - This step is used to obtain the number of reads for normalization, used to scale the pileups from the primary alignment. After normalization, sample pileups/peak may then be appropriately compared to one another assuming an equal use of spike-in material during library preparation. Note the default genome index for this step should be *E. coli* K-12 if no spike-in material was called out in the library protocol. Refer to [Step 16](https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=16#step-4A3D8C70DC3011EABA5FF3676F0827C5) of the "CUT&Tag Data Processing and Analysis Tutorial" by Zheng Y et al (2020). Protocol.io.
  7. Formatting alignment file to account for fragments based on paired-end BAM.
      - Generates a filtered and normalized bed file to be used as input for peak calling.
  8. Call enriched regions using MACS2.
      - This step called peaks (broad and narrow) using the MACS2 tool with default parameters and no normalization to a control sample.
  9. Generation and formatting of output files.
      - This step collects read, alignment, and peak statistics, as well asgenerates BigWig coverage/pileup files for display on the browser using IGV. The coverage shows the number of fragments that cover each base in the genome both normalized and unnormalized to the calculated spike-in scaling factor.


  ### __References__
    - Meers MP, Tenenbaum D, Henikoff S. (2019). Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling. Epigenetics and Chromatin 12(1):42.
    - Langmead B, Trapnell C, Pop M, Salzberg SL. Ultrafast and memory-efficient alignment of short DNA sequences to the human genome. Genome Biol 10:R25.
