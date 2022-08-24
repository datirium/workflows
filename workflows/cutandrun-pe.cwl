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
    label: "Indexed genome folder (bowtie)"
    doc: "Path to bowtie indexed genome folder for treatment organism"

  indices_folder_for_spikein:
    type: Directory
    'sd:upstreamSource': "genome_indices_spikein/bowtie_indices"
    label: "Indexed genome folder (bowtie)"
    doc: "Path to bowtie indexed genome folder for spike-in organism"

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
    label: "FASTQ 1 input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, received after paired end sequencing"

  fastq_file_downstream:
    type:
      - File
      - type: array
        items: File
    label: "FASTQ 2 input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, received after paired end sequencing"

  exp_fragment_size:
    type: int?
    default: 150
    'sd:layout':
      advanced: true
    label: "Expected fragment size"
    doc: "Expected fragment size for experimental sample"

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

  remove_duplicates:
    type: boolean?
    default: false
    'sd:layout':
      advanced: true
    label: "Remove duplicates"
    doc: "Calls samtools rmdup to remove duplicates from sortesd BAM file"

  promoter_dist:
    type: int?
    default: 1000
    'sd:layout':
      advanced: true
    label: "Max distance from gene TSS (in both direction) overlapping which the peak will be assigned to the promoter region"
    doc: "Max distance from gene TSS (in both direction) overlapping which the peak will be assigned to the promoter region"

  upstream_dist:
    type: int?
    default: 20000
    'sd:layout':
      advanced: true
    label: "Max distance from the promoter (only in upstream direction) overlapping which the peak will be assigned to the upstream region"
    doc: "Max distance from the promoter (only in upstream direction) overlapping which the peak will be assigned to the upstream region"

  threads:
    type: int?
    default: 2
    'sd:layout':
      advanced: true
    doc: "Number of threads for those steps that support multithreading"
    label: "Number of threads"

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
  
  get_stat_markdown:
    type: File?
    label: "Markdown formatted combined log"
    format: "http://edamontology.org/format_3835"
    doc: "Markdown formatted combined log"
    outputSource: get_stat/collected_statistics_md
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  get_stat_formatted_log:
    type: File?
    label: "Bowtie & Samtools Rmdup combined formatted log"
    format: "http://edamontology.org/format_3475"
    doc: "Processed and combined Bowtie aligner and Samtools rmdup formatted log"
    outputSource: get_stat/collected_statistics_tsv
    'sd:visualPlugins':
    - tableView:
        vertical: true
        tab: 'Overview'
    'sd:preview':
      'sd:visualPlugins':
      - pie:
          colors: ['#b3de69', '#99c0db', '#fb8072', '#fdc381']
          data: [$2, $3, $4, $5]

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

  bedgraph:
    type: File
    format: "http://edamontology.org/format_3583"
    label: "bedgraph file"
    doc: "Generated bedgraph file"
    outputSource: bam_to_bigwig/bedgraph_file 

  peak_tsv:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "bedgraph file"
    doc: "Bed file of enriched regions"
    outputSource: seacr_callpeak/bedgraph_file

  norm_peak_tsv:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "bedgraph file"
    doc: "Bed file of enriched regions"
    outputSource: seacr_callpeak_spikein_norm/bedgraph_file 
    

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

  preseq:
    label: "Sequencing depth estimation"
    doc: |
      Estimates the complexity of the sequencing library, evaluates how many reads can
      be expected from the additional sequencing of the same experiment.
    run: ../tools/preseq-lc-extrap.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
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
    out: [bigwig_file, bedgraph_file]

  bam_to_bigwig_scaled:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index_after_rmdup/bam_bai_pair
      chrom_length_file: chrom_length
      scale: get_scale_from_spikein/scaling_factor
      pairchip:
        default: true
    out: [bigwig_file, bedgraph_file]

  seacr_callpeak:
    run: ../tools/seacr.cwl
    in:
      treatment_bedgraph: bam_to_bigwig/bedgraph_file
      numeric_threshold:
        default: 0.01
      norm_control_to_treatment:
        default: "non"
      peakcalling_mode:
        default: "stringent"
      output_prefix:
        source: bam_to_bigwig/bedgraph_file
        valueFrom: $(get_root(self.basename))
    out: [peak_tsv_file]


  seacr_callpeak_spikein_norm:
    run: ../tools/seacr.cwl
    in:
      treatment_bedgraph: bam_to_bigwig_scaled/bedgraph_file
      numeric_threshold:
        default: 0.01
      norm_control_to_treatment:
        default: "non"
      peakcalling_mode:
        default: "stringent"
      output_prefix:
        source: bam_to_bigwig_scaled/bedgraph_file
        valueFrom: $(get_root(self.basename)+"_spikein_norm")
    out: [peak_tsv_file]
    doc: |
      input is normalized depth data using spike-in mapped reads (E. coli by default)
      Henikoff protocol, Section 16: https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=16#step-4A3D8C70DC3011EABA5FF3676F0827C5)

  get_stat_raw:
    run: ../tools/collect-statistics-cutandrun.cwl
    in:
      trimgalore_report_fastq_1: bypass_trim/selected_report_file_1
      trimgalore_report_fastq_2: bypass_trim/selected_report_file_2
      bowtie_alignment_report: bowtie_aligner/log_file
      bam_statistics_report: get_bam_statistics/log_file
      bam_statistics_after_filtering_report: get_bam_statistics_after_filtering/log_file
      peaks_raw: seacr_callpeak/peak_tsv_file
      preseq_results: preseq/estimates_file
      paired_end:
        default: True
    out: [collected_statistics_yaml, collected_statistics_tsv, mapped_reads, collected_statistics_md]

  get_stat_norm:
    run: ../tools/collect-statistics-cutandrun.cwl
    in:
      trimgalore_report_fastq_1: bypass_trim/selected_report_file_1
      trimgalore_report_fastq_2: bypass_trim/selected_report_file_2
      bowtie_alignment_report: bowtie_aligner/log_file
      bam_statistics_report: get_bam_statistics/log_file
      bam_statistics_after_filtering_report: get_bam_statistics_after_filtering/log_file
      peaks_norm: seacr_callpeak_spikein_norm/peak_tsv_file
      preseq_results: preseq/estimates_file
      paired_end:
        default: True
    out: [collected_statistics_yaml, collected_statistics_tsv, mapped_reads, collected_statistics_md]



$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "cutandrun with trimgalore pipeline paired-end"
label: "cutandrun with trimgalore pipeline paired-end"
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
  **SEACR** basic analysis workflow for a **paired-read** experiment with Trim Galore.
    
  This SEACR (Sparse Enrichment Analysis of Cut & Run data) pipeline calls
  enriched regions in target data by selecting the top 1% of regions by AUC.
  It is based on the
  [CUT-RUNTools-2.0 pipeline](https://github.com/fl-yu/CUT-RUNTools-2.0)
  pipeline, and the ChIP-Seq pipeline from [BioWardrobe](https://biowardrobe.com)
  [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465) was used as a CWL
  template.

  ### Data Analysis
  SciDAP starts from the .fastq files which most DNA cores and commercial NGS
  companies return. Starting from raw data allows us to ensure that all
  experiments have been processed in the same way and simplifies the deposition
  of data to GEO upon publication. The data can be uploaded from the user's
  computer, downloaded directly from an ftp server of the core facility by
  providing a URL, or from GEO by providing an SRA accession number.
  
  Our current pipelines include the following steps:
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
          is more MACS2.
  7.  Generation of BigWig coverage files for display on the browser. The coverage
          shows the number of fragments at each base in the genome normalized to
          the number of millions of mapped reads.


  Note, the upstream analyses should not have duplicates removed




  ### Details
  _Trim Galore_ is a wrapper around [Cutadapt](https://github.com/marcelm/cutadapt)
  and [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) to consistently
  apply adapter and quality trimming to FastQ files, with extra functionality for RRBS data.

  A [FASTQ](http://maq.sourceforge.net/fastq.shtml) input file has to be provided.


  In outputs it returns coordinate sorted BAM file alongside with index BAI file,
  quality statistics for both the input FASTQ files, reads coverage in a form of BigWig file,
  peaks calling data in a form of narrowPeak or broadPeak files, islands with the assigned nearest
  genes and region type, data for average tag density plot (on the base of BAM file).

  Workflow starts with running fastx_quality_stats (steps fastx_quality_stats_upstream and
  fastx_quality_stats_downstream) from FASTX-Toolkit to calculate quality statistics for both upstream
  and downstream input FASTQ files. At the same time Bowtie is used to align reads from input FASTQ
  files to reference genome (Step bowtie_aligner). The output of this step is unsorted SAM file which
  is being sorted and indexed by samtools sort and samtools index (Step samtools_sort_index).
  Depending on workflow’s input parameters indexed and sorted BAM file
  could be processed by samtools rmdup (Step samtools_rmdup) to remove all possible read duplicates.
  In a case when removing duplicates is not necessary the step returns original input BAM and BAI
  files without any processing. If the duplicates were removed the following step
  (Step samtools_sort_index_after_rmdup) reruns samtools sort and samtools index with BAM and BAI files,
  if not - the step returns original unchanged input files. Right after that macs2 callpeak performs
  peak calling (Step seacr_callpeak). On the base of returned outputs the next step
  (Step macs2_island_count) calculates the number of islands and estimated fragment size. If the last
  one is less than 80 (hardcoded in a workflow) macs2 callpeak is rerun again with forced fixed
  fragment size value (Step seacr_callpeak_forced). If at the very beginning it was set in workflow
  input parameters to force run peak calling with fixed fragment size, this step is skipped and the
  original peak calling results are saved. In the next step workflow again calculates the number
  of islands and estimated fragment size (Step macs2_island_count_forced) for the data obtained from
  seacr_callpeak_forced step. If the last one was skipped the results from macs2_island_count_forced step
  are equal to the ones obtained from macs2_island_count step.
  Next step (Step macs2_stat) is used to define which of the islands and estimated fragment size should be used
  in workflow output: either from macs2_island_count step or from macs2_island_count_forced step. If input
  trigger of this step is set to True it means that seacr_callpeak_forced step was run and it returned different
  from seacr_callpeak step results, so macs2_stat step should return [fragments_new, fragments_old, islands_new],
  if trigger is False the step returns [fragments_old, fragments_old, islands_old], where sufix "old" defines
  results obtained from macs2_island_count step and sufix "new" - from macs2_island_count_forced step.
  The following two steps (Step bamtools_stats and bam_to_bigwig) are used to calculate coverage on the base
  of input BAM file and save it in BigWig format. For that purpose bamtools stats returns the number of
  mapped reads number which is then used as scaling factor by bedtools genomecov when it performs coverage
  calculation and saves it in BED format. The last one is then being sorted and converted to BigWig format by
  bedGraphToBigWig tool from UCSC utilities. Step get_stat is used to return a text file with statistics
  in a form of [TOTAL, ALIGNED, SUPRESSED, USED] reads count. Step island_intersect assigns genes and regions
  to the islands obtained from seacr_callpeak_forced. Step average_tag_density is used to calculate data for
  average tag density plot on the base of BAM file.
