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
  control_file: "chipseq-pe.cwl"


inputs:

  indices_folder:
    type: Directory
    'sd:upstreamSource': "genome_indices/bowtie_indices"
    label: "Indexed genome folder (bowtie)"
    doc: "Path to indexed genome folder by **bowtie**"

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
    doc: "MACS2 effective genome size: hs, mm, ce, dm or number, for example 2.7e9"

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
    label: "Use experiment as a control"
    format: "http://edamontology.org/format_2572"
    doc: "Use experiment as a control for MACS2 peak calling"

  broad_peak:
    type: boolean?
    default: False
    label: "Callpeak broad"
    doc: "Set to call broad peak for MACS2"

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

  iaintersect_log:
    type: File
    label: "Island intersect log"
    format: "http://edamontology.org/format_3475"
    doc: "Iaintersect generated log"
    outputSource: island_intersect/log_file

  iaintersect_result:
    type: File
    label: "Island intersect results"
    format: "http://edamontology.org/format_3475"
    doc: "Iaintersect generated results"
    outputSource: island_intersect/result_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Peak Calling'
        Title: 'Islands list'

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

  macs2_fragment_stat:
    type: File?
    label: "FRAGMENT, FRAGMENTE, ISLANDS"
    format: "http://edamontology.org/format_2330"
    doc: "fragment, calculated fragment, islands count from MACS2 results"
    outputSource: macs2_callpeak/macs2_stat_file

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
  
  estimated_fragment_size:
    type: int
    label: "Estimated fragment size"
    doc: "Estimated fragment size for downstream analyses"
    outputSource: macs2_callpeak/macs2_fragments_calculated

  mapped_reads_number:
    type: int
    label: "Mapped reads number"
    doc: "Mapped reads number for downstream analyses"
    outputSource: get_stat/mapped_reads
    

steps:

  extract_fastq_upstream:
    label: "Loading unmapped sequence data for read 1"
    doc: |
      Most DNA cores and commercial NGS companies return unmapped sequence data in FASTQ format.
      The data can be uploaded from users computer, downloaded directly from an ftp server of
      the core facility by providing a URL or from GEO by providing SRA accession number.
    run: ../tools/extract-fastq.cwl
    in:
      output_prefix:  
        default: "read_1"
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
      output_prefix:  
        default: "read_2"
      compressed_file: fastq_file_downstream
    out: [fastq_file]

  fastx_quality_stats_upstream:
    label: "Quality control of unmapped sequence data for read 1"
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: extract_fastq_upstream/fastq_file
    out: [statistics_file]

  fastx_quality_stats_downstream:
    label: "Quality control of unmapped sequence data for read 2"
    doc: |
      Evaluates the quality of your sequence data. Provides per base quality scores as well as
      base frequencies along the reads. These metrics can be used to identify whether your data
      has any problems that should be taken into account in the subsequent analysis steps.
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: extract_fastq_downstream/fastq_file
    out: [statistics_file]

  bowtie_aligner:
    label: "Alignment to reference genome"
    doc: |
      Aligns reads to the reference genome keeping only uniquely mapped reads with
      less than 3 mismatches.
    run: ../tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: extract_fastq_upstream/fastq_file
      downstream_filelist: extract_fastq_downstream/fastq_file
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

  macs2_callpeak:
    label: "Peak detection"
    doc: |
      Identifies enriched with aligned reads genome areas. Those areas correspond to the
      transcription factor binding sites.
    run: ../tools/macs2-callpeak-biowardrobe-only.cwl
    in:
      treatment_file: samtools_sort_index_after_rmdup/bam_bai_pair
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
        default: BAMPE
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

  bam_to_bigwig:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index_after_rmdup/bam_bai_pair
      chrom_length_file: chrom_length
      mapped_reads_number: get_stat/mapped_reads
      pairchip:
        default: true
    out: [bigwig_file]

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

  get_stat:
    run: ../tools/collect-statistics-chip-seq.cwl
    in:
      bowtie_alignment_report: bowtie_aligner/log_file
      bam_statistics_report: get_bam_statistics/log_file
      bam_statistics_after_filtering_report: get_bam_statistics_after_filtering/log_file
      macs2_called_peaks: macs2_callpeak/peak_xls_file
      preseq_results: preseq/estimates_file
      paired_end:
        default: True
    out: [collected_statistics_yaml, collected_statistics_tsv, mapped_reads, collected_statistics_md]

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

label: "ChIP-Seq pipeline paired-end"
s:name: "ChIP-Seq pipeline paired-end"
s:alternateName: "ChIP-Seq basic analysis workflow for a paired-end experiment"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/chipseq-pe.cwl
s:codeRepository: https://github.com/datirium/workflows
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

# doc:
#   $include: ../descriptions/chipseq-pe.md


doc: |
  The original [BioWardrobe's](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465)
  **ChIP-Seq** basic analysis workflow for a **paired-end** experiment.
  A [FASTQ](http://maq.sourceforge.net/fastq.shtml) input file has to be provided.

  The pipeline produces a sorted BAM file alongside with index BAI file, quality
  statistics of the input FASTQ file, coverage by estimated fragments as a BigWig file, peaks calling
  data in a form of narrowPeak or broadPeak files, islands with the assigned nearest genes and
  region type, data for average tag density plot.

  Workflow starts with step *fastx\_quality\_stats* from FASTX-Toolkit
  to calculate quality statistics for input FASTQ file.

  At the same time `bowtie` is used to align
  reads from input FASTQ file to reference genome *bowtie\_aligner*. The output of this step
  is an unsorted SAM file which is being sorted and indexed by `samtools sort` and `samtools index`
  *samtools\_sort\_index*.

  Depending on workflow’s input parameters indexed and sorted BAM file
  can be processed by `samtools rmdup` *samtools\_rmdup* to get rid of duplicated reads.
  If removing duplicates is not required the original BAM and BAI
  files are returned. Otherwise step *samtools\_sort\_index\_after\_rmdup* repeat `samtools sort` and `samtools index` with BAM and BAI files without duplicates.

  Next `macs2 callpeak` performs peak calling *macs2\_callpeak* and the next step
  reports *macs2\_island\_count*  the number of islands and estimated fragment size. If the latter
  is less that 80bp (hardcoded in the workflow) `macs2 callpeak` is rerun again with forced fixed
  fragment size value (*macs2\_callpeak\_forced*). It is also possible to force MACS2 to use pre set  fragment size in the first place.

  Next step (*macs2\_stat*) is used to define which of the islands and estimated fragment size should be used
  in workflow output: either from *macs2\_island\_count* step or from *macs2\_island\_count\_forced* step. If input
  trigger of this step is set to True it means that *macs2\_callpeak\_forced* step was run and it returned different
  from *macs2\_callpeak* step results, so *macs2\_stat* step should return [fragments\_new, fragments\_old, islands\_new],
  if trigger is False the step returns [fragments\_old, fragments\_old, islands\_old], where sufix "old" defines
  results obtained from *macs2\_island\_count* step and sufix "new" - from *macs2\_island\_count\_forced* step.

  The following two steps (*bamtools\_stats* and *bam\_to\_bigwig*) are used to calculate coverage from BAM file and save it in BigWig format. For that purpose bamtools stats returns the number of
  mapped reads  which is then used as scaling factor by bedtools genomecov when it performs coverage
  calculation and saves it as a BEDgraph file whichis then  sorted and converted to BigWig format by
  bedGraphToBigWig tool from UCSC utilities. Step *get\_stat* is used to return a text file with statistics
  in a form of [TOTAL, ALIGNED, SUPRESSED, USED] reads count.

  Step *island\_intersect* assigns nearest genes and regions to the islands obtained from *macs2\_callpeak\_forced*.
  Step *average\_tag\_density* is used to calculate data for average tag density plot from the BAM file.