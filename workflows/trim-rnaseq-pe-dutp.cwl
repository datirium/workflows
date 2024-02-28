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
    type:
    - File
    - type: array
      items: File
    label: "FASTQ read1 input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, received after paired end sequencing"

  fastq_file_downstream:
    type:
    - File
    - type: array
      items: File
    label: "FASTQ read2 input file"
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

  minimum_rpkm:
    type: float?
    default: 1
    label: "Minimum RPKM for Gene Body Average Tag Density Plot"
    doc: "Minimum RPKM for Gene Body Average Tag Density Plot"
    'sd:layout':
      advanced: true

  max_multimap:
    type: int?
    default: 1
    label: "Maximum number of loci the read is allowed to map to"
    doc: "Maximum number of loci the read is allowed to map to"
    'sd:layout':
      advanced: true

  max_multimap_anchor:
    type: int?
    default: 50
    label: "Maximum number of loci anchors are allowed to map to"
    doc: "Maximum number of loci anchors are allowed to map to"
    'sd:layout':
      advanced: true

# System dependent

  threads:
    type: int?
    default: 2
    'sd:layout':
      advanced: true
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"


outputs:

  unaligned_fastq:
    type:
      - "null"
      - File[]
    format: "http://edamontology.org/format_1930"
    label: "Unaligned FASTQ file(s)"
    doc: "Unaligned FASTQ file(s)"
    outputSource: bowtie_aligner/unaligned_fastq

  bigwig_upstream:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated BigWig file for (+)strand reads"
    outputSource: bam_to_bigwig_upstream/bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "(+)strand BigWig"
        height: 120

  bigwig_downstream:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated BigWig file for (-)strand reads"
    outputSource: bam_to_bigwig_downstream/bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "(-)strand BigWig"
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

  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index/bam_bai_pair
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

  htseq_count_gene_expression_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "HTSeq: read counts grouped by gene_id"
    doc: "HTSeq: read counts grouped by gene_id"
    outputSource: htseq_count_gene_expression/feature_counts_report_file

  htseq_count_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "HTSeq: stdout log"
    doc: "HTSeq: stdout log"
    outputSource: htseq_count_gene_expression/stdout_log

  htseq_count_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "HTSeq: stderr log"
    doc: "HTSeq: stderr log"
    outputSource: htseq_count_gene_expression/stderr_log

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
    label: "BioWardrobe compatible log"
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
        comparable: "isdp"

  trim_report_upstream:
    type: File
    label: "TrimGalore report FASTQ 1"
    doc: "TrimGalore generated log for FASTQ 1"
    outputSource: trim_fastq/report_file

  trim_report_downstream:
    type: File
    label: "TrimGalore report FASTQ 2"
    doc: "TrimGalore generated log for FASTQ 2"
    outputSource: trim_fastq/report_file_pair

  gene_body_report:
    type: File?
    format: "http://edamontology.org/format_3475"
    label: "Gene body average tag density plot for all isoforms longer than 1000 bp"
    doc: "Gene body average tag density plot for all isoforms longer than 1000 bp in TSV format"
    outputSource: get_gene_body/gene_body_report_file
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'Gene body average tag density plot'
        xAxisTitle: "Gene body percentile (5' -> 3')"
        yAxisTitle: "Average Tag Density (per percentile)"
        colors: ["#232C15"]
        data: [$2]
        comparable: "gbatdp"

  gene_body_plot_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "Gene body average tag density plot for all isoforms longer than 1000 bp"
    doc: "Gene body average tag density plot for all isoforms longer than 1000 bp in PDF format"
    outputSource: get_gene_body/gene_body_plot_pdf

  rpkm_distribution_plot_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "RPKM distribution plot for isoforms"
    doc: "RPKM distribution plot for isoforms in PDF format"
    outputSource: get_gene_body/rpkm_distribution_plot_pdf


steps:

  extract_fastq_upstream:
    run: ../tools/extract-fastq.cwl
    in:
      output_prefix:  
        default: "read_1"
      compressed_file: fastq_file_upstream
    out: [fastq_file]

  extract_fastq_downstream:
    run: ../tools/extract-fastq.cwl
    in:
      output_prefix:  
        default: "read_2"
      compressed_file: fastq_file_downstream
    out: [fastq_file]

  trim_fastq:
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

  star_aligner:
    run: ../tools/star-alignreads.cwl
    in:
      readFilesIn: [rename_upstream/target_file, rename_downstream/target_file]
      genomeDir: star_indices_folder
      outFilterMultimapNmax: max_multimap
      winAnchorMultimapNmax: max_multimap_anchor
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

  fastx_quality_stats_upstream:
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename_upstream/target_file
    out: [statistics_file]

  fastx_quality_stats_downstream:
    run: ../tools/fastx-quality-stats.cwl
    in:
      input_file: rename_downstream/target_file
    out: [statistics_file]

  samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: star_aligner/aligned_file
      sort_output_filename:
        source: rename_upstream/target_file
        valueFrom: $(self.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+'.bam')
      threads: threads
    out: [bam_bai_pair]

  bam_to_bigwig_upstream:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
      chrom_length_file: chrom_length_file
      mapped_reads_number:
        source: star_aligner/uniquely_mapped_reads_number
        valueFrom: $(self*2)
      bigwig_filename:
        source: samtools_sort_index/bam_bai_pair
        valueFrom: |
          ${
            let root = self.basename.split('.').slice(0,-1).join('.');
            let ext = "_upstream.bigWig";
            return (root == "")?self.basename+ext:root+ext;
          }
      strand:
        default: '+'
      dutp:
        default: true
    out: [bigwig_file]

  bam_to_bigwig_downstream:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
      chrom_length_file: chrom_length_file
      mapped_reads_number:
        source: star_aligner/uniquely_mapped_reads_number
        valueFrom: $(-self*2)
      bigwig_filename:
        source: samtools_sort_index/bam_bai_pair
        valueFrom: |
          ${
            let root = self.basename.split('.').slice(0,-1).join('.');
            let ext = "_downstream.bigWig";
            return (root == "")?self.basename+ext:root+ext;
          }
      strand:
        default: '-'
      dutp:
        default: true
    out: [bigwig_file]

  bowtie_aligner:
    run: ../tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: rename_upstream/target_file
      downstream_filelist: rename_downstream/target_file
      indices_folder: bowtie_indices_folder
      clip_3p_end: clip_3p_end
      clip_5p_end: clip_5p_end
      v:
        default: 3
      m:
        default: 1
      unaligned_prefix:
        default: "unaligned_reads"
      sam:
        default: true
      threads: threads
    out: [log_file, unaligned_fastq]

  rpkm_calculation:
    run: ../tools/geep.cwl
    in:
      bam_file: samtools_sort_index/bam_bai_pair
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

  get_annotation_gtf:
    run: ../tools/ucsc-genepredtogtf.cwl
    in:
      annotation_tsv_file: annotation_file
    out:
    - annotation_gtf_file

  htseq_count_gene_expression:
    run: ../tools/htseq-count.cwl
    in:
      alignment_bam_file: samtools_sort_index/bam_bai_pair
      annotation_gtf_file: get_annotation_gtf/annotation_gtf_file
      strand_specific:
        default: "reverse"
      feature_type:
        default: "exon"
      feature_id:
        default: "gene_id"
    out:
    - feature_counts_report_file
    - stdout_log
    - stderr_log

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
        trimgalore_report_fastq_1: bypass_trim/selected_report_file_1
        trimgalore_report_fastq_2: bypass_trim/selected_report_file_2
        star_alignment_report: star_aligner/log_final
        bowtie_alignment_report: bowtie_aligner/log_file
        bam_statistics_report: get_bam_statistics/log_file
        isoforms_file: rpkm_calculation/isoforms_file
        paired_end:
          default: true
      out: [collected_statistics_yaml, collected_statistics_tsv, collected_statistics_md]

  get_gene_body:
    run: ../tools/plugin-plot-rna.cwl
    in:
      annotation_file: annotation_file
      bambai_pair: samtools_sort_index/bam_bai_pair
      isoforms_file: rpkm_calculation/isoforms_file
      mapped_reads_number: star_aligner/uniquely_mapped_reads_number
      minimum_rpkm: minimum_rpkm
      pair:
        default: true
      strand_specificity:
        default: "reverse"
      threads: threads
    out:
    - gene_body_report_file
    - gene_body_plot_pdf
    - rpkm_distribution_plot_pdf


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "Trim Galore RNA-Seq pipeline paired-end strand specific"
label: "Trim Galore RNA-Seq pipeline paired-end strand specific"
s:alternateName: "Run RNA-Seq BioWardrobe basic analysis with strand specific pair-end data file"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/trim-rnaseq-pe-dutp.cwl
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
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898
      - class: s:Person
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0001-9102-5681


# doc:
#   $include: ../descriptions/trim-rnaseq-pe-dutp.md


doc: |
  Modified original [BioWardrobe's](https://biowardrobe.com) [PubMed ID:26248465](https://www.ncbi.nlm.nih.gov/pubmed/26248465)
  **RNA-Seq** basic analysis for a **pair-end** experiment.
  A corresponded input [FASTQ](http://maq.sourceforge.net/fastq.shtml) file has to be provided.

  Current workflow should be used only with the single-end RNA-Seq data. It performs the following steps:
  1. Trim adapters from input FASTQ files
  2. Use STAR to align reads from input FASTQ files according to the predefined reference indices; generate unsorted BAM file and alignment statistics file
  3. Use fastx_quality_stats to analyze input FASTQ files and generate quality statistics files
  4. Use samtools sort to generate coordinate sorted BAM(+BAI) file pair from the unsorted BAM file obtained on the step 1 (after running STAR)
  5. Generate BigWig file on the base of sorted BAM file
  6. Map input FASTQ files to predefined rRNA reference indices using Bowtie to define the level of rRNA contamination; export resulted statistics to file
  7. Calculate isoform expression level for the sorted BAM file and GTF/TAB annotation file using GEEP reads-counting utility; export results to file
