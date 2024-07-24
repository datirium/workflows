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
    format: "http://edamontology.org/format_3475"
    'sd:upstreamSource': "genome_indices/annotation"
    doc: "GTF or TAB-separated annotation file"

  annotation_gtf_file:
    type: File
    label: "GTF annotation file"
    format: "http://edamontology.org/format_2306"
    'sd:upstreamSource': "genome_indices/annotation_gtf"
    doc: "GTF annotation file"

  fastq_file:
    type:
    - File
    - type: array
      items: File
    label: "FASTQ input file(s)"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format"

  use_umi:
    type: boolean?
    default: true
    'sd:layout':
      advanced: true
    label: "Use UMIs"
    doc: "Use UMIs (for FWD-UMI libraries)"

  strand_specificity:
    type:
    - "null"
    - type: enum
      symbols:
      - "yes"
      - "no"
      - "reverse"
    default: "yes"
    'sd:layout':
      advanced: true
    label: "Strand specificity. 'Yes' for FWD or FWD-UMI analyses, 'Reverse' for REV, 'No' to disable"
    doc: |
      Whether the data is from a strand-specific assay. For stranded=no, a read is
      considered overlapping with a feature regardless of whether it is mapped to
      the same or the opposite strand as the feature. For stranded=yes and single-end
      reads, the read has to be mapped to the same strand as the feature. For paired-end
      reads, the first read has to be on the same strand and the second read on the
      opposite strand. For stranded=reverse, these rules are reversed.

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

  threads:
    type: int?
    default: 1
    'sd:layout':
      advanced: true
    label: "Number of threads"
    doc: "Number of threads for those steps that support multi-threading"


outputs:

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

  fastx_statistics:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "FASTQ statistics"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats/statistics_file
    'sd:visualPlugins':
    - line:
        tab: 'QC Plots'
        Title: 'Base frequency plot'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Frequency'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$13, $14, $15, $16, $17]
    - boxplot:
        tab: 'QC Plots'
        Title: 'Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"]
        data: [$11, $7, $8, $9, $12]

  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index_after_dedup/bam_bai_pair
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

  gene_expression_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Gene expression"
    doc: "Gene expression"
    outputSource: group_transcript_expression/gene_expression_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Gene Expression'
        Title: 'Read counts grouped by gene'

  # common_tss_expression_file:
  #   type: File
  #   format: "http://edamontology.org/format_3475"
  #   label: "Common TSS expression"
  #   doc: "Common TSS expression"
  #   outputSource: group_transcript_expression/common_tss_expression_file

  rpkm_genes:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "GEEP: expression grouped by gene name"
    doc: "GEEP: expression grouped by gene name"
    outputSource: group_geep_transcript_expression/genes_file

  combined_gene_expression_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "HTSeq vs GEEP gene expression comparison"
    doc: |
      Merged by GeneId, Chrom, TxStart, TxEnd and Strand gene expression files
      with reported and renamed TotalReads columns.
    outputSource: feature_expression_merge/merged_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Gene Expression Comparison'
        Title: 'HTSeq vs GEEP gene expression comparison (read counts)'

  feature_expression_merge_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "HTSeq vs GEEP gene expression comparison stdout log"
    doc: "HTSeq vs GEEP gene expression comparison stdout log"
    outputSource: feature_expression_merge/stdout_log

  feature_expression_merge_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "HTSeq vs GEEP gene expression comparison stderr log"
    doc: "HTSeq vs GEEP gene expression comparison stderr log"
    outputSource: feature_expression_merge/stderr_log

  # geep_common_tss_expression_file:
  #   type: File
  #   format: "http://edamontology.org/format_3475"
  #   label: "GEEP: expression grouped by common TSS"
  #   doc: "GEEP: expression grouped by common TSS"
  #   outputSource: group_geep_transcript_expression/common_tss_file

  htseq_count_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "HTSeq: stdout log"
    doc: "HTSeq: stdout log"
    outputSource: htseq_count_transcript_expression/stdout_log

  htseq_count_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "HTSeq: stderr log"
    doc: "HTSeq: stderr log"
    outputSource: htseq_count_transcript_expression/stderr_log

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
          colors: ['#b3de69', '#99c0db', '#fdc381', '#fb8072', '#778899']
          data: [$2, $3, $4, $5, $6]

  bam_statistics_report:
    type: File
    label: "BAM statistics report"
    format: "http://edamontology.org/format_2330"
    doc: "BAM statistics report (after deduplication step)"
    outputSource: get_bam_statistics/log_file

  trim_adapters_stdout_log:
    type: File
    label: "cutadapt: stdout log"
    doc: "cutadapt: stdout log"
    outputSource: trim_adapters/stdout_log

  trim_adapters_stderr_log:
    type: File
    label: "cutadapt: stderr log"
    doc: "cutadapt: stderr log"
    outputSource: trim_adapters/stderr_log

  umi_tools_dedup_stdout_log:
    type: File
    label: "umi_tools dedup: stdout log"
    doc: "umi_tools dedup: stdout log"
    outputSource: umi_tools_dedup/stdout_log

  umi_tools_dedup_stderr_log:
    type: File
    label: "umi_tools dedup: stderr log"
    doc: "umi_tools dedup: stderr log"
    outputSource: umi_tools_dedup/stderr_log

  umi_tools_dedup_stats:
    type:
    - "null"
    - File[]
    label: "umi_tools dedup statistics"
    doc: "umi_tools dedup statistics"
    outputSource: umi_tools_dedup/output_stats

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

  extract_fastq:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file
      output_prefix:
        default: "read_1"
    out:
    - fastq_file

  move_umi_to_read_name:
    run: ../tools/custom-bash.cwl
    in:
      input_file: extract_fastq/fastq_file
      param:
        source: use_umi
        valueFrom: $(self?"true":"false")
      script:
        default: |
          #!/bin/bash
          FILE=$0
          BASENAME=$(basename "$FILE")
          if [ "$1" = "true" ]; then
            cat ${FILE} | awk '
            NR%4==1{ rd_name=$1; rd_info=$2 }
            NR%4==2{ umi=substr($1,1,10); rd_seq=substr($1,11) }
            NR%4==0{ print rd_name"_"umi" "rd_info; print rd_seq; print "+"; print substr($1,11) }' > ${BASENAME}
          else
            cp ${FILE} ${BASENAME}
          fi
    out:
    - output_file

  trim_adapters:
    in:
      fastq_file: move_umi_to_read_name/output_file
    out:
    - trimmed_file
    - stdout_log
    - stderr_log
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: scidap/trimgalore:v0.6.6
      inputs:
        bash_script:
          type: string?
          default: |
            #!/bin/bash
            FILE=$0
            BASENAME=$(basename "$FILE")
            cat ${FILE} |
            cutadapt -m 20 -O 20 -a "polyA=A{20}" -a "QUALITY=G{20}" -n 2 - |
            cutadapt -m 20 -O 3 --nextseq-trim=10  -a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000" - |
            cutadapt -m 20 -O 3 -a "r1polyA=A{18}" - |
            cutadapt -m 20 -O 20 -g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20" --discard-trimmed -o ${BASENAME} -
          inputBinding:
            position: 1
        fastq_file:
          type: File
          inputBinding:
            position: 3
      outputs:
        trimmed_file:
          type: File
          outputBinding:
            glob: $(inputs.fastq_file.basename)
        stdout_log:
          type: stdout
        stderr_log:
          type: stderr
      baseCommand: ["bash", "-c"]
      stdout: cutadapt_stdout.log
      stderr: cutadapt_stderr.log

  star_aligner:
    run: ../tools/star-alignreads.cwl
    in:
      readFilesIn: trim_adapters/trimmed_file
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
      input_file: trim_adapters/trimmed_file
    out:
    - statistics_file

  samtools_sort_index_before_dedup:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: star_aligner/aligned_file
      sort_output_filename:
        source: trim_adapters/trimmed_file
        valueFrom: $(get_root(self.basename)+".bam")
      threads: threads
    out:
    - bam_bai_pair

  umi_tools_dedup:
    run: ../tools/umi-tools-dedup.cwl
    in:
      bam_file: samtools_sort_index_before_dedup/bam_bai_pair
      multimapping_detection_method:
        default: "NH"
      trigger: use_umi
    out:
    - dedup_bam_file
    - output_stats
    - stdout_log
    - stderr_log

  samtools_sort_index_after_dedup:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: umi_tools_dedup/dedup_bam_file
      sort_output_filename:
        source: trim_adapters/trimmed_file
        valueFrom: $(get_root(self.basename)+".bam")
      threads: threads
      trigger: use_umi
    out:
    - bam_bai_pair

  htseq_count_transcript_expression:
    run: ../tools/htseq-count.cwl
    in:
      alignment_bam_file: samtools_sort_index_after_dedup/bam_bai_pair
      annotation_gtf_file: annotation_gtf_file
      strand_specific: strand_specificity
      feature_type:
        default: "exon"
      feature_id:
        default: "transcript_id"
      additional_id:
        default: "gene_id"
    out:
    - feature_counts_report_file
    - stdout_log
    - stderr_log

  group_transcript_expression:
    in:
      transcript_expression_file: htseq_count_transcript_expression/feature_counts_report_file
      annotation_file: annotation_file
      mapped_reads_number: star_aligner/uniquely_mapped_reads_number
    out:
    - transcript_expression_file
    - gene_expression_file
    - common_tss_expression_file
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap-deseq:v0.0.21
      requirements:          
      - class: InitialWorkDirRequirement
        listing:
        - entryname: process.R
          entry: |
            #!/usr/bin/env Rscript
            options(warn=-1)
            options("width"=300)
            suppressMessages(library(data.table))

            args <- commandArgs(trailingOnly = TRUE)

            transcript_counts <- read.table(args[1], sep="\t", header=FALSE, stringsAsFactors=FALSE)
            colnames(transcript_counts) <- c("RefseqId", "GeneId", "TotalReads")

            annotation <- read.table(args[2], sep="\t", header=TRUE, stringsAsFactors=FALSE)
            colnames(annotation) <- c("bin", "RefseqId", "Chrom",	"Strand",	"TxStart", "TxEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "GeneId",	"cdsStartStat",	"cdsEndStat",	"exonFrames")

            mapped_reads_number <- as.numeric(args[3])

            transcript_counts <- merge(transcript_counts, annotation, by=c("RefseqId", "GeneId"), sort = FALSE)
            transcript_counts <- transcript_counts[, c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand", "TotalReads")]

            transcript_counts_table <- setDT(transcript_counts)

            gene_counts <- as.data.frame(
              transcript_counts_table[
                ,
                .(
                  RefseqId = paste(sort(unique(RefseqId)), collapse = ","),
                  Chrom = Chrom[1],
                  TxStart = TxStart[1],
                  TxEnd = TxEnd[1],
                  Strand = Strand[1],
                  TotalReads = sum(TotalReads)
                ),
                by = GeneId
              ]
            )

            common_tss_counts_positive_strand_table <- transcript_counts_table[
              Strand=="+",
              .(
                RefseqId = paste(sort(unique(RefseqId)), collapse = ","),
                GeneId = paste(sort(unique(GeneId)), collapse = ","),
                TxEnd = max(TxEnd),
                TotalReads = sum(TotalReads)),
                by = .(Chrom, TxStart, Strand)
            ]
            
            common_tss_counts_negative_strand_table <- transcript_counts_table[
              Strand=="-",
              .(
                RefseqId = paste(sort(unique(RefseqId)), collapse = ","),
                GeneId = paste(sort(unique(GeneId)), collapse = ","),
                TxStart = min(TxStart),
                TotalReads = sum(TotalReads)),
                by = .(Chrom, TxEnd, Strand)
            ]
            
            common_tss_counts <- as.data.frame(
              rbind(common_tss_counts_positive_strand_table, common_tss_counts_negative_strand_table)
            )

            reordered_transcript_counts <- transcript_counts[, c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand", "TotalReads")]
            reordered_transcript_counts$Rpkm <- reordered_transcript_counts$TotalReads / ( (reordered_transcript_counts$TxEnd - reordered_transcript_counts$TxStart) / 1000 * mapped_reads_number / 1000000 )
            reordered_gene_counts <- gene_counts[, c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand", "TotalReads")]
            reordered_common_tss_counts <- common_tss_counts[, c("RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand", "TotalReads")]

            write.table(reordered_transcript_counts, file="transcript_expression.csv", sep=",", row.names=FALSE, col.names=TRUE, quote=FALSE)
            write.table(reordered_gene_counts, file="gene_expression.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
            write.table(reordered_common_tss_counts, file="common_tss_expression.tsv", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

      inputs:
        transcript_expression_file:
          type: File
          inputBinding:
            position: 1
        annotation_file:
          type: File
          inputBinding:
            position: 2
        mapped_reads_number:
          type: int
          inputBinding:
            position: 3
      outputs:
        transcript_expression_file:
          type: File
          outputBinding:
            glob: "transcript_expression.csv"          
        gene_expression_file:
          type: File
          outputBinding:
            glob: "gene_expression.tsv"
        common_tss_expression_file:
          type: File
          outputBinding:
            glob: "common_tss_expression.tsv"
      baseCommand: ["Rscript", "process.R"]

  bam_to_bigwig_upstream:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index_after_dedup/bam_bai_pair
      chrom_length_file: chrom_length_file
      mapped_reads_number: star_aligner/uniquely_mapped_reads_number
      bigwig_filename:
        source: samtools_sort_index_after_dedup/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_upstream.bigWig")
      strand:
        default: '+'
    out:
    - bigwig_file

  bam_to_bigwig_downstream:
    run: ../tools/bam-bedgraph-bigwig.cwl
    in:
      bam_file: samtools_sort_index_after_dedup/bam_bai_pair
      chrom_length_file: chrom_length_file
      mapped_reads_number:
        source: star_aligner/uniquely_mapped_reads_number
        valueFrom: $(-self)
      bigwig_filename:
        source: samtools_sort_index_after_dedup/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_downstream.bigWig")
      strand:
        default: '-'
    out:
    - bigwig_file

  bowtie_aligner:
    run: ../tools/bowtie-alignreads.cwl
    in:
      upstream_filelist: trim_adapters/trimmed_file
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
    out:
    - log_file

  get_bam_statistics:
    run: ../tools/samtools-stats.cwl
    in:
      bambai_pair: samtools_sort_index_after_dedup/bam_bai_pair
      output_filename:
        source: samtools_sort_index_after_dedup/bam_bai_pair
        valueFrom: $(get_root(self.basename)+"_bam_statistics_report.txt")
    out:
    - log_file

  get_stat:
    run: ../tools/collect-statistics-rna-quantseq.cwl
    in:
      star_alignment_report: star_aligner/log_final
      bowtie_alignment_report: bowtie_aligner/log_file
      bam_statistics_report: get_bam_statistics/log_file
      isoforms_file: group_transcript_expression/transcript_expression_file
    out:
    - collected_statistics_yaml
    - collected_statistics_tsv
    - collected_statistics_md

  geep_count_transcript_expression:
    run: ../tools/geep.cwl
    in:
      bam_file: samtools_sort_index_after_dedup/bam_bai_pair
      annotation_file: annotation_file
      rpkm_threshold:
        default: 0
      max_cycles:
        default: 0
      threads: threads
    out:
    - isoforms_file

  group_geep_transcript_expression:
    in:
      isoforms_file: geep_count_transcript_expression/isoforms_file
    out:
    - genes_file
    - common_tss_file
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap-deseq:v0.0.20
      inputs:
        bash_script:
          type: string?
          default: |
            #!/bin/bash
            FILE=$0
            BASENAME=$(basename "$FILE")
            get_gene_n_tss.R --isoforms "${FILE}" --gene grouped.genes.tsv --tss grouped.common_tss.tsv
            sed -ibak 's/[[:space:]]\{1,\}[^[:space:]]\{1,\}$//' grouped.genes.tsv
            sed -ibak 's/[[:space:]]\{1,\}[^[:space:]]\{1,\}$//' grouped.common_tss.tsv
            rm -f ./*bak
          inputBinding:
            position: 1
        isoforms_file:
          type: File
          inputBinding:
            position: 5
      outputs:
        genes_file:
          type: File
          outputBinding:
            glob: $(inputs.genes_filename?inputs.genes_filename:"*genes.tsv")
        common_tss_file:
          type: File
          outputBinding:
            glob: $(inputs.common_tss_file?inputs.common_tss_file:"*common_tss.tsv")
      baseCommand: ["bash", "-c"]

  feature_expression_merge:
    run: ../tools/feature-merge.cwl
    in:
      feature_files:
        source:
        - group_transcript_expression/gene_expression_file
        - group_geep_transcript_expression/genes_file
      feature_aliases:
        default:
        - "HTSeq"
        - "GEEP"
      mergeby:
        default: ["RefseqId", "GeneId", "Chrom", "TxStart", "TxEnd", "Strand"]
    out:
    - merged_file
    - stdout_log
    - stderr_log

  get_gene_body:
    run: ../tools/plugin-plot-rna.cwl
    in:
      annotation_file: annotation_file
      bambai_pair: samtools_sort_index_after_dedup/bam_bai_pair
      isoforms_file: group_transcript_expression/transcript_expression_file
      mapped_reads_number: star_aligner/uniquely_mapped_reads_number
      minimum_rpkm: minimum_rpkm
      strand_specificity: strand_specificity
      threads: threads
    out:
    - gene_body_report_file
    - gene_body_plot_pdf
    - rpkm_distribution_plot_pdf


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "QuantSeq 3' FWD, FWD-UMI or REV for single-read mRNA-Seq data"
label: "QuantSeq 3' FWD, FWD-UMI or REV for single-read mRNA-Seq data"
s:alternateName: "Runs QuantSeq 3' FWD, FWD-UMI or REV analysis for single-read mRNA-Seq data"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/trim-quantseq-mrnaseq-se-strand-specific.cwl
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
        s:name: Artem Barski
        s:email: mailto:Artem.Barski@datirum.com
      - class: s:Person
        s:name: Andrey Kartashov
        s:email: mailto:Andrey.Kartashov@datirium.com
        s:sameAs:
          - id: http://orcid.org/0000-0001-9102-5681
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  ### QuantSeq 3' FWD, FWD-UMI or REV for single-read mRNA-Seq data