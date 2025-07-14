cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: DockerRequirement
  dockerPull: biowardrobe2/diffbind:v0.0.16
inputs:
  read_files_cond_1:
    type: File[]
    inputBinding:
      prefix: -r1
    doc: Read files for condition 1. Minimim 2 files in BAM format
  read_files_cond_2:
    type: File[]
    inputBinding:
      prefix: -r2
    doc: Read files for condition 2. Minimim 2 files in BAM format
  peak_files_cond_1:
    type: File[]
    inputBinding:
      prefix: -p1
    doc: Peak files for condition 1. Minimim 2 files in format set with -pf
  peak_files_cond_2:
    type: File[]
    inputBinding:
      prefix: -p2
    doc: Peak files for condition 2. Minimim 2 files in format set with -pf
  sample_names_cond_1:
    type:
    - 'null'
    - string[]
    inputBinding:
      prefix: -n1
    doc: 'Sample names for condition 1. Default: basenames of -r1 without extensions'
  sample_names_cond_2:
    type:
    - 'null'
    - string[]
    inputBinding:
      prefix: -n2
    doc: 'Sample names for condition 2. Default: basenames of -r2 without extensions'
  blocked_attributes:
    type:
    - 'null'
    - string[]
    inputBinding:
      prefix: -bl
    doc: |
      Blocking attribute for multi-factor analysis. Minimum 2.
      Either names from --name1 or/and --name2 or array of strings that can be parsed by R to bool.
      In the later case the order and size should correspond to [--read1]+[--read2].
      Default: not applied
  blocked_file:
    type: File?
    inputBinding:
      prefix: -bf
    doc: |
      Blocking attribute metadata file for multi-factor analysis. Headerless TSV/CSV file.
      First column - names from --name1 and --name2, second column - group name. --block is ignored
  peakformat:
    type:
    - 'null'
    - type: enum
      name: peakformat
      symbols:
      - raw
      - bed
      - narrow
      - macs
      - bayes
      - tpic
      - sicer
      - fp4
      - swembl
      - csv
      - report
    inputBinding:
      prefix: -pf
    doc: 'Peak files format. One of [raw, bed, narrow, macs, bayes, tpic, sicer, fp4, swembl, csv, report]. Default: macs'
  name_cond_1:
    type: string?
    inputBinding:
      prefix: -c1
    doc: 'Condition 1 name, single word with letters and numbers only. Default: condition_1'
  name_cond_2:
    type: string?
    inputBinding:
      prefix: -c2
    doc: 'Condition 2 name, single word with letters and numbers only. Default: condition_2'
  cutoff_value:
    type: float?
    inputBinding:
      prefix: -cu
    doc: 'P-value or FDR cutoff for reported results. Default: 0.05'
  cutoff_param:
    type:
    - 'null'
    - type: enum
      name: cutoff
      symbols:
      - pvalue
      - fdr
    inputBinding:
      prefix: -cp
    doc: 'Parameter to which cutoff should be applied (fdr or pvalue). Default: fdr'
  analysis_method:
    type:
    - 'null'
    - type: enum
      name: method
      symbols:
      - deseq2
      - edger
      - all
    inputBinding:
      prefix: -me
    doc: 'Method by which to analyze differential binding affinity. Default: all'
  min_overlap:
    type: int?
    inputBinding:
      prefix: -mo
    doc: 'Min peakset overlap. Only include peaks in at least this many peaksets when generating consensus peakset. Default: 2'
  rec_summits:
    type: int?
    inputBinding:
      prefix: --summits
    doc: |
      Width in bp to extend peaks around their summits in both directions
      and replace the original ones. Set it to 100 bp for ATAC-Seq and 200
      bp for ChIP-Seq datasets. To skip peaks extension and replacement, set
      it to negative value. Default: 200 bp (results in 401 bp wide peaks)
  use_common:
    type: boolean?
    inputBinding:
      prefix: -uc
    doc: 'Derive consensus peaks only from the common peaks within each condition. Min peakset overlap is ignored. Default: false'
  output_prefix:
    type: string?
    inputBinding:
      prefix: --output
    doc: 'Output prefix. Default: diffbind'
  threads:
    type: int?
    inputBinding:
      prefix: -th
    doc: 'Threads number. Default: 1'
outputs:
  diff_filtered_report_deseq:
    type: File?
    outputBinding:
      glob: '*_filtered_report_deseq.tsv'
    doc: Differential binding analysis report for significantly differentially bound sites exported as TSV, DESeq2
  diff_filtered_report_deseq_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_report_deseq_block.tsv'
    doc: Differential binding analysis report for significantly differentially bound sites exported as TSV, DESeq2 Blocked
  diff_filtered_report_edger:
    type: File?
    outputBinding:
      glob: '*_filtered_report_edger.tsv'
    doc: Differential binding analysis report for significantly differentially bound sites exported as TSV, EdgeR
  diff_filtered_report_edger_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_report_edger_block.tsv'
    doc: Differential binding analysis report for significantly differentially bound sites exported as TSV, EdgeR Blocked
  all_report_deseq:
    type: File?
    outputBinding:
      glob: '*_all_report_deseq.tsv'
    doc: Not filtered differential binding analysis report exported as TSV, DESeq2
  all_report_deseq_blocked:
    type: File?
    outputBinding:
      glob: '*_all_report_deseq_block.tsv'
    doc: Not filtered differential binding analysis report exported as TSV, DESeq2 Blocked
  all_report_edger:
    type: File?
    outputBinding:
      glob: '*_all_report_edger.tsv'
    doc: Not filtered differential binding analysis report exported as TSV, EdgeR
  all_report_edger_blocked:
    type: File?
    outputBinding:
      glob: '*_all_report_edger_block.tsv'
    doc: Not filtered differential binding analysis report exported as TSV, EdgeR Blocked
  boxplot_deseq:
    type: File?
    outputBinding:
      glob: '*_filtered_box_plot_deseq.png'
    doc: Box plots of read distributions for significantly differentially bound sites, DESeq2
  boxplot_deseq_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_box_plot_deseq.pdf'
    doc: Box plots of read distributions for significantly differentially bound sites, DESeq2
  boxplot_deseq_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_box_plot_deseq_block.png'
    doc: Box plots of read distributions for significantly differentially bound sites, DESeq2 Blocked
  boxplot_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_box_plot_deseq_block.pdf'
    doc: Box plots of read distributions for significantly differentially bound sites, DESeq2 Blocked
  boxplot_edger:
    type: File?
    outputBinding:
      glob: '*_filtered_box_plot_edger.png'
    doc: Box plots of read distributions for significantly differentially bound sites, EdgeR
  boxplot_edger_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_box_plot_edger.pdf'
    doc: Box plots of read distributions for significantly differentially bound sites, EdgeR
  boxplot_edger_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_box_plot_edger_block.png'
    doc: Box plots of read distributions for significantly differentially bound sites, EdgeR Blocked
  boxplot_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_box_plot_edger_block.pdf'
    doc: Box plots of read distributions for significantly differentially bound sites, EdgeR Blocked
  volcano_plot_deseq:
    type: File?
    outputBinding:
      glob: '*_filtered_volcano_plot_deseq.png'
    doc: Volcano plot for significantly differentially bound sites, DESeq2
  volcano_plot_deseq_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_volcano_plot_deseq.pdf'
    doc: Volcano plot for significantly differentially bound sites, DESeq2
  volcano_plot_deseq_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_volcano_plot_deseq_block.png'
    doc: Volcano plot for for significantly differentially bound sites, DESeq2 Blocked
  volcano_plot_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_volcano_plot_deseq_block.pdf'
    doc: Volcano plot for for significantly differentially bound sites, DESeq2 Blocked
  volcano_plot_edger:
    type: File?
    outputBinding:
      glob: '*_filtered_volcano_plot_edger.png'
    doc: Volcano plot for for significantly differentially bound sites, EdgeR
  volcano_plot_edger_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_volcano_plot_edger.pdf'
    doc: Volcano plot for for significantly differentially bound sites, EdgeR
  volcano_plot_edger_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_volcano_plot_edger_block.png'
    doc: Volcano plot for for significantly differentially bound sites, EdgeR Blocked
  volcano_plot_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_volcano_plot_edger_block.pdf'
    doc: Volcano plot for for significantly differentially bound sites, EdgeR Blocked
  ma_plot_deseq:
    type: File?
    outputBinding:
      glob: '*_filtered_ma_plot_deseq.png'
    doc: MA plot for significantly differentially bound sites, DESeq2
  ma_plot_deseq_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_ma_plot_deseq.pdf'
    doc: MA plot for significantly differentially bound sites, DESeq2
  ma_plot_deseq_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_ma_plot_deseq_block.png'
    doc: MA plot for significantly differentially bound sites, DESeq2 Blocked
  ma_plot_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_ma_plot_deseq_block.pdf'
    doc: MA plot for significantly differentially bound sites, DESeq2 Blocked
  ma_plot_edger:
    type: File?
    outputBinding:
      glob: '*_filtered_ma_plot_edger.png'
    doc: MA plot for significantly differentially bound sites, EdgeR
  ma_plot_edger_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_ma_plot_edger.pdf'
    doc: MA plot for significantly differentially bound sites, EdgeR
  ma_plot_edger_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_ma_plot_edger_block.png'
    doc: MA plot for significantly differentially bound sites, EdgeR Blocked
  ma_plot_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_ma_plot_edger_block.pdf'
    doc: MA plot for significantly differentially bound sites, EdgeR Blocked
  diff_filtered_pca_plot_deseq:
    type: File?
    outputBinding:
      glob: '*_filtered_pca_plot_deseq.png'
    doc: PCA plot for significantly differentially bound sites, DESeq2
  diff_filtered_pca_plot_deseq_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_pca_plot_deseq.pdf'
    doc: PCA plot for significantly differentially bound sites, DESeq2
  diff_filtered_pca_plot_deseq_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_pca_plot_deseq_block.png'
    doc: PCA plot for significantly differentially bound sites, DESeq2 Blocked
  diff_filtered_pca_plot_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_pca_plot_deseq_block.pdf'
    doc: PCA plot for significantly differentially bound sites, DESeq2 Blocked
  diff_filtered_pca_plot_edger:
    type: File?
    outputBinding:
      glob: '*_filtered_pca_plot_edger.png'
    doc: PCA plot for significantly differentially bound sites, EdgeR
  diff_filtered_pca_plot_edger_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_pca_plot_edger.pdf'
    doc: PCA plot for significantly differentially bound sites, EdgeR
  diff_filtered_pca_plot_edger_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_pca_plot_edger_block.png'
    doc: PCA plot for significantly differentially bound sites, EdgeR Blocked
  diff_filtered_pca_plot_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_pca_plot_edger_block.pdf'
    doc: PCA plot for significantly differentially bound sites, EdgeR Blocked
  all_pca_plot_deseq:
    type: File?
    outputBinding:
      glob: '*_all_pca_plot_deseq.png'
    doc: PCA plot for not filtered bound sites, DESeq2
  all_pca_plot_deseq_pdf:
    type: File?
    outputBinding:
      glob: '*_all_pca_plot_deseq.pdf'
    doc: PCA plot for not filtered bound sites, DESeq2
  all_pca_plot_deseq_blocked:
    type: File?
    outputBinding:
      glob: '*_all_pca_plot_deseq_block.png'
    doc: PCA plot for not filtered bound sites, DESeq2 Blocked
  all_pca_plot_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_all_pca_plot_deseq_block.pdf'
    doc: PCA plot for not filtered bound sites, DESeq2 Blocked
  all_pca_plot_edger:
    type: File?
    outputBinding:
      glob: '*_all_pca_plot_edger.png'
    doc: PCA plot for not filtered bound sites, EdgeR
  all_pca_plot_edger_pdf:
    type: File?
    outputBinding:
      glob: '*_all_pca_plot_edger.pdf'
    doc: PCA plot for not filtered bound sites, EdgeR
  all_pca_plot_edger_blocked:
    type: File?
    outputBinding:
      glob: '*_all_pca_plot_edger_block.png'
    doc: PCA plot for not filtered bound sites, EdgeR Blocked
  all_pca_plot_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_all_pca_plot_edger_block.pdf'
    doc: PCA plot for not filtered bound sites, EdgeR Blocked
  binding_heatmap_deseq:
    type: File?
    outputBinding:
      glob: '*_filtered_binding_heatmap_deseq.png'
    doc: Binding heatmap for significantly differentially bound sites, DESeq2
  binding_heatmap_deseq_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_binding_heatmap_deseq.pdf'
    doc: Binding heatmap for significantly differentially bound sites, DESeq2
  binding_heatmap_deseq_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_binding_heatmap_deseq_block.png'
    doc: Binding heatmap for significantly differentially bound sites, DESeq2 Blocked
  binding_heatmap_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_binding_heatmap_deseq_block.pdf'
    doc: Binding heatmap for significantly differentially bound sites, DESeq2 Blocked
  binding_heatmap_edger:
    type: File?
    outputBinding:
      glob: '*_filtered_binding_heatmap_edger.png'
    doc: Binding heatmap for significantly differentially bound sites, EdgeR
  binding_heatmap_edger_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_binding_heatmap_edger.pdf'
    doc: Binding heatmap for significantly differentially bound sites, EdgeR
  binding_heatmap_edger_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_binding_heatmap_edger_block.png'
    doc: Binding heatmap for significantly differentially bound sites, EdgeR Blocked
  binding_heatmap_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_binding_heatmap_edger_block.pdf'
    doc: Binding heatmap for significantly differentially bound sites, EdgeR Blocked
  diff_filtered_norm_counts_corr_heatmap_deseq:
    type: File?
    outputBinding:
      glob: '*_filtered_normalized_counts_correlation_heatmap_deseq.png'
    doc: Normalized counts correlation heatmap for significantly differentially bound sites, DESeq2
  diff_filtered_norm_counts_corr_heatmap_deseq_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_normalized_counts_correlation_heatmap_deseq.pdf'
    doc: Normalized counts correlation heatmap for significantly differentially bound sites, DESeq2
  diff_filtered_norm_counts_corr_heatmap_deseq_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_normalized_counts_correlation_heatmap_deseq_block.png'
    doc: Normalized counts correlation heatmap for significantly differentially bound sites, DESeq2 Blocked
  diff_filtered_norm_counts_corr_heatmap_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_normalized_counts_correlation_heatmap_deseq_block.pdf'
    doc: Normalized counts correlation heatmap for significantly differentially bound sites, DESeq2 Blocked
  diff_filtered_norm_counts_corr_heatmap_edger:
    type: File?
    outputBinding:
      glob: '*_filtered_normalized_counts_correlation_heatmap_edger.png'
    doc: Normalized counts correlation heatmap for significantly differentially bound sites, EdgeR
  diff_filtered_norm_counts_corr_heatmap_edger_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_normalized_counts_correlation_heatmap_edger.pdf'
    doc: Normalized counts correlation heatmap for significantly differentially bound sites, EdgeR
  diff_filtered_norm_counts_corr_heatmap_edger_blocked:
    type: File?
    outputBinding:
      glob: '*_filtered_normalized_counts_correlation_heatmap_edger_block.png'
    doc: Normalized counts correlation heatmap for significantly differentially bound sites, EdgeR Blocked
  diff_filtered_norm_counts_corr_heatmap_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_filtered_normalized_counts_correlation_heatmap_edger_block.pdf'
    doc: Normalized counts correlation heatmap for significantly differentially bound sites, EdgeR Blocked
  all_norm_counts_corr_heatmap_deseq:
    type: File?
    outputBinding:
      glob: '*_all_normalized_counts_correlation_heatmap_deseq.png'
    doc: Not filtered normalized counts correlation heatmap, DESeq2
  all_norm_counts_corr_heatmap_deseq_pdf:
    type: File?
    outputBinding:
      glob: '*_all_normalized_counts_correlation_heatmap_deseq.pdf'
    doc: Not filtered normalized counts correlation heatmap, DESeq2
  all_norm_counts_corr_heatmap_deseq_blocked:
    type: File?
    outputBinding:
      glob: '*_all_normalized_counts_correlation_heatmap_deseq_block.png'
    doc: Not filtered normalized counts correlation heatmap, DESeq2 Blocked
  all_norm_counts_corr_heatmap_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_all_normalized_counts_correlation_heatmap_deseq_block.pdf'
    doc: Not filtered normalized counts correlation heatmap, DESeq2 Blocked
  all_norm_counts_corr_heatmap_edger:
    type: File?
    outputBinding:
      glob: '*_all_normalized_counts_correlation_heatmap_edger.png'
    doc: Not filtered normalized counts correlation heatmap, EdgeR
  all_norm_counts_corr_heatmap_edger_pdf:
    type: File?
    outputBinding:
      glob: '*_all_normalized_counts_correlation_heatmap_edger.pdf'
    doc: Not filtered normalized counts correlation heatmap, EdgeR
  all_norm_counts_corr_heatmap_edger_blocked:
    type: File?
    outputBinding:
      glob: '*_all_normalized_counts_correlation_heatmap_edger_block.png'
    doc: Not filtered normalized counts correlation heatmap, EdgeR Blocked
  all_norm_counts_corr_heatmap_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: '*_all_normalized_counts_correlation_heatmap_edger_block.pdf'
    doc: Not filtered normalized counts correlation heatmap, EdgeR Blocked
  consensus_peak_venn_diagram:
    type: File?
    outputBinding:
      glob: '*_consensus_peak_venn_diagram.png'
    doc: Consensus peak Venn Diagram
  consensus_peak_venn_diagram_pdf:
    type: File?
    outputBinding:
      glob: '*_consensus_peak_venn_diagram.pdf'
    doc: Consensus peak Venn Diagram
  raw_counts_corr_heatmap:
    type: File?
    outputBinding:
      glob: '*_raw_counts_correlation_heatmap.png'
    doc: Raw counts correlation heatmap
  raw_counts_corr_heatmap_pdf:
    type: File?
    outputBinding:
      glob: '*_raw_counts_correlation_heatmap.pdf'
    doc: Raw counts correlation heatmap
  peak_overlap_corr_heatmap:
    type: File?
    outputBinding:
      glob: '*_peak_overlap_correlation_heatmap.png'
    doc: Peak overlap correlation heatmap
  peak_overlap_corr_heatmap_pdf:
    type: File?
    outputBinding:
      glob: '*_peak_overlap_correlation_heatmap.pdf'
    doc: Peak overlap correlation heatmap
  peak_overlap_rate_plot_cond_1:
    type: File?
    outputBinding:
      glob: '*_condition_1_peak_overlap_rate.png'
    doc: Peak overlap rate plot, condition 1
  peak_overlap_rate_plot_cond_1_pdf:
    type: File?
    outputBinding:
      glob: '*_condition_1_peak_overlap_rate.pdf'
    doc: Peak overlap rate plot, condition 1
  peak_overlap_rate_plot_cond_2:
    type: File?
    outputBinding:
      glob: '*_condition_2_peak_overlap_rate.png'
    doc: Peak overlap rate plot, condition 2
  peak_overlap_rate_plot_cond_2_pdf:
    type: File?
    outputBinding:
      glob: '*_condition_2_peak_overlap_rate.pdf'
    doc: Peak overlap rate plot, condition 2
  all_peak_overlap_rate_plot:
    type: File?
    outputBinding:
      glob: '*_all_peak_overlap_rate.png'
    doc: All peak overlap rate plot
  all_peak_overlap_rate_plot_pdf:
    type: File?
    outputBinding:
      glob: '*_all_peak_overlap_rate.pdf'
    doc: All peak overlap rate plot
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- run_diffbind.R
stdout: diffbind_stdout.log
stderr: diffbind_stderr.log
label: DiffBind - Differential Binding Analysis of ChIP-Seq Peak Data
doc: |
  Runs R script to compute differentially bound sites from multiple ChIP-seq experiments using affinity (quantitative) and occupancy data.
