cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-diffbind:v1.0.0


inputs:

  read_files_cond_1:
    type: File[]
    inputBinding:
      prefix: "-r1"
    doc: "Read files for condition 1. Minimim 2 files in BAM format"

  read_files_cond_2:
    type: File[]
    inputBinding:
      prefix: "-r2"
    doc: "Read files for condition 2. Minimim 2 files in BAM format"

  spikein_read_files_cond_1:
    type: File[]
    inputBinding:
      prefix: "-s1"
    doc: "Spike-in read files for condition 1. Minimim 2 files in BAM format"

  spikein_read_files_cond_2:
    type: File[]
    inputBinding:
      prefix: "-s2"
    doc: "Spike-in read files for condition 2. Minimim 2 files in BAM format"

  peak_files_cond_1:
    type: File[]
    inputBinding:
      prefix: "-p1"
    doc: "Peak files for condition 1. Minimim 2 files in format set with -pf"

  peak_files_cond_2:
    type: File[]
    inputBinding:
      prefix: "-p2"
    doc: "Peak files for condition 2. Minimim 2 files in format set with -pf"

  sample_names_cond_1:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "-n1"
    doc: "Sample names for condition 1. Default: basenames of -r1 without extensions"

  sample_names_cond_2:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "-n2"
    doc: "Sample names for condition 2. Default: basenames of -r2 without extensions"

  blocked_attributes:
    type:
      - "null"
      - string[]
    inputBinding:
      prefix: "-bl"
    doc: |
      Blocking attribute for multi-factor analysis. Minimum 2.
      Either names from --name1 or/and --name2 or array of strings that can be parsed by R to bool.
      In the later case the order and size should correspond to [--read1]+[--read2].
      Default: not applied

  blocked_file:
    type: File?
    inputBinding:
      prefix: "-bf"
    doc: |
      Blocking attribute metadata file for multi-factor analysis. Headerless TSV/CSV file.
      First column - names from --name1 and --name2, second column - group name. --block is ignored

  peakformat:
    type:
      - "null"
      - type: enum
        name: "peakformat"
        symbols: ["raw","bed","narrow","macs","bayes","tpic","sicer","fp4","swembl","csv","report"]
    inputBinding:
      prefix: "-pf"
    doc: "Peak files format. One of [raw, bed, narrow, macs, bayes, tpic, sicer, fp4, swembl, csv, report]. Default: macs"

  name_cond_1:
    type: string?
    inputBinding:
      prefix: "-c1"
    doc: "Condition 1 name, single word with letters and numbers only. Default: condition_1"

  name_cond_2:
    type: string?
    inputBinding:
      prefix: "-c2"
    doc: "Condition 2 name, single word with letters and numbers only. Default: condition_2"

  cutoff_value:
    type: float?
    inputBinding:
      prefix: "-cu"
    doc: "P-value or FDR cutoff for reported results. Default: 0.05"

  cutoff_param:
    type:
      - "null"
      - type: enum
        name: "cutoff"
        symbols: ["pvalue", "fdr"]
    inputBinding:
      prefix: "-cp"
    doc: "Parameter to which cutoff should be applied (fdr or pvalue). Default: fdr"

  analysis_method:
    type:
      - "null"
      - type: enum
        name: "method"
        symbols: ["deseq2", "edger", "all"]
    inputBinding:
      prefix: "-me"
    doc: "Method by which to analyze differential binding affinity. Default: all"

  min_overlap:
    type: int?
    inputBinding:
      prefix: "-mo"
    doc: "Min peakset overlap. Only include peaks in at least this many peaksets when generating consensus peakset. Default: 2"

  rec_summits:
    type: int?
    inputBinding:
      prefix: "--summits"
    doc: |
      Width in bp to extend peaks around their summits in both directions
      and replace the original ones. Set it to 100 bp for ATAC-Seq and 200
      bp for ChIP-Seq datasets. To skip peaks extension and replacement, set
      it to negative value. Default: 200 bp (results in 401 bp wide peaks)

  use_common:
    type: boolean?
    inputBinding:
      prefix: "-uc"
    doc: "Derive consensus peaks only from the common peaks within each condition. Min peakset overlap is ignored. Default: false"

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: "Output prefix. Default: diffbind"

  threads:
    type: int?
    inputBinding:
      prefix: "-th"
    doc: "Threads number. Default: 1"


outputs:

  diff_filtered_report_deseq:
    type: File?
    outputBinding:
      glob: "*_filtered_report_deseq.tsv"
    doc: "Differential binding analysis report for significantly differentially bound sites exported as TSV, DESeq2"

  diff_filtered_report_deseq_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_report_deseq_block.tsv"
    doc: "Differential binding analysis report for significantly differentially bound sites exported as TSV, DESeq2 Blocked"

  diff_filtered_report_edger:
    type: File?
    outputBinding:
      glob: "*_filtered_report_edger.tsv"
    doc: "Differential binding analysis report for significantly differentially bound sites exported as TSV, EdgeR"

  diff_filtered_report_edger_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_report_edger_block.tsv"
    doc: "Differential binding analysis report for significantly differentially bound sites exported as TSV, EdgeR Blocked"

  all_report_deseq:
    type: File?
    outputBinding:
      glob: "*_all_report_deseq.tsv"
    doc: "Not filtered differential binding analysis report exported as TSV, DESeq2"

  all_report_deseq_blocked:
    type: File?
    outputBinding:
      glob: "*_all_report_deseq_block.tsv"
    doc: "Not filtered differential binding analysis report exported as TSV, DESeq2 Blocked"

  all_report_edger:
    type: File?
    outputBinding:
      glob: "*_all_report_edger.tsv"
    doc: "Not filtered differential binding analysis report exported as TSV, EdgeR"

  all_report_edger_blocked:
    type: File?
    outputBinding:
      glob: "*_all_report_edger_block.tsv"
    doc: "Not filtered differential binding analysis report exported as TSV, EdgeR Blocked"

  boxplot_deseq:
    type: File?
    outputBinding:
      glob: "*_filtered_box_plot_deseq.png"
    doc: "Box plots of read distributions for significantly differentially bound sites, DESeq2"

  boxplot_deseq_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_box_plot_deseq.pdf"
    doc: "Box plots of read distributions for significantly differentially bound sites, DESeq2"

  boxplot_deseq_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_box_plot_deseq_block.png"
    doc: "Box plots of read distributions for significantly differentially bound sites, DESeq2 Blocked"

  boxplot_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_box_plot_deseq_block.pdf"
    doc: "Box plots of read distributions for significantly differentially bound sites, DESeq2 Blocked"

  boxplot_edger:
    type: File?
    outputBinding:
      glob: "*_filtered_box_plot_edger.png"
    doc: "Box plots of read distributions for significantly differentially bound sites, EdgeR"
  
  boxplot_edger_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_box_plot_edger.pdf"
    doc: "Box plots of read distributions for significantly differentially bound sites, EdgeR"

  boxplot_edger_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_box_plot_edger_block.png"
    doc: "Box plots of read distributions for significantly differentially bound sites, EdgeR Blocked"

  boxplot_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_box_plot_edger_block.pdf"
    doc: "Box plots of read distributions for significantly differentially bound sites, EdgeR Blocked"

  volcano_plot_deseq:
    type: File?
    outputBinding:
      glob: "*_filtered_volcano_plot_deseq.png"
    doc: "Volcano plot for significantly differentially bound sites, DESeq2"

  volcano_plot_deseq_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_volcano_plot_deseq.pdf"
    doc: "Volcano plot for significantly differentially bound sites, DESeq2"

  volcano_plot_deseq_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_volcano_plot_deseq_block.png"
    doc: "Volcano plot for for significantly differentially bound sites, DESeq2 Blocked"

  volcano_plot_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_volcano_plot_deseq_block.pdf"
    doc: "Volcano plot for for significantly differentially bound sites, DESeq2 Blocked"

  volcano_plot_edger:
    type: File?
    outputBinding:
      glob: "*_filtered_volcano_plot_edger.png"
    doc: "Volcano plot for for significantly differentially bound sites, EdgeR"

  volcano_plot_edger_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_volcano_plot_edger.pdf"
    doc: "Volcano plot for for significantly differentially bound sites, EdgeR"

  volcano_plot_edger_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_volcano_plot_edger_block.png"
    doc: "Volcano plot for for significantly differentially bound sites, EdgeR Blocked"

  volcano_plot_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_volcano_plot_edger_block.pdf"
    doc: "Volcano plot for for significantly differentially bound sites, EdgeR Blocked"

  ma_plot_deseq:
    type: File?
    outputBinding:
      glob: "*_filtered_ma_plot_deseq.png"
    doc: "MA plot for significantly differentially bound sites, DESeq2"

  ma_plot_deseq_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_ma_plot_deseq.pdf"
    doc: "MA plot for significantly differentially bound sites, DESeq2"

  ma_plot_deseq_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_ma_plot_deseq_block.png"
    doc: "MA plot for significantly differentially bound sites, DESeq2 Blocked"

  ma_plot_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_ma_plot_deseq_block.pdf"
    doc: "MA plot for significantly differentially bound sites, DESeq2 Blocked"

  ma_plot_edger:
    type: File?
    outputBinding:
      glob: "*_filtered_ma_plot_edger.png"
    doc: "MA plot for significantly differentially bound sites, EdgeR"

  ma_plot_edger_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_ma_plot_edger.pdf"
    doc: "MA plot for significantly differentially bound sites, EdgeR"
    
  ma_plot_edger_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_ma_plot_edger_block.png"
    doc: "MA plot for significantly differentially bound sites, EdgeR Blocked"

  ma_plot_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_ma_plot_edger_block.pdf"
    doc: "MA plot for significantly differentially bound sites, EdgeR Blocked"

  diff_filtered_pca_plot_deseq:
    type: File?
    outputBinding:
      glob: "*_filtered_pca_plot_deseq.png"
    doc: "PCA plot for significantly differentially bound sites, DESeq2"

  diff_filtered_pca_plot_deseq_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_pca_plot_deseq.pdf"
    doc: "PCA plot for significantly differentially bound sites, DESeq2"

  diff_filtered_pca_plot_deseq_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_pca_plot_deseq_block.png"
    doc: "PCA plot for significantly differentially bound sites, DESeq2 Blocked"
  
  diff_filtered_pca_plot_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_pca_plot_deseq_block.pdf"
    doc: "PCA plot for significantly differentially bound sites, DESeq2 Blocked"

  diff_filtered_pca_plot_edger:
    type: File?
    outputBinding:
      glob: "*_filtered_pca_plot_edger.png"
    doc: "PCA plot for significantly differentially bound sites, EdgeR"

  diff_filtered_pca_plot_edger_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_pca_plot_edger.pdf"
    doc: "PCA plot for significantly differentially bound sites, EdgeR"

  diff_filtered_pca_plot_edger_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_pca_plot_edger_block.png"
    doc: "PCA plot for significantly differentially bound sites, EdgeR Blocked"

  diff_filtered_pca_plot_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_pca_plot_edger_block.pdf"
    doc: "PCA plot for significantly differentially bound sites, EdgeR Blocked"

  all_pca_plot_deseq:
    type: File?
    outputBinding:
      glob: "*_all_pca_plot_deseq.png"
    doc: "PCA plot for not filtered bound sites, DESeq2"

  all_pca_plot_deseq_pdf:
    type: File?
    outputBinding:
      glob: "*_all_pca_plot_deseq.pdf"
    doc: "PCA plot for not filtered bound sites, DESeq2"

  all_pca_plot_deseq_blocked:
    type: File?
    outputBinding:
      glob: "*_all_pca_plot_deseq_block.png"
    doc: "PCA plot for not filtered bound sites, DESeq2 Blocked"

  all_pca_plot_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_all_pca_plot_deseq_block.pdf"
    doc: "PCA plot for not filtered bound sites, DESeq2 Blocked"

  all_pca_plot_edger:
    type: File?
    outputBinding:
      glob: "*_all_pca_plot_edger.png"
    doc: "PCA plot for not filtered bound sites, EdgeR"

  all_pca_plot_edger_pdf:
    type: File?
    outputBinding:
      glob: "*_all_pca_plot_edger.pdf"
    doc: "PCA plot for not filtered bound sites, EdgeR"

  all_pca_plot_edger_blocked:
    type: File?
    outputBinding:
      glob: "*_all_pca_plot_edger_block.png"
    doc: "PCA plot for not filtered bound sites, EdgeR Blocked"

  all_pca_plot_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_all_pca_plot_edger_block.pdf"
    doc: "PCA plot for not filtered bound sites, EdgeR Blocked"

  binding_heatmap_deseq:
    type: File?
    outputBinding:
      glob: "*_filtered_binding_heatmap_deseq.png"
    doc: "Binding heatmap for significantly differentially bound sites, DESeq2"

  binding_heatmap_deseq_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_binding_heatmap_deseq.pdf"
    doc: "Binding heatmap for significantly differentially bound sites, DESeq2"

  binding_heatmap_deseq_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_binding_heatmap_deseq_block.png"
    doc: "Binding heatmap for significantly differentially bound sites, DESeq2 Blocked"

  binding_heatmap_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_binding_heatmap_deseq_block.pdf"
    doc: "Binding heatmap for significantly differentially bound sites, DESeq2 Blocked"

  binding_heatmap_edger:
    type: File?
    outputBinding:
      glob: "*_filtered_binding_heatmap_edger.png"
    doc: "Binding heatmap for significantly differentially bound sites, EdgeR"

  binding_heatmap_edger_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_binding_heatmap_edger.pdf"
    doc: "Binding heatmap for significantly differentially bound sites, EdgeR"

  binding_heatmap_edger_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_binding_heatmap_edger_block.png"
    doc: "Binding heatmap for significantly differentially bound sites, EdgeR Blocked"

  binding_heatmap_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_binding_heatmap_edger_block.pdf"
    doc: "Binding heatmap for significantly differentially bound sites, EdgeR Blocked"

  diff_filtered_norm_counts_corr_heatmap_deseq:
    type: File?
    outputBinding:
      glob: "*_filtered_normalized_counts_correlation_heatmap_deseq.png"
    doc: "Normalized counts correlation heatmap for significantly differentially bound sites, DESeq2"

  diff_filtered_norm_counts_corr_heatmap_deseq_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_normalized_counts_correlation_heatmap_deseq.pdf"
    doc: "Normalized counts correlation heatmap for significantly differentially bound sites, DESeq2"

  diff_filtered_norm_counts_corr_heatmap_deseq_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_normalized_counts_correlation_heatmap_deseq_block.png"
    doc: "Normalized counts correlation heatmap for significantly differentially bound sites, DESeq2 Blocked"

  diff_filtered_norm_counts_corr_heatmap_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_normalized_counts_correlation_heatmap_deseq_block.pdf"
    doc: "Normalized counts correlation heatmap for significantly differentially bound sites, DESeq2 Blocked"

  diff_filtered_norm_counts_corr_heatmap_edger:
    type: File?
    outputBinding:
      glob: "*_filtered_normalized_counts_correlation_heatmap_edger.png"
    doc: "Normalized counts correlation heatmap for significantly differentially bound sites, EdgeR"

  diff_filtered_norm_counts_corr_heatmap_edger_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_normalized_counts_correlation_heatmap_edger.pdf"
    doc: "Normalized counts correlation heatmap for significantly differentially bound sites, EdgeR"

  diff_filtered_norm_counts_corr_heatmap_edger_blocked:
    type: File?
    outputBinding:
      glob: "*_filtered_normalized_counts_correlation_heatmap_edger_block.png"
    doc: "Normalized counts correlation heatmap for significantly differentially bound sites, EdgeR Blocked"

  diff_filtered_norm_counts_corr_heatmap_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_filtered_normalized_counts_correlation_heatmap_edger_block.pdf"
    doc: "Normalized counts correlation heatmap for significantly differentially bound sites, EdgeR Blocked"

  all_norm_counts_corr_heatmap_deseq:
    type: File?
    outputBinding:
      glob: "*_all_normalized_counts_correlation_heatmap_deseq.png"
    doc: "Not filtered normalized counts correlation heatmap, DESeq2"

  all_norm_counts_corr_heatmap_deseq_pdf:
    type: File?
    outputBinding:
      glob: "*_all_normalized_counts_correlation_heatmap_deseq.pdf"
    doc: "Not filtered normalized counts correlation heatmap, DESeq2"

  all_norm_counts_corr_heatmap_deseq_blocked:
    type: File?
    outputBinding:
      glob: "*_all_normalized_counts_correlation_heatmap_deseq_block.png"
    doc: "Not filtered normalized counts correlation heatmap, DESeq2 Blocked"

  all_norm_counts_corr_heatmap_deseq_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_all_normalized_counts_correlation_heatmap_deseq_block.pdf"
    doc: "Not filtered normalized counts correlation heatmap, DESeq2 Blocked"

  all_norm_counts_corr_heatmap_edger:
    type: File?
    outputBinding:
      glob: "*_all_normalized_counts_correlation_heatmap_edger.png"
    doc: "Not filtered normalized counts correlation heatmap, EdgeR"

  all_norm_counts_corr_heatmap_edger_pdf:
    type: File?
    outputBinding:
      glob: "*_all_normalized_counts_correlation_heatmap_edger.pdf"
    doc: "Not filtered normalized counts correlation heatmap, EdgeR"

  all_norm_counts_corr_heatmap_edger_blocked:
    type: File?
    outputBinding:
      glob: "*_all_normalized_counts_correlation_heatmap_edger_block.png"
    doc: "Not filtered normalized counts correlation heatmap, EdgeR Blocked"

  all_norm_counts_corr_heatmap_edger_blocked_pdf:
    type: File?
    outputBinding:
      glob: "*_all_normalized_counts_correlation_heatmap_edger_block.pdf"
    doc: "Not filtered normalized counts correlation heatmap, EdgeR Blocked"

  consensus_peak_venn_diagram:
    type: File?
    outputBinding:
      glob: "*_consensus_peak_venn_diagram.png"
    doc: "Consensus peak Venn Diagram" 

  consensus_peak_venn_diagram_pdf:
    type: File?
    outputBinding:
      glob: "*_consensus_peak_venn_diagram.pdf"
    doc: "Consensus peak Venn Diagram" 

  raw_counts_corr_heatmap:
    type: File?
    outputBinding:
      glob: "*_raw_counts_correlation_heatmap.png"
    doc: "Raw counts correlation heatmap"

  raw_counts_corr_heatmap_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_counts_correlation_heatmap.pdf"
    doc: "Raw counts correlation heatmap"

  peak_overlap_corr_heatmap:
    type: File?
    outputBinding:
      glob: "*_peak_overlap_correlation_heatmap.png"
    doc: "Peak overlap correlation heatmap"

  peak_overlap_corr_heatmap_pdf:
    type: File?
    outputBinding:
      glob: "*_peak_overlap_correlation_heatmap.pdf"
    doc: "Peak overlap correlation heatmap"

  peak_overlap_rate_plot_cond_1:
    type: File?
    outputBinding:
      glob: "*_condition_1_peak_overlap_rate.png"
    doc: "Peak overlap rate plot, condition 1"

  peak_overlap_rate_plot_cond_1_pdf:
    type: File?
    outputBinding:
      glob: "*_condition_1_peak_overlap_rate.pdf"
    doc: "Peak overlap rate plot, condition 1"

  peak_overlap_rate_plot_cond_2:
    type: File?
    outputBinding:
      glob: "*_condition_2_peak_overlap_rate.png"
    doc: "Peak overlap rate plot, condition 2"

  peak_overlap_rate_plot_cond_2_pdf:
    type: File?
    outputBinding:
      glob: "*_condition_2_peak_overlap_rate.pdf"
    doc: "Peak overlap rate plot, condition 2"

  all_peak_overlap_rate_plot:
    type: File?
    outputBinding:
      glob: "*_all_peak_overlap_rate.png"
    doc: "All peak overlap rate plot"

  all_peak_overlap_rate_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_all_peak_overlap_rate.pdf"
    doc: "All peak overlap rate plot"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["run_diffbind-for-spikein.R"]
stdout: diffbind_stdout.log
stderr: diffbind_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "DiffBind spike-in - Differential Binding Analysis of Peak Data"
s:name: "DiffBind spike-in - Differential Binding Analysis of Peak Data"
s:alternateName: "Compute differentially bound sites from multiple ChIP-seq spike-in experiments using affinity (quantitative) data."

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/diffbind.cwl
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

doc: |
  Runs R script to compute differentially bound sites from multiple ChIP-seq spike-in experiments using affinity (quantitative) and occupancy data.

  Selected excerpts on Normalization from diffbind manual (https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf):

    A third method is also provided that avoids making any assumptions regarding the distribution
    of reads between binding sites in different samples that may be specific to RNA-seq analysis
    and inappropriate for ChIP-seq analysis. This method ("lib" or DBA_NORM_LIB) is based on
    the different library sizes for each sample, computing normalization factors to weight each
    sample so as to appear to have the same library size. For DESeq2 , this is accomplished
    by dividing the number of reads in each library by the mean library size. For edgeR , the
    normalization factors are all set to 1.0, indicating that only library-size normalization should
    occur.

    ...

    In the previous section, we saw how a library-size normalization based on the full library
    size, rather than on the reads in the enriched areas designated by the consensus peaks set,
    enabled an analysis that avoided making assumption about the changes in binding patterns
    (specifically, assumptions that binding changes are roughly balanced). Another way to do this,
    while gaining some of the benefits of the native normalization methods, is to use background
    normalization. Given the difficulty in differentiating between technical biases and biological
    signal when attempting to normalize ChIP-seq data, the idea is to normalize based on a more
    neutral sample of reads.

    The core background normalization technique is to divide the genome into large bins and
    count overlapping reads. As the enrichment expected in ChIP-seq (and ATAC-seq) is expected to
    occur over relatively narrow intervals (roughly between 100bp and 600bp), it is
    expected that there should not be systematic differences in signals over much larger intervals
    (on the order of 10,000bp and greater). Any differences seen should be technical rather than
    biological, so it is safer to normalize based these differences.

    Note also that this type of background normalization is one of the methods recommended
    for ATAC-seq analysis.

    By specifying background=TRUE in dba.normalize, all chromosomes that contains peaks in the
    consensus set are divided into non-overlapping bins (default size 15,000bp) and overlapping
    reads counted. These can then be normalized using any of the normalization methods. While
    generally doing a library size based normalization will yield the same result as on the full library
    size, the native normalization methods can be applied to the background bins.

    ...

    An alternative for avoiding the use of the consensus count matrix
            (A matrix of offsets can be supplied via the dba.normalize function,
             or an offset matrix can be calculated using a loess fit)
                                                                      when normalizing ChIP-seq
    data is to use spike-in data, where exogenous chromatin (usually from Drosophila melanogaster)
    is "spiked in" to the ChIP. If the amount of spiked-in chromatin can be precisely controlled,
    then we can use the relative amounts of reads that map to the alternative reference genome
    for each sample. Counting these reads then forms a background we can use
    in the same manner as background normalization discussed previously (either a library size
    adjustment using the total number of exogenous reads for each sample, or a native RNA-seq
    method (TMM or RLE) taking into account how those reads are distributed).

    And from one of the dev of diffbind, Rory Stark:

        RLE is a good choice when you expect most of the peaks not to be different between the
        conditions, which is exactly the assumption for spike-ins, where differences in the raw
        data should be due to IP efficiency and sequencing depth.

    So for this SciDAP workflow, we use 'spikes <- dba.normalize(spikes, normalize="RLE", background=TRUE)'


s:about: |
  usage: run_diffbind.R [-h] -r1 READ1 [READ1 ...] -r2 READ2 [READ2 ...] -p1
                        PEAK1 [PEAK1 ...] -p2 PEAK2 [PEAK2 ...]
                        [-n1 [NAME1 ...]] [-n2 [NAME2 ...]] [-bl [BLOCK ...]]
                        [-bf BLOCKFILE]
                        [-pf {raw,bed,narrow,macs,bayes,tpic,sicer,fp4,swembl,csv,report}]
                        [-c1 CONDITION1] [-c2 CONDITION2]
                        [-me {edger,deseq2,all}] [-mo MINOVERLAP] [-uc]
                        [--summits SUMMITS] [-cu CUTOFF] [-cp {pvalue,fdr}]
                        [-co {Reds,Greens,Blues,Greys,YlOrRd,Oranges}]
                        [-th THREADS] [-pa PADDING] [-o OUTPUT]

  Differential binding analysis of ChIP-Seq experiments using affinity (read
  count) data

  options:
    -h, --help            show this help message and exit
    -r1 READ1 [READ1 ...], --read1 READ1 [READ1 ...]
                          Read files for condition 1. Minimim 2 files in BAM
                          format
    -r2 READ2 [READ2 ...], --read2 READ2 [READ2 ...]
                          Read files for condition 2. Minimim 2 files in BAM
                          format
    -s1 SPIKE1 [SPIKE1 ...], --spike1 SPIKE1 [SPIKE1 ...]
                          Spike-in read files for condition 1. Minimim 2 files
                          in BAM format
    -s2 SPIKE2 [SPIKE2 ...], --spike2 SPIKE2 [SPIKE2 ...]
                          Spike-in read files for condition 2. Minimim 2 files
                          in BAM format
    -p1 PEAK1 [PEAK1 ...], --peak1 PEAK1 [PEAK1 ...]
                          Peak files for condition 1. Minimim 2 files in format
                          set with -pf
    -p2 PEAK2 [PEAK2 ...], --peak2 PEAK2 [PEAK2 ...]
                          Peak files for condition 2. Minimim 2 files in format
                          set with -pf
    -n1 [NAME1 ...], --name1 [NAME1 ...]
                          Sample names for condition 1. Default: basenames of
                          -r1 without extensions
    -n2 [NAME2 ...], --name2 [NAME2 ...]
                          Sample names for condition 2. Default: basenames of
                          -r2 without extensions
    -bl [BLOCK ...], --block [BLOCK ...]
                          Blocking attribute for multi-factor analysis. Minimum
                          2. Either names from --name1 or/and --name2 or array
                          of bool based on [read1]+[read2]. Default: not applied
    -bf BLOCKFILE, --blockfile BLOCKFILE
                          Blocking attribute metadata file for multi-factor
                          analysis. Headerless TSV/CSV file. First column -
                          names from --name1 and --name2, second column - group
                          name. --block is ignored
    -pf {raw,bed,narrow,macs,bayes,tpic,sicer,fp4,swembl,csv,report}, --peakformat {raw,bed,narrow,macs,bayes,tpic,sicer,fp4,swembl,csv,report}
                          Peak files format. One of [raw, bed, narrow, macs,
                          bayes, tpic, sicer, fp4, swembl, csv, report].
                          Default: macs
    -c1 CONDITION1, --condition1 CONDITION1
                          Condition 1 name, single word with letters and numbers
                          only. Default: condition_1
    -c2 CONDITION2, --condition2 CONDITION2
                          Condition 2 name, single word with letters and numbers
                          only. Default: condition_2
    -me {edger,deseq2,all}, --method {edger,deseq2,all}
                          Method by which to analyze differential binding
                          affinity. Default: all
    -mo MINOVERLAP, --minoverlap MINOVERLAP
                          Min peakset overlap. Only include peaks in at least
                          this many peaksets when generating consensus peakset.
                          Default: 2
    -uc, --usecommon      Derive consensus peaks only from the common peaks
                          within each condition. Min peakset overlap and min
                          read counts are ignored. Default: false
    --summits SUMMITS     Width in bp to extend peaks around their summits in
                          both directions and replace the original ones. Set it
                          to 100 bp for ATAC-Seq and 200 bp for ChIP-Seq
                          datasets. To skip peaks extension and replacement, set
                          it to negative value. Default: 200 bp (results in 401
                          bp wide peaks)
    -cu CUTOFF, --cutoff CUTOFF
                          Cutoff for reported results. Applied to the parameter
                          set with -cp. Default: 0.05
    -cp {pvalue,fdr}, --cparam {pvalue,fdr}
                          Parameter to which cutoff should be applied (fdr or
                          pvalue). Default: fdr
    -co {Reds,Greens,Blues,Greys,YlOrRd,Oranges}, --color {Reds,Greens,Blues,Greys,YlOrRd,Oranges}
                          Color scheme. Default: Greens
    -th THREADS, --threads THREADS
                          Threads to use
    -pa PADDING, --padding PADDING
                          Padding for generated heatmaps. Default: 20
    -o OUTPUT, --output OUTPUT
                          Output prefix. Default: diffbind