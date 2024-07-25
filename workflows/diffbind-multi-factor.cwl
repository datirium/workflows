cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var split_by_common_delim = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          var splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };


'sd:upstream':
  dna_experiment:
  - "chipseq-se.cwl"
  - "chipseq-pe.cwl"
  - "trim-chipseq-se.cwl"
  - "trim-chipseq-pe.cwl"
  - "trim-atacseq-se.cwl"
  - "trim-atacseq-pe.cwl"
  genome_indices:
  - "genome-indices.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name"
    sd:preview:
      position: 1

  alignment_files:
    type: File[]
    secondaryFiles:
    - .bai
    label: "ChIP-Seq/ATAC-Seq experiments"
    doc: |
      Sorted and indexed alignment files
      in BAM format
    'sd:upstreamSource': "dna_experiment/bambai_pair"
    'sd:localLabel': true

  peak_files:
    type: File[]
    label: "ChIP-Seq/ATAC-Seq experiments"
    doc:
      Peak files in the MACS2 xls format
    'sd:upstreamSource': "dna_experiment/macs2_called_peaks"
    'sd:localLabel': true

  dataset_names:
    type: string[]
    label: "ChIP-Seq/ATAC-Seq experiments"
    doc: |
      Unique names for samples
    'sd:upstreamSource': "dna_experiment/alias"
    'sd:localLabel': true

  genome_coverage_files:
    type: File[]
    label: "ChIP-Seq/ATAC-Seq experiments"
    doc: |
      Genome coverage files in bigWig format
    'sd:upstreamSource': "dna_experiment/bigwig"
    'sd:localLabel': true

  narrow_peak_files:
    type:
    - "null"
    - File[]
    label: "ChIP-Seq/ATAC-Seq experiments"
    doc: |
      Called peaks files in narrowPeak format
    'sd:upstreamSource': "dna_experiment/macs2_narrow_peaks"
    'sd:localLabel': true

  broad_peak_files:
    type:
    - "null"
    - File[]
    label: "ChIP-Seq/ATAC-Seq experiments"
    doc: |
      Called peaks files in broadPeak format
    'sd:upstreamSource': "dna_experiment/macs2_broad_peaks"
    'sd:localLabel': true

  annotation_file:
    type: File
    label: "Reference genome"
    doc: |
      Reference genome
    'sd:upstreamSource': "genome_indices/annotation"
    'sd:localLabel': true

  chrom_length_file:
    type: File
    label: "Reference genome"
    doc: |
      Reference genome
    'sd:upstreamSource': "genome_indices/chrom_length"
    'sd:localLabel': true

  metadata_file:
    type: File
    label: "Diff. analysis. Metadata file to describe samples categories"
    doc: |
      Metadata file in TSV/CSV format to describe
      input samples categories. First column should
      have the name 'sample', all other columns names
      should be selected from the following list:
      Tissue, Factor, Condition, Treatment, Replicate.
      The values from the 'sample' column should
      correspond to the names of the selected
      ChIP-Seq/ATAC-Seq experiments. Values defined in
      each metadata column should not be used in any of
      the other columns. All metadata columns are treated
      as factors (no covariates are supported).

  groupby:
    type: string?
    default: null
    label: "Diff. analysis. Metadata column(s) that should be used for samples grouping"
    doc: |
      Column(s) from the metadata table to define
      samples grouping. Minimum peakset overlap
      threshold will be applied within each group
      independently. Union of the resulted peaks
      from each of the groups will be used in the
      differential analysis. If not provided,
      minimum peakset overlap filtering threshold
      will be applied for all samples jointly.
      For grouping by multiple columns provide
      space separated values, for example,
      'Treatment Tissue'

  design_formula:
    type: string
    label: "Diff. analysis. Design formula"
    doc: |
      Design formula comprised of the metadata
      columns names. For example, to model the
      effect of Treatment, Tissue, and their
      interaction use
      ~Treatment%2BTissue%2BTreatment%3ATissue

  base_levels:
    type: string?
    default: null
    label: "Diff. analysis. Base levels (optional)"
    doc: |
      Base levels for each of the metadata columns.
      Number and order of the provided values should
      correspond to the metadata columns. If not
      provided, the defauls base levels will be
      defined alphabetically.

  contrast:
    type: string?
    default: null
    label: "Diff. analysis. Contrast for calculating log2 fold changes"
    doc: |
      Contrast applied to the analysis results when
      calculating log2 fold changes. It should be
      formatted as a mathematical formula of values
      present in the metadata table. If not provided,
      the last term from the design formula
      will be used.

  padj_threshold:
    type: float?
    default: 0.05
    label: "Diff. analysis. Maximum allowed adjusted P-value for differentially bound sites"
    doc: |
      Filtering threshold to report only differentially
      bound sites with adjusted P-value less than or
      equal to the provided value.

  scoreby:
    type:
    - "null"
    - type: enum
      symbols:
      - "pvalue"
      - "qvalue"
    default: "pvalue"
    label: "Peak selection. Score metrics to exclude low quality peaks"
    doc: |
      Score metrics to build peak overlap correlation
      heatmap and exclude low quality peaks based on
      the specific threshold value

  score_threshold:
    type: float?
    default: 0.05
    label: "Peak selection. Maximum allowed peak score (pvalue/qvalue)"
    doc: |
      Filtering threshold to keep only those peaks
      where the selected metric is less than or equal
      to the provided value

  overlap_threshold:
    type: float?
    default: 2
    label: "Peak selection. Minimum peakset overlap threshold"
    doc: |
      Filtering threshold to keep only those peaks
      that are present in at least this many samples
      when generating consensus set of peaks used in
      differential analysis. If this threshold has a
      value between zero and one, only those peaks
      will be included that are present in at least
      this proportion of samples. If input samples
      are grouped by the certain metadata columns,
      minimum peakset overlap threshold will be first
      applied per group, then union of the resulted
      peaks will be used in the differential analysis.

  rpkm_threshold:
    type: float?
    default: 1
    label: "Peak selection. Minimum allowed RPKM for consensus peaks"
    doc: |
      Filtering threshold to keep only those consensus
      peaks where the maximum RPKM for all samples is
      bigger than or equal to the provided value.

  rec_summits:
    type: int?
    default: 200
    label: "Width in bp to extend peaks around summits"
    doc: |
      Width in bp to extend peaks around their summits
      in both directions and replace the original ones.
      Set it to 100 bp for ATAC-Seq and 200 bp for
      ChIP-Seq datasets. To skip peaks extension and
      replacement, set it to negative value.
      Default: 200 bp (results in 401 bp wide peaks)
    'sd:layout':
      advanced: true

  promoter_dist:
    type: int?
    default: 1000
    label: "Peak annotation. Promoter distance, bp"
    doc: |
      Maximum distance from gene TSS (in both
      direction) overlapping which the peak will
      be assigned to the promoter region.
    'sd:layout':
      advanced: true

  upstream_dist:
    type: int?
    default: 20000
    label: "Peak annotation. Upstream distance, bp"
    doc: |
      Maximum distance from the promoter (only in
      upstream direction) overlapping which the peak
      will be assigned to the upstream region.
    'sd:layout':
      advanced: true

  cluster_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "row"
      - "column"
      - "both"
      - "none"
    default: "none"
    label: "Peak clustering. Clustering method"
    doc: |
      Hierarchical clustering method to be run
      on normalized read counts
    'sd:layout':
      advanced: true

  row_distance:
    type:
    - "null"
    - type: enum
      symbols:
      - "cosangle"
      - "abscosangle"
      - "euclid"
      - "abseuclid"
      - "cor"
      - "abscor"
    default: "cosangle"
    label: "Peak clustering. Distance metric for row clustering"
    doc: |
      Distance metric for hierarchical row clustering
    'sd:layout':
      advanced: true

  column_distance:
    type:
    - "null"
    - type: enum
      symbols:
      - "cosangle"
      - "abscosangle"
      - "euclid"
      - "abseuclid"
      - "cor"
      - "abscor"
    default: "euclid"
    label: "Peak clustering. Distance metric for column clustering"
    doc: |
      Distance metric for hierarchical
      column clustering
    'sd:layout':
      advanced: true

  threads:
    type:
    - "null"
    - type: enum
      symbols:
      - "1"
      - "2"
      - "3"
      - "4"
    default: "1"
    label: "Number of cores/cpus to use"
    doc: |
      Number of cores/cpus to use
    'sd:layout':
      advanced: true


outputs:

  gc_files:
    type: File[]
    label: "Genome coverage"
    doc: |
      Genome coverage files in bigWig format
    outputSource: pipe/gc_files
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "Genome coverage"
        height: 120

  np_files:
    type:
    - "null"
    - File[]
    label: "Called peaks (narrowPeak format)"
    doc: |
      Called peaks files in narrowPeak format
    outputSource: pipe/np_files
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Called peaks"
        displayMode: "COLLAPSE"
        height: 40

  bp_files:
    type:
    - "null"
    - File[]
    label: "Called peaks (broadPeak format)"
    doc: |
      Called peaks files in broadPeak format
    outputSource: pipe/bp_files
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Called peaks"
        displayMode: "COLLAPSE"
        height: 40

  diff_sts_bigbed:
    type: File
    label: "Differentially bound sites (bigBed format)"
    doc: |
      Differentially bound sites in bigBed format
    outputSource: bed_to_bigbed/bigbed_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        format: 'bigbed'
        name: "Differentially bound sites"
        height: 40

  pk_vrlp_s_plot_png:
    type:
    - "null"
    - type: array
      items: File
    label: "Peakset overlap rate"
    doc: |
      Peakset overlap rate
      PNG format
    outputSource: diffbind/pk_vrlp_s_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'Peakset overlap rate'

  all_pk_scr_corr_plot_png:
    type: File?
    label: "Samples correlation (all peaks)"
    doc: |
      Samples correlation (all peaks)
      PNG format
    outputSource: diffbind/all_pk_scr_corr_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'Samples correlation (all peaks)'

  cns_pk_scr_corr_plot_png:
    type: File?
    label: "Samples correlation (opt. rec. cons. peaks)"
    doc: |
      Samples correlation (optionally
      recentered consensus peaks)
      PNG format
    outputSource: diffbind/cns_pk_scr_corr_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'Samples correlation (opt. rec. cons. peaks)'

  rw_rds_corr_plot_png:
    type: File?
    label: "Samples correlation (raw reads)"
    doc: |
      Samples correlation (raw reads)
      PNG format
    outputSource: diffbind/rw_rds_corr_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'Samples correlation (raw reads)'

  nr_rds_corr_plot_png:
    type: File?
    label: "Samples correlation (normalized reads)"
    doc: |
      Samples correlation (normalized reads)
      PNG format
    outputSource: diffbind/nr_rds_corr_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'Samples correlation (normalized reads)'

  pk_prfl_plot_png:
    type: File?
    label: "Peak profiles"
    doc: |
      Peak profiles
      PNG format
    outputSource: diffbind/pk_prfl_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Differential plots'
          Caption: 'Peak profiles'

  diff_vlcn_plot_png:
    type: File?
    label: "Volcano plot for differentially bound sites"
    doc: |
      Volcano plot for differentially bound sites
      PNG format
    outputSource: diffbind/diff_vlcn_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Differential plots'
          Caption: 'Volcano plot for differentially bound sites'

  diff_ma_plot_png:
    type: File?
    label: "MA-plot for differentially bound sites"
    doc: |
      MA-plot for differentially bound sites
      PNG format
    outputSource: diffbind/diff_ma_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Differential plots'
          Caption: 'MA-plot for differentially bound sites'

  nr_rds_pca_1_2_plot_png:
    type: File?
    label: "PCA (1,2) of not filtered normalized counts"
    doc: |
      PCA (1,2) of not filtered normalized counts
      PNG format
    outputSource: diffbind/nr_rds_pca_1_2_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'PCA (1,2) of not filtered normalized counts'

  nr_rds_pca_2_3_plot_png:
    type: File?
    label: "PCA (2,3) of not filtered normalized counts"
    doc: |
      PCA (2,3) of not filtered normalized counts
      PNG format
    outputSource: diffbind/nr_rds_pca_2_3_plot_png
    'sd:visualPlugins':
      - image:
          tab: 'Exploratory plots'
          Caption: 'PCA (2,3) of not filtered normalized counts'

  nr_rds_mds_html:
    type: File?
    outputSource: diffbind/nr_rds_mds_html
    label: "MDS plot of normalized counts"
    doc: |
      MDS plot of normalized counts.
      HTML format
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  diff_sts_tsv:
    type: File
    label: "Differentially bound sites with assigned nearest genes"
    doc: |
      Differentially bound sites with assigned nearest genes
      TSV format
    outputSource: restore_columns/output_file
    'sd:visualPlugins':
      - syncfusiongrid:
          tab: 'Differentially bound sites'
          Title: 'Differentially bound sites'

  diff_sts_labeled_tsv:
    type: File
    label: "Differentially bound sites with labels"
    doc: |
      Differentially bound sites with labels
      TSV format
    outputSource: add_label_column/output_file

  volcano_plot_html_file:
    type: File
    label: "Volcano Plot"
    doc: |
      HTML index file for Volcano Plot
    outputSource: make_volcano_plot/html_file
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  volcano_plot_html_data:
    type: Directory
    label: "Directory html data for Volcano Plot"
    doc: |
      Directory html data for Volcano Plot
    outputSource: make_volcano_plot/html_data

  ma_plot_html_file:
    type: File
    label: "MA-plot"
    doc: |
      HTML index file for MA-plot
    outputSource: make_ma_plot/html_file
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  ma_plot_html_data:
    type: Directory
    label: "Directory html data for Volcano Plot"
    doc: |
      Directory html data for MA-plot
    outputSource: make_ma_plot/html_data

  heatmap_html:
    type: File
    label: "Heatmap of normalized counts"
    doc: |
      Morpheus heatmap in HTML format
    outputSource: morpheus_heatmap/heatmap_html
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  nr_rds_gct:
    type: File
    label: "GCT file with normalized read counts per peak"
    doc: |
      GCT file with normalized read counts per peak
    outputSource: extend_gct/extended_gct

  experiment_info:
    type: File
    label: "Samples order for IGV"
    doc: |
      Markdown file to explain the sample order for IGV
    outputSource: create_metadata/output_file
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  pdf_plots:
    type: File
    outputSource: compress_pdf_plots/compressed_folder
    label: "Plots in PDF format"
    doc: |
      Compressed folder with plots
      in PDF format

  diffbind_stdout_log:
    type: File
    label: "DiffBind stdout log"
    doc: |
      DiffBind stdout log
    outputSource: diffbind/stdout_log

  diffbind_stderr_log:
    type: File
    label: "DiffBind stderr log"
    doc: |
      DiffBind stderr log
    outputSource: diffbind/stderr_log

  morpheus_stdout_log:
    type: File
    label: "Morpheus heatmap stdout log"
    doc: |
      Morpheus heatmap stdout log
    outputSource: morpheus_heatmap/stdout_log

  morpheus_stderr_log:
    type: File
    label: "Morpheus heatmap stderr log"
    doc: |
      Morpheus heatmap stderr log
    outputSource: morpheus_heatmap/stderr_log

  iaintersect_result:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Differential binding analysis results formatted as chip/atac/cutandrun results"
    doc: "Differential binding analysis results  formatted as chip/atac/cutandrun results exported as TSV"
    outputSource: assign_genes/result_file


steps:

  pipe:
    run:
      cwlVersion: v1.0
      class: ExpressionTool
      inputs:
        genome_coverage_files:
          type: File[]
        narrow_peak_files:
          type:
          - "null"
          - File[]
        broad_peak_files:
          type:
          - "null"
          - File[]
      outputs:
        gc_files:
          type: File[]
        np_files:
          type:
          - "null"
          - File[]
        bp_files:
          type:
          - "null"
          - File[]
      expression: |
        ${
          var results = {};
          var output_names = [
            "gc_files",
            "np_files",
            "bp_files"
          ];
          var sources = [
            inputs.genome_coverage_files,
            inputs.narrow_peak_files,
            inputs.broad_peak_files
          ];
          for (var i = 0; i < sources.length; i++){
            var current_source = sources[i];
            var current_output_name = output_names[i];
            results[current_output_name] = null;
            if (current_source != null && current_source.length > 0){
              for (var j = 0; j < current_source.length; j++){
                    var new_item = current_source[j];
                    new_item["basename"] = "u" + "_" + i + "_" + j+ "_" + new_item.basename;
                    if (results[current_output_name] == null){
                      results[current_output_name] = [new_item];
                    } else {
                      results[current_output_name].push(new_item);
                    }
              }
            }
          }
          return results;
        }
    in:
      genome_coverage_files: genome_coverage_files
      narrow_peak_files: narrow_peak_files
      broad_peak_files: broad_peak_files
    out:
    - gc_files
    - np_files
    - bp_files

  diffbind:
    run: ../tools/diffbind-multi-factor.cwl
    in:
      alignment_files: alignment_files
      peak_files: peak_files
      dataset_names: dataset_names
      metadata_file: metadata_file
      scoreby: scoreby
      score_threshold: score_threshold
      rpkm_threshold: rpkm_threshold
      rec_summits: rec_summits
      overlap_threshold: overlap_threshold
      groupby:
        source: groupby
        valueFrom: $(split_by_common_delim(self))
      design_formula: design_formula
      contrast:
        source: contrast
        valueFrom: $(self==""?null:self)                 # safety measure
      base_levels:
        source: base_levels
        valueFrom: $(split_by_common_delim(self))
      analysis_method:
        default: "deseq2"                                # hardcoded to always use DESeq2 because EdgeR fails to run without contrast
      normalization_method:
        default: "auto"                                  # harcoded to auto as we don't allow to use EdgeR
      padj_threshold: padj_threshold
      cluster_method:
        source: cluster_method
        valueFrom: $(self=="none"?null:self)
      row_distance: row_distance
      column_distance: column_distance
      center_row:
        default: true
      export_pdf_plots:
        default: true
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - pk_vrlp_s_plot_png
    - all_pk_scr_corr_plot_png
    - cns_pk_scr_corr_plot_png
    - rw_rds_corr_plot_png
    - nr_rds_corr_plot_png
    - pk_prfl_plot_png
    - diff_vlcn_plot_png
    - diff_ma_plot_png
    - nr_rds_pca_1_2_plot_png
    - nr_rds_pca_2_3_plot_png
    - pk_vrlp_s_plot_pdf
    - all_pk_scr_corr_plot_pdf
    - cns_pk_scr_corr_plot_pdf
    - rw_rds_corr_plot_pdf
    - nr_rds_corr_plot_pdf
    - pk_prfl_plot_pdf
    - diff_vlcn_plot_pdf
    - diff_ma_plot_pdf
    - nr_rds_pca_1_2_plot_pdf
    - nr_rds_pca_2_3_plot_pdf
    - nr_rds_mds_html
    - diff_sts_tsv
    - nr_rds_gct
    - stdout_log
    - stderr_log

  folder_pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - diffbind/pk_vrlp_s_plot_pdf
        - diffbind/all_pk_scr_corr_plot_pdf
        - diffbind/cns_pk_scr_corr_plot_pdf
        - diffbind/rw_rds_corr_plot_pdf
        - diffbind/nr_rds_corr_plot_pdf
        - diffbind/pk_prfl_plot_pdf
        - diffbind/diff_vlcn_plot_pdf
        - diffbind/diff_ma_plot_pdf
        - diffbind/nr_rds_pca_1_2_plot_pdf
        - diffbind/nr_rds_pca_2_3_plot_pdf
        valueFrom: $(self.flat().filter(n => n))
      folder_basename:
        default: "pdf_plots"
    out:
    - folder

  compress_pdf_plots:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: folder_pdf_plots/folder
    out:
    - compressed_folder

  filter_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file: diffbind/diff_sts_tsv
      script:
        default: >
          cat $0 | grep -v "Start" | awk
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t"NR"\t0\t0\t0\t0"}' > `basename $0`
    out:
    - output_file

  assign_genes:
    run: ../tools/iaintersect.cwl
    in:
      input_filename: filter_columns/output_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
    out:
    - result_file

  restore_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file:
      - assign_genes/result_file
      - diffbind/diff_sts_tsv
      script:
        default: |
          cat $0 | grep -v "start" | sort -k 11n | cut -f 1-5,15 > iaintersect_result.tsv
          cat $1 | grep -v "Start" > diffbind_result.tsv
          HEADER=`head -n 1 $1`;
          echo -e "Refseq_id\tGene_id\ttxStart\ttxEnd\tStrand\tRegion\t${HEADER}" > `basename $0`;
          cat iaintersect_result.tsv | paste - diffbind_result.tsv >> `basename $0`
          rm iaintersect_result.tsv diffbind_result.tsv
    out:
    - output_file

  convert_to_bed:
    run: ../tools/custom-bash.cwl
    in:
      input_file: restore_columns/output_file
      script:
        default: |
          cat "$0" | awk -F "\t" 'NR==1 {for (i=1; i<=NF; i++) {ix[$i]=i} } NR>1 {color="255,0,0"; if ($ix["log2FoldChange"]<0) color="0,255,0"; print $ix["Chr"]"\t"$ix["Start"]"\t"$ix["End"]"\tpvalue="$ix["pvalue"]+0.0";padj="$ix["padj"]+0.0";log2FC="$ix["log2FoldChange"]"\t"1000"\t"$ix["Strand"]"\t"$ix["Start"]"\t"$ix["End"]"\t"color}' > `basename $0`
    out:
    - output_file

  sort_bed:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: convert_to_bed/output_file
      key:
        default: ["1,1","2,2n"]
    out:
    - sorted_file

  overlap_with_chr_length:
    run: ../tools/custom-bedops.cwl
    in:
      input_file:
      - chrom_length_file
      - sort_bed/sorted_file
      script:
        default: |
          cat "$0" | awk '{print $1"\t0\t"$2}' | sort-bed - > temp_chrom_length.bed
          cat "$1" | awk '$2 >= 0' > temp_sorted.bed
          bedops --element-of 100% temp_sorted.bed temp_chrom_length.bed > `basename $1`
          rm -f temp_chrom_length.bed temp_sorted.bed
    out:
    - output_file

  bed_to_bigbed:
    run: ../tools/ucsc-bedtobigbed.cwl
    in:
      input_bed: overlap_with_chr_length/output_file
      bed_type:
        default: "bed4+5"
      chrom_length_file: chrom_length_file
      output_filename:
        source: overlap_with_chr_length/output_file
        valueFrom: $(self.basename.split('.').slice(0,-1).join('.') + ".bigBed")
    out: 
    - bigbed_file

  add_label_column:
    run: ../tools/custom-bash.cwl
    in:
      input_file: diffbind/diff_sts_tsv
      script:
        default: |
          HEADER=`head -n 1 $0`;
          echo -e "label\t${HEADER}" > diff_sts_labeled.tsv;
          cat "$0" | grep -v "Start" | awk -F "\t" '{print $1":"$2"-"$3"\t"$0}' >> diff_sts_labeled.tsv
    out:
    - output_file

  make_volcano_plot:
    run: ../tools/volcano-plot.cwl
    in:
      diff_expr_file: add_label_column/output_file
      x_axis_column:
        default: "log2FoldChange"
      y_axis_column:
        default: "padj"
      label_column:
        default: "label"
    out:
    - html_data
    - html_file

  make_ma_plot:
    run: ../tools/ma-plot.cwl
    in:
      diff_expr_file: add_label_column/output_file
      x_axis_column:
        default: "baseMean"
      y_axis_column:
        default: "log2FoldChange"
      label_column:
        default: "label"
    out:
    - html_data
    - html_file

  extend_gct:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/morpheus:v0.0.2
      - class: InitialWorkDirRequirement
        listing:
        - entryname: extend.R
          entry: |
            options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})
            suppressMessages(library("cmapR"))
            suppressMessages(library("dplyr"))
            suppressMessages(library("tibble"))
            suppressMessages(library("morpheus"))
            suppressMessages(library("argparse"))
            args = commandArgs(trailingOnly=TRUE)
            gct_data <- read.gct(args[1])
            metadata <- read.table(args[2], sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE) %>%
                        mutate(id=paste(Chr, paste(Start, End, sep="-"), sep=":")) %>%
                        select(id, Gene_id, Region)
            row_metadata <- gct_data$rowAnnotations %>%
                            rownames_to_column("id") %>%
                            left_join(metadata, by="id") %>%
                            mutate_at("id", as.vector)
            col_metadata <- gct_data$columnAnnotations %>%
                            rownames_to_column("id") %>%
                            mutate_at("id", as.vector)
            gct_data <- new(
                "GCT",
                mat=gct_data$data[row_metadata$id, col_metadata$id],
                rdesc=row_metadata,
                cdesc=col_metadata
            )
            write_gct(ds=gct_data, ofile="extended.gct", appenddim=FALSE)
      inputs:
        input_files:
          type: File[]
          inputBinding:
            position: 5
      outputs:
        extended_gct:
          type: File
          outputBinding:
            glob: "extended.gct"
      baseCommand: ["Rscript", "extend.R"]
    in:
      input_files:
      - diffbind/nr_rds_gct
      - restore_columns/output_file
    out:
    - extended_gct

  morpheus_heatmap:
    run: ../tools/morpheus-heatmap.cwl
    in:
     read_counts_gct: extend_gct/extended_gct
    out:
    - heatmap_html
    - stdout_log
    - stderr_log

  create_metadata:
    run: ../tools/custom-bash.cwl
    in:
      input_file: peak_files
      param: dataset_names
      script:
        default: |
          #!/bin/bash
          set -- "$0" "$@"
          COUNT=`expr $# / 2`
          echo "| Sample | Index |" > experiment_info.md
          echo "| :-- | --: |" >> experiment_info.md
          j=1
          for i in "${@:$COUNT+1:$#}"; do
            echo "| $i | $j |" >> experiment_info.md
            (( j++ ))
          done;
    out:
    - output_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "DiffBind Multi-factor Analysis"
s:name: "DiffBind Multi-factor Analysis"
s:alternateName: "Runs DiffBind multi-factor analysis with manual control over major parameters"


s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/diffbind-multi-factor.cwl
s:codeRepository: https://github.com/Barski-lab/workflows-datirium
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


doc: |
  DiffBind Multi-factor Analysis
  ------------------------------
  
  DiffBind processes ChIP-Seq data enriched for genomic loci where specific protein/DNA binding occurs, including peak sets identified by ChIP-Seq peak callers and
  aligned sequence read datasets. It is designed to work with multiple peak sets simultaneously, representing different ChIP experiments (antibodies, transcription
  factor and/or histone marks, experimental conditions, replicates) as well as managing the results of multiple peak callers.

  For more information please refer to:
  -------------------------------------
  Ross-Innes CS, Stark R, Teschendorff AE, Holmes KA, Ali HR, Dunning MJ, Brown GD, Gojis O, Ellis IO, Green AR, Ali S, Chin S, Palmieri C, Caldas C, Carroll JS (2012).
  “Differential oestrogen receptor binding is associated with clinical outcome in breast cancer.” Nature, 481, -4.