cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var split_features = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          let splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };
    - var split_by_comma = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          var splitted_line = line?line.split(/,+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };


"sd:upstream":
  sc_tools_sample:
  - "sc-atac-reduce.cwl"
  - "sc-atac-cluster.cwl"
  - "sc-rna-reduce.cwl"
  - "sc-rna-cluster.cwl"
  - "sc-wnn-cluster.cwl"
  - "sc-ctype-assign.cwl"
  sc_atac_sample:
  - "cellranger-arc-count.cwl"
  - "cellranger-arc-aggr.cwl"
  - "cellranger-atac-count.cwl"
  - "cellranger-atac-aggr.cwl"
  genome_indices:
  - "genome-indices.cwl"


inputs:

  alias:
    type: string
    label: "Analysis name"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Single-cell Analysis with ATAC-Seq Datasets"
    doc: |
      Analysis that includes single-cell
      multiome RNA and ATAC-Seq or just
      ATAC-Seq datasets run through either
      "Single-Cell ATAC-Seq Dimensionality
      Reduction Analysis", "Single-Cell
      RNA-Seq Dimensionality Reduction
      Analysis" or "Single-Cell WNN Cluster
      Analysis" at any of the processing stages.
    "sd:upstreamSource": "sc_tools_sample/seurat_data_rds"
    "sd:localLabel": true

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    label: "Cell Ranger ATAC or RNA+ATAC Sample"
    doc: |
      Any "Cell Ranger ATAC or RNA+ATAC Sample"
      for loading chromatin accessibility data
      from. This sample can be obtained from
      one of the following pipelines: "Cell
      Ranger Count (RNA+ATAC)", "Cell Ranger
      Aggregate (RNA+ATAC)", "Cell Ranger Count
      (ATAC)", or "Cell Ranger Aggregate (ATAC)".
    "sd:upstreamSource": "sc_atac_sample/atac_fragments_file"
    "sd:localLabel": true

  genome_type:
    type: string
    label: "Genome type"
    doc: |
      Reference genome
    "sd:upstreamSource": "genome_indices/genome"
    "sd:localLabel": true

  annotation_file:
    type: File
    "sd:upstreamSource": "genome_indices/annotation"

  chrom_length_file:
    type: File
    "sd:upstreamSource": "genome_indices/chrom_length"

  blacklist_regions_file:
    type:
      type: enum
      symbols:
      - "hg19"
      - "hg38"
      - "mm10"
    "sd:upstreamSource": "genome_indices/genome"

  groupby:
    type: string?
    default: null
    label: "Subsetting category (optional)"
    doc: |
      Single cell metadata column to group
      cells for optional subsetting before
      running differential accessibility analysis.
      To group cells by dataset, use "dataset".
      Custom groups can be defined based on
      any single cell metadata added through
      the "Datasets metadata (optional)" or
      "Selected cell barcodes (optional)"
      inputs. Default: do not subset cells

  subset:
    type: string?
    default: null
    label: "Subsetting values (optional)"
    doc: |
      Comma separated list of values from
      the single cell metadata column
      selected in "Subsetting category
      (optional)" input. Ignored if grouping
      category is not provided. Default: do
      not subset cells

  splitby:
    type: string
    label: "Comparison category"
    doc: |
      Single cell metadata column to split
      cells into two comparison groups before
      running differential accessibility analysis.
      To split cells by dataset, use "dataset".
      Custom groups can be defined based on
      any single cell metadata added through
      the "Datasets metadata (optional)" or
      "Selected cell barcodes (optional)"
      inputs. The direction of comparison is
      always "Second comparison group" vs
      "First comparison group".

  first_cond:
    type: string
    label: "First comparison group"
    doc: |
      Value from the single cell metadata
      column selected in "Comparison category"
      input to define the first group of cells
      for differential accessibility analysis.

  second_cond:
    type: string
    label: "Second comparison group"
    doc: |
      Value from the single cell metadata
      column selected in "Comparison category"
      input to define the second group of cells
      for differential accessibility analysis.

  analysis_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "negative-binomial"          # (negbinom) Negative Binomial Generalized Linear Model (use FindMarkers with peaks from Seurat object)
      - "poisson"                    # (poisson) Poisson Generalized Linear Model (use FindMarkers with peaks from Seurat object)
      - "logistic-regression"        # (LR) Logistic Regression (use FindMarkers with peaks from Seurat object)
      - "mast"                       # (MAST) MAST package (use FindMarkers with peaks from Seurat object)
      - "manorm2"                    # call peaks for each group with MACS2, run MAnorm2
    default: "logistic-regression"
    label: "Statistical test"
    doc: |
      Statistical test to use in the
      differential accessibility analysis.
      Chromatin accessibility data will first
      be filtered to include fragments from
      cells retained in the loaded single-cell
      datasets after optional filtering by
      "Subsetting category/values" and/or
      "Selected cell barcodes". The resulting
      fragments will then be splitted by
      "Comparison category" and, if "manorm2"
      is selected, aggregated to the pseudo
      bulk form for peak calling with MACS2.
      Othwerwise, the analysis will be run on
      the cells level using peaks already
      available in the loaded single-cell
      datasets.
      Default: logistic-regression

  maximum_padj:
    type: float?
    default: 0.05
    label: "Maximum adjusted p-value"
    doc: |
      Maximum adjusted p-value threshold
      for selecting differentially accessible
      regions to be visualized on the heatmap.
      Default: 0.05

  minimum_logfc:
    type: float?
    default: 1
    label: "Minimum log2 fold change absolute value"
    doc: |
      Minimum log2 fold change absolute value
      for selecting differentially accessible
      regions to be visualized on the heatmap.
      Default: 1.0

  datasets_metadata:
    type: File?
    label: "Datasets metadata (optional)"
    doc: |
      If the selected single-cell analysis
      includes multiple aggregated datasets,
      each of them can be assigned to a
      separate group by one or multiple
      categories. This can be achieved by
      providing a TSV/CSV file with
      "library_id" as the first column and
      any number of additional columns with
      unique names, representing the desired
      grouping categories. To obtain a proper
      template of this file, download
      "datasets_metadata.tsv" output from
      the "Files" tab of either "Single-Cell
      Multiome ATAC-Seq and RNA-Seq Filtering
      Analysis" or "Single-Cell ATAC-Seq
      Filtering Analysis" and add extra
      columns as needed. The pipeline selection
      depends on whether the "Single-cell
      Analysis with ATAC-Seq Datasets" includes
      both RNA and ATAC-Seq or only ATAC-Seq
      datasets.

  barcodes_data:
    type: File?
    label: "Selected cell barcodes (optional)"
    doc: |
      A TSV/CSV file to optionally prefilter
      the single cell data by including only
      the cells with the selected barcodes.
      The provided file should include at
      least one column named "barcode", with
      one cell barcode per line. All other
      columns, except for "barcode", will be
      added to the single cell metadata loaded
      from "Single-cell Analysis with ATAC-Seq
      Datasets" and can be utilized in the
      current or future steps of analysis.

  minimum_qvalue:
    type: float?
    default: 0.05
    label: "Minimum peak calling FDR (for manorm2)"
    doc: |
      Minimum FDR (q-value) cutoff for MACS2
      peak detection. Ignored if "Statistical
      test" input is not set to "manorm2".
      Default: 0.05
    "sd:layout":
      advanced: true

  maximum_peaks:
    type: int?
    default: 0
    label: "Maximum number of peaks to keep, set to 0 to keep all (for manorm2)"
    doc: |
      The maximum number of the most significant
      peaks to select for each of the comparison
      groups when constructing reference genomic
      bins. The top significant peaks are selected
      based on the score column which is calculated
      by MACS2 as int(-10*log10qvalue). Ignored
      if "Statistical test" input is not set to
      "manorm2".
      Default: 0 - keep all peaks.
    "sd:layout":
      advanced: true

  minimum_peak_gap:
    type: int?
    default: 150
    label: "Minimum distance for merging peaks (for manorm2)"
    doc: |
      Peaks remained after optional filtering
      by the "Maximum number of peaks to keep"
      will be merged if the distance between
      them is smaller than the provided value.
      Merging is first applied within comparison
      group and then to all peaks together before
      splitting them into the reference genomic
      bins of the selected size. Ignored if
      "Statistical test" input is not set to
      "manorm2".
      Default: 150
    "sd:layout":
      advanced: true

  bin_size:
    type: int?
    default: 1000
    label: "Bin size (for manorm2)"
    doc: |
      The size of non-overlapping reference
      genomic bins used by MAnorm2 when
      generating a table of reads per bins
      counts. Ignored if "Statistical test"
      input is not set to "manorm2".
      Default: 1000
    "sd:layout":
      advanced: true

  promoter_dist:
    type: int?
    default: 1000
    label: "Promoter distance, bp"
    doc: |
      Max distance from the gene TSS (in both
      direction) overlapping which the
      differentially accessible site will be
      assigned to the promoter region.
      Default: 1000 bp
    "sd:layout":
      advanced: true

  upstream_dist:
    type: int?
    default: 20000
    label: "Upstream distance, bp"
    doc: |
      Max distance from the promoter (only
      in upstream direction) overlapping
      which the differentially accessible site
      will be assigned to the upstream region.
      Default: 20,000 bp
    "sd:layout":
      advanced: true

  export_html_report:
    type: boolean?
    default: true
    label: "Show HTML report"
    doc: |
      Export tehcnical report in HTML format.
      Default: true
    "sd:layout":
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
      - "5"
      - "6"
    default: "4"
    label: "Number of cores/cpus to use"
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 4
    "sd:layout":
      advanced: true


outputs:

  umap_rd_rnaumap_plot_png:
    type: File?
    outputSource: sc_atac_dbinding/umap_rd_rnaumap_plot_png
    label: "UMAP, split by comparison category, RNA"
    doc: |
      UMAP, split by the single cell metadata
      column defined in the "Comparison category",
      optionally subsetted to include only cells
      with "Subsetting values (optional)" from
      the "Subsetting category (optional)", RNA.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "UMAP, split by comparison category, RNA"

  umap_rd_atacumap_plot_png:
    type: File?
    outputSource: sc_atac_dbinding/umap_rd_atacumap_plot_png
    label: "UMAP, split by comparison category, ATAC"
    doc: |
      UMAP, split by the single cell metadata
      column defined in the "Comparison category",
      optionally subsetted to include only cells
      with "Subsetting values (optional)" from
      the "Subsetting category (optional)", ATAC.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "UMAP, split by comparison category, ATAC"

  umap_rd_wnnumap_plot_png:
    type: File?
    outputSource: sc_atac_dbinding/umap_rd_wnnumap_plot_png
    label: "UMAP, split by comparison category, WNN"
    doc: |
      UMAP, split by the single cell metadata
      column defined in the "Comparison category",
      optionally subsetted to include only cells
      with "Subsetting values (optional)" from
      the "Subsetting category (optional)", WNN.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "QC"
        Caption: "UMAP, split by comparison category, WNN"

  dbnd_vlcn_plot_png:
    type: File?
    outputSource: sc_atac_dbinding/dbnd_vlcn_plot_png
    label: "Volcano plot of differentially accessible regions"
    doc: |
      Volcano plot of differentially
      accessible regions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Volcano plot"
        Caption: "Volcano plot of differentially accessible regions"

  tag_density_heatmap:
    type: File
    outputSource: make_heatmap/heatmap_file
    label: "Tag density heatmap around centers of filtered diff. accessible regions"
    doc: |
      Tag density heatmap around centers of
      filtered by "Maximum adjusted p-value"
      and "Minimum log2 fold change absolute
      value" differentially accessible regions.
      PNG format.
    "sd:visualPlugins":
    - image:
        tab: "Heatmap"
        Caption: "Tag density heatmap around centers of filtered diff. accessible regions"

  tag_density_matrix:
    type: File
    outputSource: compute_score_matrix/scores_matrix
    label: "Score matrix for tag density heatmap around centers of filtered diff. accessible regions"
    doc: |
      Scores matrix generated by Deeptools
      with tag density information around
      centers of filrered differentially
      accessible regions.

  first_fragments_bigwig_file:
    type: File
    outputSource: sc_atac_dbinding/first_fragments_bigwig_file
    label: "Genome coverage (first comparison group)"
    doc: |
      Normalized genome coverage calculated
      for all ATAC fragments from the cells
      that belong to the "First comparison
      group".
      BigWig format.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "Genome coverage (first comparison group)"
        height: 120

  first_peaks_bed_file:
    type: File?
    outputSource: sc_atac_dbinding/first_peaks_bed_file
    label: "Called peaks (first comparison group)"
    doc: |
      Peaks called from the aggregated to
      the pseudo bulk form ATAC fragments.
      Only cells that belong to the "First
      comparison group", optionally filtered
      by the "Subsetting category/values" and/or
      "Selected cell barcodes", were selected.
      NarrowPeak format.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "Called peaks (first comparison group)"
        displayMode: "COLLAPSE"
        height: 40

  first_peaks_xls_file:
    type: File?
    outputSource: sc_atac_dbinding/first_peaks_xls_file
    label: "Called peaks (first comparison group, xls format)"
    doc: |
      Peaks called from the aggregated to
      the pseudo bulk form ATAC fragments.
      Only cells that belong to the "First
      comparison group", optionally filtered
      by the "Subsetting category/values" and/or
      "Selected cell barcodes", were selected.
      XLS format.

  first_summits_bed_file:
    type: File?
    outputSource: sc_atac_dbinding/first_summits_bed_file
    label: "Peaks summits (first comparison group)"
    doc: |
      Summits of the peaks called from
      the aggregated to the pseudo bulk
      form ATAC fragments. Only cells that
      belong to the "First comparison group",
      optionally filtered by the "Subsetting
      category/values" and/or "Selected cell
      barcodes", were selected.
      BED format.

  second_fragments_bigwig_file:
    type: File
    outputSource: sc_atac_dbinding/second_fragments_bigwig_file
    label: "Genome coverage (second comparison group)"
    doc: |
      Normalized genome coverage calculated
      for all ATAC fragments from the cells
      that belong to the "Second comparison
      group".
      BigWig format.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "Genome coverage (second comparison group)"
        height: 120

  second_peaks_bed_file:
    type: File?
    outputSource: sc_atac_dbinding/second_peaks_bed_file
    label: "Called peaks (second comparison group)"
    doc: |
      Peaks called from the aggregated to
      the pseudo bulk form ATAC fragments.
      Only cells that belong to the "Second
      comparison group", optionally filtered
      by the "Subsetting category/values" and/or
      "Selected cell barcodes", were selected.
      NarrowPeak format.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "Called peaks (second comparison group)"
        displayMode: "COLLAPSE"
        height: 40

  second_peaks_xls_file:
    type: File?
    outputSource: sc_atac_dbinding/second_peaks_xls_file
    label: "Called peaks (second comparison group, xls format)"
    doc: |
      Peaks called from the aggregated to
      the pseudo bulk form ATAC fragments.
      Only cells that belong to the "Second
      comparison group", optionally filtered
      by the "Subsetting category/values" and/or
      "Selected cell barcodes", were selected.
      XLS format.

  second_summits_bed_file:
    type: File?
    outputSource: sc_atac_dbinding/second_summits_bed_file
    label: "Peaks summits (second comparison group)"
    doc: |
      Summits of the peaks called from
      the aggregated to the pseudo bulk
      form ATAC fragments. Only cells that
      belong to the "Second comparison group",
      optionally filtered by the "Subsetting
      category/values" and/or "Selected cell
      barcodes", were selected.
      BED format.

  seurat_peaks_bigbed_file:
    type: File?
    outputSource: sc_atac_dbinding/seurat_peaks_bigbed_file
    label: "Original peaks (both comparison groups)"
    doc: |
      Original peaks from the loaded "Single-cell
      Analysis with ATAC-Seq Datasets".
      BigBed format.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        format: "bigbed"
        name: "Original peaks (both comparison groups)"
        displayMode: "COLLAPSE"
        height: 40

  diff_bound_sites:
    type: File
    outputSource: restore_columns/output_file
    label: "Differentially accessible regions (not filtered)"
    doc: |
      Not filtered by "Maximum adjusted p-value"
      and "Minimum log2 fold change absolute value"
      differentially accessible regions with the
      nearest genes assigned.
      TSV format.
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Diff. accessible regions"
        Title: "Differentially accessible regions (not filtered)"

  diff_bound_sites_with_labels:
    type: File
    outputSource: add_label_column/output_file
    label: "Differentially accessible regions with labels (not filtered)"
    doc: |
      Not filtered by "Maximum adjusted p-value"
      and "Minimum log2 fold change absolute value"
      differentially accessible regions with labels.
      TSV format.
    'sd:visualPlugins':
    - queryRedirect:
        tab: "Overview"
        label: "Volcano Plot"
        url: "https://scidap.com/vp/volcano"
        query_eval_string: "`data_file=${this.getSampleValue('outputs', 'diff_bound_sites_with_labels')}&data_col=label&x_col=log2FoldChange&y_col=padj`"

  diff_bound_sites_bigbed:
    type: File
    label: "Differentially accessible regions (not filtered)"
    doc: |
      Not filtered by "Maximum adjusted p-value"
      and "Minimum log2 fold change absolute value"
      differentially accessible regions with the
      nearest genes assigned.
      BigBed format.
    outputSource: bed_to_bigbed/bigbed_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        format: "bigbed"
        name: "Differentially accessible regions (not filtered)"
        displayMode: "COLLAPSE"
        height: 40

  pdf_plots:
    type: File
    outputSource: compress_pdf_plots/compressed_folder
    label: "Compressed folder with all PDF plots"
    doc: |
      Compressed folder with all PDF plots.

  sc_report_html_file:
    type: File?
    outputSource: sc_atac_dbinding/sc_report_html_file
    label: "Analysis log"
    doc: |
      Tehcnical report.
      HTML format.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  sc_atac_dbinding_stdout_log:
    type: File
    outputSource: sc_atac_dbinding/stdout_log
    label: "Output log"
    doc: |
      Stdout log from the sc_atac_dbinding step.

  sc_atac_dbinding_stderr_log:
    type: File
    outputSource: sc_atac_dbinding/stderr_log
    label: "Error log"
    doc: |
      Stderr log from the sc_atac_dbinding step.


steps:

  sc_atac_dbinding:
    run: ../tools/sc-atac-dbinding.cwl
    in:
      query_data_rds: query_data_rds
      atac_fragments_file: atac_fragments_file
      datasets_metadata: datasets_metadata
      barcodes_data: barcodes_data
      groupby:
        source: groupby
        valueFrom: |
          ${
            if (self == "dataset") {
              return "new.ident";
            } else if (self == "") {
              return null;
            } else {
              return self;
            }
          }
      subset:
        source: subset
        valueFrom: $(split_by_comma(self))
      splitby:
        source: splitby
        valueFrom: |
          ${
            if (self == "dataset") {
              return "new.ident";
            } else {
              return self;
            }
          }
      first_cond: first_cond
      second_cond: second_cond
      analysis_method: analysis_method
      genome_type:
        source: genome_type
        valueFrom: $(self=="mm10"?"mm":"hs")
      minimum_qvalue: minimum_qvalue
      minimum_peak_gap: minimum_peak_gap
      bin_size: bin_size
      maximum_peaks:
        source: maximum_peaks
        valueFrom: $(self==0?null:self)                 # to return null for 0
      blacklist_regions_file: blacklist_regions_file
      maximum_padj: maximum_padj
      minimum_logfc: minimum_logfc
      verbose:
        default: true
      export_pdf_plots:
        default: true
      parallel_memory_limit:
        default: 32
      vector_memory_limit:
        default: 128
      export_html_report: export_html_report
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - umap_rd_rnaumap_plot_png
    - umap_rd_atacumap_plot_png
    - umap_rd_wnnumap_plot_png
    - seurat_peaks_bigbed_file
    - first_fragments_bigwig_file
    - second_fragments_bigwig_file
    - first_peaks_xls_file
    - second_peaks_xls_file
    - first_peaks_bed_file
    - second_peaks_bed_file
    - first_summits_bed_file
    - second_summits_bed_file
    - diff_bound_sites
    - dbnd_vlcn_plot_png
    - first_enrch_bed_file
    - second_enrch_bed_file
    - all_plots_pdf
    - sc_report_html_file
    - stdout_log
    - stderr_log

  folder_pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - sc_atac_dbinding/all_plots_pdf
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
      input_file: sc_atac_dbinding/diff_bound_sites
      script:
        default: >
          cat $0 | grep -v "start" | awk
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
      - sc_atac_dbinding/diff_bound_sites
      script:
        default: |
          cat $0 | grep -v "start" | sort -k 11n | cut -f 1-5,15 > iaintersect_result.tsv
          cat $1 | grep -v "start" > sc_atac_dbinding_result.tsv
          HEADER=`head -n 1 $1`;
          echo -e "refseq_id\tgene_id\ttxStart\ttxEnd\tstrand\tregion\t${HEADER}" > `basename $0`;
          cat iaintersect_result.tsv | paste - sc_atac_dbinding_result.tsv >> `basename $0`
          rm iaintersect_result.tsv sc_atac_dbinding_result.tsv
    out:
    - output_file

  convert_to_bed:
    run: ../tools/custom-bash.cwl
    in:
      input_file: restore_columns/output_file
      script:
        default: |
          cat "$0" | awk -F "\t" 'NR==1 {for (i=1; i<=NF; i++) {ix[$i]=i} } NR>1 {color="255,0,0"; if ($ix["log2FoldChange"]<0) color="0,255,0"; print $ix["chr"]"\t"$ix["start"]"\t"$ix["end"]"\tpvalue="$ix["pvalue"]+0.0";padj="$ix["padj"]+0.0";log2FC="$ix["log2FoldChange"]"\t"1000"\t"$ix["strand"]"\t"$ix["start"]"\t"$ix["end"]"\t"color}' > `basename $0`
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
      input_file: sc_atac_dbinding/diff_bound_sites
      script:
        default: |
          HEADER=`head -n 1 $0`;
          echo -e "label\t${HEADER}" > sc_db_sites_labeled.tsv;
          cat "$0" | grep -v "start" | awk -F "\t" '{print $1":"$2"-"$3"-"$NF"\t"$0}' >> sc_db_sites_labeled.tsv
    out:
    - output_file

  recenter_first_enrch_bed:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap:v0.0.3
      requirements:
      - class: InitialWorkDirRequirement
        listing:
        - entryname: dummy.csv
          entry: |
            chr1,1,10
      inputs:
        script:
          type: string?
          default: |
            #!/bin/bash
            cat "$0" | tr -d "\r" | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n | awk '{center=$2+int(($3-$2)/2); print $1"\t"center"\t"center+1}' > first_recentered.bed
          inputBinding:
            position: 1
        regions_file:
          type:
          - "null"
          - string
          - File
          inputBinding:
            position: 5
            valueFrom: $(self==""?"dummy.csv":self)
          default: ""
      outputs:
        recentered_regions_file:
          type: File
          outputBinding:
            glob: "first_recentered.bed"
      baseCommand: [bash, '-c']
    in:
      regions_file: sc_atac_dbinding/first_enrch_bed_file
    out:
    - recentered_regions_file

  recenter_second_enrch_bed:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap:v0.0.3
      requirements:
      - class: InitialWorkDirRequirement
        listing:
        - entryname: dummy.csv
          entry: |
            chr1,1,10
      inputs:
        script:
          type: string?
          default: |
            #!/bin/bash
            cat "$0" | tr -d "\r" | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n | awk '{center=$2+int(($3-$2)/2); print $1"\t"center"\t"center+1}' > second_recentered.bed
          inputBinding:
            position: 1
        regions_file:
          type:
          - "null"
          - string
          - File
          inputBinding:
            position: 5
            valueFrom: $(self==""?"dummy.csv":self)
          default: ""
      outputs:
        recentered_regions_file:
          type: File
          outputBinding:
            glob: "second_recentered.bed"
      baseCommand: [bash, '-c']
    in:
      regions_file: sc_atac_dbinding/second_enrch_bed_file
    out:
    - recentered_regions_file

  compute_score_matrix:
    run: ../tools/deeptools-computematrix-referencepoint.cwl
    in:
      score_files:
      - sc_atac_dbinding/first_fragments_bigwig_file
      - sc_atac_dbinding/second_fragments_bigwig_file
      regions_files:
      - recenter_first_enrch_bed/recentered_regions_file
      - recenter_second_enrch_bed/recentered_regions_file
      reference_point: 
        default: "TSS"                 # doesn't matter what we set here because we centered regions ourlselves
      before_region_start_length:
        default: 5000
      after_region_start_length:
        default: 5000
      bin_size:
        default: 10
      sort_regions:
        default: "descend"
      samples_label:
        source:
        - first_cond
        - second_cond
        valueFrom: $(["Reads " + self[0], "Reads " + self[1]])
      output_filename:
        default: "score_matrix.gz"
      missing_data_as_zero:
        default: true
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - scores_matrix
    - stdout_log
    - stderr_log

  make_heatmap:
    run: ../tools/deeptools-plotheatmap.cwl
    in:
      plot_title:
        default: "Tag density around peak centers"
      scores_matrix: compute_score_matrix/scores_matrix
      output_filename:
        default: "tag_density_heatmap.png"
      plot_type:
        default: "lines"
      sort_regions:
        default: "descend"
      average_type_summary_plot:
        default: "mean"
      what_to_show:
        default: "plot, heatmap and colorbar"
      ref_point_label:
        default: "Peak Center"
      regions_label:
        source:
        - first_cond
        - second_cond
        valueFrom: $(["Peaks " + self[0], "Peaks " + self[1]])
      samples_label:
        source:
        - first_cond
        - second_cond
        valueFrom: $(["Reads " + self[0], "Reads " + self[1]])
      x_axis_label:
        default: "distance (bp)"
      y_axisLabel:
        default: "Signal mean"
      per_group:
        default: false
      plot_file_format:
        default: "png"
      legend_location:
        default: "upper-left"
    out:
    - heatmap_file
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-Cell ATAC-Seq Differential Accessibility Analysis"
s:name: "Single-Cell ATAC-Seq Differential Accessibility Analysis"
s:alternateName: "Identifies differentially accessible regions between two groups of cells"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/sc-atac-dbinding.cwl
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
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Single-Cell ATAC-Seq Differential Accessibility Analysis

  Identifies differentially accessible regions between any
  two groups of cells, optionally aggregating chromatin
  accessibility data from single-cell to pseudo bulk form.