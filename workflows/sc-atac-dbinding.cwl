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
  - "sc-atac-cluster.cwl"
  - "sc-wnn-cluster.cwl"
  - "sc-ctype-assign.cwl"
  - "sc-rna-azimuth.cwl"
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
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  query_data_rds:
    type: File
    label: "Experiment run through any pipeline related Single-cell ATAC-Seq"
    doc: |
      Path to the RDS file to load Seurat object from.
      This file should include chromatin accessibility
      information stored in the ATAC assay. Additionally
      'rnaumap', and/or 'atacumap', and/or 'wnnumap'
      dimensionality reductions should be present.
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
    label: "Genome"
    doc: |
      Reference genome
    "sd:upstreamSource": "genome_indices/genome"
    "sd:localLabel": true

  datasets_metadata:
    type: File?
    label: "Optional TSV/CSV file to extend metadata by dataset"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat
      object metadata with categorical values using samples
      identities. First column - 'library_id' should
      correspond to all unique values from the 'new.ident'
      column of the loaded Seurat object. If any of the
      provided in this file columns are already present in
      the Seurat object metadata, they will be overwritten.
      When combined with --barcodes parameter, first the
      metadata will be extended, then barcode filtering will
      be applied. Default: no extra metadata is added

  barcodes_data:
    type: File?
    label: "Optional TSV/CSV file to prefilter and extend metadata by barcodes. First column should be named as 'barcode'"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and
      extend Seurat object metadata by selected barcodes.
      First column should be named as 'barcode'. If file
      includes any other columns they will be added to the
      Seurat object metadata ovewriting the existing ones if
      those are present. Default: all cells used, no extra
      metadata is added

  groupby:
    type: string?
    default: null
    label: "Category to group cells for optional subsetting"
    doc: |
      Column from the Seurat object metadata to group cells
      for optional subsetting when combined with --subset
      parameter. May be one of the extra metadata columns
      added with --metadata or --barcodes parameters.
      Ignored if --subset is not set. Default: do not
      subset, include all cells into analysis.

  subset:
    type: string?
    default: null
    label: "List of values to subset cells from the selected category"
    doc: |
      Values from the column set with --groupby parameter to
      subset cells before running differential binding
      analysis. Ignored if --groupby is not provided.
      Default: do not subset cells, include all of them.

  splitby:
    type: string
    label: "Category to split cell into two groups"
    doc: |
      Column from the Seurat object metadata to split cells
      into two groups to run --second vs --first
      differential binding analysis. May be one of the extra
      metadata columns added with --metadata or --barcodes
      parameters.

  first_cond:
    type: string
    label: "Value from the selected category to define the first group of cells"
    doc: |
      Value from the Seurat object metadata column set with
      --splitby parameter to define the first group of cells
      for differential binding analysis.

  second_cond:
    type: string
    label: "Value from the selected category to define the second group of cells"
    doc: |
      Value from the Seurat object metadata column set with
      --splitby parameter to define the second group of
      cells for differential binding analysis.

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
    label: "Test type to use in differential binding analysis"
    doc: |
      Test type to use in differential binding analysis. For
      all tests except manorm2, peaks present in the loaded
      Seurat object will be used. If manorm2 test selected,
      peaks will be called per group defined by --splitby
      parameter. Default: logistic-regression

  maximum_padj:
    type: float?
    default: 0.05
    label: "Maximum adjusted P-value to show in IGV"
    doc: |
      In the exploratory visualization part of the analysis
      output only differentially accessible regions with adjusted
      P-value not bigger than this value. Default: 0.05

  minimum_logfc:
    type: float?
    default: 1
    label: "Maximum log2 Fold Change value to show in IGV"
    doc: |
      In the exploratory visualization part of the analysis
      output only differentially accessible regions with log2 Fold
      Change not smaller than this value. Default: 1.0

  blacklist_regions_file:
    type: File?
    label: "Optional BED file with the genomic blacklist regions (for manorm2)"
    doc: |
      Path to the optional BED file with the genomic
      blacklist regions to be filtered out before running
      differential binding analysis. Any reference genomic
      bin overlapping a blacklist region will be removed
      from the output. Ignored if --test is not set to
      manorm2.

  minimum_qvalue:
    type: float?
    default: 0.05
    label: "Minimum FDR (q-value) cutoff for MACS2 peak detection (for manorm2)"
    doc: |
      Minimum FDR (q-value) cutoff for MACS2 peak detection.
      Ignored if --test is not set to manorm2. Default: 0.05
    "sd:layout":
      advanced: true

  minimum_peak_gap:
    type: int?
    default: 150
    label: "Minimum distabce between the peaks to be merged (for manorm2)"
    doc: |
      If a distance between peaks is smaller than the
      provided value they will be merged before splitting
      them into reference genomic bins of size --binsize.
      Ignored if --test is not set to manorm2. Default: 150
    "sd:layout":
      advanced: true

  bin_size:
    type: int?
    default: 1000
    label: "The size of non-overlapping reference genomic bins (for manorm2)"
    doc: |
      The size of non-overlapping reference genomic bins
      used by MAnorm2 when generating a table of reads
      counts per peaks. Ignored if --test is not set to
      manorm2. Default: 1000
    "sd:layout":
      advanced: true

  maximum_peaks:
    type: int?
    default: 0
    label: "The maximum number of the most significant peaks to keep, 0 - keep all (for manorm2)"
    doc: |
      The maximum number of the most significant (based on
      qvalue) peaks to keep from each group of cells when
      constructing reference genomic bins. Ignored if --test
      is not set to manorm2. Default: keep all peaks
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
    label: "Cells RNA UMAP split by selected criteria"
    doc: |
      Cells UMAP split by selected criteria,
      optionally subsetted to the specific
      group (rnaumap dim. reduction).
      PNG format
    "sd:visualPlugins":
    - image:
        tab: "Overall"
        Caption: "Cells RNA UMAP split by selected criteria"

  umap_rd_atacumap_plot_png:
    type: File?
    outputSource: sc_atac_dbinding/umap_rd_atacumap_plot_png
    label: "Cells ATAC UMAP split by selected criteria"
    doc: |
      Cells UMAP split by selected criteria,
      optionally subsetted to the specific
      group (atacumap dim. reduction).
      PNG format
    "sd:visualPlugins":
    - image:
        tab: "Overall"
        Caption: "Cells ATAC UMAP split by selected criteria"

  umap_rd_wnnumap_plot_png:
    type: File?
    outputSource: sc_atac_dbinding/umap_rd_wnnumap_plot_png
    label: "Cells WNN UMAP split by selected criteria"
    doc: |
      Cells UMAP split by selected criteria,
      optionally subsetted to the specific
      group (atacumap dim. reduction).
      PNG format
    "sd:visualPlugins":
    - image:
        tab: "Overall"
        Caption: "Cells WNN UMAP split by selected criteria"

  dbnd_vlcn_plot_png:
    type: File?
    outputSource: sc_atac_dbinding/dbnd_vlcn_plot_png
    label: "Volcano plot of differentially accessible regions"
    doc: |
      Volcano plot of differentially accessible regions.
      PNG format
    "sd:visualPlugins":
    - image:
        tab: "Overall"
        Caption: "Volcano plot of differentially accessible regions"

  seurat_peaks_bigbed_file:
    type: File
    outputSource: sc_atac_dbinding/seurat_peaks_bigbed_file
    label: "Peaks from the provided Seurat object"
    doc: |
      Peaks in bigBed format extracted
      from the loaded from provided RDS
      file Seurat object. 
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        format: "bigbed"
        name: "Seurat peaks"
        height: 40

  first_fragments_bigwig_file:
    type: File
    outputSource: sc_atac_dbinding/first_fragments_bigwig_file
    label: "Genome coverage for ATAC fragments (first)"
    doc: |
      Genome coverage in bigWig format calculated
      for ATAC fragments from the cells that belong to
      the group defined by the --first and
      --groupby parameters.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "ATAC fragments coverage (first)"
        height: 120

  second_fragments_bigwig_file:
    type: File
    outputSource: sc_atac_dbinding/second_fragments_bigwig_file
    label: "Genome coverage for ATAC fragments (second)"
    doc: |
      Genome coverage in bigWig format calculated
      for ATAC fragments from the cells that belong to
      the group defined by the --second and
      --groupby parameters.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "ATAC fragments coverage (second)"
        height: 120

  first_tn5ct_bigwig_file:
    type: File?
    outputSource: sc_atac_dbinding/first_tn5ct_bigwig_file
    label: "Genome coverage for Tn5 cut sites (first)"
    doc: |
      Genome coverage in bigWig format calculated
      for Tn5 cut sites from the cells that belong
      to the group defined by the --first and
      --groupby parameters.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "Tn5 coverage (first)"
        height: 120

  second_tn5ct_bigwig_file:
    type: File?
    outputSource: sc_atac_dbinding/second_tn5ct_bigwig_file
    label: "Genome coverage for Tn5 cut sites (second)"
    doc: |
      Genome coverage in bigWig format calculated
      for Tn5 cut sites from the cells that belong
      to the group defined by the --second and
      --groupby parameters.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "Tn5 coverage (second)"
        height: 120

  first_peaks_xls_file:
    type: File?
    outputSource: sc_atac_dbinding/first_peaks_xls_file
    label: "MACS2 report in XLS format (first)"
    doc: |
      MACS2 report in XLS format for peaks
      called from the Tn5 cut sites of the
      cells that belong to the group defined
      by the --first and --groupby parameters.

  second_peaks_xls_file:
    type: File?
    outputSource: sc_atac_dbinding/second_peaks_xls_file
    label: "MACS2 report in XLS format (second)"
    doc: |
      MACS2 report in XLS format for peaks
      called from the Tn5 cut sites of the
      cells that belong to the group defined
      by the --second and --groupby parameters.

  first_peaks_bed_file:
    type: File?
    outputSource: sc_atac_dbinding/first_peaks_bed_file
    label: "MACS2 peaks in narrowPeak format (first)"
    doc: |
      MACS2 peaks in narrowPeak format called
      from the Tn5 cut sites of the cells that
      belong to the group defined by the --first
      and --groupby parameters.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "Called peaks (first)"
        displayMode: "COLLAPSE"
        height: 40

  second_peaks_bed_file:
    type: File?
    outputSource: sc_atac_dbinding/second_peaks_bed_file
    label: "MACS2 peaks in narrowPeak format (second)"
    doc: |
      MACS2 peaks in narrowPeak format called
      from the Tn5 cut sites of the cells that
      belong to the group defined by the --second
      and --groupby parameters.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "Called peaks (second)"
        displayMode: "COLLAPSE"
        height: 40

  first_summits_bed_file:
    type: File?
    outputSource: sc_atac_dbinding/first_summits_bed_file
    label: "MACS2 peaks summits in BED format (first)"
    doc: |
      MACS2 peaks summits in BED format called
      from the Tn5 cut sites of the cells that
      belong to the group defined by the --first
      and --groupby parameters.

  second_summits_bed_file:
    type: File?
    outputSource: sc_atac_dbinding/second_summits_bed_file
    label: "MACS2 peaks summits in BED format (second)"
    doc: |
      MACS2 peaks summits in BED format called
      from the Tn5 cut sites of the cells that
      belong to the group defined by the --second
      and --groupby parameters.

  diff_bound_sites:
    type: File
    outputSource: sc_atac_dbinding/diff_bound_sites
    label: "Differentially accessible regions"
    doc: |
      Not filtered differentially accessible regions
      in TSV format
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Diff. accessible regions"
        Title: "Differentially accessible regions. Not filtered"

  diff_bound_sites_with_labels:
    type: File
    outputSource: add_label_column/output_file
    label: "Differentially accessible regions with labels"
    doc: |
      Not filtered differentially accessible regions
      with labels in TSV format
    'sd:visualPlugins':
    - queryRedirect:
        tab: "Overview"
        label: "Volcano Plot"
        url: "https://scidap.com/vp/volcano"
        query_eval_string: "`data_file=${this.getSampleValue('outputs', 'diff_bound_sites_with_labels')}&data_col=label&x_col=log2FoldChange&y_col=padj`"

  first_enrch_bigbed_file:
    type: File?
    outputSource: sc_atac_dbinding/first_enrch_bigbed_file
    label: "Significant differentially accessible regions (first)"
    doc: |
      Peaks in bigBed format filtered by
      --padj and --logfc thresholds enriched
      in the group of cells defined by the
      --first and --groupby parameters.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        format: "bigbed"
        name: "Diff. accessible regions (first)"
        height: 40

  second_enrch_bigbed_file:
    type: File?
    outputSource: sc_atac_dbinding/second_enrch_bigbed_file
    label: "Significant differentially accessible regions (second)"
    doc: |
      Peaks in bigBed format filtered by
      --padj and --logfc thresholds enriched
      in the group of cells defined by the
      --second and --groupby parameters.
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        format: "bigbed"
        name: "Diff. accessible regions (second)"
        height: 40

  first_enrch_bed_file:
    type: File?
    outputSource: sc_atac_dbinding/first_enrch_bed_file
    label: "Significant differentially accessible regions (first)"
    doc: |
      Peaks in BED format filtered by
      --padj and --logfc thresholds enriched
      in the group of cells defined by the
      --first and --groupby parameters.

  second_enrch_bed_file:
    type: File?
    outputSource: sc_atac_dbinding/second_enrch_bed_file
    label: "Significant differentially accessible regions (second)"
    doc: |
      Peaks in BED format filtered by
      --padj and --logfc thresholds enriched
      in the group of cells defined by the
      --second and --groupby parameters.

  tag_density_matrix:
    type: File
    outputSource: compute_score_matrix/scores_matrix
    label: "Score matrix for tag density heatmap"
    doc: |
      Scores matrix generated by
      Deeptools with tag density
      information around centers
      of regions of interest.

  tag_density_heatmap:
    type: File
    outputSource: make_heatmap/heatmap_file
    label: "Tag density heatmap"
    doc: |
      Tag density heatmap around centers
      of differentially accessible regions in
      PNG format
    "sd:visualPlugins":
    - image:
        tab: "Overall"
        Caption: "Tag density heatmap around centers of diff. accessible regions"

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
        valueFrom: $(self==""?null:self)                # safety measure
      subset:
        source: subset
        valueFrom: $(split_by_comma(self))
      splitby: splitby
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
    - first_tn5ct_bigwig_file
    - second_tn5ct_bigwig_file
    - first_peaks_xls_file
    - second_peaks_xls_file
    - first_peaks_bed_file
    - second_peaks_bed_file
    - first_summits_bed_file
    - second_summits_bed_file
    - diff_bound_sites
    - dbnd_vlcn_plot_png
    - first_enrch_bigbed_file
    - second_enrch_bigbed_file
    - first_enrch_bed_file
    - second_enrch_bed_file
    - umap_rd_rnaumap_plot_pdf
    - umap_rd_atacumap_plot_pdf
    - umap_rd_wnnumap_plot_pdf
    - dbnd_vlcn_plot_pdf
    - sc_report_html_file
    - stdout_log
    - stderr_log

  folder_pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - sc_atac_dbinding/umap_rd_rnaumap_plot_pdf
        - sc_atac_dbinding/umap_rd_atacumap_plot_pdf
        - sc_atac_dbinding/umap_rd_wnnumap_plot_pdf
        - sc_atac_dbinding/dbnd_vlcn_plot_pdf
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

  add_label_column:
    run: ../tools/custom-bash.cwl
    in:
      input_file: sc_atac_dbinding/diff_bound_sites
      script:
        default: |
          HEADER=`head -n 1 $0`;
          echo -e "label\t${HEADER}" > diff_sts_labeled.tsv;
          cat "$0" | grep -v "start" | awk -F "\t" '{print $1":"$2"-"$3"-"$NF"\t"$0}' >> diff_sts_labeled.tsv
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
        source:
        - analysis_method
        - sc_atac_dbinding/first_fragments_bigwig_file
        - sc_atac_dbinding/second_fragments_bigwig_file
        - sc_atac_dbinding/first_tn5ct_bigwig_file
        - sc_atac_dbinding/second_tn5ct_bigwig_file
        valueFrom: |
          ${
              if (self[0].analysis_method != "manorm2" ) {
                return [self[1], self[2]];
              } else {
                return [self[3], self[4]];
              }
          }
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
  accessibility data from single-cell to pseudobulk form.