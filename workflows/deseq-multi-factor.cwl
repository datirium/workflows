cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var split_by_common_delim = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          let splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };


'sd:upstream':
  rnaseq_experiment:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"
    - "trim-quantseq-mrnaseq-se-strand-specific.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name"
    sd:preview:
      position: 1

  isoforms_expression_files:
    type:
    - "null"
    - File[]
    label: "RNA-Seq experiments"
    doc: |
      Path to the TSV/CSV files with expression data.
      All files should have the following header:
      RefseqId GeneId Chrom TxStart TxEnd Strand TotalReads Rpkm
    'sd:upstreamSource': "rnaseq_experiment/rpkm_isoforms"
    'sd:localLabel': true

  genes_expression_files:
    type:
    - "null"
    - File[]
    label: "RNA-Seq experiments"
    doc: |
      Path to the TSV/CSV files with expression data.
      All files should have the following header:
      RefseqId GeneId Chrom TxStart TxEnd Strand TotalReads Rpkm
    'sd:upstreamSource': "rnaseq_experiment/rpkm_genes"
    'sd:localLabel': true

  expression_names:
    type: string[]
    label: "RNA-Seq experiments"
    doc: |
      Unique names for files provided in --expression,
      no special characters or spaces are allowed.
      Number and order of the names should corresponds
      to values from --expression.
    'sd:upstreamSource': "rnaseq_experiment/alias"

  feature_type:
    type:
      - "null"
      - type: enum
        symbols:
        - "transcript"
        - "gene"
    default: "gene"
    label: "Feature type to use for differential expression"
    doc: |
      Feature type to use for differential expression.
      If set to 'gene', use 'GeneId' column from the provided in --expression files.
      If set to 'transcript', use 'RefseqId' from the provided in --expression files.
      Default: gene

  design_formula:
    type: string
    label: "Design formula"
    doc: |
      Design formula. Should start with ~ and include terms from
      the --metadata table.

  reduced_formula:
    type: string?
    default: "~ 1"
    label: "Reduced formula. If removed, force to use Wald instead of LRT."
    doc: |
      Reduced formula with the term(s) of interest removed.
      Should start with ~. If design formula includes only
      one criteria, use ~ 1. If removed, forces DESeq2 to run
      Wald test instead of the LRT.

  contrast:
    type: string?
    label: "Contrast. If not provided, use the last term from the design formula."
    doc: |
      Contrast to be be applied for the output, formatted as
      a mathematical formula of values from the --metadata table.
      If not provided, the last term from the design formula will
      be used.

  remove:
    type: string?
    label: "Column from the metadata file to remove batch effect"
    doc: |
      Column from the metadata file to remove batch effect
      before running differential expression analysis. If
      present, all components that include this term will be
      removed from the design and reduced formulas.
      Default: do not remove batch effect

  base:
    type: string?
    label: "Values from each column of the metadata file to be set as base levels. Order matters."
    doc: |
      Value(s) from each metadata file column(s) to be set as
      the base level(s). Number and order of provided values should
      correspond the order of columns in --metadata file. Default:
      define base levels alphabetically for each metadata column.

  metadata_file:
    type: File
    label: "Metadata file to assign categories to datasets"
    doc: |
      Path to the TSV/CSV file to provide metadata for the
      samples from --expression. First column should have
      the name 'sample', other columns may have arbitrary names.
      The values from the 'sample' column should correspond to
      the values provided in --aliases. For a proper --contrast
      intepretation, values defined in each column should not be
      used in other columns. All metadata columns are treated as
      factors (no covariates are supported).

  normalization_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "vst"
      - "rlog"
    default: "rlog"
    label: "Read counts normalization for the exploratory visualization analysis"
    doc: |
      Read counts normalization for the exploratory visualization analysis.
      Use 'vst' for medium-to-large datasets (n > 30) and 'rlog' for
      small datasets (n < 30), when there is a wide range of sequencing
      depth across samples.
      Default: rlog
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
    label: "Hopach clustering method to be run on normalized read counts"
    doc: |
      Hopach clustering method to be run on normalized read counts for the
      exploratory visualization analysis. Default: do not run clustering
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
    label: "Distance metric for HOPACH row clustering"
    doc: |
      Distance metric for HOPACH row clustering. Ignored if --cluster is not
      provided. Default: cosangle
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
    label: "Distance metric for HOPACH column clustering"
    doc: |
      Distance metric for HOPACH column clustering. Ignored if --cluster is not
      provided. Default: euclid
    'sd:layout':
      advanced: true

  selected_features:
    type: string?
    label: "Features of interest to label on the volcano plot"
    doc: |
      Features of interest to label on the generated volcanot plot. Default:
      top 10 features with the highest and the lowest log2 fold change
      expression values.
    'sd:layout':
      advanced: true

  excluded_features:
    type: string?
    label: "Features to be excluded from the differential expression analysis"
    doc: |
      Features to be excluded from the differential expression analysis.
      Default: include all features
    'sd:layout':
      advanced: true

  maximum_padj:
    type: float?
    default: 0.05
    label: "Maximum P-adjusted to show features in the exploratory visualization analysis"
    doc: |
      In the exploratory visualization analysis output only features with
      adjusted P-value not bigger than this value. Default: 0.05
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
    label: "Cores/CPUs"
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 1
    "sd:layout":
      advanced: true


outputs:

  diff_expr_file:
    type: File
    outputSource: deseq_multi_factor/diff_expr_features
    label: "TSV file with not filtered differentially expressed features"
    doc: |
      TSV file with not filtered differentially
      expressed features
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'DE features'
        Title: 'Differentially expressed features'

  read_counts_gct:
    type: File
    outputSource: deseq_multi_factor/read_counts_gct
    label: "GCT file with normalized, optionally batch corrected, read counts"
    doc: |
      GCT file with normalized, optionally batch
      corrected, read counts

  mds_plot_html:
    type: File?
    outputSource: deseq_multi_factor/mds_plot_html
    label: "MDS plot of normalized, optionally batch corrected, read counts"
    doc: |
      MDS plot of normalized, optionally batch
      corrected, read counts.
      HTML format
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  volcano_plot_png:
    type: File?
    outputSource: deseq_multi_factor/volcano_plot_png
    label: "Volcano plot of differentially expressed features"
    doc: |
      Volcano plot of differentially expressed
      features.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'Volcano plot of differentially expressed features'

  pca_plot_png:
    type: File?
    outputSource: deseq_multi_factor/pca_plot_png
    label: "PCA plot of normalized, optionally batch corrected, read counts"
    doc: |
      PCA plot of normalized, optionally batch
      corrected, read counts based on the top
      500 features selected by the highest row
      variance
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'PCA plot of normalized, optionally batch corrected, read counts'

  volcano_plot_html_file:
    type: File
    outputSource: make_volcano_plot/html_file
    label: "Volcano Plot"
    doc: |
      HTML index file for Volcano Plot
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  volcano_plot_html_data:
    type: Directory
    outputSource: make_volcano_plot/html_data
    label: "Directory html data for Volcano Plot"
    doc: |
      Directory html data for Volcano Plot

  ma_plot_html_file:
    type: File
    outputSource: make_ma_plot/html_file
    label: "MA-plot"
    doc: |
      HTML index file for MA-plot
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  ma_plot_html_data:
    type: Directory
    outputSource: make_ma_plot/html_data
    label: "Directory html data for Volcano Plot"
    doc: |
      Directory html data for MA-plot

  heatmap_html:
    type: File
    outputSource: morpheus_heatmap/heatmap_html
    label: "Heatmap of normalized counts"
    doc: |
      Morpheus heatmap in HTML format
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  deseq_stdout_log:
    type: File
    outputSource: deseq_multi_factor/stdout_log
    label: "DESeq stdout log"
    doc: "DESeq stdout log"

  deseq_stderr_log:
    type: File
    outputSource: deseq_multi_factor/stderr_log
    label: "DESeq stderr log"
    doc: "DESeq stderr log"

  morpheus_stdout_log:
    type: File
    outputSource: morpheus_heatmap/stdout_log
    label: "Morpheus heatmap stdout log"
    doc: "Morpheus heatmap stdout log"

  morpheus_stderr_log:
    type: File
    outputSource: morpheus_heatmap/stderr_log
    label: "Morpheus heatmap stderr log"
    doc: "Morpheus heatmap stderr log"


steps:

  deseq_multi_factor:
    run: ../tools/deseq-multi-factor.cwl
    in:
      expression_files:
        source: [feature_type, isoforms_expression_files, genes_expression_files]
        valueFrom: |
          ${
              if (self[0] == "transcript") {
                return self[1];
              } else {
                return self[2];
              }
          }
      expression_names: expression_names
      metadata_file: metadata_file
      design_formula: design_formula
      reduced_formula:
        source: reduced_formula
        valueFrom: $(self==""?null:self)            # safety measure
      contrast: 
        source: contrast
        valueFrom: $(self==""?null:self)            # safety measure
      base:
        source: base
        valueFrom: $(split_by_common_delim(self))
      feature_type: feature_type
      excluded_features:
        source: excluded_features
        valueFrom: $(split_by_common_delim(self))
      normalization_method: normalization_method
      remove:
        source: remove
        valueFrom: $(self==""?null:self)            # safety measure
      cluster_method:
        source: cluster_method
        valueFrom: $(self=="none"?null:self)
      row_distance: row_distance
      column_distance: column_distance
      center_row:
        default: true
      selected_features:
        source: selected_features
        valueFrom: $(split_by_common_delim(self))
      maximum_padj: maximum_padj
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - diff_expr_features
    - read_counts_gct
    - volcano_plot_png
    - pca_plot_png
    - mds_plot_html
    - stdout_log
    - stderr_log

  make_volcano_plot:
    run: ../tools/volcano-plot.cwl
    in:
      diff_expr_file: deseq_multi_factor/diff_expr_features
      x_axis_column:
        default: "log2FoldChange"
      y_axis_column:
        default: "padj"
      label_column:
        default: "feature"
    out:
    - html_data
    - html_file

  make_ma_plot:
    run: ../tools/ma-plot.cwl
    in:
      diff_expr_file: deseq_multi_factor/diff_expr_features
      x_axis_column:
        default: "baseMean"
      y_axis_column:
        default: "log2FoldChange"
      label_column:
        default: "feature"
    out:
    - html_data
    - html_file

  morpheus_heatmap:
    run: ../tools/morpheus-heatmap.cwl
    in:
     read_counts_gct: deseq_multi_factor/read_counts_gct
    out:
    - heatmap_html
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "DESeq2 Multi-factor Analysis"
label: "DESeq2 Multi-factor Analysis"
s:alternateName: "Runs DeSeq2 multi-factor analysis with manual control over major parameters"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-scidap/master/workflows/deseq-multi-factor.cwl
s:codeRepository: https://github.com/Barski-lab/workflows-scidap
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
  DESeq2 Multi-factor Analysis

  Runs DeSeq2 multi-factor analysis with manual control over major parameters