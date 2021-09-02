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
    - var parse_splitby = function(line) {
          return (line == "dataset")?"new.ident":line;
      };
    - var parse_resolution = function(line) {
          return "integrated_snn_res."+line;
      };


'sd:upstream':
  seurat_cluster_sample:
  - "seurat-cluster.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  seurat_data_rds:
    type: File
    label: "Seurat Cluster Experiment"
    doc: |
      Path to the RDS file to load Seurat object from.
      RDS file can be produced by run_seurat.R script.
    'sd:upstreamSource': "seurat_cluster_sample/seurat_clst_data_rds"
    'sd:localLabel': true

  splitby:
    type:
    - "null"
    - type: enum
      symbols:
      - "condition"
      - "dataset"
    default: "condition"
    label: "Divide cell based on"
    doc: |
      Field from the Seurat object metadata to split cells into groups
      for differential expression analysis.

  first_cond:
    type: string
    label: "First group of cells"
    doc: |
      Value from the column set with --splitby to define a first group of cells.

  second_cond:
    type: string
    label: "Second group of cells"
    doc: |
      Value from the column set with --splitby to define a second group of cells.

  resolution:
    type: string
    label: "Clustering resolution to subset cells"
    doc: |
      Clustering resolution to subset cells. Will be used to define a dield from
      the Seurat object metadata to group cells for subsetting.

  selected_clusters:
    type: string
    label: "Comma or space separated list of clusters to subset cells"
    doc: |
      Value(s) from the column set with --groupby (inferred from resolution) to
      subset cells before running differential expression analysis.

  selected_features:
    type: string?
    default: null
    label: "Comma or space separated list of genes of interest"
    doc: |
      Genes of interest to label on the generated plots.
      Default: --topn N genes with the highest and the
      lowest log2 fold change expression values.
    'sd:layout':
      advanced: true

  excluded_features:
    type: string?
    default: null
    label: "Comma or space separated list of genes to be excluded"
    doc: |
      Genes to be excluded from the differential expression analysis.
      Excluded genes will be still present in the dataset, but they won't
      be used in the FindMarkers function.
      Default: include all genes
    'sd:layout':
      advanced: true

  minimum_logfc:
    type: float?
    default: 0.25
    label: "Include only those genes that on average have the absolute value of log2 fold change expression difference not lower than this value"
    doc: |
      Include only those genes that on average have the absolute value of log2
      fold change expression difference not lower than this value.
    'sd:layout':
      advanced: true

  minimum_pct:
    type: float?
    default: 0.1
    label: "Include only those genes that are detected in not lower than this fraction of cells in either of the two tested groups"
    doc: |
      Include only those genes that are detected in not lower than this
      fraction of cells in either of the two tested groups.
    'sd:layout':
      advanced: true

  maximum_pvadj:
    type: float?
    default: 0.05
    label: "Include only those genes for which adjusted P-val is not bigger that this value"
    doc: |
      Include only those genes for which adjusted P-val is not bigger that this value.
    'sd:layout':
      advanced: true

  test_use:
    type:
    - "null"
    - type: enum
      symbols:
      - "wilcox"
      - "bimod"
      - "roc"
      - "t"
      - "negbinom"
      - "poisson"
      - "LR"
      - "MAST"
      - "DESeq2"
    default: "wilcox"
    label: "Statistical test to use for differential gene expression analysis"
    doc: |
      Statistical test to use for differential gene expression analysis.
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 2
    label: "Threads number to use"
    doc: |
      Threads number
    'sd:layout':
      advanced: true


outputs:

  avg_gene_expr_plot_png:
    type: File?
    outputSource: sc_diff_expr/avg_gene_expr_plot_png
    label: "Log normalized average gene expression"
    doc: |
      Log normalized average gene expression for first_cond vs second_cond cells split
      by criteria set in splitby (a.k.a condition). Cells are optionally subsetted by
      selected_groups (a.k.a clusters) from the groups defined in groupby.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression plots'
        Caption: 'Log normalized average gene expression'

  avg_gene_expr_plot_pdf:
    type: File?
    outputSource: sc_diff_expr/avg_gene_expr_plot_pdf
    label: "Log normalized average gene expression"
    doc: |
      Log normalized average gene expression for first_cond vs second_cond cells split
      by criteria set in splitby (a.k.a condition). Cells are optionally subsetted by
      selected_groups (a.k.a clusters) from the groups defined in groupby.
      PDF format

  diff_expr_genes_plot_png:
    type: File?
    outputSource: sc_diff_expr/diff_expr_genes_plot_png
    label: "Differentially expressed genes"
    doc: |
      Volcano plot of differentially expressed genes for first_cond vs second_cond cells
      split by criteria set in splitby (a.k.a condition). Cells are optionally subsetted
      by selected_groups (a.k.a clusters) from the groups defined in groupby.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression plots'
        Caption: 'Differentially expressed genes'

  diff_expr_genes_plot_pdf:
    type: File?
    outputSource: sc_diff_expr/diff_expr_genes_plot_pdf
    label: "Differentially expressed genes"
    doc: |
      Volcano plot of differentially expressed genes for first_cond vs second_cond cells
      split by criteria set in splitby (a.k.a condition). Cells are optionally subsetted
      by selected_groups (a.k.a clusters) from the groups defined in groupby.
      PDF format

  diff_expr_genes:
    type: File
    outputSource: sc_diff_expr/diff_expr_genes
    label: "Differentially expressed genes"
    doc: |
      Differentially expressed genes for first_cond vs second_cond cells split by criteria
      set in splitby (a.k.a condition). Cells are optionally subsetted by selected_groups
      (a.k.a clusters) from the groups defined in groupby.
      TSV format
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Diff expressed genes'
        Title: 'Differentially expressed genes'

  sc_diff_expr_stdout_log:
    type: File
    outputSource: sc_diff_expr/stdout_log
    label: stdout log generated by Seurat Differential Expression Analysis
    doc: |
      stdout log generated by Seurat Differential Expression Analysis

  sc_diff_expr_stderr_log:
    type: File
    outputSource: sc_diff_expr/stderr_log
    label: stderr log generated by Seurat Differential Expression Analysis
    doc: |
      stderr log generated by Seurat Differential Expression Analysis


steps:

  sc_diff_expr:
    run: ../tools/sc_diff_expr.cwl
    in:
      seurat_data_rds: seurat_data_rds
      splitby:
        source: splitby
        valueFrom: $(parse_splitby(self))
      first_cond: first_cond
      second_cond: second_cond
      groupby:
        source: resolution
        valueFrom: $(parse_resolution(self))
      selected_groups:
        source: selected_clusters
        valueFrom: $(split_features(self))
      topn_genes_count:
        default: 10
      selected_features:
        source: selected_features
        valueFrom: $(split_features(self))
      excluded_features:
        source: excluded_features
        valueFrom: $(split_features(self))
      minimum_logfc: minimum_logfc
      minimum_pct: minimum_pct
      maximum_pvadj: maximum_pvadj
      test_use: test_use
      export_pdf_plots:
        default: true
      threads: threads
    out:
    - avg_gene_expr_plot_png
    - avg_gene_expr_plot_pdf
    - diff_expr_genes_plot_png
    - diff_expr_genes_plot_pdf
    - diff_expr_genes
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Seurat Differential Expression"
s:name: "Seurat Differential Expression"
s:alternateName: "Runs differential expression analysis between two biological conditions for a group of cells"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/sc_diff_expr.cwl
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


doc: |
  Seurat Differential Expression
  ==============================

  Runs differential expression analysis between two biological conditions for a group of cells.