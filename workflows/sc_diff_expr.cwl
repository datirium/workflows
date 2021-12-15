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
  - "https://github.com/datirium/workflows/workflows/seurat-cluster.cwl"
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
      Column from the Seurat object metadata to split cells into two groups
      to run second_cond vs first_cond differential expression analysis. May include
      columns from the metadata fields added with conditions_data.

  first_cond:
    type: string
    label: "First group of cells"
    doc: |
      Value from the Seurat object metadata column set with splitby to define the
      first group of cells or pseudobulk RNA-Seq samples (when using pseudo).

  second_cond:
    type: string
    label: "Second group of cells"
    doc: |
      Value from the Seurat object metadata column set with splitby to define the
      the second group of cells or pseudobulk RNA-Seq samples (when using pseudo).

  resolution:
    type: string
    label: "Clustering resolution to subset cells"
    doc: |
      Clustering resolution to subset cells. Will be used to define a field from
      the Seurat object metadata to group cells for subsetting.

  selected_clusters:
    type: string
    label: "Comma or space separated list of clusters to subset cells"
    doc: |
      Value(s) from the column set with groupby (inferred from resolution) to
      subset cells before running differential expression analysis.

  selected_features:
    type: string?
    default: null
    label: "Comma or space separated list of genes of interest"
    doc: |
      Genes of interest to label on the generated plots.
      Default: 10 genes with the highest and the
      lowest log2 fold change expression values.
    'sd:layout':
      advanced: true

  excluded_features:
    type: string?
    default: null
    label: "Comma or space separated list of genes to be excluded"
    doc: |
      Genes to be excluded from the differential expression analysis.
    'sd:layout':
      advanced: true

  minimum_logfc:
    type: float?
    default: 0.25
    label: "Include only those genes that on average have the absolute value of log2 fold change expression difference not lower than this value"
    doc: |
      Include only those genes that on average have the absolute value of log2
      fold change expression difference not lower than this value. Increasing
      minimum_logfc speeds up calculations, but can cause missing weaker signals.
      Ignored with pseudo.
    'sd:layout':
      advanced: true

  minimum_pct:
    type: float?
    default: 0.1
    label: "Include only those genes that are detected in not lower than this fraction of cells in either of the two tested groups"
    doc: |
      Include only those genes that are detected in not lower than this fraction of cells
      in either of the two tested groups. Increasing minimum_pct speeds up calculations by not
      testing genes that are very infrequently expressed. Ignored with pseudo.
    'sd:layout':
      advanced: true

  maximum_pvadj:
    type: float?
    default: 0.1
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
      Ignored with pseudo.
    'sd:layout':
      advanced: true

  batchby:
    type: string?
    default: null
    label: "Column from the metadata to define the variable that should be modelled as a batch effect"
    doc: |
      Column from the Seurat object metadata to define the variable that should
      be modelled as a batch effect when running differential expression analysis.
      Applied only when test_use is one of 'LR', 'negbinom', 'poisson', or 'MAST',
      or when using pseudo. May include columns from the metadata fields added
      with conditions_data. Values selected from the column set with batchby should
      establish 1:1 relation with the 'new.ident' column of the Seurat object loaded
      from seurat_data_rds.
    'sd:layout':
      advanced: true

  pseudo:
    type: boolean?
    default: false
    label: "Aggregate gene expression of the cells from the same dataset into a pseudobulk RNA-Seq sample"
    doc: |
      Aggregate gene expression of the cells from the same dataset into a pseudobulk
      RNA-Seq sample before running differential expression analysis with DESeq2.
      The following parameters will be ignored: test_use, minimum_pct, minimum_logfc.
    'sd:layout':
      advanced: true

  lrt:
    type: boolean?
    default: false
    label: "Use LRT instead of the pair-wise Wald test"
    doc: |
      Use LRT instead of the pair-wise Wald test. Shows any differences across the variable
      set with batchby whith the log2 fold changes calculated as the average expression
      changes due to criteria set with splitby. Ignored when pseudo or batchby
      parameters are not provided.
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

  conditions_data:
    type: File?
    label: "TSV/CSV file to optionally extend metadata"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat object metadata. First
      column 'library_id' should include all unique values from the 'new.ident'
      column of the loaded from seurat_data_rds object metadata. All other columns will
      be added to the Seurat object metadata. If any of the provided in this file
      columns were already present in the Seurat object metadata, they will be
      overwritten.
    'sd:layout':
      advanced: true


outputs:

  cell_abundance_plot_png:
    type: File?
    outputSource: sc_diff_expr/cell_abundance_plot_png
    label: "Cell abundance"
    doc: |
      Cell abundance plot split by criteria set in splitby (a.k.a condition) and optionally
      subsetted by selected_groups (a.k.a clusters) from the groups defined in groupby.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression plots'
        Caption: 'Cell abundance plot'

  cell_abundance_plot_pdf:
    type: File?
    outputSource: sc_diff_expr/cell_abundance_plot_pdf
    label: "Cell abundance"
    doc: |
      Cell abundance plot split by criteria set in splitby (a.k.a condition) and optionally
      subsetted by selected_groups (a.k.a clusters) from the groups defined in groupby.
      PDF format

  aggr_gene_expr_plot_png:
    type: File?
    outputSource: sc_diff_expr/aggr_gene_expr_plot_png
    label: "Log normalized aggregated gene expression"
    doc: |
      Log normalized aggregated gene expression split by criteria set in splitby
      (a.k.a condition).
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Gene expression plots'
        Caption: 'Log normalized aggregated gene expression'

  aggr_gene_expr_plot_pdf:
    type: File?
    outputSource: sc_diff_expr/aggr_gene_expr_plot_pdf
    label: "Log normalized aggregated gene expression"
    doc: |
      Log normalized aggregated gene expression split by criteria set in splitby
      (a.k.a condition).
      PDF format

  diff_expr_genes_plot_png:
    type: File?
    outputSource: sc_diff_expr/diff_expr_genes_plot_png
    label: "Differentially expressed genes"
    doc: |
      Volcano plot of differentially expressed genes for second_cond vs first_cond cells
      or pseudobulk RNA-Seq samples split by criteria set in splitby (a.k.a condition)
      and optionally subsetted by selected_clusters from the groups defined in groupby.
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
      Volcano plot of differentially expressed genes for second_cond vs first_cond cells
      or pseudobulk RNA-Seq samples split by criteria set in splitby (a.k.a condition)
      and optionally subsetted by selected_clusters from the groups defined in groupby.
      PDF format

  diff_expr_genes:
    type: File
    outputSource: sc_diff_expr/diff_expr_genes
    label: "Differentially expressed genes"
    doc: |
      Differentially expressed genes for second_cond vs first_cond cells or pseudobulk
      RNA-Seq samples split by criteria set in splitby (a.k.a condition) and optionally
      subsetted by selected_clusters from the groups defined in groupby.
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
      conditions_data: conditions_data
      splitby:
        source: splitby
        valueFrom: $(parse_splitby(self))
      first_cond: first_cond
      second_cond: second_cond
      batchby: batchby
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
      pseudo: pseudo
      lrt: lrt
      export_pdf_plots:
        default: true
      threads: threads
    out:
    - cell_abundance_plot_png
    - cell_abundance_plot_pdf
    - aggr_gene_expr_plot_png
    - aggr_gene_expr_plot_pdf
    - diff_expr_genes_plot_png
    - diff_expr_genes_plot_pdf
    - diff_expr_genes
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Single-cell Differential Expression"
s:name: "Single-cell Differential Expression"
s:alternateName: "Runs differential expression analysis for a subset of cells between two selected conditions"

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
  Single-cell Differential Expression
  ===================================

  Runs differential expression analysis for a subset of cells between two selected conditions.