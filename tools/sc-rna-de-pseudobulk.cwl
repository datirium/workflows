cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.15


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from. This file should include genes
      expression information stored in the RNA assay. Additionally, 'rnaumap', and/or
      'atacumap', and/or 'wnnumap' dimensionality reductions should be present.

  datasets_metadata:
    type: File?
    inputBinding:
      prefix: "--metadata"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat object metadata with
      categorical values using samples identities. First column - 'library_id'
      should correspond to all unique values from the 'new.ident' column of the
      loaded Seurat object. If any of the provided in this file columns are already
      present in the Seurat object metadata, they will be overwritten. Default: no
      extra metadata is added

  splitby:
    type: string
    inputBinding:
      prefix: "--splitby"
    doc: |
      Column from the Seurat object metadata to split datasets into two groups
      to run --second vs --first pseudobulk DE analysis, i.e., calculate log2FC.
      May be one of the columns from the extra metadata added with --metadata
      parameter. Provided value should group the datasets, not cells, therefore
      do not use a column with clustering results.

  first_cond:
    type: string
    inputBinding:
      prefix: "--first"
    doc: |
      Value from the Seurat object metadata column set with --splitby to define the
      first group of datasets for pseudobulk DE analysis.

  second_cond:
    type: string
    inputBinding:
      prefix: "--second"
    doc: |
      Value from the Seurat object metadata column set with --splitby to define the
      second group of datasets for pseudobulk DE analysis.

  batchby:
    type: string?
    inputBinding:
      prefix: "--batchby"
    doc: |
      Column from the Seurat object metadata to group datasets into batches. It will be used
      as a factor variable to model batch effect when running pseudobulk DE analysis (makes
      design formula look like ~splitby+batchby). May be one of the columns from the extra
      metadata added with --metadata parameter. Provided value should batch the datasets, not
      cells, therefore do not use a column with clustering results. Default: do not model
      batch effect.

  groupby:
    type: string?
    inputBinding:
      prefix: "--groupby"
    doc: |
      Column from the Seurat object metadata to group cells for optional subsetting
      when combined with --subset parameter. May be one of the columns from the extra
      metadata added with --metadata parameter. Ignored if --subset is not set. Provided
      value defines the groups of cells, therefore any metadata column, including the
      clustering results, may be used. Default: do not subset, run pseudobulk DE analysis
      for all cells jointly

  subset:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--subset"
    doc: |
      Value(s) from the column set with --groupby parameter to subset cells
      before running pseudobulk DE analysis. If multiple values are provided
      run analysis jointly for selected groups of cells. Ignored if --groupby
      is not set. Default: do not subset, run pseudobulk DE analysis for all
      cells jointly

  lrt:
    type: boolean?
    inputBinding:
      prefix: "--lrt"
    doc: |
      Use LRT instead of the pair-wise Wald test. If --batchby is not provided
      use ~1 as a reduced formula, otherwise ~batchby. Default: use Wald test

  maximum_padj:
    type: float?
    inputBinding:
      prefix: "--padj"
    doc: |
      In the exploratory visualization part of the analysis output only features
      with adjusted P-value not bigger than this value. Default: 0.05

  genes_of_interest:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--genes"
    doc: |
      Genes of interest to label on the generated plots. Default: top 10 genes
      with the highest and the lowest log2FC expression values.

  exclude_pattern:
    type: string?
    inputBinding:
      prefix: "--exclude"
    doc: |
      Regex pattern to identify and exclude non-coding RNA genes from the pseudobulk
      DE analysis (not case-sensitive). If any of such genes were provided in the --genes
      parameter, they will be excluded from there as well.
      Default: use all genes

  normalization_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "vst"
      - "rlog"
    inputBinding:
      prefix: "--norm"
    doc: |
      Read counts normalization for the exploratory visualization part of the analysis.
      Use 'vst' for medium-to-large datasets (n > 30) and 'rlog' for small datasets
      (n < 30), when there is a wide range of sequencing depth across samples.
      Default: rlog

  remove:
    type: boolean?
    inputBinding:
      prefix: "--remove"
    doc: |
      Remove batch effect when generating normalized read counts for the exploratory
      visualization part of the analysis. Ignored if --batchby is not provided.
      Default: do not remove batch effect from normalized read counts.

  cluster_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "row"
      - "column"
      - "both"
    inputBinding:
      prefix: "--cluster"
    doc: |
      Hopach clustering method to be run on normalized read counts for the
      exploratory visualization part of the analysis. Default: do not run
      clustering

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
    inputBinding:
      prefix: "--rowdist"
    doc: |
      Distance metric for HOPACH row clustering. Ignored if --cluster is set
      to column or not provided. Default: cosangle

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
    inputBinding:
      prefix: "--columndist"
    doc: |
      Distance metric for HOPACH column clustering. Ignored if --cluster is set
      to row or not provided. Default: euclid

  center_row:
    type: boolean?
    inputBinding:
      prefix: "--center"
    doc: |
      Apply mean centering for gene expression prior to running
      clustering by row. Ignored if --cluster is set to column or
      not provided. Default: do not centered

  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: "--pdf"
    doc: |
      Export plots in PDF.
      Default: false

  color_theme:
    type:
    - "null"
    - type: enum
      symbols:
      - "gray"
      - "bw"
      - "linedraw"
      - "light"
      - "dark"
      - "minimal"
      - "classic"
      - "void"
    inputBinding:
      prefix: "--theme"
    doc: |
      Color theme for all generated plots. One of gray, bw, linedraw, light,
      dark, minimal, classic, void.
      Default: classic

  verbose:
    type: boolean?
    inputBinding:
      prefix: "--verbose"
    doc: |
      Print debug information.
      Default: false

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix.
      Default: ./sc

  parallel_memory_limit:
    type: int?
    inputBinding:
      prefix: "--memory"
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Default: 32

  vector_memory_limit:
    type: int?
    default: 128
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Default: 128

  threads:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of cores/cpus to use.
      Default: 1


outputs:

  umap_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_rd_rnaumap.png"
    doc: |
      Cells UMAP split by selected biological condition, optionally
      subsetted to the specific cluster or cell type (rnaumap dim.
      reduction).
      PNG format

  umap_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_rd_rnaumap.pdf"
    doc: |
      Cells UMAP split by selected biological condition, optionally
      subsetted to the specific cluster or cell type (rnaumap dim.
      reduction).
      PDF format

  umap_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_rd_atacumap.png"
    doc: |
      Cells UMAP split by selected biological condition, optionally
      subsetted to the specific cluster or cell type (atacumap dim.
      reduction).
      PNG format

  umap_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_rd_atacumap.pdf"
    doc: |
      Cells UMAP split by selected biological condition, optionally
      subsetted to the specific cluster or cell type (atacumap dim.
      reduction).
      PDF format

  umap_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_rd_wnnumap.png"
    doc: |
      Cells UMAP split by selected biological condition, optionally
      subsetted to the specific cluster or cell type (wnnumap dim.
      reduction).
      PNG format

  umap_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_rd_wnnumap.pdf"
    doc: |
      Cells UMAP split by selected biological condition, optionally
      subsetted to the specific cluster or cell type (wnnumap dim.
      reduction).
      PDF format

  mds_plot_html:
    type: File?
    outputBinding:
      glob: "*_mds_plot.html"
    doc: |
      MDS plot of normalized counts. Optionally batch corrected
      if --remove was set to True.
      HTML format

  pca_1_2_plot_png:
    type: File?
    outputBinding:
      glob: "*_pca_1_2.png"
    doc: |
      Normalized counts PCA (PC1 and PC2) subsetted to all DE genes regardless
      of Padj, optionally batch corrected by the selected criteria.
      PNG format

  pca_1_2_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_pca_1_2.pdf"
    doc: |
      Normalized counts PCA (PC1 and PC2) subsetted to all DE genes regardless
      of Padj, optionally batch corrected by the selected criteria.
      PDF format

  pca_2_3_plot_png:
    type: File?
    outputBinding:
      glob: "*_pca_2_3.png"
    doc: |
      Normalized counts PCA (PC2 and PC3) subsetted to all DE genes regardless
      of Padj, optionally batch corrected by the selected criteria.
      PNG format

  pca_2_3_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_pca_2_3.pdf"
    doc: |
      Normalized counts PCA (PC2 and PC3) subsetted to all DE genes regardless
      of Padj, optionally batch corrected by the selected criteria.
      PDF format

  dxpr_vlcn_plot_png:
    type: File?
    outputBinding:
      glob: "*_dxpr_vlcn.png"
    doc: |
      Volcano plot of differentially expressed genes. Highlighed genes are either
      provided by user or top 10 genes with the highest log2FC values. The direction
      of comparison is defined by --second vs --first groups of cells optionally
      subsetted to the specific cluster or cell type and coerced to the pseudobulk
      RNA-Seq samples.
      PNG format

  dxpr_vlcn_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_dxpr_vlcn.pdf"
    doc: |
      Volcano plot of differentially expressed genes. Highlighed genes are either
      provided by user or top 10 genes with the highest log2FC values. The direction
      of comparison is defined by --second vs --first groups of cells optionally
      subsetted to the specific cluster or cell type and coerced to the pseudobulk
      RNA-Seq samples.
      PDF format

  xpr_dnst_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_dnst_*.png"
    doc: |
      Log normalized gene expression density per dataset optionally subsetted to the
      specific cluster or cell type.
      PNG format

  xpr_dnst_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_dnst_*.pdf"
    doc: |
      Log normalized gene expression density per dataset optionally subsetted to the
      specific cluster or cell type.
      PDF format

  xpr_per_cell_rd_rnaumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_rd_rnaumap_*.png"
    doc: |
      Log normalized gene expression on cells UMAP per dataset optionally subsetted
      to the specific cluster or cell type (rnaumap dim. reduction).
      PNG format

  xpr_per_cell_rd_rnaumap_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_rd_rnaumap_*.pdf"
    doc: |
      Log normalized gene expression on cells UMAP per dataset optionally subsetted
      to the specific cluster or cell type (rnaumap dim. reduction).
      PDF format

  xpr_per_cell_rd_atacumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_rd_atacumap_*.png"
    doc: |
      Log normalized gene expression on cells UMAP per dataset optionally subsetted
      to the specific cluster or cell type (atacumap dim. reduction).
      PNG format

  xpr_per_cell_rd_atacumap_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_rd_atacumap_*.pdf"
    doc: |
      Log normalized gene expression on cells UMAP per dataset optionally subsetted
      to the specific cluster or cell type (atacumap dim. reduction).
      PDF format

  xpr_per_cell_rd_wnnumap_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_rd_wnnumap_*.png"
    doc: |
      Log normalized gene expression on cells UMAP per dataset optionally subsetted
      to the specific cluster or cell type (wnnumap dim. reduction).
      PNG format

  xpr_per_cell_rd_wnnumap_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_xpr_per_cell_rd_wnnumap_*.pdf"
    doc: |
      Log normalized gene expression on cells UMAP per dataset optionally subsetted
      to the specific cluster or cell type (wnnumap dim. reduction).
      PDF format

  xpr_htmp_plot_png:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.png"
    doc: |
      Normalized gene expression heatmap optionally subsetted
      to the specific cluster or cell type.
      PNG format

  xpr_htmp_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_xpr_htmp.pdf"
    doc: |
      Normalized gene expression heatmap optionally subsetted
      to the specific cluster or cell type.
      PDF format

  diff_expr_genes:
    type: File?
    outputBinding:
      glob: "*_de_genes.tsv"
    doc: |
      Differentially expressed genes.
      TSV format

  read_counts_gct:
    type: File?
    outputBinding:
      glob: "*_norm_read_counts.gct"
    doc: |
      GSEA compatible normalized counts, optionally, batch corrected.
      GCT format

  phenotypes_cls:
    type: File?
    outputBinding:
      glob: "*_phenotypes.cls"
    doc: |
      GSEA compatible phenotypes file defined based on --splitby, --first,
      and --second parameters.
      CLS format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_rna_de_pseudobulk.R"]

stdout: sc_rna_de_pseudobulk_stdout.log
stderr: sc_rna_de_pseudobulk_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell Pseudobulk Differential Expression Analysis Between Datasets"
s:name: "Single-cell Pseudobulk Differential Expression Analysis Between Datasets"
s:alternateName: "Identifies differentially expressed genes between groups of cells coerced to pseudobulk datasets"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-rna-de-pseudobulk.cwl
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
  Single-cell Pseudobulk Differential Expression Analysis Between Datasets

  Identifies differentially expressed genes between groups
  of cells coerced to pseudobulk datasets.


s:about: |
  usage: sc_rna_de_pseudobulk.R
        [-h] --query QUERY [--metadata METADATA] --splitby SPLITBY --first
        FIRST --second SECOND [--batchby BATCHBY] [--groupby GROUPBY]
        [--subset [SUBSET ...]] [--lrt] [--padj PADJ] [--genes [GENES ...]]
        [--exclude EXCLUDE] [--norm {vst,rlog}] [--remove]
        [--cluster {row,column,both}]
        [--rowdist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}]
        [--columndist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}]
        [--center] [--pdf] [--verbose] [--output OUTPUT]
        [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
        [--cpus CPUS] [--memory MEMORY]

  Single-cell Pseudobulk Differential Expression Analysis Between Datasets

  options:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include genes expression information
                          stored in the RNA assay. Additionally, 'rnaumap',
                          and/or 'atacumap', and/or 'wnnumap' dimensionality
                          reductions should be present.
    --metadata METADATA   Path to the TSV/CSV file to optionally extend Seurat
                          object metadata with categorical values using samples
                          identities. First column - 'library_id' should
                          correspond to all unique values from the 'new.ident'
                          column of the loaded Seurat object. If any of the
                          provided in this file columns are already present in
                          the Seurat object metadata, they will be overwritten.
                          Default: no extra metadata is added
    --splitby SPLITBY     Column from the Seurat object metadata to split
                          datasets into two groups to run --second vs --first
                          pseudobulk DE analysis, i.e., calculate log2FC. May be
                          one of the columns from the extra metadata added with
                          --metadata parameter. Provided value should group the
                          datasets, not cells, therefore do not use a column
                          with clustering results.
    --first FIRST         Value from the Seurat object metadata column set with
                          --splitby to define the first group of datasets for
                          pseudobulk DE analysis.
    --second SECOND       Value from the Seurat object metadata column set with
                          --splitby to define the second group of datasets for
                          pseudobulk DE analysis.
    --batchby BATCHBY     Column from the Seurat object metadata to group
                          datasets into batches. It will be used as a factor
                          variable to model batch effect when running pseudobulk
                          DE analysis (makes design formula look like
                          ~splitby+batchby). May be one of the columns from the
                          extra metadata added with --metadata parameter.
                          Provided value should batch the datasets, not cells,
                          therefore do not use a column with clustering results.
                          Default: do not model batch effect.
    --groupby GROUPBY     Column from the Seurat object metadata to group cells
                          for optional subsetting when combined with --subset
                          parameter. May be one of the columns from the extra
                          metadata added with --metadata parameter. Ignored if
                          --subset is not set. Provided value defines the groups
                          of cells, therefore any metadata column, including the
                          clustering results, may be used. Default: do not
                          subset, run pseudobulk DE analysis for all cells
                          jointly
    --subset [SUBSET ...]
                          Value(s) from the column set with --groupby parameter
                          to subset cells before running pseudobulk DE analysis.
                          If multiple values are provided run analysis jointly
                          for selected groups of cells. Ignored if --groupby is
                          not set. Default: do not subset, run pseudobulk DE
                          analysis for all cells jointly
    --lrt                 Use LRT instead of the pair-wise Wald test. If
                          --batchby is not provided use ~1 as a reduced formula,
                          otherwise ~batchby. Default: use Wald test
    --padj PADJ           In the exploratory visualization part of the analysis
                          output only features with adjusted P-value not bigger
                          than this value. Default: 0.05
    --genes [GENES ...]   Genes of interest to label on the generated plots.
                          Default: top 10 genes with the highest and the lowest
                          log2FC expression values.
    --exclude EXCLUDE     Regex pattern to identify and exclude non-coding RNA
                          genes from the pseudobulk DE analysis (not case-
                          sensitive). If any of such genes were provided in the
                          --genes parameter, they will be excluded from there as
                          well. Default: use all genes
    --norm {vst,rlog}     Read counts normalization for the exploratory
                          visualization part of the analysis. Use 'vst' for
                          medium-to-large datasets (n > 30) and 'rlog' for small
                          datasets (n < 30), when there is a wide range of
                          sequencing depth across samples. Default: rlog
    --remove              Remove batch effect when generating normalized read
                          counts for the exploratory visualization part of the
                          analysis. Ignored if --batchby is not provided.
                          Default: do not remove batch effect from normalized
                          read counts.
    --cluster {row,column,both}
                          Hopach clustering method to be run on normalized read
                          counts for the exploratory visualization part of the
                          analysis. Default: do not run clustering
    --rowdist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}
                          Distance metric for HOPACH row clustering. Ignored if
                          --cluster is set to column or not provided. Default:
                          cosangle
    --columndist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}
                          Distance metric for HOPACH column clustering. Ignored
                          if --cluster is set to row or not provided. Default:
                          euclid
    --center              Apply mean centering for gene expression prior to
                          running clustering by row. Ignored if --cluster is set
                          to column or not provided. Default: do not centered
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32