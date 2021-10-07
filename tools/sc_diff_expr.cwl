cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/seurat:v0.0.12


inputs:

  seurat_data_rds:
    type: File
    inputBinding:
      prefix: "--rds"
    doc: |
      Path to the RDS file to load Seurat object from.
      RDS file can be produced by run_seurat.R script.

  splitby:
    type: string?
    inputBinding:
      prefix: "--splitby"
    doc: |
      Field from the Seurat object metadata to split cells into groups
      for differential expression analysis.
      Default: condition

  first_cond:
    type: string
    inputBinding:
      prefix: "--first"
    doc: |
      Value from the column set with --splitby to define a first group of cells.

  second_cond:
    type: string
    inputBinding:
      prefix: "--second"
    doc: |
      Value from the column set with --splitby to define a second group of cells.

  groupby:
    type: string?
    inputBinding:
      prefix: "--groupby"
    doc: |
      Field from the Seurat object metadata to group cells for optional subsetting
      (for example, clustering or predicted cell type field). Ignored if select is
      not set.
      Default: do not subset, use all cells.

  selected_groups:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--select"
    doc: |
      Value(s) from the column set with --groupby to optionally subset cells before
      running differential expression analysis. Ignored if groupby is not set.
      Default: do not subset, use all cells.

  selected_features:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--genes"
    doc: |
      Genes of interest to label on the generated plots.
      Default: --topn N genes with the highest and the
      lowest log2 fold change expression values.

  excluded_features:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--exgenes"
    doc: |
      Genes to be excluded from the differential expression analysis.
      Excluded genes will be still present in the dataset, but they won't
      be used in the FindMarkers function.
      Default: include all genes

  topn_genes_count:
    type: int?
    inputBinding:
      prefix: "--topn"
    doc: |
      Show N genes with the highest and N genes with the lowest log2 fold
      change expression values. Ignored when --genes is provided.
      Default: 10

  minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--minlogfc"
    doc: |
      Include only those genes that on average have the absolute value of log2
      fold change expression difference not lower than this value.
      Default: 0.25

  minimum_pct:
    type: float?
    inputBinding:
      prefix: "--minpct"
    doc: |
      Include only those genes that are detected in not lower than this
      fraction of cells in either of the two tested groups.
      Default: 0.1

  maximum_pvadj:
    type: float?
    inputBinding:
      prefix: "--maxpvadj"
    doc: |
      Include only those genes for which adjusted P-val is not bigger that this value.
      Default: 0.05

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
    inputBinding:
      prefix: "--testuse"
    doc: |
      Statistical test to use for differential gene expression analysis.
      Default: wilcox

  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: "--pdf"
    doc: |
      Export plots in PDF.
      Default: false

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix.
      Default: ./seurat

  threads:
    type: int?
    inputBinding:
      prefix: "--threads"
    doc: |
      Threads number
      Default: 1


outputs:

  avg_gene_expr_plot_png:
    type: File?
    outputBinding:
      glob: "*_avg_gene_expr.png"
    doc: |
      Log normalized average gene expression for first_cond vs second_cond cells split
      by criteria set in splitby (a.k.a condition). Cells are optionally subsetted by
      selected_groups (a.k.a clusters) from the groups defined in groupby.
      PNG format

  avg_gene_expr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_avg_gene_expr.pdf"
    doc: |
      Log normalized average gene expression for first_cond vs second_cond cells split
      by criteria set in splitby (a.k.a condition). Cells are optionally subsetted by
      selected_groups (a.k.a clusters) from the groups defined in groupby.
      PDF format

  diff_expr_genes_plot_png:
    type: File?
    outputBinding:
      glob: "*_diff_expr_genes.png"
    doc: |
      Volcano plot of differentially expressed genes for first_cond vs second_cond cells
      split by criteria set in splitby (a.k.a condition). Cells are optionally subsetted
      by selected_groups (a.k.a clusters) from the groups defined in groupby.
      PNG format

  diff_expr_genes_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_diff_expr_genes.pdf"
    doc: |
      Volcano plot of differentially expressed genes for first_cond vs second_cond cells
      split by criteria set in splitby (a.k.a condition). Cells are optionally subsetted
      by selected_groups (a.k.a clusters) from the groups defined in groupby.
      PDF format

  diff_expr_genes:
    type: File
    outputBinding:
      glob: "*_diff_expr_genes.tsv"
    doc: |
      Differentially expressed genes for first_cond vs second_cond cells split by criteria
      set in splitby (a.k.a condition). Cells are optionally subsetted by selected_groups
      (a.k.a clusters) from the groups defined in groupby.
      TSV format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_diff_expr.R"]


stdout: seurat_diff_expr_stdout.log
stderr: seurat_diff_expr_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Seurat Differential Expression"
s:name: "Seurat Differential Expression"
s:alternateName: "Runs differential expression analysis between two biological conditions for a group of cells"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc_diff_expr.cwl
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
  Seurat Differential Expression
  ==============================

  Runs differential expression analysis between two biological conditions for a group of cells


s:about: |
  usage: sc_diff_expr.R [-h] --rds RDS [--splitby SPLITBY] --first FIRST
                        --second SECOND [--groupby GROUPBY]
                        [--select [SELECT [SELECT ...]]]
                        [--genes [GENES [GENES ...]]]
                        [--exgenes [EXGENES [EXGENES ...]]] [--topn TOPN]
                        [--minlogfc MINLOGFC] [--minpct MINPCT]
                        [--maxpvadj MAXPVADJ]
                        [--testuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
                        [--pdf] [--output OUTPUT] [--threads THREADS]

  Runs differential expression analysis for a subset of cells between two
  selected groups

  optional arguments:
    -h, --help            show this help message and exit
    --rds RDS             Path to the RDS file to load Seurat object from. RDS
                          file can be produced by run_seurat.R script.
    --splitby SPLITBY     Field from the Seurat object metadata to split cells
                          into groups for differential expression analysis.
                          Default: condition
    --first FIRST         Value from the column set with --splitby to define a
                          first group of cells.
    --second SECOND       Value from the column set with --splitby to define a
                          second group of cells.
    --groupby GROUPBY     Field from the Seurat object metadata to group cells
                          for optional subsetting (for example, clustering or
                          predicted cell type field).
    --select [SELECT [SELECT ...]]
                          Value(s) from the column set with --groupby to
                          optionally subset cells before running differential
                          expression analysis. Default: do not subset, use all
                          cells.
    --genes [GENES [GENES ...]]
                          Genes of interest to label on the generated plots.
                          Default: --topn N genes with the highest and the
                          lowest log2 fold change expression values.
    --exgenes [EXGENES [EXGENES ...]]
                          Genes to be excluded from the differential expression
                          analysis. Excluded genes will be still present in the
                          dataset, but they won't be used in the FindMarkers
                          function. Default: include all genes
    --topn TOPN           Show N genes with the highest and N genes with the
                          lowest log2 fold change expression values. Ignored
                          when --genes is provided. Default: 10
    --minlogfc MINLOGFC   Include only those genes that on average have the
                          absolute value of log2 fold change expression
                          difference not lower than this value. Default: 0.25
    --minpct MINPCT       Include only those genes that are detected in not
                          lower than this fraction of cells in either of the two
                          tested groups. Default: 0.1
    --maxpvadj MAXPVADJ   Include only those genes for which adjusted P-val is
                          not bigger that this value. Default: 0.05
    --testuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}
                          Statistical test to use for differential gene
                          expression analysis. Default: wilcox
    --pdf                 Export plots in PDF. Default: false
    --output OUTPUT       Output prefix. Default: ./seurat
    --threads THREADS     Threads. Default: 1