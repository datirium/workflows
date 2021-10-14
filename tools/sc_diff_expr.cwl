cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/seurat:v0.0.13


inputs:

  seurat_data_rds:
    type: File
    inputBinding:
      prefix: "--rds"
    doc: |
      Path to the RDS file to load Seurat object from.
      RDS file can be produced by run_seurat.R script.

  conditions_data:
    type: File?
    inputBinding:
      prefix: "--condition"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat object metadata. First
      column 'library_id' should include all unique values from the 'new.ident'
      column of the loaded from --rds Seurat object metadata. All other columns will
      be added to the Seurat object metadata. If any of the provided in this file
      columns were already present in the Seurat object metadata, they will be
      overwritten.
      Default: no metadata columns will be added or overwritten

  splitby:
    type: string
    inputBinding:
      prefix: "--splitby"
    doc: |
      Column from the Seurat object metadata to split cells into two groups
      to run --second vs --first differential expression analysis. May include
      columns from the metadata fields added with --condition.

  first_cond:
    type: string
    inputBinding:
      prefix: "--first"
    doc: |
      Value from the Seurat object metadata column set with --splitby to define the
      first group of cells or pseudobulk RNA-Seq samples (when using --pseudo).

  second_cond:
    type: string
    inputBinding:
      prefix: "--second"
    doc: |
      Value from the Seurat object metadata column set with --splitby to define the
      the second group of cells or pseudobulk RNA-Seq samples (when using --pseudo).

  batchby:
    type: string?
    inputBinding:
      prefix: "--batchby"
    doc: |
      Column from the Seurat object metadata to define the variable that should
      be modelled as a batch effect when running differential expression analysis.
      Applied only when --testuse is one of 'LR', 'negbinom', 'poisson', or 'MAST',
      or when using --pseudo. May include columns from the metadata fields added
      with --condition. Values selected from the column set with --batchby should
      establish 1:1 relation with the 'new.ident' column of the Seurat object loaded
      from --rds.
      Default: do not model batch effect.

  groupby:
    type: string?
    inputBinding:
      prefix: "--groupby"
    doc: |
      Column from the Seurat object metadata to group cells for optional
      subsetting (for example, subset to the specific cluster or predicted
      cell type). May include columns from the metadata fields added with
      --condition.

  selected_groups:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--select"
    doc: |
      Value(s) from the column set with --groupby to optionally subset cells
      before running differential expression analysis.
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
      Default: include all genes

  topn_genes_count:
    type: int?
    inputBinding:
      prefix: "--topn"
    doc: |
      Show N genes with the highest and N genes with the lowest log2 fold
      change expression values. Ignored with --genes.
      Default: 10

  minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--minlogfc"
    doc: |
      Include only those genes that on average have the absolute value of log2
      fold change expression difference not lower than this value. Increasing
      --minlogfc speeds up calculations, but can cause missing weaker signals.
      Ignored with --pseudo.
      Default: 0.25

  minimum_pct:
    type: float?
    inputBinding:
      prefix: "--minpct"
    doc: |
      Include only those genes that are detected in not lower than this fraction of cells
      in either of the two tested groups. Increasing --minpct speeds up calculations by not
      testing genes that are very infrequently expressed. Ignored with --pseudo.
      Default: 0.1

  maximum_pvadj:
    type: float?
    inputBinding:
      prefix: "--maxpvadj"
    doc: |
      Include only those genes for which adjusted P-val is not bigger that this value.
      Default: 0.1

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
      Ignored with --pseudo.
      Default: wilcox

  pseudo:
    type: boolean?
    inputBinding:
      prefix: "--pseudo"
    doc: |
      Aggregate gene expression of the cells from the same dataset into a pseudobulk
      RNA-Seq sample before running differential expression analysis with DESeq2.
      The following parameters will be ignored: --testuse, --minpct, --minlogfc.
      Default: false

  lrt:
    type: boolean?
    inputBinding:
      prefix: "--lrt"
    doc: |
      Use LRT instead of the pair-wise Wald test. Shows any differences across the variable
      set with --batchby whith the log2 fold changes calculated as the average expression
      changes due to criteria set with --splitby. Ignored when --pseudo or --batchby
      parameters are not provided.
      Default: use Wald test

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

  cell_abundance_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap.png"
    doc: |
      Cell abundance plot split by criteria set in splitby (a.k.a condition) and optionally
      subsetted by selected_groups (a.k.a clusters) from the groups defined in groupby.
      PNG format

  cell_abundance_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap.pdf"
    doc: |
      Cell abundance plot split by criteria set in splitby (a.k.a condition) and optionally
      subsetted by selected_groups (a.k.a clusters) from the groups defined in groupby.
      PDF format

  aggr_gene_expr_plot_png:
    type: File?
    outputBinding:
      glob: "*_counts.png"
    doc: |
      Log normalized aggregated gene expression split by criteria set in splitby
      (a.k.a condition).
      PNG format

  aggr_gene_expr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_counts.pdf"
    doc: |
      Log normalized aggregated gene expression split by criteria set in splitby
      (a.k.a condition).
      PDF format

  diff_expr_genes_plot_png:
    type: File?
    outputBinding:
      glob: "*_diff_expr_genes.png"
    doc: |
      Volcano plot of differentially expressed genes for second_cond vs first_cond cells
      or pseudobulk RNA-Seq samples split by criteria set in splitby (a.k.a condition)
      and optionally subsetted by selected_groups (a.k.a clusters) from the groups defined
      in groupby.
      PNG format

  diff_expr_genes_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_diff_expr_genes.pdf"
    doc: |
      Volcano plot of differentially expressed genes for second_cond vs first_cond cells
      or pseudobulk RNA-Seq samples split by criteria set in splitby (a.k.a condition)
      and optionally subsetted by selected_groups (a.k.a clusters) from the groups defined
      in groupby.
      PDF format

  diff_expr_genes:
    type: File
    outputBinding:
      glob: "*_diff_expr_genes.tsv"
    doc: |
      Differentially expressed genes for second_cond vs first_cond cells or pseudobulk
      RNA-Seq samples split by criteria set in splitby (a.k.a condition) and optionally
      subsetted by selected_groups (a.k.a clusters) from the groups defined in groupby.
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


label: "Single-cell Differential Expression Analysis"
s:name: "Single-cell Differential Expression Analysis"
s:alternateName: "Runs differential expression analysis for a subset of cells between two selected conditions"

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
  Single-cell Differential Expression Analysis
  =============================================

  Runs differential expression analysis for a subset of cells between two selected conditions


s:about: |
  usage: /Users/kot4or/workspaces/cwl_ws/workflows/tools/dockerfiles/scripts/sc_diff_expr.R
        [-h] --rds RDS [--condition CONDITION] --splitby SPLITBY --first FIRST
        --second SECOND [--batchby BATCHBY] [--groupby GROUPBY]
        [--select [SELECT ...]] [--genes [GENES ...]] [--exgenes [EXGENES ...]]
        [--topn TOPN] [--minlogfc MINLOGFC] [--minpct MINPCT]
        [--maxpvadj MAXPVADJ]
        [--testuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}]
        [--pseudo] [--lrt] [--pdf] [--output OUTPUT] [--threads THREADS]

  Differential expression analysis for a subset of cells between two selected
  conditions

  optional arguments:
    -h, --help            show this help message and exit
    --rds RDS             Path to the RDS file to load Seurat object from. RDS
                          file can be produced by run_seurat.R script.
    --condition CONDITION
                          Path to the TSV/CSV file to optionally extend Seurat
                          object metadata. First column 'library_id' should
                          include all unique values from the 'new.ident' column
                          of the loaded from --rds Seurat object metadata. All
                          other columns will be added to the Seurat object
                          metadata. If any of the provided in this file columns
                          were already present in the Seurat object metadata,
                          they will be overwritten. Default: no metadata columns
                          will be added or overwritten
    --splitby SPLITBY     Column from the Seurat object metadata to split cells
                          into two groups to run --second vs --first
                          differential expression analysis. May include columns
                          from the metadata fields added with --condition.
    --first FIRST         Value from the Seurat object metadata column set with
                          --splitby to define the first group of cells or
                          pseudobulk RNA-Seq samples (when using --pseudo).
    --second SECOND       Value from the Seurat object metadata column set with
                          --splitby to define the the second group of cells or
                          pseudobulk RNA-Seq samples (when using --pseudo)
    --batchby BATCHBY     Column from the Seurat object metadata to define the
                          variable that should be modelled as a batch effect
                          when running differential expression analysis. Applied
                          only when --testuse is one of 'LR', 'negbinom',
                          'poisson', or 'MAST', or when using --pseudo. May
                          include columns from the metadata fields added with
                          --condition. Values selected from the column set with
                          --batchby should establish 1:1 relation with the
                          'new.ident' column of the Seurat object loaded from
                          --rds. Default: do not model batch effect.
    --groupby GROUPBY     Column from the Seurat object metadata to group cells
                          for optional subsetting (for example, subset to the
                          specific cluster or predicted cell type). May include
                          columns from the metadata fields added with
                          --condition.
    --select [SELECT ...]
                          Value(s) from the column set with --groupby to
                          optionally subset cells before running differential
                          expression analysis. Default: do not subset, use all
                          cells.
    --genes [GENES ...]   Genes of interest to label on the generated plots.
                          Default: --topn N genes with the highest and the
                          lowest log2 fold change expression values.
    --exgenes [EXGENES ...]
                          Genes to be excluded from the differential expression
                          analysis. Default: include all genes
    --topn TOPN           Show N genes with the highest and N genes with the
                          lowest log2 fold change expression values. Ignored
                          with --genes. Default: 10
    --minlogfc MINLOGFC   Include only those genes that on average have the
                          absolute value of log2 fold change expression
                          difference not lower than this value. Increasing
                          --minlogfc speeds up calculations, but can cause
                          missing weaker signals. Ignored with --pseudo.
                          Default: 0.25
    --minpct MINPCT       Include only those genes that are detected in not
                          lower than this fraction of cells in either of the two
                          tested groups. Increasing --minpct speeds up
                          calculations by not testing genes that are very
                          infrequently expressed. Ignored with --pseudo.
                          Default: 0.1
    --maxpvadj MAXPVADJ   Include only those genes for which adjusted P-val is
                          not bigger that this value. Default: 0.1
    --testuse {wilcox,bimod,roc,t,negbinom,poisson,LR,MAST,DESeq2}
                          Statistical test to use for differential gene
                          expression analysis. Ignored with --pseudo. Default:
                          wilcox
    --pseudo              Aggregate gene expression of the cells from the same
                          dataset into a pseudobulk RNA-Seq sample before
                          running differential expression analysis with DESeq2.
                          The following parameters will be ignored: --testuse,
                          --minpct, --minlogfc. Default: false
    --lrt                 Use LRT instead of the pair-wise Wald test. Shows any
                          differences across the variable set with --batchby
                          whith the log2 fold changes calculated as the average
                          expression changes due to criteria set with --splitby.
                          Ignored when --pseudo or --batchby parameters are not
                          provided. Default: use Wald test
    --pdf                 Export plots in PDF. Default: false
    --output OUTPUT       Output prefix. Default: ./seurat
    --threads THREADS     Threads. Default: 1