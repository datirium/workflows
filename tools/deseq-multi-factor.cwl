cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/deseq:v0.0.5


inputs:

  expression_files:
    type: File[]
    inputBinding:
      prefix: "--expression"
    doc: |
      Path to the TSV/CSV files with expression data.
      All files should have the following header:
      RefseqId GeneId Chrom TxStart TxEnd Strand TotalReads Rpkm

  expression_names:
    type: string[]
    inputBinding:
      prefix: "--aliases"
    doc: |
      Unique names for files provided in --expression,
      no special characters or spaces are allowed.
      Number and order of the names should corresponds
      to values from --expression.

  metadata_file:
    type: File
    inputBinding:
      prefix: "--metadata"
    doc: |
      Path to the TSV/CSV file to provide metadata for the
      samples from --expression. First column should have
      the name 'sample', other columns may have arbitrary names.
      The values from the 'sample' column should correspond to
      the values provided in --aliases. For a proper --contrast
      intepretation, values defined in each column should not be
      used in other columns. All metadata columns are treated as
      factors (no covariates are supported).

  design_formula:
    type: string
    inputBinding:
      prefix: "--design"
    doc: |
      Design formula. Should start with ~ and include terms from
      the --metadata table.

  reduced_formula:
    type: string?
    inputBinding:
      prefix: "--reduced"
    doc: |
      Reduced formula with the term(s) of interest removed.
      Should start with ~. If provided, force DESeq2 to run
      LRT test instead of the Wald.

  contrast:
    type: string?
    inputBinding:
      prefix: "--contrast"
    doc: |
      Contrast to be be applied for the output, formatted as
      a mathematical formula of values from the --metadata table.
      If not provided, the last term from the design formula will
      be used.

  base:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--base"
    doc: |
      Value(s) from each metadata file column(s) to be set as
      the base level(s). Number and order of provided values should
      correspond the order of columns in --metadata file. Default:
      define base levels alphabetically for each metadata column.

  feature_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "gene"
      - "transcript"
    inputBinding:
      prefix: "--type"
    doc: |
      Feature type to use for differential expression.
      If set to 'gene', use 'GeneId' column from the provided in --expression files.
      If set to 'transcript', use 'RefseqId' from the provided in --expression files.
      Default: gene

  excluded_features:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--exclude"
    doc: |
      Features to be excluded from the differential expression analysis.
      Default: include all features

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
      Read counts normalization for the exploratory visualization analysis.
      Use 'vst' for medium-to-large datasets (n > 30) and 'rlog' for
      small datasets (n < 30), when there is a wide range of sequencing
      depth across samples.
      Default: vst

  remove:
    type: string?
    inputBinding:
      prefix: "--remove"
    doc: |
      Column from the metadata file to remove batch effect
      before running differential expression analysis. If
      present, all components that include this term will be
      removed from the design and reduced formulas.
      Default: do not remove batch effect

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
      exploratory visualization analysis. Default: do not run clustering

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
      Distance metric for HOPACH row clustering. Ignored if --cluster is not
      provided. Default: cosangle

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
      Distance metric for HOPACH column clustering. Ignored if --cluster is not
      provided. Default: euclid

  center_row:
    type: boolean?
    inputBinding:
      prefix: "--center"
    doc: |
      Apply mean centering for feature expression prior to running
      clustering by row. Ignored when --cluster is not row or both.
      Default: do not centered

  selected_features:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--label"
    doc: |
      Features of interest to label on the generated volcanot plot. Default:
      top 10 features with the highest and the lowest log2 fold change
      expression values.

  maximum_padj:
    type: float?
    inputBinding:
      prefix: "--padj"
    doc: |
      In the exploratory visualization analysis output only features with
      adjusted P-value not bigger than this value. Default: 0.05

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
      Output prefix for generated files

  threads:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of cores/cpus to use. Default: 1


outputs:

  diff_expr_features:
    type: File
    outputBinding:
      glob: "*_diff_expr_features.tsv"
    doc: |
      TSV file with not filtered differentially expressed features

  read_counts_gct:
    type: File
    outputBinding:
      glob: "*_norm_read_counts.gct"
    doc: |
      GCT file with normalized, optionally batch corrected, read counts

  volcano_plot_png:
    type: File?
    outputBinding:
      glob: "*_volcano_plot.png"
    doc: |
      Volcano plot of differentially expressed features.
      PNG format

  volcano_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_volcano_plot.pdf"
    doc: |
      Volcano plot of differentially expressed features.
      PDF format

  pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_pca_plot.png"
    doc: |
      PCA plot of normalized, optionally batch corrected,
      read counts based on the top 500 features selected
      by the highest row variance
      PNG format

  pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_pca_plot.pdf"
    doc: |
      PCA plot of normalized, optionally batch corrected,
      read counts based on the top 500 features selected
      by the highest row variance
      PDF format

  mds_plot_html:
    type: File?
    outputBinding:
      glob: "*_mds_plot.html"
    doc: |
      MDS plot of normalized, optionally batch corrected,
      read counts
      HTML format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: [run_deseq_manual.R]
stdout: deseq_manual_factor_stdout.log
stderr: deseq_manual_factor_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


s:name: "DESeq2 Multi-factor Analysis"
label:  "DESeq2 Multi-factor Analysis"
s:alternateName: "Runs DeSeq2 multi-factor analysis with manual control over major parameters"


s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/deseq-multi-factor.cwl
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
  DESeq2 Multi-factor Analysis

  Runs DeSeq2 multi-factor analysis with manual control over major parameters


s:about: |
  usage: run_deseq_manual.R
        [-h] --expression EXPRESSION [EXPRESSION ...] --aliases ALIASES
        [ALIASES ...] --metadata METADATA --design DESIGN [--reduced REDUCED]
        [--contrast CONTRAST] [--base [BASE ...]] [--type {gene,transcript}]
        [--exclude [EXCLUDE ...]] [--norm {vst,rlog}] [--remove REMOVE]
        [--cluster {row,column,both}] [--center] [--label [LABEL ...]]
        [--padj PADJ] [--pdf] [--output OUTPUT] [--cpus CPUS]

  DESeq2 Multi-factor Analysis

  options:
    -h, --help            show this help message and exit
    --expression EXPRESSION [EXPRESSION ...]
                          Path to the TSV/CSV files with expression data. All
                          files should have the following header: RefseqId
                          GeneId Chrom TxStart TxEnd Strand TotalReads Rpkm
    --aliases ALIASES [ALIASES ...]
                          Unique names for files provided in --expression, no
                          special characters or spaces are allowed. Number and
                          order of the names should corresponds to values from
                          --expression.
    --metadata METADATA   Path to the TSV/CSV file to provide metadata for the
                          samples from --expression. First column should have
                          the name 'sample', other columns may have arbitrary
                          names. The values from the 'sample' column should
                          correspond to the values provided in --aliases. For a
                          proper --contrast intepretation, values defined in
                          each column should not be used in other columns. All
                          metadata columns are treated as factors (no covariates
                          are supported).
    --design DESIGN       Design formula. Should start with ~ and include terms
                          from the --metadata table.
    --reduced REDUCED     Reduced formula with the term(s) of interest removed.
                          Should start with ~. If provided, force DESeq2 to run
                          LRT test instead of the Wald.
    --contrast CONTRAST   Contrast to be be applied for the output, formatted as
                          a mathematical formula of values from the --metadata
                          table. If not provided, the last term from the design
                          formula will be used.
    --base [BASE ...]     Value(s) from each metadata file column(s) to be set
                          as the base level(s). Number and order of provided
                          values should correspond the order of columns in
                          --metadata file. Default: define base levels
                          alphabetically for each metadata column.
    --type {gene,transcript}
                          Feature type to use for differential expression. If
                          set to 'gene', use 'GeneId' column from the provided
                          in --expression files. If set to 'transcript', use
                          'RefseqId' from the provided in --expression files.
                          Default: gene
    --exclude [EXCLUDE ...]
                          Features to be excluded from the differential
                          expression analysis. Default: include all features
    --norm {vst,rlog}     Read counts normalization for the exploratory
                          visualization analysis. Use 'vst' for medium-to-large
                          datasets (n > 30) and 'rlog' for small datasets (n <
                          30), when there is a wide range of sequencing depth
                          across samples. Default: vst
    --remove REMOVE       Column from the metadata file to remove batch effect
                          before running differential expression analysis. If
                          present, all components that include this term will be
                          removed from the design and reduced formulas.
                          Default: do not remove batch effect
    --cluster {row,column,both}
                          Hopach clustering method to be run on normalized read
                          counts for the exploratory visualization analysis.
                          Default: do not run clustering
    --center              Apply mean centering for feature expression prior to
                          running clustering by row. Ignored when --cluster is
                          not row or both. Default: do not centered
    --label [LABEL ...]   Features of interest to label on the generated
                          volcanot plot. Default: top 10 features with the
                          highest and the lowest log2 fold change expression
                          values.
    --padj PADJ           In the exploratory visualization analysis output only
                          features with adjusted P-value not bigger than this
                          value. Default: 0.05
    --pdf                 Export plots in PDF. Default: false
    --output OUTPUT       Output prefix for generated files
    --cpus CPUS           Number of cores/cpus to use. Default: 1