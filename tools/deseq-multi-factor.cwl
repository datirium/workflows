cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap-deseq:v0.0.25


inputs:

  expression_files:
    type: File[]
    inputBinding:
      prefix: "--expression"
    doc: |
      Path to the TSV/CSV files with expression data. All files should have the following header
      RefseqId GeneId Chrom TxStart TxEnd Strand TotalReads Rpkm

  expression_names:
    type: string[]
    inputBinding:
      prefix: "--aliases"
    doc: |
      Unique names for files provided in --expression, no special characters or spaces are allowed.
      Number and order of the names should corresponds to values from --expression

  metadata_file:
    type: File
    inputBinding:
      prefix: "--metadata"
    doc: |
      Path to the TSV/CSV file to provide metadata for the samples from --expression.
      First column should have the name 'sample', other columns may have arbitrary names.
      The values from the 'sample' column should correspond to the values provided in --aliases.
      For a proper --contrast intepretation, values defined in each column should not be used in others.

  design_formula:
    type: string
    inputBinding:
      prefix: "--design"
    doc: |
      Design formula. Should start with ~ and include terms from the --metadata table

  reduced_formula:
    type: string?
    inputBinding:
      prefix: "--reduced"
    doc: |
      Reduced formula to compare against with the term(s) of interest removed.
      Should start with ~. Ignored when run with --wald

  contrast:
    type: string
    inputBinding:
      prefix: "--contrast"
    doc: |
      Contrast to be be applied for the output, formatted as a mathematical formula
      of values from the --metadata table

  base:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--base"
    doc: |
      Value from each column of metadata file to be set as base levels.
      Number and order of provided values should correspond the order of columns
      in --metadata file.
      Default: define base levels alphabetically for each columns of metadata table

  feature_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "gene"
      - "transcript"
    inputBinding:
      prefix: "--ftype"
    doc: |
      Feature type to use for differential expression.
      If set to 'gene', use 'GeneId' column from the provided in --expression files.
      If set to 'transcript', use 'RefseqId' from the provided in --expression files.
      Default: gene

  minimum_counts:
    type: int?
    inputBinding:
      prefix: "--mincounts"
    doc: |
      Keep only those features where the total number of counts for all samples
      is bigger than this value.
      Default: 0

  use_wald:
    type: boolean?
    inputBinding:
      prefix: "--wald"
    doc: |
      Use pair-wise Wald test instead of LRT. --reduced parameter will be ignored
      Default: use LRT test

  splitby:
    type: string?
    inputBinding:
      prefix: "--splitby"
    doc: |
      Used only in plots. Column from the metadata file to split samples into categories.
      Default: the first after the 'sample' column from the metadata file

  groupby:
    type: string?
    inputBinding:
      prefix: "--groupby"
    doc: |
      Used only in plots. Column from the metadata file to combine samples into groups.
      Default: the last column from the metadata file

  selected_features:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--features"
    doc: |
      Used only in plots. Features of interest to label on the generated plots.
      Default: --topn N features with the highest and the lowest log2 fold change
      expression values.

  excluded_features:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--exfeatures"
    doc: |
      Used only in plots. Features to be excluded from the differential expression analysis.
      Default: include all features

  topn_features:
    type: int?
    inputBinding:
      prefix: "--topn"
    doc: |
      Used only in plots. Show N features with the highest and N features with the lowest log2 fold
      change expression values. Ignored with --features.
      Default: 10

  maximum_padj:
    type: float?
    inputBinding:
      prefix: "--padj"
    doc: |
      Used only in plots. Output only features with adjusted P-value not bigger than this treshold.
      Default: 0.05

  use_pvalue:
    type: boolean?
    inputBinding:
      prefix: "--usepvalue"
    doc: |
      Used only in plots. Treat --padj as a theshold for P-value
      Default: --padj defines the treshold for adjusted P-value

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
      position: 9
      prefix: "--output"
    doc: |
      Output prefix for generated files

  threads:
    type: int?
    inputBinding:
      position: 10
      prefix: "--threads"
    doc: |
      Threads number


outputs:

  diff_expr_features:
    type: File
    outputBinding:
      glob: "*_diff_expr_features.tsv"
    doc: |
      TSV file with not filtered differentially expressed features

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
      PCA plot of rlog-normalized counts based on the top 500
      features selected by the highest row variance
      PNG format

  pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_pca_plot.pdf"
    doc: |
      PCA plot of rlog-normalized counts based on the top 500
      features selected by the highest row variance
      PDF format

  counts_plot_png:
    type: File?
    outputBinding:
      glob: "*_counts_plot.png"
    doc: |
      rlog-normalized counts plots
      PNG format

  counts_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_counts_plot.pdf"
    doc: |
      rlog-normalized counts plots
      PDF format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: [run_deseq_manual.R]
stdout: deseq_multi_factor_stdout.log
stderr: deseq_multi_factor_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


s:name: "DESeq2 Multi-factor Analysis"
label:  "DESeq2 Multi-factor Analysis"
s:alternateName: "Runs DeSeq2 multi-factor analysis with manual control over major parameters"


s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/deseq-lrt.cwl
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
  ============================
  Runs DeSeq2 multi-factor analysis with manual control over major parameters


s:about: |
  usage: run_deseq_manual.R
        [-h] --expression EXPRESSION [EXPRESSION ...] --aliases ALIASES
        [ALIASES ...] --metadata METADATA --design DESIGN [--reduced REDUCED]
        --contrast CONTRAST [--base [BASE ...]] [--ftype {gene,transcript}]
        [--mincounts MINCOUNTS] [--wald] [--splitby SPLITBY]
        [--groupby GROUPBY] [--features [FEATURES ...]]
        [--exfeatures [EXFEATURES ...]] [--topn TOPN] [--padj PADJ]
        [--usepvalue] [--pdf] [--output OUTPUT] [--threads THREADS]

  Run DeSeq2 with manual control over major parameters

  optional arguments:
    -h, --help            show this help message and exit
    --expression EXPRESSION [EXPRESSION ...]
                          Path to the TSV/CSV files with expression data. All
                          files should have the following header: RefseqId
                          GeneId Chrom TxStart TxEnd Strand TotalReads Rpkm
    --aliases ALIASES [ALIASES ...]
                          Unique names for files provided in --expression, no
                          special characters or spaces are allowed. Number and
                          order of the names should corresponds to values from
                          --expression
    --metadata METADATA   Path to the TSV/CSV file to provide metadata for the
                          samples from --expression. First column should have
                          the name 'sample', other columns may have arbitrary
                          names. The values from the 'sample' column should
                          correspond to the values provided in --aliases. For a
                          proper --contrast intepretation, values defined in
                          each column should not be used in others.
    --design DESIGN       Design formula. Should start with ~ and include terms
                          from the --metadata table
    --reduced REDUCED     Reduced formula to compare against with the term(s) of
                          interest removed. Should start with ~. Ignored when
                          run with --wald
    --contrast CONTRAST   Contrast to be be applied for the output, formatted as
                          a mathematical formula of values from the --metadata
                          table
    --base [BASE ...]     Value from each column of metadata file to be set as
                          base levels. Number and order of provided values
                          should correspond the order of columns in --metadata
                          file. Default: define base levels alphabetically for
                          each columns of metadata table
    --ftype {gene,transcript}
                          Feature type to use for differential expression. If
                          set to 'gene', use 'GeneId' column from the provided
                          in --expression files. If set to 'transcript', use
                          'RefseqId' from the provided in --expression files.
                          Default: gene
    --mincounts MINCOUNTS
                          Keep only those features where the total number of
                          counts for all samples is bigger than this value.
                          Default: 0
    --wald                Use pair-wise Wald test instead of LRT. --reduced
                          parameter will be ignored Default: use LRT test
    --splitby SPLITBY     Used only in plots. Column from the metadata file to
                          split samples into categories. Default: the first
                          after the 'sample' column from the metadata file
    --groupby GROUPBY     Used only in plots. Column from the metadata file to
                          combine samples into groups. Default: the last column
                          from the metadata file
    --features [FEATURES ...]
                          Used only in plots. Features of interest to label on
                          the generated plots. Default: --topn N features with
                          the highest and the lowest log2 fold change expression
                          values.
    --exfeatures [EXFEATURES ...]
                          Used only in plots. Features to be excluded from the
                          differential expression analysis. Default: include all
                          features
    --topn TOPN           Used only in plots. Show N features with the highest
                          and N features with the lowest log2 fold change
                          expression values. Ignored with --features. Default:
                          10
    --padj PADJ           Used only in plots. Output only features with adjusted
                          P-value not bigger than this treshold. Default: 0.05
    --usepvalue           Used only in plots. Treat --padj as a theshold for
                          P-value Default: --padj defines the treshold for
                          adjusted P-value
    --pdf                 Export plots in PDF. Default: false
    --output OUTPUT       Output prefix for generated files
    --threads THREADS     Threads number