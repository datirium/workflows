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


'sd:upstream':
  rnaseq_experiment:
    - "https://github.com/datirium/workflows/workflows/rnaseq-se.cwl"
    - "https://github.com/datirium/workflows/workflows/rnaseq-pe.cwl"
    - "https://github.com/datirium/workflows/workflows/rnaseq-se-dutp.cwl"
    - "https://github.com/datirium/workflows/workflows/rnaseq-pe-dutp.cwl"
    - "https://github.com/datirium/workflows/workflows/trim-rnaseq-pe.cwl"
    - "https://github.com/datirium/workflows/workflows/trim-rnaseq-se.cwl"
    - "https://github.com/datirium/workflows/workflows/trim-rnaseq-pe-dutp.cwl"
    - "https://github.com/datirium/workflows/workflows/trim-rnaseq-se-dutp.cwl"


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
    default: null
    label: "RNA-Seq experiments"
    doc: |
      TSV/CSV files with expression data grouped by isoforms.
      The following header is required:
      RefseqId GeneId Chrom TxStart TxEnd Strand TotalReads Rpkm
    'sd:upstreamSource': "rnaseq_experiment/rpkm_isoforms"
    'sd:localLabel': true

  genes_expression_files:
    type:
    - "null"
    - File[]
    default: null
    label: "RNA-Seq experiments"
    doc: |
      TSV/CSV files with expression data grouped by genes.
      The following header is required:
      RefseqId GeneId Chrom TxStart TxEnd Strand TotalReads Rpkm
    'sd:upstreamSource': "rnaseq_experiment/rpkm_genes"
    'sd:localLabel': true

  expression_names:
    type: string[]
    label: "RNA-Seq experiments"
    doc: |
      Unique names for files provided in 'isoforms_expression_files' or
      'genes_expression_files' inputs. No special characters or spaces are allowed.
      Number and order of the names should corresponds to order of files.
    'sd:upstreamSource': "rnaseq_experiment/alias"

  feature_type:
    type:
      - "null"
      - type: enum
        symbols: ["transcript", "gene"]
    default: "gene"
    label: "Feature type to use for differential expression"
    doc: |
      Feature type to use for differential expression.
      If set to 'gene', the 'GeneId' column will be renamed to 'feature' and used
      for differential expression. If set to 'transcript', the 'RefseqId' column
      will be used instead.

  design_formula:
    type: string
    label: "Design formula"
    doc: |
      Design formula should start with ~ and include terms from the 'metadata_file'

  reduced_formula:
    type: string?
    label: "Reduced formula"
    doc: |
      Reduced formula to compare against with the term(s) of interest removed.
      Should start with ~. Ignored when 'use_wald' is set to true.

  use_wald:
    type: boolean?
    default: false
    label: "Use pair-wise Wald test instead of LRT"
    doc: |
      Use pair-wise Wald test instead of LRT. 'reduced_formula' parameter will be ignored
      Default: use LRT test

  contrast:
    type: string
    label: "Contrast to be be applied for the output"
    doc: |
      Contrast to be be applied for the output, formatted as a mathematical formula
      of values from the 'metadata_file'

  base:
    type:
    - "null"
    - string
    - string[]
    label: "Values from each column of the metadata file to be set as base levels"
    doc: |
      Value from each column of 'metadata_file' to be set as base levels.
      Number and order of provided values should correspond the order of columns
      in the 'metadata_file'.
      Default: define base levels alphabetically for each columns of 'metadata_file'

  metadata_file:
    type: File
    label: "Metadata file to assign categories to datasets"
    doc: |
      TSV/CSV file to provide metadata for the samples from 'isoforms_expression_files'
      or 'genes_expression_files' inputs. First column should have the name 'sample',
      other columns may have arbitrary names. The values from the 'sample' column should
      correspond to the values provided in 'expression_names' input. For a proper 'contrast'
      intepretation, values defined in each column should not be used in others.

  minimum_counts:
    type: int?
    default: 0
    label: "Minimum number of counts among all samples for feature to be included in the analysis"
    doc: |
      Keep only those features where the total number of counts for all samples
      is bigger than this value.
    'sd:layout':
      advanced: true

  splitby:
    type: string?
    label: "Column from the metadata file to split samples into categories (plots only)"
    doc: |
      Used only in plots. Column from the metadata file to split samples into categories.
      Default: the first after the 'sample' column from the metadata file
    'sd:layout':
      advanced: true

  groupby:
    type: string?
    label: "Column from the metadata file to combine samples into groups (plots only)"
    doc: |
      Used only in plots. Column from the metadata file to combine samples into groups.
      Default: the last column from the metadata file
    'sd:layout':
      advanced: true

  selected_features:
    type: string?
    label: "Features of interest to label on the generated plots (plots only)"
    doc: |
      Used only in plots. Features of interest to label on the generated plots.
      Default: 'topn_features' features with the highest and the lowest log2 fold change
      expression values.
    'sd:layout':
      advanced: true

  excluded_features:
    type: string?
    label: "Features to be excluded from the differential expression analysis (plots only)"
    doc: |
      Used only in plots. Features to be excluded from the differential expression analysis.
      Default: include all features
    'sd:layout':
      advanced: true

  topn_features:
    type: int?
    default: 10
    label: "Top 2 x N features with the highest absolute log2 fold change values (plots only)"
    doc: |
      Used only in plots. Show N features with the highest and N features with the lowest log2 fold
      change expression values. Ignored with 'selected_features'.
      Default: 10
    'sd:layout':
      advanced: true

  maximum_padj:
    type: float?
    default: 0.05
    label: "Use only features with the adjusted P-value not bigger than this treshold (plots only)"
    doc: |
      Used only in plots. Output only features with adjusted P-value not bigger than this treshold.
    'sd:layout':
      advanced: true

  use_pvalue:
    type: boolean?
    default: false
    label: "Treat 'maximum_padj' as a theshold for P-value (plots only)"
    doc: |
      Used only in plots. Treat --padj as a theshold for P-value
      Default: 'maximum_padj' defines the treshold for adjusted P-value
    'sd:layout':
      advanced: true

  threads:
    type: int?
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    default: 1
    'sd:layout':
      advanced: true


outputs:

  diff_expr_features:
    type: File
    outputSource: deseq_multi_factor/diff_expr_features
    label: "TSV file with not filtered differentially expressed features"
    doc: |
      TSV file with not filtered differentially expressed features
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'DE Analysis'
        Title: 'Differentially expressed features'

  volcano_plot_png:
    type: File?
    outputSource: deseq_multi_factor/volcano_plot_png
    label: "Volcano plot of differentially expressed features"
    doc: |
      Volcano plot of differentially expressed features.
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'Volcano plot of differentially expressed features'

  volcano_plot_pdf:
    type: File?
    outputSource: deseq_multi_factor/volcano_plot_pdf
    label: "Volcano plot of differentially expressed features"
    doc: |
      Volcano plot of differentially expressed features.
      PDF format

  pca_plot_png:
    type: File?
    outputSource: deseq_multi_factor/pca_plot_png
    label: "PCA plot of rlog-normalized counts based on the top 500 features with the highest row variance"
    doc: |
      PCA plot of rlog-normalized counts based on the top 500
      features selected by the highest row variance
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'PCA plot of rlog-normalized counts based on the top 500 features with the highest row variance'

  pca_plot_pdf:
    type: File?
    outputSource: deseq_multi_factor/pca_plot_pdf
    label: "PCA plot of rlog-normalized counts based on the top 500 features with the highest row variance"
    doc: |
      PCA plot of rlog-normalized counts based on the top 500
      features selected by the highest row variance
      PDF format

  counts_plot_png:
    type: File?
    outputSource: deseq_multi_factor/counts_plot_png
    label: "rlog-normalized counts plots"
    doc: |
      rlog-normalized counts plots
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'rlog-normalized counts plots'

  counts_plot_pdf:
    type: File?
    outputSource: deseq_multi_factor/counts_plot_pdf
    label: "rlog-normalized counts plots"
    doc: |
      rlog-normalized counts plots
      PDF format

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
      reduced_formula: reduced_formula
      contrast: contrast
      base: base
      feature_type: feature_type
      minimum_counts: minimum_counts
      use_wald: use_wald
      splitby: splitby
      groupby: groupby
      selected_features:
        source: selected_features
        valueFrom: $(split_features(self))
      excluded_features:
        source: excluded_features
        valueFrom: $(split_features(self))
      topn_features: topn_features
      maximum_padj: maximum_padj
      use_pvalue: use_pvalue
      export_pdf_plots:
        default:  true
      threads: threads
    out:
    - diff_expr_features
    - volcano_plot_png
    - volcano_plot_pdf
    - pca_plot_png
    - pca_plot_pdf
    - counts_plot_png
    - counts_plot_pdf
    - stdout_log
    - stderr_log

