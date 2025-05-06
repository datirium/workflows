cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement

'sd:upstream':
  rnaseq_experiment:
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-pe-smarter-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"
    - "https://github.com/datirium/workflows/workflows/trim-rnaseq-pe.cwl"
    - "https://github.com/datirium/workflows/workflows/trim-rnaseq-se.cwl"
    - "https://github.com/datirium/workflows/workflows/trim-rnaseq-pe-dutp.cwl"
    - "https://github.com/datirium/workflows/workflows/trim-rnaseq-pe-smarter-dutp.cwl"
    - "https://github.com/datirium/workflows/workflows/trim-rnaseq-se-dutp.cwl"

inputs:

  alias_trigger:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1
  
  expression_files:
    type: File[]
    format: "http://edamontology.org/format_3752"
    label: "RNA-Seq experiments"
    doc: "CSV/TSV input files grouped by isoforms"
    'sd:upstreamSource': "rnaseq_experiment/rpkm_isoforms"
    'sd:localLabel': true

  expression_file_names:
    type: string[]
    label: "RNA-Seq experiments"
    doc: "Aliases for RNA-Seq experiments. The same aliases should be used in metadata file"
    'sd:upstreamSource': "rnaseq_experiment/alias"

  group_by:
    type:
      - "null"
      - type: enum
        symbols: [ "isoforms", "genes", "common tss" ]
    default: "genes"
    label: "Group by"
    doc: "Grouping method for features: isoforms, genes or common tss"

  metadata_file:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Metadata file to describe categories. See workflow description for details"
    doc: "Metadata file to describe relation between samples, formatted as CSV/TSV"

  design_formula:
    type: string
    label: "Design formula. See workflow description for details"
    doc: "Design formula. Should start with ~. See DESeq2 manual for details"

  reduced_formula:
    type: string
    label: "Reduced formula to compare against with the term(s) of interest removed. See workflow description for details"
    doc: "Reduced formula to compare against with the term(s) of interest removed. Should start with ~. See DESeq2 manual for details"

  batchcorrection:
    type:
      - "null"
      - type: enum
        symbols:
          - "none"
          - "combatseq"
          - "model"
    default: "none"
    label: "Batch Correction Method"
    doc: |
      Specifies the batch correction method to be applied.
      'combatseq' applies ComBat_seq at the beginning of the analysis.
      'limmaremovebatcheffect' notes the batch correction to be applied in step 2.
      Default: none
    'sd:layout':
      advanced: true

  scaling_type:
    type:
      - "null"
      - type: enum
        symbols:
          - "minmax"
          - "zscore"
    default: "zscore"
    label: "Expression Data Scaling Method"
    doc: |
      Specifies the type of scaling to be applied to the expression data.
      'minmax' applies Min-Max scaling, normalizing values to a range of [-2, 2].
      'zscore' applies Z-score standardization, centering data to mean = 0 and standard deviation = 1.
      Default: zscore.
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
          - "cor"
          - "abscor"
    default: "euclid"
    label: "Distance metric for HOPACH column clustering"
    doc: |
      Distance metric for HOPACH column clustering. Ignored if --cluster is not
      provided. Default: euclid
    'sd:layout':
      advanced: true

  fdr:
    type: float?
    default: 0.1
    label: "Maximum P-adjusted to show features in the exploratory visualization analysis"
    doc: |
      In the exploratory visualization part of the analysis, output only features
      with adjusted p-value (FDR) not bigger than this value. Also, the significance
      cutoff used for optimizing the independent filtering. Default: 0.1.
    'sd:layout':
      advanced: true

  lfcthreshold:
    type: float?
    default: 0.59
    label: "Log2 Fold Change Threshold"
    doc: |
      Log2 fold change threshold for determining significant differential expression.
      Genes with absolute log2 fold change greater than this threshold will be considered.
      Default: 0.59 (about 1.5 fold change)
    'sd:layout':
      advanced: true

  use_lfc_thresh:
    type: boolean
    default: false
    label: "Use lfcthreshold as the null hypothesis value in the results function call"
    doc: "Use lfcthreshold as the null hypothesis value in the results function call. Default: TRUE"
    'sd:layout':
      advanced: true

  rpkm_cutoff:
    type: int?
    default: null
    label: "RPKM cutoff for filtering expression data"
    doc: |
      Integer cutoff for filtering rows in the expression data.
      Rows will be retained if any column with "Rpkm" in its name exceeds this cutoff.
      If not provided (i.e. remains null), no filtering is applied.
      Recommended values are: 3, 5.
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

  k_hopach:
    type: int?
    default: 3
    label: "Number of levels for HOPACH clustering"
    doc: "Number of levels (depth) for Hopach clustering: min - 1, max - 15. Default: 3."
    'sd:layout':
      advanced: true

  kmax_hopach:
    type: int?
    default: 5
    label: "Maximum number of clusters at each level for HOPACH clustering"
    doc: "Maximum number of clusters at each level for Hopach clustering: min - 2, max - 9. Default: 5."
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 6
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true

  lrt_only_mode:
    type: boolean
    default: false
    label: "Run the LRT only"
    doc: "LRT only mode. Contrasts are not run in this mode"
    'sd:layout':
      advanced: true

  test_mode:
    type: boolean
    default: false
    label: "Run only 100 genes for testing purposes to speed up DESeq2"
    'sd:layout':
      advanced: true

outputs:

  lrt_diff_expr:
    type: File
    label: "Differentially expressed features grouped by isoforms, genes or common TSS"
    format: "http://edamontology.org/format_3475"
    doc: "DESeq2 generated file of differentially expressed features grouped by isoforms, genes or common TSS in TSV format"
    outputSource: deseq/lrt_diff_expr
    'sd:visualPlugins':
      - syncfusiongrid:
          tab: 'Differential Expression Analysis'
          Title: 'Combined DESeq2 results'

  contrasts_table:
    type: File?
    label: "Comprehensive List of DESeq2 Contrasts"
    format: "http://edamontology.org/format_3475"
    doc: "This file contains all possible contrasts (main and interaction effects) generated by DESeq2 Wald test for the complex interaction design formula. The contrasts are grouped by isoforms, genes, or common TSS in TSV format."
    outputSource: deseq/contrasts_table
    'sd:visualPlugins':
      - syncfusiongrid:
          tab: 'Complex Interaction Analysis'
          Title: 'DESeq2 Wald Test Contrasts'

  dsq_obj_data:
    type: File?
    label: "Comprehensive List of DESeq2 Contrasts in RDS format"
    doc: "This RDS file contains all possible contrasts (main and interaction effects) generated by DESeq2 Wald test for the complex interaction design formula."
    outputSource: deseq/dsq_obj_data

  lrt_summary_md:
    type: File
    label: "DESeq2 Results Summary"
    doc: |
      Markdown file that includes a warning message if batch_file is not provided
      but batchcorrection is set to "combatseq" or "limmaremovebatcheffect". Additionally,
      it contains a detailed summary of the DESeq2 analysis results, including total genes
      with non-zero read count, log fold changes (LFC), outliers, and low count genes.
    outputSource: deseq/lrt_summary_md
    "sd:visualPlugins":
      - markdownView:
          tab: "Overview"

  read_counts_file_all:
    type: File
    label: "Normalized read counts in GCT format without padj filtering"
    format: "http://edamontology.org/format_3709"
    doc: "DESeq generated files of all normalized read counts in GCT format. Compatible with GSEA"
    outputSource: deseq/counts_all_gct

  read_counts_file_filtered:
    type: File?
    label: "Normalized read counts in GCT format filtered by padj"
    format: "http://edamontology.org/format_3709"
    doc: "DESeq generated files of padj-filtered normalized read counts in GCT format. Compatible with Morpheus heatmap"
    outputSource: deseq/counts_filtered_gct

  mds_plots_html:
    type: File
    outputSource: deseq/mds_plots_html
    label: "MDS plots of normalized counts"
    doc: |
      MDS plots of normalized counts for each contrast
      HTML format
    'sd:visualPlugins':
      - linkList:
          tab: 'Overview'
          target: "_blank"

  mds_plots_corrected_html:
    type: File?
    outputSource: deseq/mds_plots_corrected_html
    label: "MDS plots of normalized counts corrected after batch effect removal"
    doc: |
      MDS plots of normalized counts for each contrast
      HTML format
    'sd:visualPlugins':
      - linkList:
          tab: 'Overview'
          target: "_blank"

  heatmap_html:
    type: File?
    outputSource: morpheus_heatmap/heatmap_html
    label: "Combined Heatmap of normalized counts"
    doc: |
      Morpheus heatmap in HTML format combining all contrasts
    'sd:visualPlugins':
      - linkList:
          tab: 'Overview'
          target: "_blank"

  alignment_stats_barchart:
    type: File
    label: "Alignment statistics bar chart"
    doc: "Alignment statistics bar chart"
    outputSource: deseq/alignment_stats_barchart

  deseq_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "DESeq2 LRT Step 1 stdout log"
    doc: "DESeq2 LRT Step 1 stdout log"
    outputSource: deseq/stdout_log

  deseq_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "DESeq2 LRT Step 1 stderr log"
    doc: "DESeq2 LRT Step 1 stderr log"
    outputSource: deseq/stderr_log

steps:

  group_isoforms:
    run: ../tools/group-isoforms-batch.cwl
    in:
      isoforms_file: expression_files
    out:
      - genes_file
      - common_tss_file

  deseq:
    run: ../tools/deseq-lrt-step-1.cwl
    in:
      expression_files:
        source: [ group_by, expression_files, group_isoforms/genes_file, group_isoforms/common_tss_file ]
        valueFrom: |
          ${
            if (self[0] == "isoforms") {
              return self[1];
            } else if (self[0] == "genes") {
              return self[2];
            } else {
              return self[3];
            }
          }
      expression_file_names: expression_file_names
      metadata_file: metadata_file
      design_formula: design_formula
      reduced_formula: reduced_formula
      batchcorrection: batchcorrection
      fdr: fdr
      lfcthreshold: lfcthreshold
      use_lfc_thresh: use_lfc_thresh
      row_distance: row_distance
      column_distance: column_distance
      rpkm_cutoff: rpkm_cutoff
      cluster_method: cluster_method
      scaling_type: scaling_type
      k_hopach: k_hopach
      kmax_hopach: kmax_hopach
      output_prefix: alias_trigger
      threads: threads
      lrt_only_mode: lrt_only_mode
      test_mode: test_mode
    out:
      - lrt_diff_expr
      - contrasts_table
      - dsq_obj_data
      - mds_plots_html
      - mds_plots_corrected_html
      - counts_all_gct
      - counts_filtered_gct
      - lrt_summary_md
      - alignment_stats_barchart
      - stdout_log
      - stderr_log

  morpheus_heatmap:
    run: ../tools/morpheus-heatmap.cwl
    in:
      read_counts_gct: deseq/counts_filtered_gct
    out:
      - heatmap_html
      - stdout_log
      - stderr_log

$namespaces:
  s: http://schema.org/

$schemas:
  - https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "DESeq2 (LRT, step 1) - differential gene expression analysis using likelihood ratio test"
label: "DESeq2 (LRT, step 1) - differential gene expression analysis using likelihood ratio test"
s:alternateName: "Differential gene expression analysis based on the LRT (likelihood ratio test)"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/deseq-lrt-step-1.cwl
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
  Runs DESeq2 using LRT (Likelihood Ratio Test)
  =============================================

  The LRT examines two models for the counts: a full model with a certain number of terms and a reduced model, in which some of the terms of the full model are removed. The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero.

  The LRT is useful for testing multiple terms at once, for example, testing 3 or more levels of a factor at once or all interactions between two variables. The LRT for count data is conceptually similar to an analysis of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a reduced model.

  When performing a likelihood ratio test, the p-values and the test statistic (the 'stat' column) are values for the test that removes all of the variables which are present in the full design and not in the reduced design. This tests the null hypothesis that all the coefficients from these variables and levels of these factors are equal to zero.

  The likelihood ratio test p-values therefore represent a test of all the variables and all the levels of factors which are among these variables. However, the results table only has space for one column of log fold change, so a single variable and a single comparison is shown (among the potentially multiple log fold changes which were tested in the likelihood ratio test). This indicates that the p-value is for the likelihood ratio test of all the variables and all the levels, while the log fold change is a single comparison from among those variables and levels.

  **Technical notes**

  1. **Biological Replicates:** At least two biological replicates are required for every compared category.
  2. **Metadata File:** The metadata file describes relations between compared experiments. For example:

     ```csv
     ,time,condition
     DH1,day5,WT
     DH2,day5,KO
     DH3,day7,WT
     DH4,day7,KO
     DH5,day7,KO
     ```
     where `time`, `condition`, `day5`, `day7`, `WT`, `KO` should be single words (without spaces), and `DH1`, `DH2`, `DH3`, `DH4`, `DH5` correspond to the experiment aliases set in **RNA-Seq experiments** input.
  3. **Design and Reduced Formulas:** Design and reduced formulas should start with `~` and include categories or, optionally, their interactions from the metadata file header. See details in the DESeq2 manual [here](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions) and [here](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#likelihood-ratio-test).
  4. **Batch Correction:** If batch correction is required, provide the `batch_file` input. This file should be a headerless TSV/CSV file where the first column contains sample names matching `expression_file_names`, and the second column contains the batch group name.
