cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  rnaseq_experiment:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-pe-smarter-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"


inputs:

  alias:
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
        symbols: ["isoforms", "genes", "common tss"]
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
    doc: "Design formula. Should start with ~. See DeSeq2 manual for details"

  reduced_formula:
    type: string
    label: "Reduced formula to compare against with the term(s) of interest removed. See workflow description for details"
    doc: "Reduced formula to compare against with the term(s) of interest removed. Should start with ~. See DeSeq2 manual for details"

  contrast:
    type: string
    label: "Contrast to be be applied for output, formatted as Factor Numerator Denominator"
    doc: "Contrast to be be applied for output, formatted as Factor Numerator Denominator or 'Factor Numerator Denominator'"

  threads:
    type: int?
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    default: 1
    'sd:layout':
      advanced: true


outputs:

  diff_expr_file:
    type: File
    label: "Differentially expressed features grouped by isoforms, genes or common TSS"
    format: "http://edamontology.org/format_3475"
    doc: "DESeq2 generated file of differentially expressed features grouped by isoforms, genes or common TSS in TSV format"
    outputSource: deseq/diff_expr_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Differential Expression Analysis'
        Title: 'Combined DESeq2 results'

  volcano_plot:
    type: File
    label: "Volcano plot"
    format: "http://edamontology.org/format_3475"
    doc: "TSV file with input data to build volcano plot - log2FoldChange vs -LOG10(padj)"
    outputSource: make_volcano_plot/output_file
    'sd:visualPlugins':
    - scatter:
        tab: 'Volcano Plot'
        Title: 'Volcano'
        xAxisTitle: 'log fold change'
        yAxisTitle: '-log10(pAdj)'
        colors: ["#b3de69"]
        height: 600
        data: [$1, $2]

  ma_plot:
    type: File
    label: "Plot of normalised mean versus log2 fold change"
    format: "http://edamontology.org/format_3603"
    doc: "Plot of the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples"
    outputSource: deseq/ma_plot
    'sd:visualPlugins':
    - image:
        tab: 'Other Plots'
        Caption: 'LFC vs mean'

  deseq_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "DeSeq2 stdout log"
    doc: "DeSeq2 stdout log"
    outputSource: deseq/stdout_log

  deseq_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "DeSeq2 stderr log"
    doc: "DeSeq2 stderr log"
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
    run: ../tools/deseq-lrt.cwl
    in:
      expression_files:
        source: [group_by, expression_files, group_isoforms/genes_file, group_isoforms/common_tss_file]
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
      contrast: contrast
      threads: threads
    out: [diff_expr_file, ma_plot, stdout_log, stderr_log]

  make_volcano_plot:
    run: ../tools/custom-bash.cwl
    in:
      input_file: deseq/diff_expr_file
      script:
        default: |
          cat "$0" | awk -F "\t" 'BEGIN {xl="log2FoldChange"; yl="-LOG10(padj)"} NR==1 {for (i=1; i<=NF; i++) {ix[$i]=i} print xl"\t"yl } NR>1 {print $ix[xl]"\t"$ix[yl]}' > volcano_plot_data.tsv
    out: [output_file]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "DESeq2 (LRT) - differential gene expression analysis using likelihood ratio test"
label:  "DESeq2 (LRT) - differential gene expression analysis using likelihood ratio test"
s:alternateName: "Differential gene expression analysis based on the LRT (likelihood ratio test)"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/deseq-lrt.cwl
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


# doc:
#   $include: ../descriptions/deseq-lrt.md


doc: |
  Runs DESeq2 using LRT (Likelihood Ratio Test)
  =============================================

  The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model, in which some of the terms of the full model are removed. The test determines if the increased likelihood of the data using the extra terms in the full model is more than expected if those extra terms are truly zero.

  The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor at once, or all interactions between two variables. The LRT for count data is conceptually similar to an analysis of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a reduced model.

  When one performs a likelihood ratio test, the p values and the test statistic (the stat column) are values for the test that removes all of the variables which are present in the full design and not in the reduced design. This tests the null hypothesis that all the coefficients from these variables and levels of these factors are equal to zero.

  The likelihood ratio test p values therefore represent a test of all the variables and all the levels of factors which are among these variables. However, the results table only has space for one column of log fold change, so a single variable and a single comparison is shown (among the potentially multiple log fold changes which were tested in the likelihood ratio test). This indicates that the p value is for the likelihood ratio test of all the variables and all the levels, while the log fold change is a single comparison from among those variables and levels.

  **Technical notes**

  1. At least two biological replicates are required for every compared category
  2. Metadata file describes relations between compared experiments, for example

     ```
      ,time,condition
      DH1,day5,WT
      DH2,day5,KO
      DH3,day7,WT
      DH4,day7,KO
      DH5,day7,KO
     ```
     where `time, condition, day5, day7, WT, KO` should be a single words (without spaces) and `DH1, DH2, DH3, DH4, DH5` correspond to the experiment aliases set in **RNA-Seq experiments** input.
  3. Design and reduced formulas should start with **~** and include categories or, optionally, their interactions from the metadata file header. See details in DESeq2 manual [here](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#interactions) and [here](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#likelihood-ratio-test)
  4. Contrast should be set based on your metadata file header and available categories in a form of `Factor Numerator Denominator`, where `Factor` - column name from metadata file, `Numerator`  - category from metadata file to be used as numerator in fold change calculation, `Denominator` - category from metadata file to be used as denominator in fold change calculation. For example `condition WT KO`.