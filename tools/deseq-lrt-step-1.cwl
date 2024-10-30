cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: biowardrobe2/scidap-deseq:v0.0.30

inputs:

  expression_files:
    type: File[]
    inputBinding:
      position: 1
      prefix: "--input"
    doc: "Grouped by gene / TSS/ isoform expression files, formatted as CSV/TSV"

  expression_file_names:
    type: string[]
    inputBinding:
      position: 2
      prefix: "--name"
    doc: "Unique names for input files, no special characters, spaces are allowed. Number and order corresponds to --input"

  metadata_file:
    type: File
    inputBinding:
      position: 3
      prefix: "--meta"
    doc: "Metadata file to describe relation between samples, where first column corresponds to --name, formatted as CSV/TSV"

  design_formula:
    type: string
    inputBinding:
      position: 4
      prefix: "--design"
    doc: "Design formula. Should start with ~. See DESeq2 manual for details"

  reduced_formula:
    type: string
    inputBinding:
      position: 5
      prefix: "--reduced"
    doc: "Reduced formula to compare against with the term(s) of interest removed. Should start with ~. See DESeq2 manual for details"

  batchcorrection:
    type:
      - "null"
      - type: enum
        symbols:
          - "none"
          - "combatseq"
          - "limmaremovebatcheffect"
    inputBinding:
      position: 6
      prefix: "--batchcorrection"
    default: "none"
    doc: |
      Specifies the batch correction method to be applied.
      - 'combatseq' applies ComBat_seq at the beginning of the analysis.
      - 'limmaremovebatcheffect' notes the batch correction to be applied in step 2.
      - Default: none

  batchfile:
    type: File?
    inputBinding:
      position: 7
      prefix: "-bf"
    doc: |
      Metadata file for batch correction. Headerless TSV/CSV file.
      First column - sample names matching --name, second column - batch group name.

  fdr:
    type: float?
    inputBinding:
      position: 8
      prefix: "--fdr"
    default: 0.1
    doc: |
      In the exploratory visualization part of the analysis, output only features with adjusted p-value (FDR) not bigger than this value.
      Also, the significance cutoff used for optimizing the independent filtering. Default: 0.1.

  lfcthreshold:
    type: float?
    inputBinding:
      position: 9
      prefix: "--lfcthreshold"
    default: 0.59
    doc: |
      Log2 fold change threshold for determining significant differential expression.
      Genes with absolute log2 fold change greater than this threshold will be considered.
      Default: 0.59 (about 1.5 fold change)

  use_lfc_thresh:
    type: boolean
    inputBinding:
      position: 10
      prefix: "--use_lfc_thresh"
    default: true
    doc: "Use lfcthreshold as the null hypothesis value in the results function call. Default: TRUE"

  output_prefix:
    type: string?
    inputBinding:
      position: 11
      prefix: "--output"
    default: "./deseq"
    doc: "Output prefix for generated files"

  threads:
    type: int?
    inputBinding:
      position: 12
      prefix: "--threads"
    default: 1
    doc: "Number of threads"

  test_mode:
    type: boolean
    inputBinding:
      position: 13
      prefix: "--test_mode"
    default: false
    doc: "Run for test, only first 100 rows"

outputs:

  contrasts_file:
    type: File
    outputBinding:
      glob: "*_contrasts_table.tsv"

  contrasts_rds:
    type: File
    outputBinding:
      glob: "*_contrasts.rds"

  diff_expr_file:
    type: File
    outputBinding:
      glob: "*_gene_exp_table.tsv"

  expression_data_rds:
    type: File
    outputBinding:
      glob: "*_expression_data.rds"

  metadata_rds:
    type: File
    outputBinding:
      glob: "*_metadata.rds"

  read_counts_rds:
    type: File
    outputBinding:
      glob: "*_read_counts.rds"

  dsq_wald_rds:
    type: File
    outputBinding:
      glob: "*_dsq_wald.rds"

  dsq_lrt_rds:
    type: File
    outputBinding:
      glob: "*_dsq_lrt.rds"

  dsq_lrt_res_rds:
    type: File
    outputBinding:
      glob: "*_dsq_lrt_res.rds"

  batch_correction_method_rds:
    type: File
    outputBinding:
      glob: "*_batch_correction_method.rds"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr

baseCommand: [ run_deseq_lrt_step_1.R ]
stdout: deseq_stdout.log
stderr: deseq_stderr.log

$namespaces:
  s: http://schema.org/

$schemas:
  - https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "DESeq2 (LRT) - differential gene expression analysis using likelihood ratio test"
label: "DESeq2 (LRT) - differential gene expression analysis using likelihood ratio test"
s:alternateName: "Differential gene expression analysis based on the LRT (likelihood ratio test)"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/deseq-lrt-step-1.cwl
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
  Runs DESeq2 using LRT (Likelihood Ratio Test)

  The LRT examines two models for the counts, a full model with a certain number of terms and a reduced model,
  in which some of the terms of the full model are removed. The test determines if the increased likelihood of
  the data using the extra terms in the full model is more than expected if those extra terms are truly zero.

  The LRT is therefore useful for testing multiple terms at once, for example testing 3 or more levels of a factor
  at once, or all interactions between two variables. The LRT for count data is conceptually similar to an analysis
  of variance (ANOVA) calculation in linear regression, except that in the case of the Negative Binomial GLM, we use
  an analysis of deviance (ANODEV), where the deviance captures the difference in likelihood between a full and a
  reduced model.

  When one performs a likelihood ratio test, the p values and the test statistic (the stat column) are values for
  the test that removes all of the variables which are present in the full design and not in the reduced design.
  This tests the null hypothesis that all the coefficients from these variables and levels of these factors are
  equal to zero.

  The likelihood ratio test p values therefore represent a test of all the variables and all the levels of factors
  which are among these variables. However, the results table only has space for one column of log fold change, so
  a single variable and a single comparison is shown (among the potentially multiple log fold changes which were
  tested in the likelihood ratio test). This indicates that the p value is for the likelihood ratio test of all the
  variables and all the levels, while the log fold change is a single comparison from among those variables and levels.

  **Note:** At least two biological replicates are required for every compared category.

  **Input Files:**

  All input CSV/TSV files should have the following header (case-sensitive):

  - CSV: `RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm`
  - TSV: `RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm`

  The format of the input files is identified based on the file's extension:

  - `*.csv` - CSV
  - `*.tsv` - TSV
  - Otherwise, CSV format is assumed by default.

  **Metadata File Example:**

  ```csv
  ,time,condition
  DH1,day5,WT
  DH2,day5,KO
  DH3,day7,WT
  DH4,day7,KO
  DH5,day7,KO