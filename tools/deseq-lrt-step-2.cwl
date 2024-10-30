cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: biowardrobe2/scidap-deseq:v0.0.30

inputs:

  expression_data_rds:
    type: File
    inputBinding:
      position: 1
      prefix: "--expression_data_rds"
    doc: "RDS file containing the expression data from step 1"

  contrasts_rds:
    type: File
    inputBinding:
      position: 2
      prefix: "--contrasts_rds"
    doc: "RDS file containing the contrasts list from step 1"

  dsq_wald_rds:
    type: File
    inputBinding:
      position: 3
      prefix: "--dsq_wald_rds"
    doc: "RDS file containing the DESeq2 object from the Wald test in step 1"

  metadata_rds:
    type: File
    inputBinding:
      position: 4
      prefix: "--metadata_rds"
    doc: "RDS file containing the metadata from step 1"

  batch_correction_method_rds:
    type: File
    inputBinding:
      position: 5
      prefix: "--batch_correction_method_rds"
    doc: "RDS file containing the batch correction method used in step 1"

  contrast_indices:
    type: string
    inputBinding:
      position: 6
      prefix: "--contrast_indices"
    doc: "Comma-separated list of integers representing contrast indices (e.g., 1,2,3)"

  fdr:
    type: float?
    inputBinding:
      position: 7
      prefix: "--fdr"
    default: 0.1
    doc: |
      In the exploratory visualization part of the analysis, output only features with adjusted p-value (FDR) not bigger than this value.
      Default: 0.1.

  lfcthreshold:
    type: float?
    inputBinding:
      position: 8
      prefix: "--lfcthreshold"
    default: 0.59
    doc: |
      Log2 fold change threshold for determining significant differential expression.
      Genes with absolute log2 fold change greater than this threshold will be considered.
      Default: 0.59 (about 1.5 fold change)

  regulation:
    type:
      - "null"
      - type: enum
        symbols:
          - "both"
          - "up"
          - "down"
    inputBinding:
      position: 9
      prefix: "--regulation"
    default: "both"
    doc: |
      Direction of differential expression comparison. β is the log2 fold change.
      - 'both' for both up and downregulated genes (|β| > lfcThreshold);
      - 'up' for upregulated genes (β > lfcThreshold);
      - 'down' for downregulated genes (β < -lfcThreshold).
      Default: both

  output_prefix:
    type: string?
    inputBinding:
      position: 10
      prefix: "--output"
    default: "./deseq"
    doc: "Output prefix"

  threads:
    type: int?
    inputBinding:
      position: 11
      prefix: "--threads"
    default: 1
    doc: "Number of threads"

  test_mode:
    type: boolean
    inputBinding:
      position: 12
      prefix: "--test_mode"
    default: false
    doc: "Run for test, only first 100 rows"

outputs:

  diff_expr_files:
    type: File[]
    outputBinding:
      glob: "*_contrast_*_gene_exp_table.tsv"

  ma_plots_png:
    type: File[]
    outputBinding:
      glob: "*_contrast_*_ma_plot.png"

  ma_plots_pdf:
    type: File[]
    outputBinding:
      glob: "*_contrast_*_ma_plot.pdf"

  heatmaps_png:
    type: File[]
    outputBinding:
      glob: "*_contrast_*_heatmap.png"

  heatmaps_pdf:
    type: File[]
    outputBinding:
      glob: "*_contrast_*_heatmap.pdf"

  mds_plots_html:
    type: File[]
    outputBinding:
      glob: "*_contrast_*_mds_plot.html"

  counts_all_gct:
    type: File[]
    outputBinding:
      glob: "*_contrast_*_counts_all.gct"

  counts_filtered_gct:
    type: File[]
    outputBinding:
      glob: "*_contrast_*_counts_filtered.gct"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr

baseCommand: [ run_deseq_lrt_step_2.R ]
stdout: deseq_step2_stdout.log
stderr: deseq_step2_stderr.log

$namespaces:
  s: http://schema.org/

$schemas:
  - https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "DESeq2 (LRT Step 2) - Differential gene expression analysis using contrasts"
label: "DESeq2 (LRT Step 2) - Differential gene expression analysis using contrasts"
s:alternateName: "Differential gene expression analysis using DESeq2 contrasts"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/deseq-lrt-step-2.cwl
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
  Runs DESeq2 analysis using contrasts from previous LRT step.

  This tool takes the outputs from DESeq2 LRT step 1 and performs differential expression analysis using specified contrasts.

  **Inputs:**

  - `expression_data_rds`: RDS file containing the expression data from step 1.
  - `contrasts_rds`: RDS file containing the contrasts list from step 1.
  - `dsq_wald_rds`: RDS file containing the DESeq2 object from the Wald test in step 1.
  - `metadata_rds`: RDS file containing the metadata from step 1.
  - `batch_correction_method_rds`: RDS file containing the batch correction method used in step 1.
  - `contrast_indices`: Comma-separated list of integers representing contrast indices (e.g., 1,2,3).
  - `fdr`: Adjusted p-value (FDR) threshold for significance. Default: 0.1.
  - `lfcthreshold`: Log2 fold change threshold for significance. Default: 0.59.
  - `regulation`: Direction of differential expression comparison ('both', 'up', 'down'). Default: 'both'.

  **Outputs:**

  - Differential expression reports for each contrast.
  - MA plots (PNG and PDF) for each contrast.
  - Heatmaps (PNG and PDF) for each contrast.
  - MDS plots (HTML) for each contrast.
  - GCT files (counts_all.gct and counts_filtered.gct) for each contrast.
