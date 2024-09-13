cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement

hints:
  - class: DockerRequirement
    dockerPull: biowardrobe2/scidap-deseq:v0.0.28

inputs:

  expression_files:
    type: File[]
    inputBinding:
      position: 5
      prefix: "--input"
    doc: "Grouped by gene / TSS/ isoform expression files, formatted as CSV/TSV"

  expression_file_names:
    type: string[]
    inputBinding:
      position: 6
      prefix: "--name"
    doc: "Unique names for input files, no special characters, spaces are allowed. Number and order corresponds to --input"

  metadata_file:
    type: File
    inputBinding:
      position: 7
      prefix: "--meta"
    doc: "Metadata file to describe relation between samples, where first column corresponds to --name, formatted as CSV/TSV"

  design_formula:
    type: string
    inputBinding:
      prefix: "--design_formula"
    doc: |
      Design formula for DESeq2. Must be provided as a string, e.g., "~ condition + batch".

  contrast_indices:
    type:
      - string
    inputBinding:
      prefix: "--contrast_indices"
    doc: |
      Comma-separated list of integers representing contrast indices (e.g., 1,2,3).

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
      exploratory visualization part of the analysis. Default: do not run clustering

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
      Distance metric for HOPACH row clustering. Ignored if --cluster is not provided. Default: cosangle

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
      Distance metric for HOPACH column clustering. Ignored if --cluster is not provided. Default: euclid

  test_mode:
    type: boolean
    inputBinding:
      position: 15
      prefix: "--test_mode"
    default: false
    doc: "Run for test, only first 100 rows"

  fdr:
    type: float?
    inputBinding:
      prefix: "--fdr"
    doc: |
      In the exploratory visualization part of the analysis, output only features with adjusted p-value (FDR) not bigger than this value. Default: 0.1.

  lfcthreshold:
    type: float?
    inputBinding:
      prefix: "--lfcthreshold"
    doc: |
      Log2 fold change threshold for determining significant differential expression. Default: 0.59 (about 1.5 fold change)

  use_lfc_thresh:
    type: boolean
    inputBinding:
      prefix: "--use_lfc_thresh"
    default: true
    doc: |
      Use lfcthreshold as the null hypothesis value in the results function call. Default: TRUE

  batchcorrection:
    type:
      - "null"
      - type: enum
        symbols:
          - "none"
          - "combatseq"
          - "limmaremovebatcheffect"
    inputBinding:
      prefix: "--batchcorrection"
    doc: |
      Specifies the batch correction method to be applied. Default: none
    default: "none"

  batch_file:
    type: File?
    inputBinding:
      prefix: "-bf"
    doc: |
      Metadata file for multi-factor analysis. Headerless TSV/CSV file.
      First column - names from --ualias and --talias, second column - batch name.

  output_prefix:
    type: string?
    inputBinding:
      prefix: "-o"
    doc: |
      Output prefix. Default: deseq

  threads:
    type: int?
    inputBinding:
      prefix: '-p'
    doc: |
      Number of threads to run the analysis.

outputs:

  diff_expr_file:
    type: File[]
    outputBinding:
      glob: "*_gene_exp_table.tsv"

  read_counts_file_all:
    type: File[]
    outputBinding:
      glob: "*counts_all.gct"

  read_counts_file_filtered:
    type: File[]
    outputBinding:
      glob: "*counts_filtered.gct"

  plot_lfc_vs_mean:
    type: File[]
    outputBinding:
      glob: "*_ma_plot.png"

  gene_expr_heatmap:
    type: File[]
    outputBinding:
      glob: "*_expression_heatmap.png"

  plot_lfc_vs_mean_pdf:
    type: File[]
    outputBinding:
      glob: "*_ma_plot.pdf"

  gene_expr_heatmap_pdf:
    type: File[]
    outputBinding:
      glob: "*_expression_heatmap.pdf"

  mds_plot_html:
    type: File[]
    outputBinding:
      glob: "*_mds_plot.html"
    doc: |
      MDS plot of normalized counts, optionally batch-corrected.

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr

baseCommand: [ run_deseq_lrt_step_2.R ]
stdout: deseq_stdout.log
stderr: deseq_stderr.log

$namespaces:
  s: http://schema.org/

$schemas:
  - https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/deseq-metadata.yaml

s:name: "deseq-lrt-step-2"
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
  Tool runs DESeq/DESeq2 analysis, originally inspired by the BioWardrobe workflow. It supports complex experimental designs, batch correction, and multiple contrast generation. The tool takes untreated and treated CSV/TSV input files, applies differential gene expression analysis, and produces visualizations such as MA plots and heatmaps, as well as GCT/CLS files for downstream GSEA.

  Input files (untreated_files and treated_files) should have the following header (case-sensitive):
  <RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand, TotalReads, Rpkm>   (for CSV)
  <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>   (for TSV)

  The format of the input files is identified based on the file extension:
  - *.csv for CSV files
  - *.tsv for TSV files
  By default, CSV format is assumed if the file extension is unrecognized.

  The output file order corresponds to the rows in the first untreated group CSV/TSV file, and all output is saved in TSV format. Rows included in the output represent the intersection of input files based on the following fields:
  RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand

  DESeq2 compares untreated and treated groups and performs differential gene expression analysis. Normalized read counts and phenotype tables are exported in GCT and CLS formats for downstream GSEA.

  You can specify the design formula for DESeq2 with `--design_formula` and generate multiple contrasts using the `--contrast_indices` argument.

  The tool supports batch correction through two options:
  - `combatseq`: applies batch correction before the differential analysis.
  - `limmaremovebatcheffect`: applies batch correction after differential analysis.

  Visual outputs include MA plots and heatmaps, and clustering options are available through HOPACH for rows, columns, or both.

s:about: |
  usage: run_deseq_lrt_step_2.R
        [-h] -u UNTREATED [UNTREATED ...] -t TREATED [TREATED ...]
        [-ua [UALIAS ...]] [-ta [TALIAS ...]] [-un UNAME] [-tn TNAME]
        [--design_formula DESIGN_FORMULA] [--contrast_indices CONTRAST_INDICES]
        [-bf BATCHFILE] [-cu CUTOFF] [--fdr FDR]
        [--lfcthreshold LFCTHRESHOLD]
        [--batchcorrection {none, combatseq, limmaremovebatcheffect}]
        [--use_lfc_thresh]

  Run DESeq/DESeq2 for untreated-vs-treated groups (condition-1-vs-condition-2) using complex experimental designs and contrasts.

  Available options:
    -h, --help                Show this help message and exit
    -u UNTREATED [UNTREATED ...], --untreated UNTREATED [UNTREATED ...]
                              Untreated (condition 1) CSV/TSV isoform expression files
    -t TREATED [TREATED ...], --treated TREATED [TREATED ...]
                              Treated (condition 2) CSV/TSV isoform expression files
    -ua [UALIAS ...], --ualias [UALIAS ...]
                              Unique aliases for untreated (condition 1) expression files.
                              Default: basenames of -u without extensions
    -ta [TALIAS ...], --talias [TALIAS ...]
                              Unique aliases for treated (condition 2) expression files.
                              Default: basenames of -t without extensions
    -un UNAME, --uname UNAME   Name for untreated (condition 1), only letters and numbers
    -tn TNAME, --tname TNAME   Name for treated (condition 2), only letters and numbers
    --design_formula DESIGN_FORMULA
                              Design formula for DESeq2 (e.g., "~ condition + batch").
    --contrast_indices CONTRAST_INDICES
                              Comma-separated list of integers representing contrast indices (e.g., 1,2,3).
    -bf BATCHFILE, --batchfile BATCHFILE
                              Metadata file for multi-factor analysis. Headerless TSV/CSV file.
                              First column - names from --ualias and --talias, second column - batch group name. Default: None
    --fdr FDR                 Filter based on adjusted p-value (FDR), default: 0.1
    --lfcthreshold LFCTHRESHOLD
                              Log2 fold change threshold for significant differential expression. Default: 0.59
    --batchcorrection {none, combatseq, limmaremovebatcheffect}
                              Batch correction method, default: none
    --use_lfc_thresh           Use lfcthreshold in the hypothesis test. Default: TRUE
    --cluster {row,column,both}
                              Hopach clustering method for rows, columns, or both (optional)
    --rowdist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}
                              Distance metric for HOPACH row clustering. Default: cosangle
    --columndist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}
                              Distance metric for HOPACH column clustering. Default: euclid
    -o OUTPUT, --output OUTPUT Output file prefix. Default: deseq
    -p THREADS, --threads THREADS Number of threads for parallel processing