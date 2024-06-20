cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap-deseq:v0.0.28


inputs:

  untreated_files:
    type:
      - File
      - File[]
    inputBinding:
      prefix: "-u"
    doc: |
      Untreated input CSV/TSV files

  treated_files:
    type:
      - File
      - File[]
    inputBinding:
      prefix: "-t"
    doc: |
      Treated input CSV/TSV files

  untreated_name:
    type: string?
    inputBinding:
      prefix: "-un"
    doc: |
      Name for untreated condition, use only letters and numbers

  treated_name:
    type: string?
    inputBinding:
      prefix: "-tn"
    doc: |
      Name for treated condition, use only letters and numbers

  untreated_sample_names:
    type:
      - "null"
      - string
      - string[]
    inputBinding:
      prefix: "-ua"
    doc: |
      Unique aliases for untreated expression files. Default: basenames of -u without extensions

  treated_sample_names:
    type:
      - "null"
      - string
      - string[]
    inputBinding:
      prefix: "-ta"
    doc: |
      Unique aliases for treated expression files. Default: basenames of -t without extensions

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
      exploratory visualization part of the analysis. Default: do not run
      clustering

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

  fdr:
    type: float?
    inputBinding:
      prefix: "--fdr"
    doc: |
      In the exploratory visualization part of the analysis output only features,
      with adjusted p-value (FDR) not bigger than this value. Also the significance,
      cutoff used for optimizing the independent filtering. Default: 0.1.
  
  lfcthreshold:
    type: float?
    inputBinding:
      prefix: "--lfcthreshold"
    doc: |
      Log2 fold change threshold for determining significant differential expression.
      Genes with absolute log2 fold change greater than this threshold will be considered.
      Default: 0.59 (about 1.5 fold change)
      
  regulation:
    type:
      - "null"
      - type: enum
        symbols:
          - "up"
          - "down"
          - "both"
    inputBinding:
      prefix: "--regulation"
    doc: |
      Direction of differential expression comparison.
      'up' for upregulated genes, 'down' for downregulated genes,
      'both' for both up and downregulated genes. Default: both
    default: "both"
    
  batchcorrection:
    type:
      - "null"
      - type: enum
        symbols:
          - "combatseq"
          - "limma"
    inputBinding:
      prefix: "--batchcorrection"
    doc: |
      Batch correction method to be applied.
      'combatseq' applies ComBat_seq at the beginning of the analysis.
      'limma' applies removeBatchEffect after differential expression analysis (DEA). Default: combatseq
    default: "combatseq"
    
  batch_file:
    type: File?
    inputBinding:
      prefix: "-bf"
    doc: |
      Metadata file for multi-factor analysis. Headerless TSV/CSV file.
      First column - names from --ua and --ta, second column - batch name.
      Default: None

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
      Run script using multiple threads


outputs:

  diff_expr_file:
    type: File
    outputBinding:
      glob: "*report.tsv"

  read_counts_file:
    type: File
    outputBinding:
      glob: "*counts.gct"

  phenotypes_file:
    type: File
    outputBinding:
      glob: "*phenotypes.cls"

  plot_lfc_vs_mean:
    type: File?
    outputBinding:
      glob: "*_ma_plot.png"

  gene_expr_heatmap:
    type: File?
    outputBinding:
      glob: "*_expression_heatmap.png"

  plot_pca:
    type: File?
    outputBinding:
      glob: "*_pca_plot.png"

  plot_lfc_vs_mean_pdf:
    type: File?
    outputBinding:
      glob: "*_ma_plot.pdf"

  gene_expr_heatmap_pdf:
    type: File?
    outputBinding:
      glob: "*_expression_heatmap.pdf"

  plot_pca_pdf:
    type: File?
    outputBinding:
      glob: "*_pca_plot.pdf"

  mds_plot_html:
    type: File?
    outputBinding:
      glob: "*_mds_plot.html"
    doc: |
      MDS plot of normalized counts. Optionally batch corrected
      based on the --remove value.
      HTML format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: [run_deseq.R]
stdout: deseq_stdout.log
stderr: deseq_stderr.log

$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/deseq-metadata.yaml

s:name: "deseq-advanced"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/deseq-advanced.cwl
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
  Tool runs DESeq/DESeq2 script similar to the original one from BioWArdrobe.
  untreated_files and treated_files input files should have the following header (case-sensitive)
  <RefseqId,GeneId,Chrom,TxStart,TxEnd,Strand,TotalReads,Rpkm>         - CSV
  <RefseqId\tGeneId\tChrom\tTxStart\tTxEnd\tStrand\tTotalReads\tRpkm>  - TSV

  Format of the input files is identified based on file's extension
  *.csv - CSV
  *.tsv - TSV
  Otherwise used CSV by default

  The output file's rows order corresponds to the rows order of the first CSV/TSV file in
  the untreated group. Output is always saved in TSV format

  Output file includes only intersected rows from all input files. Intersected by
  RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand

  DESeq/DESeq2 always compares untreated_vs_treated groups.
  Normalized read counts and phenotype table are exported as GCT and CLS files for GSEA downstream analysis.


s:about: |
  usage: /Users/kot4or/workspaces/cwl_ws/workflows/tools/dockerfiles/scripts/run_deseq.R
        [-h] -u UNTREATED [UNTREATED ...] -t TREATED [TREATED ...]
        [-ua [UALIAS ...]] [-ta [TALIAS ...]] [-un UNAME] [-tn TNAME]
        [-bf BATCHFILE] [-cu CUTOFF] [--fdr FDR]
        [--regulation {up,down,both}]
        [--lfcthreshold LFCTHRESHOLD]
        [--batchcorrection {combatseq,limma}]
        [--cluster {row,column,both}]
        [--rowdist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}]
        [--columndist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}]
        [-o OUTPUT] [-d DIGITS] [-p THREADS]

  Run BioWardrobe DESeq/DESeq2 for untreated-vs-treated groups (condition-1-vs-
  condition-2)

  options:
    -h, --help            show this help message and exit
    -u UNTREATED [UNTREATED ...], --untreated UNTREATED [UNTREATED ...]
                          Untreated (condition 1) CSV/TSV isoforms expression
                          files
    -t TREATED [TREATED ...], --treated TREATED [TREATED ...]
                          Treated (condition 2) CSV/TSV isoforms expression
                          files
    -ua [UALIAS ...], --ualias [UALIAS ...]
                          Unique aliases for untreated (condition 1) expression
                          files. Default: basenames of -u without extensions
    -ta [TALIAS ...], --talias [TALIAS ...]
                          Unique aliases for treated (condition 2) expression
                          files. Default: basenames of -t without extensions
    -un UNAME, --uname UNAME
                          Name for untreated (condition 1), use only letters and
                          numbers
    -tn TNAME, --tname TNAME
                          Name for treated (condition 2), use only letters and
                          numbers
    -bf BATCHFILE, --batchfile BATCHFILE
                          Metadata file for multi-factor analysis. Headerless
                          TSV/CSV file. First column - names from --ualias and
                          --talias, second column - batch group name. Default:
                          None
    --fdr FDR             In the exploratory visualization part of the analysis output only features,
                          with adjusted p-value (FDR) not bigger than this value. Also the significance,
                          cutoff used for optimizing the independent filtering. Default: 0.1.
    --cluster {row,column,both}
                          Hopach clustering method to be run on normalized read
                          counts for the exploratory visualization part of the
                          analysis. Default: do not run clustering
    --rowdist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}
                          Distance metric for HOPACH row clustering. Ignored if
                          --cluster is not provided. Default: cosangle
    --columndist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}
                          Distance metric for HOPACH column clustering. Ignored
                          if --cluster is not provided. Default: euclid
    -o OUTPUT, --output OUTPUT
                          Output prefix. Default: deseq
    -d DIGITS, --digits DIGITS
                          Precision, number of digits to print. Default: 3
    -p THREADS, --threads THREADS
                          Threads
    --lfcthreshold LFCTHRESHOLD
                          Log2 fold change threshold for determining significant differential expression.
                          Genes with absolute log2 fold change greater than this threshold will be considered.
                          Default: 0.59 (about 1.5 fold change)
    --regulation {up,down,both}
                          Direction of differential expression comparison.
                          'up' for upregulated genes, 'down' for downregulated genes,
                          'both' for both up and downregulated genes. Default: both
    --batchcorrection {combatseq,limma}
                          Batch correction method to be applied.
                          'combatseq' applies ComBat_seq at the beginning of the analysis.
                          'limma' applies removeBatchEffect after differential expression analysis (DEA). Default: combatseq