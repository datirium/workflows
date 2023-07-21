cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.26


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object from.
      This file should include chromatin accessibility
      information stored in the ATAC assay. Additionally
      'rnaumap', and/or 'atacumap', and/or 'wnnumap'
      dimensionality reductions should be present.

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    inputBinding:
      prefix: "--fragments"
    doc: |
      Count and barcode information for every ATAC fragment
      used in the loaded Seurat object. File should be saved
      in TSV format with tbi-index file.

  datasets_metadata:
    type: File?
    inputBinding:
      prefix: "--metadata"
    doc: |
      Path to the TSV/CSV file to optionally extend Seurat
      object metadata with categorical values using samples
      identities. First column - 'library_id' should
      correspond to all unique values from the 'new.ident'
      column of the loaded Seurat object. If any of the
      provided in this file columns are already present in
      the Seurat object metadata, they will be overwritten.
      When combined with --barcodes parameter, first the
      metadata will be extended, then barcode filtering will
      be applied. Default: no extra metadata is added

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and
      extend Seurat object metadata by selected barcodes.
      First column should be named as 'barcode'. If file
      includes any other columns they will be added to the
      Seurat object metadata ovewriting the existing ones if
      those are present. Default: all cells used, no extra
      metadata is added

  groupby:
    type: string?
    inputBinding:
      prefix: "--groupby"
    doc: |
      Column from the Seurat object metadata to group cells
      for optional subsetting when combined with --subset
      parameter. May be one of the extra metadata columns
      added with --metadata or --barcodes parameters.
      Ignored if --subset is not set. Default: do not
      subset, include all cells into analysis.

  subset:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--subset"
    doc: |
      Values from the column set with --groupby parameter to
      subset cells before running differential binding
      analysis. Ignored if --groupby is not provided.
      Default: do not subset cells, include all of them.

  splitby:
    type: string
    inputBinding:
      prefix: "--splitby"
    doc: |
      Column from the Seurat object metadata to split cells
      into two groups to run --second vs --first
      differential binding analysis. May be one of the extra
      metadata columns added with --metadata or --barcodes
      parameters.

  first_cond:
    type: string
    inputBinding:
      prefix: "--first"
    doc: |
      Value from the Seurat object metadata column set with
      --splitby parameter to define the first group of cells
      for differential binding analysis.

  second_cond:
    type: string
    inputBinding:
      prefix: "--second"
    doc: |
      Value from the Seurat object metadata column set with
      --splitby parameter to define the second group of
      cells for differential binding analysis.

  analysis_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "negative-binomial"          # (negbinom) Negative Binomial Generalized Linear Model (use FindMarkers with peaks from Seurat object)
      - "poisson"                    # (poisson) Poisson Generalized Linear Model (use FindMarkers with peaks from Seurat object)
      - "logistic-regression"        # (LR) Logistic Regression (use FindMarkers with peaks from Seurat object)
      - "mast"                       # (MAST) MAST package (use FindMarkers with peaks from Seurat object)
      - "manorm2"                    # call peaks for each group with MACS2, run MAnorm2
    inputBinding:
      prefix: "--test"
    doc: |
      Test type to use in differential binding analysis. For
      all tests except manorm2, peaks present in the loaded
      Seurat object will be used. If manorm2 test selected,
      peaks will be called per group defined by --splitby
      parameter. Default: logistic-regression

  genome_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "hs"
      - "mm"
    inputBinding:
      prefix: "--genome"
    doc: |
      Genome type of the sequencing data loaded from the
      Seurat object. It will be used for effective genome
      size selection when calling peaks with MACS2. Ignored
      if --test is not set to manorm2. Default: hs (2.7e9)

  minimum_qvalue:
    type: float?
    inputBinding:
      prefix: "--qvalue"
    doc: |
      Minimum FDR (q-value) cutoff for MACS2 peak detection.
      Ignored if --test is not set to manorm2. Default: 0.05

  minimum_peak_gap:
    type: int?
    inputBinding:
      prefix: "--minpeakgap"
    doc: |
      If a distance between peaks is smaller than the
      provided value they will be merged before splitting
      them into reference genomic bins of size --binsize.
      Ignored if --test is not set to manorm2. Default: 150

  bin_size:
    type: int?
    inputBinding:
      prefix: "--binsize"
    doc: |
      The size of non-overlapping reference genomic bins
      used by MAnorm2 when generating a table of reads
      counts per peaks. Ignored if --test is not set to
      manorm2. Default: 1000

  maximum_peaks:
    type: int?
    inputBinding:
      prefix: "--maxpeaks"
    doc: |
      The maximum number of the most significant (based on
      qvalue) peaks to keep from each group of cells when
      constructing reference genomic bins. Ignored if --test
      is not set to manorm2. Default: keep all peaks

  blacklist_regions_file:
    type: File?
    inputBinding:
      prefix: "--blacklist"
    doc: |
      Path to the optional BED file with the genomic
      blacklist regions to be filtered out before running
      differential binding analysis. Any reference genomic
      bin overlapping a blacklist region will be removed
      from the output. Ignored if --test is not set to
      manorm2.

  maximum_padj:
    type: float?
    inputBinding:
      prefix: "--padj"
    doc: |
      In the exploratory visualization part of the analysis
      output only differentially bound peaks with adjusted
      P-value not bigger than this value. Default: 0.05

  minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--logfc"
    doc: |
      In the exploratory visualization part of the analysis
      output only differentially bound peaks with log2 Fold
      Change not smaller than this value. Default: 1.0

  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: "--pdf"
    doc: |
      Export plots in PDF.
      Default: false

  color_theme:
    type:
    - "null"
    - type: enum
      symbols:
      - "gray"
      - "bw"
      - "linedraw"
      - "light"
      - "dark"
      - "minimal"
      - "classic"
      - "void"
    inputBinding:
      prefix: "--theme"
    doc: |
      Color theme for all generated plots. One of gray, bw,
      linedraw, light, dark, minimal, classic, void.
      Default: classic

  verbose:
    type: boolean?
    inputBinding:
      prefix: "--verbose"
    doc: |
      Print debug information.
      Default: false

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix.
      Default: ./sc

  parallel_memory_limit:
    type: int?
    inputBinding:
      prefix: "--memory"
    doc: |
      Maximum memory in GB allowed to be shared between the workers
      when using multiple --cpus.
      Default: 32

  vector_memory_limit:
    type: int?
    default: 128
    doc: |
      Maximum vector memory in GB allowed to be used by R.
      Default: 128

  threads:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of cores/cpus to use.
      Default: 1


outputs:

  umap_rd_rnaumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_rd_rnaumap.png"
    doc: |
      Cells UMAP split by selected criteria,
      optionally subsetted to the specific
      group (rnaumap dim. reduction).
      PNG format

  umap_rd_rnaumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_rd_rnaumap.pdf"
    doc: |
      Cells UMAP split by selected criteria,
      optionally subsetted to the specific
      group (rnaumap dim. reduction).
      PDF format

  umap_rd_atacumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_rd_atacumap.png"
    doc: |
      Cells UMAP split by selected criteria,
      optionally subsetted to the specific
      group (atacumap dim. reduction).
      PNG format

  umap_rd_atacumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_rd_atacumap.pdf"
    doc: |
      Cells UMAP split by selected criteria,
      optionally subsetted to the specific
      group (atacumap dim. reduction).
      PDF format

  umap_rd_wnnumap_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_rd_wnnumap.png"
    doc: |
      Cells UMAP split by selected criteria,
      optionally subsetted to the specific
      group (wnnumap dim. reduction).
      PNG format

  umap_rd_wnnumap_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_umap_rd_wnnumap.pdf"
    doc: |
      Cells UMAP split by selected criteria,
      optionally subsetted to the specific
      group (wnnumap dim. reduction).
      PDF format

  seurat_peaks_bigbed_file:
    type: File?
    outputBinding:
      glob: "*_seurat_peaks.bigBed"
    doc: |
      Peaks in bigBed format extracted
      from the loaded from provided RDS
      file Seurat object. 

  first_fragments_bigwig_file:
    type: File
    outputBinding:
      glob: "*_first.bigWig"
    doc: |
      Genome coverage in bigWig format calculated
      for fragments from the cells that belong to
      the group defined by the --first and
      --groupby parameters.
      
  second_fragments_bigwig_file:
    type: File
    outputBinding:
      glob: "*_second.bigWig"
    doc: |
      Genome coverage in bigWig format calculated
      for fragments from the cells that belong to
      the group defined by the --second and
      --groupby parameters.

  first_tn5ct_bigwig_file:
    type: File?
    outputBinding:
      glob: "*_first_tn5ct.bigWig"
    doc: |
      Genome coverage in bigWig format calculated
      for Tn5 cut sites from the cells that belong
      to the group defined by the --first and
      --groupby parameters.

  second_tn5ct_bigwig_file:
    type: File?
    outputBinding:
      glob: "*_second_tn5ct.bigWig"
    doc: |
      Genome coverage in bigWig format calculated
      for Tn5 cut sites from the cells that belong
      to the group defined by the --second and
      --groupby parameters.

  first_peaks_xls_file:
    type: File?
    outputBinding:
      glob: "*_first_peaks.xls"
    doc: |
      MACS2 report in XLS format for peaks
      called from the Tn5 cut sites of the
      cells that belong to the group defined
      by the --first and --groupby parameters.

  second_peaks_xls_file:
    type: File?
    outputBinding:
      glob: "*_second_peaks.xls"
    doc: |
      MACS2 report in XLS format for peaks
      called from the Tn5 cut sites of the
      cells that belong to the group defined
      by the --second and --groupby parameters.

  first_peaks_bed_file:
    type: File?
    outputBinding:
      glob: "*_first_peaks.narrowPeak"
    doc: |
      MACS2 peaks in narrowPeak format called
      from the Tn5 cut sites of the cells that
      belong to the group defined by the --first
      and --groupby parameters.

  second_peaks_bed_file:
    type: File?
    outputBinding:
      glob: "*_second_peaks.narrowPeak"
    doc: |
      MACS2 peaks in narrowPeak format called
      from the Tn5 cut sites of the cells that
      belong to the group defined by the --second
      and --groupby parameters.

  first_summits_bed_file:
    type: File?
    outputBinding:
      glob: "*_first_summits.bed"
    doc: |
      MACS2 peaks summits in BED format called
      from the Tn5 cut sites of the cells that
      belong to the group defined by the --first
      and --groupby parameters.

  second_summits_bed_file:
    type: File?
    outputBinding:
      glob: "*_second_summits.bed"
    doc: |
      MACS2 peaks summits in BED format called
      from the Tn5 cut sites of the cells that
      belong to the group defined by the --second
      and --groupby parameters.

  diff_bound_sites:
    type: File
    outputBinding:
      glob: "*_db_sites.tsv"
    doc: |
      Not filtered differentially bound sites
      in TSV format

  dbnd_vlcn_plot_png:
    type: File?
    outputBinding:
      glob: "*_dbnd_vlcn.png"
    doc: |
      Volcano plot of differentially bound sites.
      PNG format

  dbnd_vlcn_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_dbnd_vlcn.pdf"
    doc: |
      Volcano plot of differentially bound sites.
      PDF format

  first_enrch_bigbed_file:
    type: File?
    outputBinding:
      glob: "*_first_enrch.bigBed"
    doc: |
      Peaks in bigBed format filtered by
      --padj and --logfc thresholds enriched
      in the group of cells defined by the
      --first and --groupby parameters.
      
  second_enrch_bigbed_file:
    type: File?
    outputBinding:
      glob: "*_second_enrch.bigBed"
    doc: |
      Peaks in bigBed format filtered by
      --padj and --logfc thresholds enriched
      in the group of cells defined by the
      --second and --groupby parameters.

  first_enrch_bed_file:
    type: File?
    outputBinding:
      glob: "*_first_enrch.bed"
    doc: |
      Peaks in BED format filtered by
      --padj and --logfc thresholds enriched
      in the group of cells defined by the
      --first and --groupby parameters.
      
  second_enrch_bed_file:
    type: File?
    outputBinding:
      glob: "*_second_enrch.bed"
    doc: |
      Peaks in BED format filtered by
      --padj and --logfc thresholds enriched
      in the group of cells defined by the
      --second and --groupby parameters.

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_atac_dbinding.R"]

stdout: sc_atac_dbinding_stdout.log
stderr: sc_atac_dbinding_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell ATAC-Seq Differential Binding Analysis"
s:name: "Single-cell ATAC-Seq Differential Binding Analysis"
s:alternateName: "Identifies differential bound sites between two groups of cells"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-atac-dbinding.cwl
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
  Single-cell ATAC-Seq Differential Binding Analysis

  Identifies differential bound sites between two groups of cells
  --tmpdir parameter is not exposed as input.


s:about: |
  usage: sc_atac_dbinding.R [-h] --query QUERY --fragments FRAGMENTS
                            [--metadata METADATA] [--barcodes BARCODES]
                            [--groupby GROUPBY] [--subset [SUBSET [SUBSET ...]]]
                            --splitby SPLITBY --first FIRST --second SECOND
                            [--test {negative-binomial,poisson,logistic-regression,mast,manorm2}]
                            [--genome {hs,mm}] [--qvalue QVALUE]
                            [--minpeakgap MINPEAKGAP] [--binsize BINSIZE]
                            [--maxpeaks MAXPEAKS] [--blacklist BLACKLIST]
                            [--padj PADJ] [--logfc LOGFC] [--pdf] [--verbose]
                            [--tmpdir TMPDIR] [--output OUTPUT]
                            [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                            [--cpus CPUS] [--memory MEMORY]

  Single-cell ATAC-Seq Differential Binding Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include chromatin accessibility
                          information stored in the ATAC assay. Additionally
                          'rnaumap', and/or 'atacumap', and/or 'wnnumap'
                          dimensionality reductions should be present.
    --fragments FRAGMENTS
                          Count and barcode information for every ATAC fragment
                          used in the loaded Seurat object. File should be saved
                          in TSV format with tbi-index file.
    --metadata METADATA   Path to the TSV/CSV file to optionally extend Seurat
                          object metadata with categorical values using samples
                          identities. First column - 'library_id' should
                          correspond to all unique values from the 'new.ident'
                          column of the loaded Seurat object. If any of the
                          provided in this file columns are already present in
                          the Seurat object metadata, they will be overwritten.
                          When combined with --barcodes parameter, first the
                          metadata will be extended, then barcode filtering will
                          be applied. Default: no extra metadata is added
    --barcodes BARCODES   Path to the TSV/CSV file to optionally prefilter and
                          extend Seurat object metadata by selected barcodes.
                          First column should be named as 'barcode'. If file
                          includes any other columns they will be added to the
                          Seurat object metadata ovewriting the existing ones if
                          those are present. Default: all cells used, no extra
                          metadata is added
    --groupby GROUPBY     Column from the Seurat object metadata to group cells
                          for optional subsetting when combined with --subset
                          parameter. May be one of the extra metadata columns
                          added with --metadata or --barcodes parameters.
                          Ignored if --subset is not set. Default: do not
                          subset, include all cells into analysis.
    --subset [SUBSET [SUBSET ...]]
                          Values from the column set with --groupby parameter to
                          subset cells before running differential binding
                          analysis. Ignored if --groupby is not provided.
                          Default: do not subset cells, include all of them.
    --splitby SPLITBY     Column from the Seurat object metadata to split cells
                          into two groups to run --second vs --first
                          differential binding analysis. May be one of the extra
                          metadata columns added with --metadata or --barcodes
                          parameters.
    --first FIRST         Value from the Seurat object metadata column set with
                          --splitby parameter to define the first group of cells
                          for differential binding analysis.
    --second SECOND       Value from the Seurat object metadata column set with
                          --splitby parameter to define the second group of
                          cells for differential binding analysis.
    --test {negative-binomial,poisson,logistic-regression,mast,manorm2}
                          Test type to use in differential binding analysis. For
                          all tests except manorm2, peaks present in the loaded
                          Seurat object will be used. If manorm2 test selected,
                          peaks will be called per group defined by --splitby
                          parameter. Default: logistic-regression
    --genome {hs,mm}      Genome type of the sequencing data loaded from the
                          Seurat object. It will be used for effective genome
                          size selection when calling peaks with MACS2. Ignored
                          if --test is not set to manorm2. Default: hs (2.7e9)
    --qvalue QVALUE       Minimum FDR (q-value) cutoff for MACS2 peak detection.
                          Ignored if --test is not set to manorm2. Default: 0.05
    --minpeakgap MINPEAKGAP
                          If a distance between peaks is smaller than the
                          provided value they will be merged before splitting
                          them into reference genomic bins of size --binsize.
                          Ignored if --test is not set to manorm2. Default: 150
    --binsize BINSIZE     The size of non-overlapping reference genomic bins
                          used by MAnorm2 when generating a table of reads
                          counts per peaks. Ignored if --test is not set to
                          manorm2. Default: 1000
    --maxpeaks MAXPEAKS   The maximum number of the most significant (based on
                          qvalue) peaks to keep from each group of cells when
                          constructing reference genomic bins. Ignored if --test
                          is not set to manorm2. Default: keep all peaks
    --blacklist BLACKLIST
                          Path to the optional BED file with the genomic
                          blacklist regions to be filtered out before running
                          differential binding analysis. Any reference genomic
                          bin overlapping a blacklist region will be removed
                          from the output. Ignored if --test is not set to
                          manorm2.
    --padj PADJ           In the exploratory visualization part of the analysis
                          output only differentially bound peaks with adjusted
                          P-value not bigger than this value. Default: 0.05
    --logfc LOGFC         In the exploratory visualization part of the analysis
                          output only differentially bound peaks with log2 Fold
                          Change not smaller than this value. Default: 1.0
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --tmpdir TMPDIR       Directory to keep temporary files. Default: either
                          /tmp or defined by environment variables TMPDIR, TMP,
                          TEMP.
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32