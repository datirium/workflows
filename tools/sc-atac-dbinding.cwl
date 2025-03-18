cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.42


inputs:

  query_data_rds:
    type: File
    inputBinding:
      prefix: "--query"
    doc: |
      Path to the RDS file to load Seurat object
      from. This file should include chromatin
      accessibility information stored in the
      ATAC assay. The dimensionality reductions
      selected in the --reduction parameter should
      be present in the loaded Seurat object.

  reduction:
    type: string
    inputBinding:
      prefix: "--reduction"
    doc: |
      Dimensionality reduction to be
      used for generating UMAP plots.

  atac_fragments_file:
    type: File
    secondaryFiles:
    - .tbi
    inputBinding:
      prefix: "--fragments"
    doc: |
      Count and barcode information for every
      ATAC fragment used in the loaded Seurat
      object. File should be saved in TSV
      format with tbi-index file.

  datasets_metadata:
    type: File?
    inputBinding:
      prefix: "--metadata"
    doc: |
      Path to the TSV/CSV file to optionally extend
      Seurat object metadata with categorical values
      using samples identities. First column -
      library_id should correspond to all unique
      values from the new.ident column of the loaded
      Seurat object. If any of the provided in this
      file columns are already present in the Seurat
      object metadata, they will be overwritten. When
      combined with --barcodes parameter, first the
      metadata will be extended, then barcode filtering
      will be applied.
      Default: no extra metadata is added.

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the TSV/CSV file to optionally prefilter
      and extend Seurat object metadata by selected
      barcodes. First column should be named as barcode.
      If file includes any other columns they will be
      added to the Seurat object metadata ovewriting
      the existing ones if those are present.
      Default: all cells used, no extra metadata is added.

  groupby:
    type: string?
    inputBinding:
      prefix: "--groupby"
    doc: |
      Column from the Seurat object metadata to
      group cells for optional subsetting when
      combined with --subset parameter. May be
      one of the extra metadata columns added
      with --metadata or --barcodes parameters.
      Ignored if --subset is not set.
      Default: do not subset, include all
      cells into analysis.

  subset:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--subset"
    doc: |
      Values from the column set with --groupby
      parameter to subset cells before running
      differential accessibility analysis.
      Ignored if --groupby is not provided.
      Default: do not subset cells, include
      all of them.

  splitby:
    type: string
    inputBinding:
      prefix: "--splitby"
    doc: |
      Column from the Seurat object metadata to split cells
      into two groups to run --second vs --first
      differential accessibility analysis. If --test
      parameter is set to manorm2-full or manorm2-half, the
      --splitby shouldn't put cells from the same dataset
      into the different comparison groups. May be one of
      the extra metadata columns added with --metadata or
      --barcodes parameters.

  first_cond:
    type: string
    inputBinding:
      prefix: "--first"
    doc: |
      Value from the Seurat object metadata column set with
      --splitby parameter to define the first group of cells
      for differential accessibility analysis.

  second_cond:
    type: string
    inputBinding:
      prefix: "--second"
    doc: |
      Value from the Seurat object metadata column set with
      --splitby parameter to define the second group of
      cells for differential accessibility analysis.

  analysis_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "negative-binomial"          # (negbinom) Negative Binomial Generalized Linear Model (use FindMarkers with peaks from Seurat object)
      - "poisson"                    # (poisson) Poisson Generalized Linear Model (use FindMarkers with peaks from Seurat object)
      - "logistic-regression"        # (LR) Logistic Regression (use FindMarkers with peaks from Seurat object)
      - "mast"                       # (MAST) MAST package (use FindMarkers with peaks from Seurat object)
      - "manorm2-full"               # call peaks for each dataset with MACS2, then run MAnorm2 with datasets
      - "manorm2-half"               # call peaks for each comparison group with MACS2, then run MAnorm2 with datasets
    inputBinding:
      prefix: "--test"
    doc: |
      Test type to use in the differential accessibility
      analysis. For all tests except manorm2-full and
      manorm2-half, peaks already present in the loaded
      Seurat object will be used. If manorm2-full or
      manorm2-half test is selected, reads will be
      aggregated to pseudo bulk form either by dataset or
      comparison group and then peaks will be called with
      MACS2 per dataset. Default: logistic-regression

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
      if --test is not set to either manorm2-full or
      manorm2-half. Default: hs (2.7e9)

  minimum_qvalue:
    type: float?
    inputBinding:
      prefix: "--qvalue"
    doc: |
      Minimum FDR (q-value) cutoff for MACS2 peak detection.
      Ignored if --test is not set to either manorm2-full or
      manorm2-half. Default: 0.05

  minimum_peak_gap:
    type: int?
    inputBinding:
      prefix: "--minpeakgap"
    doc: |
      If a distance between peaks is smaller than the
      provided value they will be merged before splitting
      them into reference genomic bins of size --binsize.
      Ignored if --test is not set to either manorm2-full or
      manorm2-half. Default: 150

  bin_size:
    type: int?
    inputBinding:
      prefix: "--binsize"
    doc: |
      The size of non-overlapping reference genomic bins
      used by MAnorm2 when generating a table of reads
      counts per peaks. Ignored if --test is not set to
      either manorm2-full or manorm2-half. Default: 1000

  minimum_overlap:
    type: float?
    inputBinding:
      prefix: "--minoverlap"
    doc: |
      Keep only those reference genomic bins that are
      present in at least this fraction of datasets within
      each of the comparison groups. Used only when --test
      is set to manorm2-full. For manorm2-half this
      parameter will be automatically set to 1. Default: 0.5

  maximum_peaks:
    type: int?
    inputBinding:
      prefix: "--maxpeaks"
    doc: |
      The maximum number of the most significant (based on
      qvalue) peaks to keep from each group of cells when
      constructing reference genomic bins. Ignored if --test
      is not set to either manorm2-full or manorm2-half.
      Default: keep all peaks

  blacklist_regions_file:
    type:
    - "null"
    - File
    - type: enum
      symbols:
      - "hg19"
      - "hg38"
      - "mm10"
    inputBinding:
      prefix: "--blacklist"
      valueFrom: |
        ${
          if (self.class && self.class == "File"){
            return self;
          } else if (self == "hg19") {
            return "/opt/sc_tools/hg19-blacklist.v2.bed";
          } else if (self == "hg38") {
            return "/opt/sc_tools/hg38-blacklist.v2.bed";
          } else if (self == "mm10") {
            return "/opt/sc_tools/mm10-blacklist.v2.bed";
          } else {
            return null;
          }
        }
    doc: |
      Path to the optional BED file with the genomic
      blacklist regions to be filtered out before running
      differential accessibility analysis. Any reference
      genomic bin overlapping a blacklist region will be
      removed from the output. Ignored if --test is not set
      to either manorm2-full or manorm2-half.

  maximum_padj:
    type: float?
    inputBinding:
      prefix: "--padj"
    doc: |
      In the exploratory visualization part of
      the analysis output only differentially
      accessible regions with adjusted P-value
      not bigger than this value. Default: 0.05

  minimum_logfc:
    type: float?
    inputBinding:
      prefix: "--logfc"
    doc: |
      In the exploratory visualization part of
      the analysis output only differentially
      accessible regions with log2 Fold Change
      not smaller than this value. Default: 1.0

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
      Color theme for all generated plots.
      One of gray, bw, linedraw, light, dark,
      minimal, classic, void. Default: classic

  verbose:
    type: boolean?
    inputBinding:
      prefix: "--verbose"
    doc: |
      Print debug information.
      Default: false

  export_html_report:
    type: boolean?
    default: false
    doc: |
      Export tehcnical report. HTML format.
      Note, stdout will be less informative.
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
      Maximum memory in GB allowed to
      be shared between the workers
      when using multiple --cpus.
      Default: 32

  vector_memory_limit:
    type: int?
    default: 128
    doc: |
      Maximum vector memory in GB
      allowed to be used by R.
      Default: 128

  threads:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of cores/cpus to use.
      Default: 1

  seed:
    type: int?
    inputBinding:
      prefix: "--seed"
    doc: |
      Seed number for random values.
      Default: 42


outputs:

  umap_spl_tst_plot_png:
    type: File?
    outputBinding:
      glob: "*_umap_spl_tst.png"
    doc: |
      UMAP colored by selected for analysis cells.
      Split by tested condition. Optionally subsetted
      to the specific group.
      PNG format.

  cell_cnts_plot_png:
    type: File?
    outputBinding:
      glob: "*_cell_cnts.png"
    doc: |
      Number of cells per dataset or tested
      condition. Colored by tested condition.
      Optionally subsetted to the specific group.
      PNG format.

  vlcn_plot_png:
    type: File?
    outputBinding:
      glob: "*_vlcn.png"
    doc: |
      Volcano plot of differentially
      accessible regions.
      PNG format.

  tag_dnst_htmp_plot_png:
    type: File?
    outputBinding:
      glob: "*_tag_dnst_htmp.png"
    doc: |
      Tag density around the centers of
      differentially accessible regions,
      sorted in descending order by the
      mean value of each region. Optionally
      subsetted to the specific group.
      PNG format.

  tag_dnst_htmp_gct:
    type: File?
    outputBinding:
      glob: "*_tag_dnst_htmp.gct"
    doc: |
      Tag density around the centers of
      differentially accessible regions,
      sorted in descending order by the
      mean value of each region. Optionally
      subsetted to the specific group.
      GCT format.

  tag_dnst_htmp_html:
    type: File?
    outputBinding:
      glob: "*_tag_dnst_htmp.html"
    doc: |
      Tag density around the centers of
      differentially accessible regions,
      sorted in descending order by the
      mean value of each region. Optionally
      subsetted to the specific group.
      HTML format.

  tag_dnst_htmp_tsv:
    type: File?
    outputBinding:
      glob: "*_tag_dnst_htmp.tsv"
    doc: |
      Tag density around the centers of
      differentially accessible regions,
      sorted in descending order by the
      mean value of each region. Optionally
      subsetted to the specific group.
      TSV format.

  diff_bound_sites:
    type: File
    outputBinding:
      glob: "*_db_sites.tsv"
    doc: |
      Not filtered differentially
      accessible regions.
      TSV format.

  first_enrch_bed_file:
    type: File?
    outputBinding:
      glob: "*_first_enrch.bed"
    doc: |
      Differentially accessible regions
      enriched in the cells from the first
      comparison group. Filtered by adjusted
      p-value and log2FoldChange thresholds.
      Optionally subsetted to the specific
      group. Sorted by log2FoldChange in the
      ascendant order to correspond to the
      tag density heatmap.
      BED format.

  second_enrch_bed_file:
    type: File?
    outputBinding:
      glob: "*_second_enrch.bed"
    doc: |
      Differentially accessible regions
      enriched in the cells from the second
      comparison group. Filtered by adjusted
      p-value and log2FoldChange thresholds.
      Optionally subsetted to the specific
      group. Sorted by log2FoldChange in the
      descendant order to correspond to the
      tag density heatmap.
      BED format.

  fragments_bigwig_file:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*.bigWig"
    doc: |
      Normalized genome coverage calculated
      from the ATAC fragments split either
      by dataset or tested condition.
      BigWig format.

  peaks_bed_file:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_peaks.narrowPeak"
    doc: |
      Peaks called by MACS2 from the Tn5
      cut sites split either by dataset
      or tested condition.
      NarrowPeak format.

  peaks_xls_file:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_peaks.xls"
    doc: |
      Peaks called by MACS2 from the Tn5
      cut sites split either by dataset
      or tested condition.
      XLS format.

  summits_bed_file:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_summits.bed"
    doc: |
      Summits of the peaks called by MACS2
      from the Tn5 cut sites split either
      by dataset or tested condition.
      BED format.

  dflt_peaks_bigbed_file:
    type: File?
    outputBinding:
      glob: "*_dflt_peaks.bigBed"
    doc: |
      Peaks extracted from the
      loaded Seurat object.
      BigBed format.

  all_plots_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*.pdf"
    doc: |
      All generated plots.
      PDF format.

  sc_report_html_file:
    type: File?
    outputBinding:
      glob: "sc_report.html"
    doc: |
      Tehcnical report.
      HTML format.

  human_log:
    type: File?
    outputBinding:
      glob: "error_report.txt"
    doc: |
      Human readable error log.
      TXT format.

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["Rscript"]
arguments:
- valueFrom: $(inputs.export_html_report?["/usr/local/bin/sc_report_wrapper.R", "/usr/local/bin/sc_atac_dbinding.R"]:"/usr/local/bin/sc_atac_dbinding.R")


stdout: sc_atac_dbinding_stdout.log
stderr: sc_atac_dbinding_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-Cell ATAC-Seq Differential Accessibility Analysis"
s:name: "Single-Cell ATAC-Seq Differential Accessibility Analysis"
s:alternateName: "Identifies differentially accessible regions between two groups of cells"

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
  Single-Cell ATAC-Seq Differential Accessibility Analysis

  Identifies differentially accessible regions between two
  groups of cells --tmpdir parameter is not exposed as input.


s:about: |
  usage: /usr/local/bin/sc_atac_dbinding.R [-h] --query QUERY --reduction
                                          REDUCTION --fragments FRAGMENTS
                                          [--metadata METADATA]
                                          [--barcodes BARCODES]
                                          [--groupby GROUPBY]
                                          [--subset [SUBSET [SUBSET ...]]]
                                          --splitby SPLITBY --first FIRST
                                          --second SECOND
                                          [--test {negative-binomial,poisson,logistic-regression,mast,manorm2-full,manorm2-half}]
                                          [--genome {hs,mm}] [--qvalue QVALUE]
                                          [--minpeakgap MINPEAKGAP]
                                          [--binsize BINSIZE]
                                          [--minoverlap MINOVERLAP]
                                          [--maxpeaks MAXPEAKS]
                                          [--blacklist BLACKLIST] [--padj PADJ]
                                          [--logfc LOGFC] [--pdf] [--verbose]
                                          [--tmpdir TMPDIR] [--output OUTPUT]
                                          [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                                          [--cpus CPUS] [--memory MEMORY]
                                          [--seed SEED]

  Single-Cell ATAC-Seq Differential Accessibility Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --query QUERY         Path to the RDS file to load Seurat object from. This
                          file should include chromatin accessibility
                          information stored in the ATAC assay. The
                          dimensionality reductions selected in the --reduction
                          parameter should be present in the loaded Seurat
                          object.
    --reduction REDUCTION
                          Dimensionality reduction to be used for generating
                          UMAP plots.
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
                          subset cells before running differential accessibility
                          analysis. Ignored if --groupby is not provided.
                          Default: do not subset cells, include all of them.
    --splitby SPLITBY     Column from the Seurat object metadata to split cells
                          into two groups to run --second vs --first
                          differential accessibility analysis. If --test
                          parameter is set to manorm2-full or manorm2-half, the
                          --splitby shouldn't put cells from the same dataset
                          into the different comparison groups. May be one of
                          the extra metadata columns added with --metadata or
                          --barcodes parameters.
    --first FIRST         Value from the Seurat object metadata column set with
                          --splitby parameter to define the first group of cells
                          for differential accessibility analysis.
    --second SECOND       Value from the Seurat object metadata column set with
                          --splitby parameter to define the second group of
                          cells for differential accessibility analysis.
    --test {negative-binomial,poisson,logistic-regression,mast,manorm2-full,manorm2-half}
                          Test type to use in the differential accessibility
                          analysis. For all tests except manorm2-full and
                          manorm2-half, peaks already present in the loaded
                          Seurat object will be used. If manorm2-full or
                          manorm2-half test is selected, reads will be
                          aggregated to pseudo bulk form either by dataset or
                          comparison group and then peaks will be called with
                          MACS2 per dataset. Default: logistic-regression
    --genome {hs,mm}      Genome type of the sequencing data loaded from the
                          Seurat object. It will be used for effective genome
                          size selection when calling peaks with MACS2. Ignored
                          if --test is not set to either manorm2-full or
                          manorm2-half. Default: hs (2.7e9)
    --qvalue QVALUE       Minimum FDR (q-value) cutoff for MACS2 peak detection.
                          Ignored if --test is not set to either manorm2-full or
                          manorm2-half. Default: 0.05
    --minpeakgap MINPEAKGAP
                          If a distance between peaks is smaller than the
                          provided value they will be merged before splitting
                          them into reference genomic bins of size --binsize.
                          Ignored if --test is not set to either manorm2-full or
                          manorm2-half. Default: 150
    --binsize BINSIZE     The size of non-overlapping reference genomic bins
                          used by MAnorm2 when generating a table of reads
                          counts per peaks. Ignored if --test is not set to
                          either manorm2-full or manorm2-half. Default: 1000
    --minoverlap MINOVERLAP
                          Keep only those reference genomic bins that are
                          present in at least this fraction of datasets within
                          each of the comparison groups. Used only when --test
                          is set to manorm2-full. For manorm2-half this
                          parameter will be automatically set to 1. Default: 0.5
    --maxpeaks MAXPEAKS   The maximum number of the most significant (based on
                          qvalue) peaks to keep from each group of cells when
                          constructing reference genomic bins. Ignored if --test
                          is not set to either manorm2-full or manorm2-half.
                          Default: keep all peaks
    --blacklist BLACKLIST
                          Path to the optional BED file with the genomic
                          blacklist regions to be filtered out before running
                          differential accessibility analysis. Any reference
                          genomic bin overlapping a blacklist region will be
                          removed from the output. Ignored if --test is not set
                          to either manorm2-full or manorm2-half.
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
    --seed SEED           Seed number for random values. Default: 42