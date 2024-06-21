cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.38


inputs:

  feature_bc_matrices_folder:
    type: Directory
    inputBinding:
      prefix: "--mex"
    doc: |
      Path to the folder with feature-barcode matrix from Cell Ranger Count/Aggregate
      experiment in MEX format.

  aggregation_metadata:
    type: File?
    inputBinding:
      prefix: "--identity"
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities, if '--mex' points to
      the Cell Ranger Aggregate outputs. The aggregation.csv file can be used. In case of
      using feature-barcode matrices the file with identities should include at least one
      column - 'library_id', and a row with aliases per each experiment from the '--mex'
      input. The order of rows should correspond to the order of feature-barcode matrices
      provided in the '--mex' parameter.

  grouping_data:
    type: File?
    inputBinding:
      prefix: "--grouping"
    doc: |
      Path to the TSV/CSV file to define datasets grouping.
      First column - 'library_id' with the values and order
      that correspond to the 'library_id' column from the '
      --identity' file, second column 'condition'.
      Default: each dataset is assigned to its own group.

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the TSV/CSV file to optionally prefilter and
      extend Seurat object metadata be selected barcodes.
      First column should be named as 'barcode'. If file
      includes any other columns they will be added to the
      Seurat object metadata ovewriting the existing ones if
      those are present.
      Default: all cells used, no extra metadata is added

  rna_minimum_cells:
    type: int?
    inputBinding:
      prefix: "--rnamincells"
    doc: |
      Include only genes detected in at least this many cells.
      Default: 5 (applied to all datasets)

  minimum_genes:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--mingenes"
    doc: |
      Include cells where at least this many genes are detected. If multiple values
      provided, each of them will be applied to the correspondent dataset from the
      '--mex' input based on the '--identity' file. Any 0 will be replaced with the
      auto-estimated threshold (median - 2.5 * MAD) calculated per dataset.
      Default: 250 (applied to all datasets)

  maximum_genes:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--maxgenes"
    doc: |
      Include cells with the number of genes not bigger than this value. If multiple
      values provided, each of them will be applied to the correspondent dataset from
      the '--mex' input based on the '--identity' file. Any 0 will be replaced with the
      auto-estimated threshold (median + 5 * MAD) calculated per dataset.
      Default: 5000 (applied to all datasets)

  minimum_umis:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--minumis"
    doc: |
      Include cells where at least this many RNA reads are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the '--mex' input based on the '--identity' file. Any 0 will be
      replaced with the auto-estimated threshold (median - 2.5 * MAD) calculated
      per dataset.
      Default: 500 (applied to all datasets)

  minimum_novelty_score:
    type:
    - "null"
    - float
    - float[]
    inputBinding:
      prefix: "--minnovelty"
    doc: |
      Include cells with the novelty score not lower than this value, calculated for
      as log10(genes)/log10(UMI). If multiple values provided, each of them will
      be applied to the correspondent dataset from the '--mex' input based on the
      '--identity' file.
      Default: 0.8 (applied to all datasets)

  mito_pattern:
    type: string?
    inputBinding:
      prefix: "--mitopattern"
    doc: |
      Regex pattern to identify mitochondrial genes.
      Default: '^mt-|^MT-'

  maximum_mito_perc:
    type: float?
    inputBinding:
      prefix: "--maxmt"
    doc: |
      Include cells with the percentage of RNA reads mapped to mitochondrial
      genes not bigger than this value. Set to 0 for using an auto-estimated
      threshold equal to the maximum among (median + 2 * MAD) values
      calculated per dataset.
      Default: 5 (applied to all datasets)

  remove_doublets:
    type: boolean?
    inputBinding:
      prefix: "--removedoublets"
    doc: |
      Remove cells that were identified as doublets. Cells with
      RNA UMI < 200 will not be evaluated. Default: do not remove
      doublets

  rna_doublet_rate:
    type: float?
    inputBinding:
      prefix: "--rnadbr"
    doc: |
      Expected RNA doublet rate. Default: 1 percent per
      thousand cells captured with 10x genomics

  rna_doublet_rate_sd:
    type: float?
    inputBinding:
      prefix: "--rnadbrsd"
    doc: |
      Uncertainty range in the RNA doublet rate, interpreted as
      a +/- around the value provided in --rnadbr. Set to 0 to
      disable. Set to 1 to make the threshold depend entirely
      on the misclassification rate. Default: 40 percents of the
      value provided in --rnadbr

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
      Color theme for all generated plots. One of gray, bw, linedraw, light,
      dark, minimal, classic, void.
      Default: classic

  verbose:
    type: boolean?
    inputBinding:
      prefix: "--verbose"
    doc: |
      Print debug information.
      Default: false

  export_h5seurat_data:
    type: boolean?
    inputBinding:
      prefix: "--h5seurat"
    doc: |
      Save Seurat data to h5seurat file.
      Default: false

  export_h5ad_data:
    type: boolean?
    inputBinding:
      prefix: "--h5ad"
    doc: |
      Save raw counts from the RNA assay to h5ad file.
      Default: false

  export_loupe_data:
    type: boolean?
    inputBinding:
      prefix: "--loupe"
    doc: |
      Save raw counts from the RNA assay to Loupe file. By
      enabling this feature you accept the End-User License
      Agreement available at https://10xgen.com/EULA.
      Default: false

  export_ucsc_cb:
    type: boolean?
    inputBinding:
      prefix: "--cbbuild"
    doc: |
      Export results to UCSC Cell Browser. Default: false

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

  seed:
    type: int?
    inputBinding:
      prefix: "--seed"
    doc: |
      Seed number for random values.
      Default: 42


outputs:

  raw_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_1_2_qc_mtrcs_pca.png"
    doc: |
      QC metrics PCA.
      Unfiltered; PC1/PC2.
      PNG format.

  raw_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_2_3_qc_mtrcs_pca.png"
    doc: |
      QC metrics PCA.
      Unfiltered; PC2/PC3.
      PNG format.

  raw_cell_cnts_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_cell_cnts.png"
    doc: |
      Number of cells per dataset.
      Unfiltered.
      PNG format.

  raw_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_umi_dnst.png"
    doc: |
      Distribution of RNA reads per cell.
      Unfiltered.
      PNG format.

  raw_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst.png"
    doc: |
      Distribution of genes per cell.
      Unfiltered.
      PNG format.

  raw_gene_umi_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi.png"
    doc: |
      Genes vs RNA reads per cell.
      Unfiltered.
      PNG format.

  raw_umi_mito_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_umi_mito.png"
    doc: |
      RNA reads vs mitochondrial percentage
      per cell.
      Unfiltered.
      PNG format.

  raw_mito_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_dnst.png"
    doc: |
      Distribution of RNA reads mapped
      to mitochondrial genes per cell.
      Unfiltered.
      PNG format.

  raw_nvlt_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_dnst.png"
    doc: |
      Distribution of novelty score per cell.
      Unfiltered.
      PNG format.

  raw_qc_mtrcs_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs_dnst.png"
    doc: |
      Distribution of QC metrics per cell.
      Unfiltered.
      PNG format.

  raw_rnadbl_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_rnadbl.png"
    doc: |
      Percentage of RNA doublets.
      Unfiltered.
      PNG format.

  raw_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_umi_dnst_spl_cnd.png"
    doc: |
      Distribution of RNA reads per cell.
      Unfiltered; split by grouping condition.
      PNG format.

  raw_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst_spl_cnd.png"
    doc: |
      Distribution of genes per cell.
      Unfiltered; split by grouping condition.
      PNG format.

  raw_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_dnst_spl_cnd.png"
    doc: |
      Distribution of RNA reads mapped
      to mitochondrial genes per cell.
      Unfiltered; split by grouping condition.
      PNG format.

  raw_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_dnst_spl_cnd.png"
    doc: |
      Distribution of novelty score per cell.
      Unfiltered; split by grouping condition.
      PNG format.

  fltr_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_1_2_qc_mtrcs_pca.png"
    doc: |
      QC metrics PCA.
      Filtered; PC1/PC2.
      PNG format.

  fltr_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_2_3_qc_mtrcs_pca.png"
    doc: |
      QC metrics PCA.
      Filtered; PC2/PC3.
      PNG format.

  fltr_cell_cnts_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_cell_cnts.png"
    doc: |
      Number of cells per dataset.
      Filtered.
      PNG format.

  fltr_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_umi_dnst.png"
    doc: |
      Distribution of RNA reads per cell.
      Filtered.
      PNG format.

  fltr_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_dnst.png"
    doc: |
      Distribution of genes per cell.
      Filtered.
      PNG format.

  fltr_gene_umi_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_umi.png"
    doc: |
      Genes vs RNA reads per cell.
      Filtered.
      PNG format.

  fltr_umi_mito_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_umi_mito.png"
    doc: |
      RNA reads vs mitochondrial percentage
      per cell.
      Filtered.
      PNG format.

  fltr_mito_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_dnst.png"
    doc: |
      Distribution of RNA reads mapped
      to mitochondrial genes per cell.
      Filtered.
      PNG format.

  fltr_nvlt_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_nvlt_dnst.png"
    doc: |
      Distribution of novelty score per cell.
      Filtered.
      PNG format.

  fltr_qc_mtrcs_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_qc_mtrcs_dnst.png"
    doc: |
      Distribution of QC metrics per cell.
      Filtered.
      PNG format.

  fltr_rnadbl_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_rnadbl.png"
    doc: |
      Percentage of RNA doublets.
      Filtered.
      PNG format.

  fltr_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_umi_dnst_spl_cnd.png"
    doc: |
      Distribution of RNA reads per cell.
      Filtered; split by grouping condition.
      PNG format.

  fltr_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_dnst_spl_cnd.png"
    doc: |
      Distribution of genes per cell.
      Filtered; split by grouping condition.
      PNG format.

  fltr_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_dnst_spl_cnd.png"
    doc: |
      Distribution of RNA reads mapped
      to mitochondrial genes per cell.
      Filtered; split by grouping condition.
      PNG format.

  fltr_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_nvlt_dnst_spl_cnd.png"
    doc: |
      Distribution of novelty score per cell.
      Filtered; split by grouping condition.
      PNG format.

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

  ucsc_cb_config_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser"
    doc: |
      UCSC Cell Browser configuration data.

  ucsc_cb_html_data:
    type: Directory?
    outputBinding:
      glob: "*_cellbrowser/html_data"
    doc: |
      UCSC Cell Browser html data.

  ucsc_cb_html_file:
    type: File?
    outputBinding:
      glob: "*_cellbrowser/html_data/index.html"
    doc: |
      UCSC Cell Browser html index.

  seurat_data_rds:
    type: File
    outputBinding:
      glob: "*_data.rds"
    doc: |
      Seurat object.
      RDS format

  datasets_metadata:
    type: File
    outputBinding:
      glob: "*_meta.tsv"
    doc: |
      Example of datasets metadata file
      in TSV format

  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_data.h5seurat"
    doc: |
      Seurat object.
      h5Seurat format

  seurat_data_h5ad:
    type: File?
    outputBinding:
      glob: "*_counts.h5ad"
    doc: |
      Seurat object.
      H5AD format

  seurat_data_cloupe:
    type: File?
    outputBinding:
      glob: "*_counts.cloupe"
    doc: |
      Seurat object.
      Loupe format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_rna_filter.R"]

stdout: sc_rna_filter_stdout.log
stderr: sc_rna_filter_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-Cell RNA-Seq Filtering Analysis"
s:name: "Single-Cell RNA-Seq Filtering Analysis"
s:alternateName: "Single-Cell RNA-Seq Filtering Analysis"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-rna-filter.cwl
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
  Single-Cell RNA-Seq Filtering Analysis

  Removes low-quality cells from the outputs of the “Cell
  Ranger Count (RNA)”, “Cell Ranger Count (RNA+VDJ)”, and
  “Cell Ranger Aggregate (RNA, RNA+VDJ)” pipelines. The
  results of this workflow are used in the “Single-Cell
  RNA-Seq Dimensionality Reduction Analysis” pipeline.


s:about: |
  usage: sc_rna_filter.R [-h] --mex MEX [MEX ...]
                                        [--identity IDENTITY]
                                        [--grouping GROUPING]
                                        [--barcodes BARCODES]
                                        [--rnamincells RNAMINCELLS]
                                        [--mingenes [MINGENES [MINGENES ...]]]
                                        [--maxgenes [MAXGENES [MAXGENES ...]]]
                                        [--minumis [MINUMIS [MINUMIS ...]]]
                                        [--minnovelty [MINNOVELTY [MINNOVELTY ...]]]
                                        [--mitopattern MITOPATTERN]
                                        [--maxmt MAXMT] [--removedoublets]
                                        [--rnadbr RNADBR] [--rnadbrsd RNADBRSD]
                                        [--pdf] [--verbose] [--h5seurat]
                                        [--h5ad] [--loupe] [--cbbuild]
                                        [--output OUTPUT]
                                        [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                                        [--cpus CPUS] [--memory MEMORY]
                                        [--seed SEED]

  Single-Cell RNA-Seq Filtering Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --mex MEX [MEX ...]   Path to the folder with feature-barcode matrix from
                          Cell Ranger Count/Aggregate experiment in MEX format.
                          If multiple locations provided data is assumed to be
                          not aggregated (outputs from the multiple Cell Ranger
                          Count experiments) and will be merged before the
                          analysis.
    --identity IDENTITY   Path to the metadata TSV/CSV file to set the datasets
                          identities, if '--mex' points to the Cell Ranger
                          Aggregate outputs. The aggregation.csv file can be
                          used. In case of using feature-barcode matrices the
                          file with identities should include at least one
                          column - 'library_id', and a row with aliases per each
                          experiment from the '--mex' input. The order of rows
                          should correspond to the order of feature-barcode
                          matrices provided in the '--mex' parameter.
    --grouping GROUPING   Path to the TSV/CSV file to define datasets grouping.
                          First column - 'library_id' with the values and order
                          that correspond to the 'library_id' column from the '
                          --identity' file, second column 'condition'. Default:
                          each dataset is assigned to its own group.
    --barcodes BARCODES   Path to the TSV/CSV file to optionally prefilter and
                          extend Seurat object metadata be selected barcodes.
                          First column should be named as 'barcode'. If file
                          includes any other columns they will be added to the
                          Seurat object metadata ovewriting the existing ones if
                          those are present. Default: all cells used, no extra
                          metadata is added
    --rnamincells RNAMINCELLS
                          Include only genes detected in at least this many
                          cells. Ignored when '--mex' points to the feature-
                          barcode matrices from the multiple Cell Ranger Count
                          experiments. Default: 5 (applied to all datasets)
    --mingenes [MINGENES [MINGENES ...]]
                          Include cells where at least this many genes are
                          detected. If multiple values provided, each of them
                          will be applied to the correspondent dataset from the
                          '--mex' input based on the '--identity' file. Any 0
                          will be replaced with the auto-estimated threshold
                          (median - 2.5 * MAD) calculated per dataset. Default:
                          250 (applied to all datasets)
    --maxgenes [MAXGENES [MAXGENES ...]]
                          Include cells with the number of genes not bigger than
                          this value. If multiple values provided, each of them
                          will be applied to the correspondent dataset from the
                          '--mex' input based on the '--identity' file. Any 0
                          will be replaced with the auto-estimated threshold
                          (median + 5 * MAD) calculated per dataset. Default:
                          5000 (applied to all datasets)
    --minumis [MINUMIS [MINUMIS ...]]
                          Include cells where at least this many RNA reads are
                          detected. If multiple values provided, each of them
                          will be applied to the correspondent dataset from the
                          '--mex' input based on the '--identity' file. Any 0
                          will be replaced with the auto-estimated threshold
                          (median - 2.5 * MAD) calculated per dataset. Default:
                          500 (applied to all datasets)
    --minnovelty [MINNOVELTY [MINNOVELTY ...]]
                          Include cells with the novelty score not lower than
                          this value, calculated for as log10(genes)/log10(UMI).
                          If multiple values provided, each of them will be
                          applied to the correspondent dataset from the '--mex'
                          input based on the '--identity' file. Default: 0.8
                          (applied to all datasets)
    --mitopattern MITOPATTERN
                          Regex pattern to identify mitochondrial genes.
                          Default: '^mt-|^MT-'
    --maxmt MAXMT         Include cells with the percentage of RNA reads mapped
                          to mitochondrial genes not bigger than this value. Set
                          to 0 for using an auto-estimated threshold equal to
                          the maximum among (median + 2 * MAD) values calculated
                          per dataset. Default: 5 (applied to all datasets)
    --removedoublets      Remove cells that were identified as doublets. Cells
                          with RNA UMI < 200 will not be evaluated. Default: do
                          not remove doublets
    --rnadbr RNADBR       Expected RNA doublet rate. Default: 1 percent per
                          thousand cells captured with 10x genomics
    --rnadbrsd RNADBRSD   Uncertainty range in the RNA doublet rate, interpreted
                          as a +/- around the value provided in --rnadbr. Set to
                          0 to disable. Set to 1 to make the threshold depend
                          entirely on the misclassification rate. Default: 40
                          percents of the value provided in --rnadbr
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --h5ad                Save raw counts from the RNA assay to h5ad file.
                          Default: false
    --loupe               Save raw counts from the RNA assay to Loupe file. By
                          enabling this feature you accept the End-User License
                          Agreement available at https://10xgen.com/EULA.
                          Default: false
    --cbbuild             Export results to UCSC Cell Browser. Default: false
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple '--cpus'. Default: 32
    --seed SEED           Seed number for random values. Default: 42