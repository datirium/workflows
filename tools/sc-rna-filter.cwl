cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entryname: dummy_metadata.csv
    entry: |
      library_id
      Experiment
- class: EnvVarRequirement
  envDef:
    R_MAX_VSIZE: $((inputs.vector_memory_limit * 1000000000).toString())


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.12


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
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities. If '--mex' points to
      the Cell Ranger Aggregate outputs, the aggregation.csv file can be used. If input is not
      provided, the default dummy_metadata.csv will be used instead.

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
      '--mex' input based on the '--identity' file.
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
      the '--mex' input based on the '--identity' file.
      Default: 5000 (applied to all datasets)

  rna_minimum_umi:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--rnaminumi"
    doc: |
      Include cells where at least this many UMI (transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the '--mex' input based on the '--identity' file.
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
      Default: '^Mt-'

  maximum_mito_perc:
    type: float?
    inputBinding:
      prefix: "--maxmt"
    doc: |
      Include cells with the percentage of transcripts mapped to mitochondrial
      genes not bigger than this value.
      Default: 5 (applied to all datasets)

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
      Save Seurat data to h5ad file.
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

  raw_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_1_2_qc_mtrcs_pca.png"
    doc: |
      PC1 and PC2 from the QC metrics PCA (not filtered).
      PNG format

  raw_1_2_qc_mtrcs_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_1_2_qc_mtrcs_pca.pdf"
    doc: |
      PC1 and PC2 from the QC metrics PCA (not filtered).
      PDF format

  raw_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_2_3_qc_mtrcs_pca.png"
    doc: |
      PC2 and PC3 from the QC metrics PCA (not filtered).
      PNG format

  raw_2_3_qc_mtrcs_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_2_3_qc_mtrcs_pca.pdf"
    doc: |
      PC2 and PC3 from the QC metrics PCA (not filtered).
      PDF format

  raw_cells_count_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_cells_count.png"
    doc: |
      Number of cells per dataset (not filtered).
      PNG format

  raw_cells_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_cells_count.pdf"
    doc: |
      Number of cells per dataset (not filtered).
      PDF format

  raw_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_umi_dnst.png"
    doc: |
      UMI per cell density (not filtered).
      PNG format

  raw_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_umi_dnst.pdf"
    doc: |
      UMI per cell density (not filtered).
      PDF format

  raw_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst.png"
    doc: |
      Genes per cell density (not filtered).
      PNG format

  raw_gene_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst.pdf"
    doc: |
      Genes per cell density (not filtered).
      PDF format

  raw_gene_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_corr.png"
    doc: |
      Genes vs UMI per cell correlation (not filtered).
      PNG format

  raw_gene_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_corr.pdf"
    doc: |
      Genes vs UMI per cell correlation (not filtered).
      PDF format

  raw_mito_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_dnst.png"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (not filtered).
      PNG format

  raw_mito_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_mito_dnst.pdf"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (not filtered).
      PDF format

  raw_nvlt_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_dnst.png"
    doc: |
      Novelty score per cell density (not filtered).
      PNG format

  raw_nvlt_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_dnst.pdf"
    doc: |
      Novelty score per cell density (not filtered).
      PDF format

  raw_qc_mtrcs_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs_dnst.png"
    doc: |
      QC metrics per cell density (not filtered).
      PNG format

  raw_qc_mtrcs_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs_dnst.pdf"
    doc: |
      QC metrics per cell density (not filtered).
      PDF format

  raw_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_umi_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition UMI per cell density (not filtered).
      PNG format

  raw_umi_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_umi_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition UMI per cell density (not filtered).
      PDF format

  raw_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition genes per cell density (not filtered).
      PNG format

  raw_gene_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition genes per cell density (not filtered).
      PDF format

  raw_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (not filtered).
      PNG format

  raw_mito_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_mito_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (not filtered).
      PDF format

  raw_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the novelty score per cell density (not filtered).
      PNG format

  raw_nvlt_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the novelty score per cell density (not filtered).
      PDF format

  fltr_1_2_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_1_2_qc_mtrcs_pca.png"
    doc: |
      PC1 and PC2 from the QC metrics PCA (filtered).
      PNG format

  fltr_1_2_qc_mtrcs_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_1_2_qc_mtrcs_pca.pdf"
    doc: |
      PC1 and PC2 from the QC metrics PCA (filtered).
      PDF format

  fltr_2_3_qc_mtrcs_pca_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_2_3_qc_mtrcs_pca.png"
    doc: |
      PC2 and PC3 from the QC metrics PCA (filtered).
      PNG format

  fltr_2_3_qc_mtrcs_pca_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_2_3_qc_mtrcs_pca.pdf"
    doc: |
      PC2 and PC3 from the QC metrics PCA (filtered).
      PDF format

  fltr_cells_count_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_cells_count.png"
    doc: |
      Number of cells per dataset (filtered).
      PNG format

  fltr_cells_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_cells_count.pdf"
    doc: |
      Number of cells per dataset (filtered).
      PDF format

  fltr_umi_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_umi_dnst.png"
    doc: |
      UMI per cell density (filtered).
      PNG format

  fltr_umi_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_umi_dnst.pdf"
    doc: |
      UMI per cell density (filtered).
      PDF format

  fltr_gene_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_dnst.png"
    doc: |
      Genes per cell density (filtered).
      PNG format

  fltr_gene_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_dnst.pdf"
    doc: |
      Genes per cell density (filtered).
      PDF format

  fltr_gene_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_umi_corr.png"
    doc: |
      Genes vs UMI per cell correlation (filtered).
      PNG format

  fltr_gene_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_umi_corr.pdf"
    doc: |
      Genes vs UMI per cell correlation (filtered).
      PDF format

  fltr_mito_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_dnst.png"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (filtered).
      PNG format

  fltr_mito_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_dnst.pdf"
    doc: |
      Percentage of transcripts mapped to mitochondrial genes per cell density (filtered).
      PDF format

  fltr_nvlt_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_nvlt_dnst.png"
    doc: |
      Novelty score per cell density (filtered).
      PNG format

  fltr_nvlt_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_nvlt_dnst.pdf"
    doc: |
      Novelty score per cell density (filtered).
      PDF format

  fltr_qc_mtrcs_dnst_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_qc_mtrcs_dnst.png"
    doc: |
      QC metrics per cell density (filtered).
      PNG format

  fltr_qc_mtrcs_dnst_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_qc_mtrcs_dnst.pdf"
    doc: |
      QC metrics per cell density (filtered).
      PDF format

  fltr_umi_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_umi_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition UMI per cell density (filtered).
      PNG format

  fltr_umi_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_umi_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition UMI per cell density (filtered).
      PDF format

  fltr_gene_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition genes per cell density (filtered).
      PNG format

  fltr_gene_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition genes per cell density (filtered).
      PDF format

  fltr_mito_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (filtered).
      PNG format

  fltr_mito_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the percentage of transcripts mapped
      to mitochondrial genes per cell density (filtered).
      PDF format

  fltr_nvlt_dnst_spl_cnd_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_nvlt_dnst_spl_cnd.png"
    doc: |
      Split by grouping condition the novelty score per cell density (filtered).
      PNG format

  fltr_nvlt_dnst_spl_cnd_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_nvlt_dnst_spl_cnd.pdf"
    doc: |
      Split by grouping condition the novelty score per cell density (filtered).
      PDF format

  seurat_data_rds:
    type: File
    outputBinding:
      glob: "*_data.rds"
    doc: |
      Filtered Seurat data in RDS format

  seurat_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_data.h5seurat"
    doc: |
      Filtered Seurat data in h5seurat format

  seurat_data_h5ad:
    type: File?
    outputBinding:
      glob: "*_data.h5ad"
    doc: |
      Reduced Seurat data in h5ad format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_rna_filter.R"]
arguments:
- valueFrom: |
    ${
      if (inputs.aggregation_metadata) {
        return inputs.aggregation_metadata;
      } else {
        return runtime.outdir + "/dummy_metadata.csv"
      }
    }
  prefix: "--identity"


stdout: sc_rna_filter_stdout.log
stderr: sc_rna_filter_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell RNA-Seq Filtering Analysis"
s:name: "Single-cell RNA-Seq Filtering Analysis"
s:alternateName: "Filters single-cell RNA-Seq datasets based on the common QC metrics"

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
  Single-cell RNA-Seq Filtering Analysis

  Filters single-cell RNA-Seq datasets based on the common QC metrics.


s:about: |
  usage: sc_rna_filter.R [-h] --mex MEX [MEX ...] --identity
                                        IDENTITY [--grouping GROUPING]
                                        [--barcodes BARCODES]
                                        [--rnamincells RNAMINCELLS]
                                        [--mingenes [MINGENES [MINGENES ...]]]
                                        [--maxgenes [MAXGENES [MAXGENES ...]]]
                                        [--rnaminumi [RNAMINUMI [RNAMINUMI ...]]]
                                        [--minnovelty [MINNOVELTY [MINNOVELTY ...]]]
                                        [--mitopattern MITOPATTERN]
                                        [--maxmt MAXMT] [--pdf] [--verbose]
                                        [--h5seurat] [--h5ad] [--output OUTPUT]
                                        [--theme {gray,bw,linedraw,light,dark,minimal,classic,void}]
                                        [--cpus CPUS] [--memory MEMORY]

  Single-cell RNA-Seq Filtering Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --mex MEX [MEX ...]   Path to the folder with feature-barcode matrix from
                          Cell Ranger Count/Aggregate experiment in MEX format.
                          If multiple locations provided data is assumed to be
                          not aggregated (outputs from the multiple Cell Ranger
                          Count experiments) and will be merged before the
                          analysis.
    --identity IDENTITY   Path to the metadata TSV/CSV file to set the datasets
                          identities. If '--mex' points to the Cell Ranger
                          Aggregate outputs, the aggregation.csv file can be
                          used. In case of using feature-barcode matrices from a
                          single or multiple Cell Ranger Count experiments the
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
                          '--mex' input based on the '--identity' file. Default:
                          250 (applied to all datasets)
    --maxgenes [MAXGENES [MAXGENES ...]]
                          Include cells with the number of genes not bigger than
                          this value. If multiple values provided, each of them
                          will be applied to the correspondent dataset from the
                          '--mex' input based on the '--identity' file. Default:
                          5000 (applied to all datasets)
    --rnaminumi [RNAMINUMI [RNAMINUMI ...]]
                          Include cells where at least this many UMI
                          (transcripts) are detected. If multiple values
                          provided, each of them will be applied to the
                          correspondent dataset from the '--mex' input based on
                          the '--identity' file. Default: 500 (applied to all
                          datasets)
    --minnovelty [MINNOVELTY [MINNOVELTY ...]]
                          Include cells with the novelty score not lower than
                          this value, calculated for as log10(genes)/log10(UMI).
                          If multiple values provided, each of them will be
                          applied to the correspondent dataset from the '--mex'
                          input based on the '--identity' file. Default: 0.8
                          (applied to all datasets)
    --mitopattern MITOPATTERN
                          Regex pattern to identify mitochondrial genes.
                          Default: '^Mt-'
    --maxmt MAXMT         Include cells with the percentage of transcripts
                          mapped to mitochondrial genes not bigger than this
                          value. Default: 5 (applied to all datasets)
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --h5ad                Save Seurat data to h5ad file. Default: false
    --output OUTPUT       Output prefix. Default: ./sc
    --theme {gray,bw,linedraw,light,dark,minimal,classic,void}
                          Color theme for all generated plots. Default: classic
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple '--cpus'. Default: 32