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
    R_MAX_VSIZE: $(inputs.vector_memory_limit * 1000000000)


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.2


inputs:

  feature_bc_matrices_folder:
    type: Directory
    inputBinding:
      prefix: "--mex"
    doc: |
      Path to the folder with feature-barcode matrix from Cell Ranger Count/Aggregate
      experiment in MEX format. If multiple locations provided data is assumed to be not
      aggregated (outputs from multiple Cell Ranger Count experiments) and will be merged
      before the analysis.

  aggregation_metadata:
    type: File?
    doc: |
      Path to the metadata TSV/CSV file to set the datasets identities. If --mex points to
      the Cell Ranger Aggregate outputs, the aggregation.csv file can be used. In case of
      using feature-barcode matrices from a single or multiple Cell Ranger Count experiments
      the file with identities should include at least one column - 'library_id', and a row
      with aliases per each experiment from the --mex input. The order of rows should correspond
      to the order of feature-barcode matrices provided in the --mex parameter.

  grouping_data:
    type: File?
    inputBinding:
      prefix: "--grouping"
    doc: |
      Path to the TSV/CSV file to define datasets grouping. First column -
      'library_id' with the values provided in the same order as in the
      correspondent column of the --identity file, second column 'condition'.
      Default: each dataset is assigned to a separate group.

  barcodes_data:
    type: File?
    inputBinding:
      prefix: "--barcodes"
    doc: |
      Path to the headerless TSV/CSV file with the list of barcodes to select
      cells of interest (one barcode per line). Prefilters input feature-barcode
      matrix to include only selected cells.
      Default: use all cells.

  gex_minimum_cells:
    type: int?
    inputBinding:
      prefix: "--gexmincells"
    doc: |
      Include only GEX features detected in at least this many cells. Ignored when
      --mex points to the feature-barcode matrices from the multiple Cell Ranger
      Count experiments.
      Default: 5 (applied to all datasets)

  gex_minimum_features:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--mingenes"
    doc: |
      Include cells where at least this many GEX features are detected.
      If multiple values provided, each of them will be applied to the
      correspondent dataset from the --mex input based on the --identity
      file.
      Default: 250 (applied to all datasets)

  gex_maximum_features:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--maxgenes"
    doc: |
      Include cells with the number of GEX features not bigger than this value.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the --mex input based on the --identity file.
      Default: 5000 (applied to all datasets)

  gex_minimum_umis:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      prefix: "--gexminumi"
    doc: |
      Include cells where at least this many GEX UMIs (transcripts) are detected.
      If multiple values provided, each of them will be applied to the correspondent
      dataset from the --mex input based on the --identity file.
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
      GEX as log10(genes)/log10(UMIs). If multiple values provided, each of them will
      be applied to the correspondent dataset from the --mex input based on the
      --identity file.
      Default: 0.8 (applied to all datasets)

  mito_pattern:
    type: string?
    inputBinding:
      prefix: "--mitopattern"
    doc: |
      Regex pattern to identify mitochondrial GEX features.
      Default: '^Mt-'

  maximum_mito_perc:
    type: float?
    inputBinding:
      prefix: "--maxmt"
    doc: |
      Include cells with the percentage of GEX transcripts mapped to mitochondrial
      genes not bigger than this value.
      Default: 5 (applied to all datasets)

  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: "--pdf"
    doc: |
      Export plots in PDF.
      Default: false

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

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix.
      Default: ./seurat

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

  raw_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_pca_1_2_qc_mtrcs.png"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (not filtered).
      PNG format

  raw_pca_1_2_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_pca_1_2_qc_mtrcs.pdf"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (not filtered).
      PDF format

  raw_pca_2_3_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_pca_2_3_qc_mtrcs.png"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (not filtered).
      PNG format

  raw_pca_2_3_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_pca_2_3_qc_mtrcs.pdf"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (not filtered).
      PDF format

  raw_cell_count_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_cell_count.png"
    doc: |
      Number of cells per dataset (not filtered).
      PNG format

  raw_cell_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_cell_count.pdf"
    doc: |
      Number of cells per dataset (not filtered).
      PDF format

  raw_gex_umi_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gex_umi_dnst_spl_by_cond.png"
    doc: |
      Split by condition GEX UMI density per cell (not filtered).
      PNG format

  raw_gex_umi_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gex_umi_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition GEX UMI density per cell (not filtered).
      PDF format

  raw_gene_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst_spl_by_cond.png"
    doc: |
      Split by condition gene density per cell (not filtered).
      PNG format

  raw_gene_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition gene density per cell (not filtered).
      PDF format

  raw_gene_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_corr.png"
    doc: |
      Genes vs GEX UMIs per cell correlation (not filtered).
      PNG format

  raw_gene_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_gene_umi_corr.pdf"
    doc: |
      Genes vs GEX UMIs per cell correlation (not filtered).
      PDF format

  raw_mito_perc_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_mito_perc_dnst_spl_by_cond.png"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PNG format

  raw_mito_perc_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_mito_perc_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (not filtered).
      PDF format

  raw_nvlt_score_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_score_dnst_spl_by_cond.png"
    doc: |
      Split by condition novelty score density per cell (not filtered).
      PNG format

  raw_nvlt_score_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_nvlt_score_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition novelty score density per cell (not filtered).
      PDF format

  raw_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs.png"
    doc: |
      QC metrics densities per cell (not filtered).
      PNG format

  raw_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_raw_qc_mtrcs.pdf"
    doc: |
      QC metrics densities per cell (not filtered).
      PDF format

  fltr_pca_1_2_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_pca_1_2_qc_mtrcs.png"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (filtered).
      PNG format

  fltr_pca_1_2_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_pca_1_2_qc_mtrcs.pdf"
    doc: |
      PC1 and PC2 of ORQ-transformed QC metrics PCA (filtered).
      PDF format

  fltr_pca_2_3_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_pca_2_3_qc_mtrcs.png"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (filtered).
      PNG format

  fltr_pca_2_3_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_pca_2_3_qc_mtrcs.pdf"
    doc: |
      PC2 and PC3 of ORQ-transformed QC metrics PCA (filtered).
      PDF format

  fltr_cell_count_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_cell_count.png"
    doc: |
      Number of cells per dataset (filtered).
      PNG format

  fltr_cell_count_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_cell_count.pdf"
    doc: |
      Number of cells per dataset (filtered).
      PDF format

  fltr_gex_umi_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gex_umi_dnst_spl_by_cond.png"
    doc: |
      Split by condition GEX UMI density per cell (filtered).
      PNG format

  fltr_gex_umi_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_gex_umi_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition GEX UMI density per cell (filtered).
      PDF format

  fltr_gene_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_dnst_spl_by_cond.png"
    doc: |
      Split by condition gene density per cell (filtered).
      PNG format

  fltr_gene_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition gene density per cell (filtered).
      PDF format

  fltr_gene_umi_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_umi_corr.png"
    doc: |
      Genes vs GEX UMIs per cell correlation (filtered).
      PNG format

  fltr_gene_umi_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_gene_umi_corr.pdf"
    doc: |
      Genes vs GEX UMIs per cell correlation (filtered).
      PDF format

  fltr_mito_perc_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_perc_dnst_spl_by_cond.png"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered).
      PNG format

  fltr_mito_perc_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_mito_perc_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition density of transcripts mapped to mitochondrial genes per cell (filtered).
      PDF format

  fltr_nvlt_score_dnst_spl_by_cond_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_nvlt_score_dnst_spl_by_cond.png"
    doc: |
      Split by condition novelty score density per cell (filtered).
      PNG format

  fltr_nvlt_score_dnst_spl_by_cond_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_nvlt_score_dnst_spl_by_cond.pdf"
    doc: |
      Split by condition novelty score density per cell (filtered).
      PDF format

  fltr_qc_mtrcs_plot_png:
    type: File?
    outputBinding:
      glob: "*_fltr_qc_mtrcs.png"
    doc: |
      QC metrics densities per cell (filtered).
      PNG format

  fltr_qc_mtrcs_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_fltr_qc_mtrcs.pdf"
    doc: |
      QC metrics densities per cell (filtered).
      PDF format

  seurat_filtered_data_rds:
    type: File
    outputBinding:
      glob: "*_filtered_data.rds"
    doc: |
      Filtered Seurat data in RDS format

  seurat_filtered_data_h5seurat:
    type: File?
    outputBinding:
      glob: "*_filtered_data.h5seurat"
    doc: |
      Filtered Seurat data in h5seurat format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sc_gex_filter.R"]
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


stdout: sc_gex_filter_stdout.log
stderr: sc_gex_filter_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Single-cell GEX Filter"
s:name: "Single-cell GEX Filter"
s:alternateName: "Filters single-cell GEX datasets based on the common QC metrics"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/sc-gex-filter.cwl
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
  Single-cell GEX Filter
  ================================================================
  Filters single-cell GEX datasets based on the common QC metrics.


s:about: |
  usage: sc_gex_filter.R
        [-h] --mex MEX [MEX ...] --identity IDENTITY [--grouping GROUPING]
        [--barcodes BARCODES] [--gexmincells GEXMINCELLS]
        [--mingenes [MINGENES ...]] [--maxgenes [MAXGENES ...]]
        [--gexminumi [GEXMINUMI ...]] [--minnovelty [MINNOVELTY ...]]
        [--mitopattern MITOPATTERN] [--maxmt MAXMT] [--pdf] [--verbose]
        [--h5seurat] [--output OUTPUT] [--cpus CPUS] [--memory MEMORY]

  Seurat GEX Filtering Analysis

  optional arguments:
    -h, --help            show this help message and exit
    --mex MEX [MEX ...]   Path to the folder with feature-barcode matrix from
                          Cell Ranger Count/Aggregate experiment in MEX format.
                          If multiple locations provided data is assumed to be
                          not aggregated (outputs from multiple Cell Ranger
                          Count experiments) and will be merged before the
                          analysis.
    --identity IDENTITY   Path to the metadata TSV/CSV file to set the datasets
                          identities. If --mex points to the Cell Ranger
                          Aggregate outputs, the aggregation.csv file can be
                          used. In case of using feature-barcode matrices from a
                          single or multiple Cell Ranger Count experiments the
                          file with identities should include at least one
                          column - 'library_id', and a row with aliases per each
                          experiment from the --mex input. The order of rows
                          should correspond to the order of feature-barcode
                          matrices provided in the --mex parameter.
    --grouping GROUPING   Path to the TSV/CSV file to define datasets grouping.
                          First column - 'library_id' with the values provided
                          in the same order as in the correspondent column of
                          the --identity file, second column 'condition'.
                          Default: each dataset is assigned to a separate group.
    --barcodes BARCODES   Path to the headerless TSV/CSV file with the list of
                          barcodes to select cells of interest (one barcode per
                          line). Prefilters input feature-barcode matrix to
                          include only selected cells. Default: use all cells.
    --gexmincells GEXMINCELLS
                          Include only GEX features detected in at least this
                          many cells. Ignored when --mex points to the feature-
                          barcode matrices from the multiple Cell Ranger Count
                          experiments. Default: 5 (applied to all datasets)
    --mingenes [MINGENES ...]
                          Include cells where at least this many GEX features
                          are detected. If multiple values provided, each of
                          them will be applied to the correspondent dataset from
                          the --mex input based on the --identity file. Default:
                          250 (applied to all datasets)
    --maxgenes [MAXGENES ...]
                          Include cells with the number of GEX features not
                          bigger than this value. If multiple values provided,
                          each of them will be applied to the correspondent
                          dataset from the --mex input based on the --identity
                          file. Default: 5000 (applied to all datasets)
    --gexminumi [GEXMINUMI ...]
                          Include cells where at least this many GEX UMIs
                          (transcripts) are detected. If multiple values
                          provided, each of them will be applied to the
                          correspondent dataset from the --mex input based on
                          the --identity file. Default: 500 (applied to all
                          datasets)
    --minnovelty [MINNOVELTY ...]
                          Include cells with the novelty score not lower than
                          this value, calculated for GEX as
                          log10(genes)/log10(UMIs). If multiple values provided,
                          each of them will be applied to the correspondent
                          dataset from the --mex input based on the --identity
                          file. Default: 0.8 (applied to all datasets)
    --mitopattern MITOPATTERN
                          Regex pattern to identify mitochondrial GEX features.
                          Default: '^Mt-'
    --maxmt MAXMT         Include cells with the percentage of GEX transcripts
                          mapped to mitochondrial genes not bigger than this
                          value. Default: 5 (applied to all datasets)
    --pdf                 Export plots in PDF. Default: false
    --verbose             Print debug information. Default: false
    --h5seurat            Save Seurat data to h5seurat file. Default: false
    --output OUTPUT       Output prefix. Default: ./seurat
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32