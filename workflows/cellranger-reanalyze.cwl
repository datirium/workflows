cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement


"sd:upstream":
  sc_experiment:
  - "single-cell-preprocess-cellranger.cwl"
  - "cellranger-multi.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/alias"
    sd:preview:
      position: 1

  filtered_feature_bc_matrix_h5:
    type: File
    label: "Single-cell Experiment"
    doc: "Filtered feature-barcode matrices in HDF5 format from cellranger count/multi"
    "sd:upstreamSource": "sc_experiment/filtered_feature_bc_matrix_h5"
    "sd:localLabel": true

  selected_barcodes:
    type: File?
    label: "CSV file containing a list of cell barcodes to use for reanalysis"
    doc: |
      A CSV file containing a list of cell barcodes to use for reanalysis,
      e.g. barcodes exported from Loupe Browser. All barcodes must be present
      in the matrix.

  selected_genes:
    type: File?
    label: "CSV file containing a list of gene IDs to use for reanalysis"
    doc: |
      A CSV file containing a list of gene IDs to use for reanalysis (corresponding
      to the gene_id field of the reference GTF). All gene IDs must be present in
      the matrix. Note that only gene features are used in secondary analysis.
      Feature Barcode features are ignored.

  excluded_genes:
    type: File?
    label: "CSV file containing a list of gene IDs to exclude for reanalysis. Applied after setting selected genes"
    doc: |
      A CSV file containing a list of gene IDs to exclude for reanalysis (corresponding
      to the gene_id field of the reference GTF). All gene IDs must be present in
      the matrix. The exclusion is applied after setting the gene list with --genes.
      Note that only gene features are used in secondary analysis. Feature Barcode features
      are ignored.

  force_cells:
    type: int?
    default: null
    label: "Force pipeline to use this number of cells, bypassing the cell detection algorithm"
    doc: |
      Force pipeline to use this number of cells, bypassing the cell detection algorithm.
      Use this if the number of cells estimated by Cell Ranger is not consistent with the
      barcode rank plot. If specifying a value that exceeds the original cell count, you
      must use the raw_gene_bc_matrices_h5.h5
    "sd:layout":
      advanced: true

  num_analysis_bcs:
    type: int?
    default: null
    label: "Randomly subset data to N barcodes for all analysis. Reduce this parameter if you want to improve performance or simulate results from lower cell counts"
    doc: |
      Randomly subset data to N barcodes for all analysis. Reduce this parameter if you
      want to improve performance or simulate results from lower cell counts. Cannot be
      set higher than the available number of cells.
      Default: null
    "sd:layout":
      advanced: true

  num_pca_bcs:
    type: int?
    default: null
    label: "Randomly subset data to N barcodes when computing PCA projection. Try reducing this parameter if your analysis is running out of memory"
    doc: |
      Randomly subset data to N barcodes when computing PCA projection (the most memory-intensive
      step). The PCA projection will still be applied to the full dataset, i.e. your final results
      will still reflect all the data. Try reducing this parameter if your analysis is running out
      of memory. Cannot be set higher than the available number of cells.
      Default: null
    "sd:layout":
      advanced: true

  num_pca_genes:
    type: int?
    default: null
    label:  "Subset data to the top N genes when computing PCA. Try reducing this parameter if your analysis is running out of memory"
    doc: |
      Subset data to the top N genes (ranked by normalized dispersion) when computing PCA.
      Differential expression will still reflect all genes. Try reducing this parameter if
      your analysis is running out of memory. Cannot be set higher than the number of genes
      in the reference transcriptome.
      Default: null
    "sd:layout":
      advanced: true

  num_principal_comps:
    type: int?
    default: 10
    label: "Compute N principal components for PCA. Setting this too high may cause spurious clusters to be called"
    doc: |
      Compute N principal components for PCA. Setting this too high may cause spurious clusters
      to be called. The default value is 100 when the chemistry batch correction is enabled.
      Set from 10 to 100, depending on the number of cell populations/clusters you expect to see.
      Default: 10
    "sd:layout":
      advanced: true

  cbc_knn:
    type: int?
    default: 10
    label: "Specify the number of nearest neighbors used to identify mutual nearest neighbors. Setting this too high will increase runtime"
    doc: |
      Specify the number of nearest neighbors used to identify mutual nearest neighbors.
      Setting this too high will increase runtime and may cause out of memory error.
      See Chemistry Batch Correction page for more details. Ranges from 5 to 20.
      Default: 10
    "sd:layout":
      advanced: true

  cbc_alpha:
    type: float?
    default: 0.1
    label: "Specify the threshold of the percentage of matched cells between two batches, which is used to determine if the batch pair will be merged"
    doc: |
      Specify the threshold of the percentage of matched cells between two batches,
      which is used to determine if the batch pair will be merged. See Chemistry
      Batch Correction page for more details. Ranges from 0.05 to 0.5.
      Default: 0.1
    "sd:layout":
      advanced: true
      
  cbc_sigma:
    type: float?
    default: 150
    label: "Specify the bandwidth of the Gaussian smoothing kernel used to compute the correction vector for each cell"
    doc: |
      Specify the bandwidth of the Gaussian smoothing kernel used to compute the correction
      vector for each cell. See Chemistry Batch Correction page for more details. Ranges
      from 10 to 500.
      Default: 150
    "sd:layout":
      advanced: true

  cbc_realign_panorama:
    type: boolean?
    default: false
    label: "Specify if two batches will be merged if they are already in the same panorama. Setting this to True will usually improve the performance"
    doc: |
      Specify if two batches will be merged if they are already in the same panorama. Setting
      this to True will usually improve the performance, but will also increase runtime and
      memory usage. See Chemistry Batch Correction page for more details. One of true or false.
      Default: false
    "sd:layout":
      advanced: true

  graphclust_neighbors:
    type: int?
    default: 0
    label: "Number of nearest-neighbors to use in the graph-based clustering. Lower values result in higher-granularity clustering"
    doc: |
      Number of nearest-neighbors to use in the graph-based clustering. Lower values result in
      higher-granularity clustering. The actual number of neighbors used is the maximum of this
      value and that determined by neighbor_a and neighbor_b. Set this value to zero to use those
      values instead. Ranged from 10 to 500, depending on desired granularity.
      Default: 0
    "sd:layout":
      advanced: true

  neighbor_a:
    type: float?
    default: -230.0
    label: "neighbor_a parameter for the number of nearest neighbors k = neighbor_a + neighbor_b * log10(n_cells)"
    doc: |
      The number of nearest neighbors, k, used in the graph-based clustering is computed as follows:
      k = neighbor_a + neighbor_b * log10(n_cells). The actual number of neighbors used is the maximum
      of this value and graphclust_neighbors. Determines how clustering granularity scales with cell count.
      Default: -230.0
    "sd:layout":
      advanced: true

  neighbor_b:
    type: float?
    default: 120.0
    label: "neighbor_b parameter for the number of nearest neighbors k = neighbor_a + neighbor_b * log10(n_cells)"
    doc: |
      The number of nearest neighbors, k, used in the graph-based clustering is computed as follows:
      k = neighbor_a + neighbor_b * log10(n_cells). The actual number of neighbors used is the maximum of
      this value and graphclust_neighbors. Determines how clustering granularity scales with cell count.
      Default: 120.0
    "sd:layout":
      advanced: true

  max_clusters:
    type: int?
    default: 10
    label: "Compute K-means clustering using K values of 2 to N. Setting this too high may cause spurious clusters to be called"
    doc: |
      Compute K-means clustering using K values of 2 to N. Setting this too high may cause spurious clusters
      to be called. Ranges from 10 to 50, depending on the number of cell populations / clusters you expect to see.
      Default: 10
    "sd:layout":
      advanced: true

  tsne_input_pcs:
    type: int?
    default: null
    label: "Subset to top N principal components for TSNE. Change this parameter if you want to see how the TSNE plot changes when using fewer PCs"
    doc: |
      Subset to top N principal components for TSNE. Change this parameter if you want to see how the TSNE plot
      changes when using fewer PCs, independent of the clustering / differential expression. You may find that TSNE
      is faster and/or the output looks better when using fewer PCs. Cannot be set higher than
      the num_principal_comps parameter.
      Default: null
    "sd:layout":
      advanced: true

  tsne_perplexity:
    type: int?
    default: 30
    label: "TSNE perplexity parameter. When analyzing 100k+ cells, increasing this parameter may improve TSNE results"
    doc: |
      TSNE perplexity parameter (see the TSNE FAQ for more details). When analyzing 100k+ cells, increasing this
      parameter may improve TSNE results, but the algorithm will be slower. Ranges from 30 to 50.
      Default: 30
    "sd:layout":
      advanced: true

  tsne_theta:
    type: float?
    default: 0.5
    label: "TSNE theta parameter. Higher values yield faster, more approximate results (and vice versa)"
    doc: |
      TSNE theta parameter (see the TSNE FAQ for more details). Higher values yield faster, more approximate results
      (and vice versa). The runtime and memory performance of TSNE will increase dramatically if you set this below 0.25.
      Ranges from 0 to 1.
      Default: 0.5
    "sd:layout":
      advanced: true

  tsne_max_dims:
    type: int?
    default: 2
    label: "Maximum number of TSNE output dimensions. Set this to 3 to produce both 2D and 3D TSNE projections"
    doc: |
      Maximum number of TSNE output dimensions. Set this to 3 to produce both 2D and 3D TSNE projections
      (note: runtime will increase significantly). Ranges from 2 to 3.
      Default: 2
    "sd:layout":
      advanced: true

  tsne_max_iter:
    type: int?
    default: 1000
    label: "Number of total TSNE iterations. Try increasing this if TSNE results do not look good on larger numbers of cells"
    doc: |
      Number of total TSNE iterations. Try increasing this if TSNE results do not look good on larger numbers
      of cells. Runtime increases linearly with number of iterations. Ranges from 1000 to 10000.
      Default: 1000
    "sd:layout":
      advanced: true

  tsne_stop_lying_iter:
    type: int?
    default: 250
    label: "Iteration at which TSNE learning rate is reduced. Try increasing this if TSNE results do not look good on larger numbers of cells"
    doc: |
      Iteration at which TSNE learning rate is reduced. Try increasing this if TSNE results do not look good
      on larger numbers of cells. Cannot be set higher than tsne_max_iter.
      Default: 250
    "sd:layout":
      advanced: true

  tsne_mom_switch_iter:
    type: int?
    default: 250
    label: "Iteration at which TSNE momentum is reduced. Try increasing this if TSNE results do not look good on larger numbers of cells"
    doc: |
      Iteration at which TSNE momentum is reduced. Try increasing this if TSNE results do not look good on
      larger numbers of cells. Cannot be set higher than tsne_max_iter. Cannot be set higher than tsne_max_iter.
      Default: 250
    "sd:layout":
      advanced: true

  umap_input_pcs:
    type: int?
    default: null
    label: "Subset to top N principal components for UMAP. Change this parameter if you want to see how the UMAP plot changes when using fewer PCs"
    doc: |
      Subset to top N principal components for UMAP. Change this parameter if you want to see how the UMAP plot
      changes when using fewer PCs, independent of the clustering / differential expression. You may find that
      UMAP is faster and/or the output looks better when using fewer PCs. Cannot be set higher than the
      num_principal_comps parameter.
      Default: null
    "sd:layout":
      advanced: true

  umap_n_neighbors:
    type: int?
    default: 30
    label: "Determines the number of neighboring points used in local approximations of manifold structure"
    doc: |
      Determines the number of neighboring points used in local approximations of manifold structure.
      Larger values will usually result in more global structure at the loss of detailed local structure.
      Ranges from 5 to 50.
      Default: 30
    "sd:layout":
      advanced: true

  umap_max_dims:
    type: int?
    default: 2
    label: "Maximum number of UMAP output dimensions. Set this to 3 to produce both 2D and 3D UMAP projections"
    doc: |
      Maximum number of UMAP output dimensions. Set this to 3 to produce both 2D and 3D UMAP projections.
      Ranges from 2 to 3.
      Default: 2
    "sd:layout":
      advanced: true

  umap_min_dist:
    type: float?
    default: 0.3
    label: "Controls how tightly the embedding is allowed to pack points together"
    doc: |
      Controls how tightly the embedding is allowed to pack points together. Larger values make embedded
      points are more evenly distributed, while smaller values make the embedding more accurately with
      regard to the local structure. Ranges from 0.001 to 0.5.
      Default: 0.3
    "sd:layout":
      advanced: true

  umap_metric:
    type:
    - "null"
    - type: enum
      symbols:
      - euclidean
      - manhattan
      - chebyshev
      - minkowski
      - canberra
      - braycurtis
      - haversine
      - mahalanobis
      - wminkowski
      - seuclidean
      - cosine
      - correlation
      - hamming
      - jaccard
      - dice
      - russellrao
      - kulsinski
      - rogerstanimoto
      - sokalmichener
      - sokalsneath
      - yule
    default: "correlation"
    label: "Determines how the distance is computed in the input space"
    doc: |
      Determines how the distance is computed in the input space.
      Default: "correlation"
    "sd:layout":
      advanced: true

  random_seed:
    type: int?
    default: 0
    label: "Random seed"
    doc: |
      Random seed. Due to the randomized nature of the algorithms, changing this will produce slightly
      different results. If the TSNE or UMAP results don't look good, try running multiple times with
      different seeds and pick the TSNE or UMAP that looks best.
      Default: 0
    "sd:layout":
      advanced: true

  threads:
    type: int?
    default: 4
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    "sd:layout":
      advanced: true

  memory_limit:
    type: int?
    default: 30
    label: "Maximum memory used (GB)"
    doc: "Maximum memory used (GB). The same will be applied to virtual memory"
    "sd:layout":
      advanced: true


outputs:

  secondary_analysis_report_folder:
    type: File
    outputSource: compress_secondary_analysis_report_folder/compressed_folder
    label: "Compressed folder with reanalyzed secondary analysis results"
    doc: |
      Compressed folder with secondary analysis results including dimensionality reduction,
      cell clustering, and differential expression of reanalyzed results

  web_summary_report:
    type: File
    outputSource: reanalyze/web_summary_report
    label: "Reanalyzed run summary metrics and charts in HTML format"
    doc: |
      Reanalyzed run summary metrics and charts in HTML format
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  filtered_feature_bc_matrix_folder:
    type: File
    outputSource: compress_filtered_feature_bc_matrix_folder/compressed_folder
    label: "Compressed folder with filtered feature-barcode matrices"
    doc: |
      Compressed folder with filtered feature-barcode matrices containing only cellular barcodes in MEX format.
      When implemented, in Targeted Gene Expression samples, the non-targeted genes won't be present.

  reanalyze_params:
    type: File
    outputSource: reanalyze/reanalyze_params
    label: "Reanalyze params in CSV format"
    doc: |
      Reanalyze params in CSV format

  loupe_browser_track:
    type: File
    outputSource: reanalyze/loupe_browser_track
    label: "Loupe Browser visualization and analysis file for reanalyzed results"
    doc: |
      Loupe Browser visualization and analysis file for reanalyzed results

  compressed_html_data_folder:
    type: File
    outputSource: compress_html_data_folder/compressed_folder
    label: "Compressed folder with CellBrowser formatted results"
    doc: |
      Compressed folder with CellBrowser formatted results

  html_data_folder:
    type: Directory
    outputSource: cellbrowser_build/html_data
    label: "Folder with not compressed CellBrowser formatted results"
    doc: |
      Folder with not compressed CellBrowser formatted results

  cellbrowser_report:
    type: File
    outputSource: cellbrowser_build/index_html_file
    label: "CellBrowser formatted Cellranger report"
    doc: |
      CellBrowser formatted Cellranger report
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  reanalyze_stdout_log:
    type: File
    outputSource: reanalyze/stdout_log
    label: "stdout log generated by cellranger reanalyze"
    doc: |
      stdout log generated by cellranger agreanalyzegr

  reanalyze_stderr_log:
    type: File
    outputSource: reanalyze/stderr_log
    label: "stderr log generated by cellranger reanalyze"
    doc: |
      stderr log generated by cellranger reanalyze


steps:

  reanalyze:
    run: ../tools/cellranger-reanalyze.cwl
    in:
      feature_bc_matrix_h5: filtered_feature_bc_matrix_h5
      selected_barcodes: selected_barcodes
      selected_genes: selected_genes
      excluded_genes: excluded_genes
      force_cells: force_cells
      threads: threads
      memory_limit: memory_limit
      virt_memory_limit: memory_limit
      num_analysis_bcs: num_analysis_bcs
      num_pca_bcs: num_pca_bcs
      num_pca_genes: num_pca_genes
      num_principal_comps: num_principal_comps
      cbc_knn: cbc_knn
      cbc_alpha: cbc_alpha
      cbc_sigma: cbc_sigma
      cbc_realign_panorama: cbc_realign_panorama
      graphclust_neighbors: graphclust_neighbors
      neighbor_a: neighbor_a
      neighbor_b: neighbor_b
      max_clusters: max_clusters
      tsne_input_pcs: tsne_input_pcs
      tsne_perplexity: tsne_perplexity
      tsne_theta: tsne_theta
      tsne_max_dims: tsne_max_dims
      tsne_max_iter: tsne_max_iter
      tsne_stop_lying_iter: tsne_stop_lying_iter
      tsne_mom_switch_iter: tsne_mom_switch_iter
      umap_input_pcs: umap_input_pcs
      umap_n_neighbors: umap_n_neighbors
      umap_max_dims: umap_max_dims
      umap_min_dist: umap_min_dist
      umap_metric: umap_metric
      random_seed: random_seed
    out:
    - secondary_analysis_report_folder
    - web_summary_report
    - filtered_feature_bc_matrix_folder
    - reanalyze_params
    - loupe_browser_track
    - stdout_log
    - stderr_log

  compress_filtered_feature_bc_matrix_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: reanalyze/filtered_feature_bc_matrix_folder
    out:
    - compressed_folder

  compress_secondary_analysis_report_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: reanalyze/secondary_analysis_report_folder
    out:
    - compressed_folder

  cellbrowser_build:
    run: ../tools/cellbrowser-build-cellranger.cwl
    in:
      secondary_analysis_report_folder: reanalyze/secondary_analysis_report_folder
      filtered_feature_bc_matrix_folder: reanalyze/filtered_feature_bc_matrix_folder
    out:
    - html_data
    - index_html_file

  compress_html_data_folder:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: cellbrowser_build/html_data
    out:
    - compressed_folder


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Cellranger Reanalyze"
s:name: "Cellranger Reanalyze"
s:alternateName: "Reruns secondary analysis for Cell Ranger Count Gene Expression or Cell Ranger Multi experiments"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/cellranger-reanalyze.cwl
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


doc: |
  Cellranger Reanalyze
  ====================