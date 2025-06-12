cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: cumulusprod/cellranger:8.0.1
hints:
- class: InitialWorkDirRequirement
  listing: |
    ${
        const skipped_ids = [
          "feature_bc_matrix_h5",
          "selected_barcodes",
          "selected_genes",
          "excluded_genes",
          "force_cells",
          "threads",
          "memory_limit",
          "virt_memory_limit"
        ]
        var entry = "";
        for (const id in inputs){
          if (skipped_ids.includes(id) || !inputs[id]){
            continue;
          }
          entry += id + "," + inputs[id] + "\n";
        }
        return [{
          "entry": entry,
          "entryname": "params.csv"
        }];
    }
inputs:
  feature_bc_matrix_h5:
    type: File
    inputBinding:
      position: 5
      prefix: --matrix
    doc: |
      A feature-barcode matrix containing
      data for one genome. Should be the
      filtered version, unless using
      --force-cells
  selected_barcodes:
    type: File?
    inputBinding:
      position: 6
      prefix: --barcodes
    doc: |
      A CSV file containing a list of cell
      barcodes to use for reanalysis, e.g.
      barcodes exported from Loupe Browser.
      All barcodes must be present in the
      matrix.
  selected_genes:
    type: File?
    inputBinding:
      position: 7
      prefix: --genes
    doc: |
      A CSV file containing a list of gene
      IDs to use for reanalysis (corresponding
      to the gene_id field of the reference
      GTF). All gene IDs must be present in
      the matrix.
  excluded_genes:
    type: File?
    inputBinding:
      position: 8
      prefix: --exclude-genes
    doc: |
      A CSV file containing a list of gene IDs
      to exclude for reanalysis (corresponding
      to the gene_id field of the reference GTF).
      All gene IDs must be present in the matrix.
      The exclusion is applied after setting the
      gene list with --genes.
  force_cells:
    type: int?
    inputBinding:
      position: 9
      prefix: --force-cells
    doc: |
      Force pipeline to use this number of cells,
      bypassing the cell detection algorithm.
      Use this if the number of cells estimated
      by Cell Ranger is not consistent with the
      barcode rank plot. If specifying a value
      that exceeds the original cell count, you
      must use the raw_gene_bc_matrices_h5.h5
  threads:
    type: int?
    inputBinding:
      position: 10
      prefix: --localcores
    doc: |
      Set max cores the pipeline may request
      at one time.
      Default: all available
  memory_limit:
    type: int?
    inputBinding:
      position: 11
      prefix: --localmem
    doc: |
      Set max GB the pipeline may request
      at one time
      Default: all available
  virt_memory_limit:
    type: int?
    inputBinding:
      position: 12
      prefix: --localvmem
    doc: |
      Set max virtual address space in
      GB for the pipeline
      Default: all available
  num_analysis_bcs:
    type: int?
    doc: "Randomly subset data to N barcodes for all analysis. \nReduce this parameter if you want to improve \nperformance or simulate results from lower cell counts. \nCannot be set higher than the available number of cells. \nDefault: null\n"
  num_pca_bcs:
    type: int?
    doc: "Randomly subset data to N barcodes when computing \nPCA projection (the most memory-intensive step). \nThe PCA projection will still be applied to the full \ndataset, i.e. your final results will still reflect all \nthe data. Try reducing this parameter if your \nanalysis is running out of memory. Cannot be set \nhigher than the available number of cells.\nDefault: null\n"
  num_pca_genes:
    type: int?
    doc: "Subset data to the top N genes (ranked by normalized \ndispersion) when computing PCA. Differential \nexpression will still reflect all genes. Try reducing \nthis parameter if your analysis is running out of memory. \nCannot be set higher than the number of genes in the \nreference transcriptome.\nDefault: null\n"
  num_principal_comps:
    type: int?
    doc: "Compute N principal components for PCA. Setting this \ntoo high may cause spurious clusters to be called. The \ndefault value is 100 when the chemistry batch correction \nis enabled. Set from 10 to 100, depending on the number \nof cell populations/clusters you expect to see. \nDefault: 10\n"
  cbc_knn:
    type: int?
    doc: "Specify the number of nearest neighbors used to identify \nmutual nearest neighbors. Setting this too high will \nincrease runtime and may cause out of memory error. See \nChemistry Batch Correction page for more details. Ranges \nfrom 5 to 20.\nDefault: 10\n"
  cbc_alpha:
    type: float?
    doc: "Specify the threshold of the percentage of matched cells \nbetween two batches, which is used to determine if the \nbatch pair will be merged. See Chemistry Batch Correction \npage for more details. Ranges from 0.05 to 0.5.\nDefault: 0.1\n"
  cbc_sigma:
    type: float?
    doc: "Specify the bandwidth of the Gaussian smoothing kernel \nused to compute the correction vector for each cell. See \nChemistry Batch Correction page for more details. Ranges \nfrom 10 to 500.\nDefault: 150\n"
  cbc_realign_panorama:
    type: boolean?
    doc: "Specify if two batches will be merged if they are already \nin the same panorama. Setting this to True will usually \nimprove performance, but will also increase runtime and \nmemory usage. See Chemistry Batch Correction page for \nmore details. One of true or false.\nDefault: false\n"
  graphclust_neighbors:
    type: int?
    doc: "Number of nearest-neighbors to use in the graph-based \nclustering. Lower values result in higher-granularity \nclustering. The actual number of neighbors used is the \nmaximum of this value and that determined by neighbor_a \nand neighbor_b. Set this value to zero to use those \nvalues instead. Ranged from 10 to 500, depending on \ndesired granularity.\nDefault: 0\n"
  neighbor_a:
    type: float?
    doc: "The number of nearest neighbors, k, used in the graph-based \nclustering is computed as follows: k = neighbor_a + neighbor_b * \nlog10(n_cells). The actual number of neighbors used is the maximum \nof this value and graphclust_neighbors. Determines how clustering \ngranularity scales with cell count.\nDefault: -230.0\n"
  neighbor_b:
    type: float?
    doc: "The number of nearest neighbors, k, used in the graph-based \nclustering is computed as follows: k = neighbor_a + neighbor_b * \nlog10(n_cells). The actual number of neighbors used is the maximum \nof this value and graphclust_neighbors. Determines how clustering \ngranularity scales with cell count.\nDefault: 120.0\n"
  max_clusters:
    type: int?
    doc: "Compute K-means clustering using K values of 2 to N. \nSetting this too high may cause spurious clusters to be \ncalled. Ranges from 10 to 50, depending on the number \nof cell populations/clusters you expect to see.\nDefault: 10\n"
  tsne_input_pcs:
    type: int?
    doc: "Subset to top N principal components for TSNE. Change \nthis parameter if you want to see how the TSNE plot \nchanges when using fewer PCs, independent of the \nclustering/differential expression. You may find that \nTSNE is faster and/or the output looks better when using \nfewer PCs. Cannot be set higher than the num_principal_comps \nparameter.\nDefault: null\n"
  tsne_perplexity:
    type: int?
    doc: "TSNE perplexity parameter (see the TSNE FAQ for more details). \nWhen analyzing 100k+ cells, increasing this parameter may \nimprove TSNE results, but the algorithm will be slower. \nRanges from 30 to 50.\nDefault: 30\n"
  tsne_theta:
    type: float?
    doc: "TSNE theta parameter (see the TSNE FAQ for more details). \nHigher values yield faster, more approximate results (and \nvice versa). The runtime and memory performance of TSNE \nwill increase dramatically if you set this below 0.25. \nRanges from 0 to 1.\nDefault: 0.5\n"
  tsne_max_dims:
    type: int?
    doc: "Maximum number of TSNE output dimensions. Set this to 3 to \nproduce both 2D and 3D TSNE projections (note: runtime will \nincrease significantly). Ranges from 2 to 3.\nDefault: 2\n"
  tsne_max_iter:
    type: int?
    doc: "Number of total TSNE iterations. Try increasing this if \nTSNE results do not look good on larger numbers of cells. \nRuntime increases linearly with the number of iterations. \nRanges from 1000 to 10000.\nDefault: 1000\n"
  tsne_stop_lying_iter:
    type: int?
    doc: "Iteration at which TSNE learning rate is reduced. Try \nincreasing this if TSNE results do not look good on larger \nnumbers of cells. Cannot be set higher than tsne_max_iter. \nDefault: 250\n"
  tsne_mom_switch_iter:
    type: int?
    doc: "Iteration at which TSNE momentum is reduced. Try \nincreasing this if TSNE results do not look good on \nlarger numbers of cells. Cannot be set higher than \ntsne_max_iter.\nDefault: 250\n"
  umap_input_pcs:
    type: int?
    doc: "Subset to top N principal components for UMAP. Change \nthis parameter if you want to see how the UMAP plot \nchanges when using fewer PCs, independent of the \nclustering/differential expression. You may find that \nUMAP is faster and/or the output looks better when \nusing fewer PCs. Cannot be set higher than the \nnum_principal_comps parameter.\nDefault: null\n"
  umap_n_neighbors:
    type: int?
    doc: "Determines the number of neighboring points used in \nlocal approximations of manifold structure. Larger values \nwill usually result in more global structure at the loss \nof detailed local structure. Ranges from 5 to 50. \nDefault: 30\n"
  umap_max_dims:
    type: int?
    doc: "Maximum number of UMAP output dimensions. Set this to 3 \nto produce both 2D and 3D UMAP projections. Ranges from 2 \nto 3.\nDefault: 2\n"
  umap_min_dist:
    type: float?
    doc: "Controls how tightly the embedding is allowed to pack points \ntogether. Larger values make embedded points more evenly \ndistributed, while smaller values make the embedding more \naccurate with regard to the local structure. Ranges from \n0.001 to 0.5.\nDefault: 0.3\n"
  umap_metric:
    type:
    - 'null'
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
    doc: |
      Determines how the distance is computed in the input space.
      Default: "correlation"
  random_seed:
    type: int?
    doc: "Random seed. Due to the randomized nature of the algorithms, \nchanging this will produce slightly different results. If \nthe TSNE or UMAP results don't look good, try running \nmultiple times with different seeds and pick the TSNE or \nUMAP that looks best.\nDefault: 0\n"
outputs:
  secondary_analysis_report_folder:
    type: Directory
    outputBinding:
      glob: reanalyzed/outs/analysis
    doc: |
      Folder with secondary analysis results including
      dimensionality reduction, cell clustering, and
      differential expression for reanalyzed results.
  web_summary_report:
    type: File
    outputBinding:
      glob: reanalyzed/outs/web_summary.html
    doc: |
      Reanalyzed run summary metrics and charts
      in HTML format.
  filtered_feature_bc_matrix_folder:
    type: Directory
    outputBinding:
      glob: reanalyzed/outs/filtered_feature_bc_matrix
    doc: |
      Folder with filtered feature-barcode matrices
      containing only cellular barcodes in MEX format.
  reanalyze_params:
    type: File
    outputBinding:
      glob: reanalyzed/outs/params.csv
    doc: |
      Copy of the input params CSV file.
  loupe_browser_track:
    type: File
    outputBinding:
      glob: reanalyzed/outs/cloupe.cloupe
    doc: |
      Loupe Browser visualization and analysis
      file for reanalyzed results.
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- cellranger
- reanalyze
- --disable-ui
- --id
- reanalyzed
arguments:
- valueFrom: params.csv
  prefix: --params
  position: 15
stdout: cellranger_reanalyze_stdout.log
stderr: cellranger_reanalyze_stderr.log
label: Cellranger Reanalyze
doc: |
  Cell Ranger Reanalyze

  Reruns secondary analysis for Cell Ranger Count
  Gene Expression or Cell Ranger Multi experiments

  Rerunning the analysis for aggregated experiments
  is not currently supported.

  Parameters set by default:
  --disable-ui - no need in any UI when running in
                 Docker container
  --id         - hardcoded to `reanalyzed` as we want
                 to return the content of the
                 output folder as separate outputs

  Skipped outputs as they are identical to inputs:
  - Filtered feature-barcode matrices HDF5

  Not implemented parameters:
  --description           - not needed for now
  --agg                   - we don't support reruning
                            secondary analysis from
                            the aggregated samples
  --dry                   - not applicable to our use
                            case
  --jobmode               - we use default local mode
  --mempercore            - not used for local mode
  --maxjobs               - not used for local mode
  --jobinterval           - not used for local mode
  --overrides             - not needed for now
  --uiport                - we disabled UI
  --noexit                - we disabled UI
  --output-dir            - not needed for now
  --nopreflight           - no reason to skip preflight
                            checks
