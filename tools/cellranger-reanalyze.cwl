cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: cumulusprod/cellranger:7.0.0

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
      prefix: "--matrix"
    doc: |
      A feature-barcode matrix containing data for one genome.
      Should be the filtered version, unless using --force-cells

  selected_barcodes:
    type: File?
    inputBinding:
      position: 6
      prefix: "--barcodes"
    doc: |
      A CSV file containing a list of cell barcodes to use for reanalysis,
      e.g. barcodes exported from Loupe Browser. All barcodes must be present
      in the matrix.

  selected_genes:
    type: File?
    inputBinding:
      position: 7
      prefix: "--genes"
    doc: |
      A CSV file containing a list of gene IDs to use for reanalysis (corresponding
      to the gene_id field of the reference GTF). All gene IDs must be present in
      the matrix. Note that only gene features are used in secondary analysis.

  excluded_genes:
    type: File?
    inputBinding:
      position: 8
      prefix: "--exclude-genes"
    doc: |
      A CSV file containing a list of gene IDs to exclude for reanalysis (corresponding
      to the gene_id field of the reference GTF). All gene IDs must be present in
      the matrix. The exclusion is applied after setting the gene list with --genes.
      Note that only gene features are used in secondary analysis.

  force_cells:
    type: int?
    inputBinding:
      position: 9
      prefix: "--force-cells"
    doc: |
      Force pipeline to use this number of cells, bypassing the cell detection algorithm.
      Use this if the number of cells estimated by Cell Ranger is not consistent with the
      barcode rank plot. If specifying a value that exceeds the original cell count, you
      must use the raw_gene_bc_matrices_h5.h5

  threads:
    type: int?
    inputBinding:
      position: 10
      prefix: "--localcores"
    doc: |
      Set max cores the pipeline may request at one time.
      Default: all available

  memory_limit:
    type: int?
    inputBinding:
      position: 11
      prefix: "--localmem"
    doc: |
      Set max GB the pipeline may request at one time
      Default: all available

  virt_memory_limit:
    type: int?
    inputBinding:
      position: 12
      prefix: "--localvmem"
    doc: |
      Set max virtual address space in GB for the pipeline
      Default: all available

  num_analysis_bcs:
    type: int?
    doc: |
      Randomly subset data to N barcodes for all analysis. Reduce this parameter if you
      want to improve performance or simulate results from lower cell counts. Cannot be
      set higher than the available number of cells.
      Default: null

  num_pca_bcs:
    type: int?
    doc: |
      Randomly subset data to N barcodes when computing PCA projection (the most memory-intensive
      step). The PCA projection will still be applied to the full dataset, i.e. your final results
      will still reflect all the data. Try reducing this parameter if your analysis is running out
      of memory. Cannot be set higher than the available number of cells.
      Default: null

  num_pca_genes:
    type: int?
    doc: |
      Subset data to the top N genes (ranked by normalized dispersion) when computing PCA.
      Differential expression will still reflect all genes. Try reducing this parameter if
      your analysis is running out of memory. Cannot be set higher than the number of genes
      in the reference transcriptome.
      Default: null

  num_principal_comps:
    type: int?
    doc: |
      Compute N principal components for PCA. Setting this too high may cause spurious clusters
      to be called. The default value is 100 when the chemistry batch correction is enabled.
      Set from 10 to 100, depending on the number of cell populations/clusters you expect to see.
      Default: 10

  cbc_knn:
    type: int?
    doc: |
      Specify the number of nearest neighbors used to identify mutual nearest neighbors.
      Setting this too high will increase runtime and may cause out of memory error.
      See Chemistry Batch Correction page for more details. Ranges from 5 to 20.
      Default: 10

  cbc_alpha:
    type: float?
    doc: |
      Specify the threshold of the percentage of matched cells between two batches,
      which is used to determine if the batch pair will be merged. See Chemistry
      Batch Correction page for more details. Ranges from 0.05 to 0.5.
      Default: 0.1
      
  cbc_sigma:
    type: float?
    doc: |
      Specify the bandwidth of the Gaussian smoothing kernel used to compute the correction
      vector for each cell. See Chemistry Batch Correction page for more details. Ranges
      from 10 to 500.
      Default: 150

  cbc_realign_panorama:
    type: boolean?
    doc: |
      Specify if two batches will be merged if they are already in the same panorama. Setting
      this to True will usually improve the performance, but will also increase runtime and
      memory usage. See Chemistry Batch Correction page for more details. One of true or false.
      Default: false

  graphclust_neighbors:
    type: int?
    doc: |
      Number of nearest-neighbors to use in the graph-based clustering. Lower values result in
      higher-granularity clustering. The actual number of neighbors used is the maximum of this
      value and that determined by neighbor_a and neighbor_b. Set this value to zero to use those
      values instead. Ranged from 10 to 500, depending on desired granularity.
      Default: 0

  neighbor_a:
    type: float?
    doc: |
      The number of nearest neighbors, k, used in the graph-based clustering is computed as follows:
      k = neighbor_a + neighbor_b * log10(n_cells). The actual number of neighbors used is the maximum
      of this value and graphclust_neighbors. Determines how clustering granularity scales with cell count.
      Default: -230.0

  neighbor_b:
    type: float?
    doc: |
      The number of nearest neighbors, k, used in the graph-based clustering is computed as follows:
      k = neighbor_a + neighbor_b * log10(n_cells). The actual number of neighbors used is the maximum of
      this value and graphclust_neighbors. Determines how clustering granularity scales with cell count.
      Default: 120.0

  max_clusters:
    type: int?
    doc: |
      Compute K-means clustering using K values of 2 to N. Setting this too high may cause spurious clusters
      to be called. Ranges from 10 to 50, depending on the number of cell populations / clusters you expect to see.
      Default: 10

  tsne_input_pcs:
    type: int?
    doc: |
      Subset to top N principal components for TSNE. Change this parameter if you want to see how the TSNE plot
      changes when using fewer PCs, independent of the clustering / differential expression. You may find that TSNE
      is faster and/or the output looks better when using fewer PCs. Cannot be set higher than
      the num_principal_comps parameter.
      Default: null

  tsne_perplexity:
    type: int?
    doc: |
      TSNE perplexity parameter (see the TSNE FAQ for more details). When analyzing 100k+ cells, increasing this
      parameter may improve TSNE results, but the algorithm will be slower. Ranges from 30 to 50.
      Default: 30

  tsne_theta:
    type: float?
    doc: |
      TSNE theta parameter (see the TSNE FAQ for more details). Higher values yield faster, more approximate results
      (and vice versa). The runtime and memory performance of TSNE will increase dramatically if you set this below 0.25.
      Ranges from 0 to 1.
      Default: 0.5

  tsne_max_dims:
    type: int?
    doc: |
      Maximum number of TSNE output dimensions. Set this to 3 to produce both 2D and 3D TSNE projections
      (note: runtime will increase significantly). Ranges from 2 to 3.
      Default: 2

  tsne_max_iter:
    type: int?
    doc: |
      Number of total TSNE iterations. Try increasing this if TSNE results do not look good on larger numbers
      of cells. Runtime increases linearly with number of iterations. Ranges from 1000 to 10000.
      Default: 1000

  tsne_stop_lying_iter:
    type: int?
    doc: |
      Iteration at which TSNE learning rate is reduced. Try increasing this if TSNE results do not look good
      on larger numbers of cells. Cannot be set higher than tsne_max_iter.
      Default: 250

  tsne_mom_switch_iter:
    type: int?
    doc: |
      Iteration at which TSNE momentum is reduced. Try increasing this if TSNE results do not look good on
      larger numbers of cells. Cannot be set higher than tsne_max_iter. Cannot be set higher than tsne_max_iter.
      Default: 250

  umap_input_pcs:
    type: int?
    doc: |
      Subset to top N principal components for UMAP. Change this parameter if you want to see how the UMAP plot
      changes when using fewer PCs, independent of the clustering / differential expression. You may find that
      UMAP is faster and/or the output looks better when using fewer PCs. Cannot be set higher than the
      num_principal_comps parameter.
      Default: null

  umap_n_neighbors:
    type: int?
    doc: |
      Determines the number of neighboring points used in local approximations of manifold structure.
      Larger values will usually result in more global structure at the loss of detailed local structure.
      Ranges from 5 to 50.
      Default: 30

  umap_max_dims:
    type: int?
    doc: |
      Maximum number of UMAP output dimensions. Set this to 3 to produce both 2D and 3D UMAP projections.
      Ranges from 2 to 3.
      Default: 2

  umap_min_dist:
    type: float?
    doc: |
      Controls how tightly the embedding is allowed to pack points together. Larger values make embedded
      points are more evenly distributed, while smaller values make the embedding more accurately with
      regard to the local structure. Ranges from 0.001 to 0.5.
      Default: 0.3

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
    doc: |
      Determines how the distance is computed in the input space.
      Default: "correlation"

  random_seed:
    type: int?
    doc: |
      Random seed. Due to the randomized nature of the algorithms, changing this will produce slightly
      different results. If the TSNE or UMAP results don't look good, try running multiple times with
      different seeds and pick the TSNE or UMAP that looks best.
      Default: 0


outputs:

  secondary_analysis_report_folder:
    type: Directory
    outputBinding:
      glob: "reanalyzed/outs/analysis"
    doc: |
      Folder with secondary analysis results including dimensionality reduction,
      cell clustering, and differential expression for reanalyzed results

  web_summary_report:
    type: File
    outputBinding:
      glob: "reanalyzed/outs/web_summary.html"
    doc: |
      Reanalyzed run summary metrics and charts in HTML format

  filtered_feature_bc_matrix_folder:
    type: Directory
    outputBinding:
      glob: "reanalyzed/outs/filtered_feature_bc_matrix"
    doc: |
      Folder with filtered feature-barcode matrices containing only cellular
      barcodes in MEX format.

  reanalyze_params:
    type: File
    outputBinding:
      glob: "reanalyzed/outs/params.csv"
    doc: |
      Copy of the input params CSV file

  loupe_browser_track:
    type: File
    outputBinding:
      glob: "reanalyzed/outs/cloupe.cloupe"
    doc: |
      Loupe Browser visualization and analysis file for reanalyzed results

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["cellranger", "reanalyze", "--disable-ui",  "--id", "reanalyzed"]
arguments:
- valueFrom: "params.csv"
  prefix: "--params"
  position: 15


stdout: cellranger_reanalyze_stdout.log
stderr: cellranger_reanalyze_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Cell Ranger Reanalyze"
s:name: "Cell Ranger Reanalyze"
s:alternateName: |
  Reruns secondary analysis performed on the GEX feature-barcode matrix (dimensionality reduction,
  clustering and visualization) using different parameter settings

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/cellranger-reanalyze.cwl
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
  Cell Ranger Reanalyze

  Runs cellranger reanalyze command to rerun secondary analysis performed on the
  GEX feature-barcode matrix (dimensionality reduction, clustering and visualization)
  using different parameter settings.

  Rerunning the analysis for aggregated experiments is not currently supported.

  Parameters set by default:
  --disable-ui - no need in any UI when running in Docker container
  --id         - hardcoded to `reanalyzed` as we want to return the content of the
                 output folder as separate outputs

  Skipped outputs as they are identical to inputs:
  - Filtered feature-barcode matrices HDF5

  Not implemented parameters:
  --description           - not needed for now
  --agg                   - we don't support reruning secondary analysis from aggregated samples
  --dry                   - not applicable to our use case
  --jobmode               - we use default local mode
  --mempercore            - not used for local mode
  --maxjobs               - not used for local mode
  --jobinterval           - not used for local mode
  --overrides             - not needed for now
  --uiport                - we disabled UI
  --noexit                - we disabled UI
  --nopreflight           - no reason to skip preflight checks


s:about: |
  Re-run secondary analysis (dimensionality reduction, clustering, etc)

  USAGE:
      cellranger reanalyze [OPTIONS] --id <ID> --matrix <MATRIX_H5>

  OPTIONS:
      --id <ID>                      A unique run id and output folder name [a-zA-Z0-9_-]+
      --description <TEXT>           Sample description to embed in output files [default: ]
      --matrix <MATRIX_H5>           A feature-barcode matrix containing data for one genome. Should be the filtered version, unless using --force-cells
      --params <PARAMS_CSV>          A CSV file specifying analysis parameters. Optional
      --barcodes <BARCODES_CSV>      A CSV file containing a list of cell barcodes to use for reanalysis, e.g. barcodes exported from Loupe Browser. Optional
      --genes <GENES_CSV>            A CSV file containing a list of feature IDs to use for reanalysis. For gene expression, this should correspond to the gene_id field in the
                                    reference GTF should be \(e.g. ENSG... for ENSEMBL-based references\). Optional
      --exclude-genes <GENES_CSV>    A CSV file containing a list of feature IDs to exclude from reanalysis. For gene expression, this should correspond to the gene_id field in
                                    the reference GTF \(e.g., ENSG... for ENSEMBL-based references\). The exclusion is applied after --genes. Optional
      --agg <AGGREGATION_CSV>        If the input matrix was produced by 'aggr', you may pass the same aggregation CSV in order to retain per-library tag information in the
                                    resulting .cloupe file.  This argument is required to enable chemistry batch correction. Optional
      --force-cells <NUM>            Force pipeline to use this number of cells, bypassing cell calling algorithm. [MINIMUM: 10]
      --dry                          Do not execute the pipeline. Generate a pipeline invocation (.mro) file and stop
      --jobmode <MODE>               Job manager to use. Valid options: local (default), sge, lsf, slurm or path to a .template file. Search for help on "Cluster Mode" at
                                    support.10xgenomics.com for more details on configuring the pipeline to use a compute cluster [default: local]
      --localcores <NUM>             Set max cores the pipeline may request at one time. Only applies to local jobs
      --localmem <NUM>               Set max GB the pipeline may request at one time. Only applies to local jobs
      --localvmem <NUM>              Set max virtual address space in GB for the pipeline. Only applies to local jobs
      --mempercore <NUM>             Reserve enough threads for each job to ensure enough memory will be available, assuming each core on your cluster has at least this much
                                    memory available. Only applies to cluster jobmodes
      --maxjobs <NUM>                Set max jobs submitted to cluster at one time. Only applies to cluster jobmodes
      --jobinterval <NUM>            Set delay between submitting jobs to cluster, in ms. Only applies to cluster jobmodes
      --overrides <PATH>             The path to a JSON file that specifies stage-level overrides for cores and memory. Finer-grained than --localcores, --mempercore and
                                    --localmem. Consult https://support.10xgenomics.com/ for an example override file
      --uiport <PORT>                Serve web UI at http://localhost:PORT
      --disable-ui                   Do not serve the web UI
      --noexit                       Keep web UI running after pipestance completes or fails
      --nopreflight                  Skip preflight checks
      -h, --help                         Print help information