cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_root = function(basename) {
          return basename.split('.').slice(0,1).join('.');
      };


'sd:upstream':
  epi_sample:
  - "chipseq-se.cwl"
  - "chipseq-pe.cwl"
  - "trim-chipseq-se.cwl"
  - "trim-chipseq-pe.cwl"
  - "trim-atacseq-se.cwl"
  - "trim-atacseq-pe.cwl"
  - "cutandrun-macs2-pe.cwl"
  - "cutandrun-seacr-pe.cwl"
  filtered_experiment:
  - "filter-peaks-for-heatmap.cwl"
  - "filter-deseq-for-heatmap.cwl"
  - "filter-diffbind-for-heatmap.cwl"
  - "filter-peaks-by-overlap.cwl"
  - "genelists-sets.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  regions_files:
    type: File[]
    format: "http://edamontology.org/format_3003"
    label: "Filtered Peaks or differential genes sample"
    doc: |
      "Regions of interest from a filtered epigenomic sample or filtered genes from a DESeq or DiffBind experiment. Formatted as a headerless BED file with [chrom start end name score strand] for gene list and [chrom start end name] for peak file."
    'sd:upstreamSource': "filtered_experiment/filtered_file"
    'sd:localLabel': true

  regions_names:
    type: string[]
    label: "regions sample aliases"
    doc: "Sample names for regions samples selected by user for regions_files. Order corresponds to the regions_files"
    'sd:upstreamSource': "filtered_experiment/alias"

  score_files:
    type: File[]
    format: "http://edamontology.org/format_2572"
    label: "Epigenomic sample(s)"
    doc: "bigWig file(s) containing the scores to be plotted. From ChIP/ATAC/C&R workflows."
    'sd:upstreamSource': "epi_sample/bigwig"
    'sd:localLabel': true

  score_names:
    type: string[]
    label: "Epigenomic sample(s)"
    doc: "Sample names for epigenomic samples selected by user for score_files. Order corresponds to the score_files"
    'sd:upstreamSource': "epi_sample/alias"

  subcommand:
    type:
    - "null"
    - type: enum
      symbols: ["reference-point", "scale-regions"]
    default: "reference-point"
    label: "Sets deeptools computeMatrix subcommand for processing the bed matrix."
    doc: "In reference-point mode, only those genomic positions before (upstream) and/or after (downstream) the center of each peak will be plotted. In scale-regions mode, all regions in the BED file are stretched or shrunken to the length (in bases) indicated by the user."
    'sd:localLabel': true

  beforeRegionStartLength:
    type: int?
    default: 3000
    label: "Distance upstream of the start site of the given regions"
    doc: "Default: 3000 bp"
    'sd:localLabel': true

  afterRegionStartLength:
    type: int?
    default: 3000
    label: "Distance downstream of the end site of the given regions"
    doc: "Default: 3000 bp"
    'sd:localLabel': true

  binSize:
    type: int?
    default: 10
    label: "Bin Size. Default: 10 bp"
    doc: "Length, in bases, of the non-overlapping bins for averaging the score over the regions length."
    'sd:layout':
      advanced: true

  regionBodyLength:
    type: int?
    default: 1000
    label: "Region Body Length. Default: 1000 bp"
    doc: "Only used in scale-regions mode. Distance between x and y (could be TSS and TES, or peak start and peak end, respectively), set to 0 for point centering of plot."
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 2
    label: "Number of Threads. Default: 2"
    doc: "Number of threads for steps that support multithreading."
    'sd:layout':
      advanced: true

  sortRegions:
    type:
    - "null"
    - type: enum
      symbols: ["descend", "ascend", "no"]
    default: "descend"
    label: "Row sorting type for heatmap. Default: descend"
    doc: "Whether the heatmap should present the regions sorted."
    'sd:layout':
      advanced: true

  sortUsing:
    type:
    - "null"
    - type: enum
      symbols: ["mean", "median", "max", "min", "sum", "region_length"]
    default: "mean"
    label: "Sorting Method. Default: mean"
    doc: "Indicate which method should be used for sorting. For each row the method is computed."
    'sd:layout':
      advanced: true

  colorMap:
    type:
    - "null"
    - type: enum
      symbols: ["RdBu", "Set1", "Set2", "Set3", "winter", "Dark2", "cool", "coolwarm", "rainbow"]
    default: "RdBu"
    label: "Color Map. Default: RdBu"
    doc: "Color map to use for the heatmap."
    'sd:layout':
      advanced: true

  kmeans:
    type: int?
    default: 0
    label: "Number of clusters to compute. Default: 0"
    doc: "Group rows by cluster instead of region set. When this option is set greater than 0, the matrix is split into clusters using the k-means algorithm."
    'sd:layout':
      advanced: true


outputs:

  matrix_file:
    type: File
    label: "gzipped matrix file"
    doc: "File for gzipped matrix needed by the plotHeatmap' and 'plotProfile' tools."
    outputSource: make_plots/matrix_file

  heatmap_file:
    type: File
    label: "Profile and heatmap plot"
    doc: "Profile and heatmap plot for scores over sets of genomic regions made by the 'plotHeatmap' tools."
    outputSource: make_plots/heatmap_file
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'Profile and heatmap plot'

  log_stdout_file:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stdout logfile"
    outputSource: make_plots/log_file_stdout
    'sd:visualPlugins':
    - markdownView:
        tab: 'Overview'

  log_stderr_file:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "stderr logfile"
    outputSource: make_plots/log_file_stderr     


steps:

  make_plots:
    run: ../tools/deeptools-plot.cwl
    in:
      regions_files: regions_files
      regions_names: regions_names
      score_files: score_files
      score_names: score_names
      beforeRegionStartLength: beforeRegionStartLength
      afterRegionStartLength: afterRegionStartLength
      regionBodyLength: regionBodyLength
      threads: threads
      sortRegions: sortRegions
      sortUsing: sortUsing
      colorMap: colorMap
      kmeans: kmeans
      subcommand: subcommand
      binSize: binSize
    out: [matrix_file, heatmap_file, log_file_stdout, log_file_stderr]
  

$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "deeptools - Tag enrichment heatmap and density profile for filtered regions"
label: "deeptools - Tag enrichment heatmap and density profile for filtered regions"
s:alternateName: "Generate tag enrichment heatmap and density profile plot around gene TSS or peak centers using deeptools"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/heatmap.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:support@datirium.com
        s:sameAs:
          - id: https://orcid.org/0000-0001-5872-259X


doc: |
  Generates tag density heatmap and histogram for the list of features in a headerless regions file.

  Inputs used are the bigWig file(s) of one or more ChIP/ATAC/C&R samples, and one or more filtered feature file(s) from the filtering and/or set operation workflows.

  The latter format contains `chrom start end name score strand`, only the first 3 columns are used in deeptools computeMatrix tool. The matrix is then used as input to plotHeatmap to generate the tag density plot and tag enrichment heatmap.



  computeMatrix paramters:
  --regionsFileName, -R
    File name, in BED format, containing the regions to plot. If multiple bed files are given, each one is considered a group that can be plotted separately. Also, adding a “#” symbol in the bed file causes all the regions until the previous “#” to be considered one group.
  --scoreFileName, -S
    bigWig file(s) containing the scores to be plotted. BigWig files can be obtained by using the bamCoverage or bamCompare tools. More information about the bigWig file format can be found at http://genome.ucsc.edu/goldenPath/help/bigWig.html
  --outFileName, -o
    File name to save the gzipped matrix file needed by the “plotHeatmap” and “plotProfile” tools.
  --beforeRegionStartLength=0, -b=0, --upstream=0
    Distance upstream of the start site of the regions defined in the region file. If the regions are genes, this would be the distance upstream of the transcription start site.
  --regionBodyLength=1000, -m=1000
    Distance in bases to which all regions will be fit.
  --afterRegionStartLength=0, -a=0, --downstream=0
    Distance downstream of the end site of the given regions. If the regions are genes, this would be the distance downstream of the transcription end site.
  --numberOfProcessors=max/2, -p=max/2
    Number of processors to use. Type “max/2” to use half the maximum number of processors or “max” to use all available processors.


  plotHeatmap parameters:
  --matrixFile, -m
    Matrix file from the computeMatrix tool.
  --outFileName, -out
    File name to save the image to. The file ending will be used to determine the image format. The available options are: “png”, “eps”, “pdf” and “svg”, e.g., MyHeatmap.png.
  --sortRegions=descend
    Whether the heatmap should present the regions sorted. The default is to sort in descending order based on the mean value per region.
    Possible choices: descend, ascend, no
  --sortUsing=mean
    Indicate which method should be used for sorting. For each row the method is computed.
    Possible choices: mean, median, max, min, sum, region_length
  --colorMap=RdYlBu
    Color map to use for the heatmap. Available values can be seen here: http://matplotlib.org/users/colormaps.html The available options are: ‘Spectral’, ‘summer’, ‘coolwarm’, ‘Set1’, ‘Set2’, ‘Set3’, ‘Dark2’, ‘hot’, ‘RdPu’, ‘YlGnBu’, ‘RdYlBu’, ‘gist_stern’, ‘cool’, ‘gray’, ‘GnBu’, ‘gist_ncar’, ‘gist_rainbow’, ‘CMRmap’, ‘bone’, ‘RdYlGn’, ‘spring’, ‘terrain’, ‘PuBu’, ‘spectral’, ‘gist_yarg’, ‘BuGn’, ‘bwr’, ‘cubehelix’, ‘YlOrRd’, ‘Greens’, ‘PRGn’, ‘gist_heat’, ‘Paired’, ‘hsv’, ‘Pastel2’, ‘Pastel1’, ‘BuPu’, ‘copper’, ‘OrRd’, ‘brg’, ‘gnuplot2’, ‘jet’, ‘gist_earth’, ‘Oranges’, ‘PiYG’, ‘YlGn’, ‘Accent’, ‘gist_gray’, ‘flag’, ‘BrBG’, ‘Reds’, ‘RdGy’, ‘PuRd’, ‘Blues’, ‘Greys’, ‘autumn’, ‘pink’, ‘binary’, ‘winter’, ‘gnuplot’, ‘RdBu’, ‘prism’, ‘YlOrBr’, ‘rainbow’, ‘seismic’, ‘Purples’, ‘ocean’, ‘PuOr’, ‘PuBuGn’, ‘nipy_spectral’, ‘afmhot’
  --kmeans	Number of clusters to compute. When this option is set, the matrix is split into clusters using the k-means algorithm. Only works for data that is not grouped, otherwise only the first group will be clustered. If more specific clustering methods are required, then save the underlying matrix and run the clustering using other software. The plotting of the clustering may fail with an error if a cluster has very few members compared to the total number or regions.
