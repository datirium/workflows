cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-deeptools:v1.0.0


inputs:

  script_command:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\n\nStdout log file for deeptools-plot.cwl tool:\n\n"
      # inputs
      regions_files=${0}
      regions_names=${1}
      score_files=${2}
      score_names=${3}
      beforeRegionStartLength=${4}
      afterRegionStartLength=${5}
      regionBodyLength=${6}
      threads=${7}
      sortRegions=${8}
      sortUsing=${9}
      colorMap=${10}
      kmeans=${11}
      subcommand=${12}
      binSize=${13}

      printf "INPUTS:\n"
      printf "\$0 - $regions_files\n"
      printf "\$1 - $regions_names\n"
      printf "\$2 - $score_files\n"
      printf "\$3 - $score_names\n"
      printf "\$4 - $beforeRegionStartLength\n"
      printf "\$5 - $afterRegionStartLength\n"
      printf "\$6 - $regionBodyLength\n"
      printf "\$7 - $threads\n"
      printf "\$8 - $sortRegions\n"
      printf "\$9 - $sortUsing\n"
      printf "\$10 - $colorMap\n"
      printf "\$11 - $kmeans\n"
      printf "\$12 - $subcommand\n"
      printf "\$13 - $binSize\n"

      # commands start
      # check that files and names for regions/scores are equal (in case there could be a comma in a name)
      count_rf=$(echo $regions_files | awk -F',' '{print NF}')
      count_rn=$(echo $regions_names | awk -F',' '{print NF}')
      if [[ $count_rf != $count_rn ]]; then echo "Error: number of regions files ($count_rf) does not match number of regions names ($count_rn). Check that there are no commas in the names strings. Exiting."; exit; fi
      count_sf=$(echo $score_files | awk -F',' '{print NF}')
      count_sn=$(echo $score_names | awk -F',' '{print NF}')
      if [[ $count_sf != $count_sn ]]; then echo "Error: number of score files ($count_sf) does not match number of score names ($count_sn). Check that there are no commas in the names strings. Exiting."; exit; fi

      # pair regions files and names
      awk -F',' '{if(NR==FNR){for(i=1;i<=NF;i++){x[i]=$i}}else{for(j in x){printf("%s\t%s\n",x[j],$j)}}}' <(echo $regions_files) <(echo $regions_names) > regions-files-names.pair
      #   replace single spaces with underscore
      sed -i 's/ /_/g' regions-files-names.pair
      printf "\nContents of 'regions-files-names.pair' file:\n"
      cat regions-files-names.pair
      printf "\n\n"
      # rename regions files using names (these correspond to the row section labels in the deeptools heatmap)
      while read pair; do cp $(echo $pair | awk -F' ' '{print $1}') $(echo $pair | awk -F' ' '{print $2}')".bed"; done < regions-files-names.pair

      # pair scores files and names
      awk -F',' '{if(NR==FNR){for(i=1;i<=NF;i++){x[i]=$i}}else{for(j in x){printf("%s\t%s\n",x[j],$j)}}}' <(echo $score_files) <(echo $score_names) > score-files-names.pair
      #   replace single spaces with underscore
      sed -i 's/ /_/g' score-files-names.pair
      printf "\nContents of 'score-files-names.pair' file:\n"
      cat score-files-names.pair
      printf "\n\n"
      # rename scores files using names (these correspond to the column section labels [sample names basically] in the deeptools heatmap)
      while read pair; do cp $(echo $pair | awk -F' ' '{print $1}') $(echo $pair | awk -F' ' '{print $2}')".bigWig"; done < score-files-names.pair



      # run deeptools compute matrix
      if [[ "$subcommand" == "reference-point" ]]; then
        printf "Running deeptools 'computeMatrix reference-point' command...\n\n"
        computeMatrix reference-point --referencePoint center \
                  -S $(find ./ -maxdepth 1 -mindepth 1 -name "*.bigWig" | sed $'$!N;s/\\\n/\t/')  \
                  -R $(find ./ -maxdepth 1 -mindepth 1 -name "*.bed" | sed $'$!N;s/\\\n/\t/') \
                  --beforeRegionStartLength $beforeRegionStartLength \
                  --afterRegionStartLength $afterRegionStartLength \
                  --binSize $binSize \
                  --skipZeros -o matrix.mat.gz \
                  --numberOfProcessors $threads
      elif [[ "$subcommand" == "scale-regions" ]]; then
        printf "Running deeptools 'computeMatrix scale-regions' command...\n\n"
        computeMatrix scale-regions \
                  -S $(find ./ -maxdepth 1 -mindepth 1 -name "*.bigWig" | sed $'$!N;s/\\\n/\t/')  \
                  -R $(find ./ -maxdepth 1 -mindepth 1 -name "*.bed" | sed $'$!N;s/\\\n/\t/') \
                  --beforeRegionStartLength $beforeRegionStartLength \
                  --regionBodyLength $regionBodyLength \
                  --afterRegionStartLength $afterRegionStartLength \
                  --binSize $binSize \
                  --skipZeros -o matrix.mat.gz \
                  --numberOfProcessors $threads
      fi

      # make plot
      #   set plot height based on 4 + number of lists * 5)
      plotheight=$(awk -F'\t' 'END{print(4+(NR*10))}' regions-files-names.pair)
      #   set plot width based on 1 + (number of samples * 3)
      plotwidth=$(awk -F'\t' 'END{print(1+(NR*3))}' score-files-names.pair)
      if [[ $kmeans -gt 0 ]]; then
        printf "Running deeptools 'plotHeatmap' command with $kmeans kmeans cluster(s), plot height = $plotheight, and width = $plotwidth . . .\n\n"
        plotHeatmap -m matrix.mat.gz \
            -out heatmap.svg \
            --sortRegions $sortRegions \
            --sortUsing $sortUsing \
            --colorMap $colorMap \
            --kmeans $kmeans \
            --legendLocation "upper-right" \
            --heatmapHeight $plotheight \
            --heatmapWidth $plotwidth \
            --plotFileFormat "svg"
      else
        printf "Running deeptools 'plotHeatmap' command without kmeans clustering, plot height = $plotheight, and width = $plotwidth . . .\n\n"
        plotHeatmap -m matrix.mat.gz \
            -out heatmap.svg \
            --sortRegions $sortRegions \
            --sortUsing $sortUsing \
            --colorMap $colorMap \
            --legendLocation "upper-right" \
            --heatmapHeight $plotheight \
            --heatmapWidth $plotwidth \
            --plotFileFormat "svg"
      fi
      printf "\nscript complete\n"
    inputBinding:
      position: 1

  regions_files:
    type: File[]
    inputBinding:
      position: 2
      itemSeparator: ","
    doc: |
      Regions of interest from a filtered epigenomic sample or filtered genes from a DESeq or
      DiffBind experiment. Formatted as a headerless BED file with [chrom start end name score
      strand] for gene list and [chrom start end name] for peak file.

  regions_names:
    type: string[]
    inputBinding:
      position: 3
      itemSeparator: ","
    doc: |
      Sample names for regions samples selected by user for regions_files. Order corresponds
      to the regions_files

  score_files:
    type: File[]
    inputBinding:
      position: 4
      itemSeparator: ","
    doc: |
      bigWig file(s) containing the scores to be plotted. From ChIP/ATAC/C&R workflows.

  score_names:
    type: string[]
    inputBinding:
      position: 5
      itemSeparator: ","
    doc: |
      Sample names for epigenomic samples selected by user for score_files. Order
      corresponds to the score_files

  beforeRegionStartLength:
    type: int
    inputBinding:
      position: 10
    doc: |
      Distance upstream of the start site of the given regions

  afterRegionStartLength:
    type: int
    inputBinding:
      position: 11
    doc: |
      Distance upstream of the start site of the given regions

  regionBodyLength:
    type: int
    inputBinding:
      position: 12
    doc: |
      Distance upstream of the start site of the given regions

  threads:
    type: int
    inputBinding:
      position: 13
    doc: |
      Number of threads for steps that support multithreading

  sortRegions:
    type:
    - "null"
    - type: enum
      symbols: ["descend", "ascend", "no"]
    inputBinding:
      position: 14
    doc: |
      Whether the heatmap should present the regions sorted. Default: descend

  sortUsing:
    type:
    - "null"
    - type: enum
      symbols: ["mean", "median", "max", "min", "sum", "region_length"]
    inputBinding:
      position: 15
    doc: |
      Indicate which method should be used for sorting. For each row the method is
      computed. Default: mean

  colorMap:
    type:
    - "null"
    - type: enum
      symbols: ["RdBu", "Set1", "Set2", "Set3", "winter", "Dark2", "cool", "coolwarm", "rainbow"]
    inputBinding:
      position: 16
    doc: |
      Color map to use for the heatmap

  kmeans:
    type: int
    inputBinding:
      position: 17
    doc: |
      Group rows by cluster instead of region set. When this option is set greater than 0, the
      matrix is split into clusters using the k-means algorithm

  subcommand:
    type:
    - "null"
    - type: enum
      symbols: ["reference-point", "scale-regions"]
    inputBinding:
      position: 18
    doc: |
      Sets deeptools computeMatrix subcommand for processing the bed matrix.
      In scale-regions mode, all regions in the BED file are stretched or shrunken to the length (in bases) indicated by the user.
      In reference-point mode, only those genomic positions before (upstream) and/or after (downstream) the center of each peak will be plotted.

  binSize:
    type: int
    inputBinding:
      position: 19
    doc: |
      Length, in bases, of the non-overlapping bins for averaging the score over the regions length.


outputs:

  matrix_file:
    type: File
    outputBinding:
      glob: matrix.mat.gz
    doc: |
      Matrix file from the computeMatrix tool

  heatmap_file:
    type: File
    outputBinding:
      glob: heatmap.svg
    doc: |
      Profile and heatmap plot for scores over sets of genomic regions made by the 'plotHeatmap' tools

  log_file_stdout:
    type: stdout

  log_file_stderr:
    type: stderr


baseCommand: ["bash", "-c"]
stdout: deeptools-plot_stdout.log
stderr: deeptools-plot_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "deeptools tag enrichment heatmap and density plot"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/deeptools-plot.cwl
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
  This tool takes as input multiple samples bigwig files from ChIP/ATAC/C&R samples, and peak/gene list
  TSV files from the filtering or set operations workflows and performs deeptools 'computeMatrix' and
  'plotHeatmap' with user defined parameters. Outputs include the computed scores matrix and heatmap png
  files.
      

  computeMatrix paramters:

  Sub-commands:
  scale-regions
  In the scale-regions mode, all regions in the BED file are stretched or shrunken to the length (in bases) indicated by the user.

  reference-point
  Reference-point refers to a position within a BED region (e.g., the starting point). In this mode, only those genomicpositions before (upstream) and/or after (downstream) of the reference point will be plotted.

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
  --binSize=10, -bs=10
    Length, in bases, of the non-overlapping bins for averaging the score over the regions length.


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
