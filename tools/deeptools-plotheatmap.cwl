cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/deeptools:v0.0.1


inputs:

  scores_matrix:
    type: File
    inputBinding:
      position: 5
      prefix: "--matrixFile"
    doc: |
      Matrix file from the computeMatrix tool

  output_filename:
    type: string
    inputBinding:
      position: 6
      prefix: "--outFileName"
    doc: |
      File name to save the image to. The file ending will be used to determine the image format.
      The available options are: “png”, “eps”, “pdf” and “svg”, e.g., MyHeatmap.png
  
  interpolation_method:
    type:
    - "null"
    - type: enum
      name: "interpolation_method"
      symbols: ["auto", "nearest", "bilinear", "bicubic", "gaussian"]
    inputBinding:
      position: 7
      prefix: "--interpolationMethod"
    doc: |
      If the heatmap image contains a large number of columns is usually better to use an
      interpolation method to produce better results. By default, plotHeatmap uses the method
      nearest if the number of columns is 1000 or less. Otherwise it uses the bilinear method.
      This default behaviour can be changed by using any of the following options: “nearest”,
      “bilinear”, “bicubic”, “gaussian”

  dpi:
    type: int?
    inputBinding:
      position: 8
      prefix: "--dpi"
    doc: "Set the DPI to save the figure."

  plot_type:
    type:
    - "null"
    - type: enum
      name: "plot_type"
      symbols: ["lines", "fill", "se", "std"]
    inputBinding:
      position: 9
      prefix: "--plotType"
    doc: |
      “lines” will plot the profile line based on the average type selected.
      “fill” fills the region between zero and the profile curve. The fill in
      color is semi transparent to distinguish different profiles.
      “se” and “std” color the region between the profile and the standard error
      or standard deviation of the data.

  sort_regions:
    type:
    - "null"
    - type: enum
      name: "sort_regions"
      symbols: ["descend", "ascend", "no", "keep"]
    inputBinding:
      position: 10
      prefix: "--sortRegions"
    doc: |
      Whether the heatmap should present the regions sorted. The default is to sort in
      descending order based on the mean value per region. Note that “keep” and “no” are
      the same thing.

  sort_using:
    type:
    - "null"
    - type: enum
      name: "sort_using"
      symbols: ["mean", "median", "max", "min", "sum", "region_length"]
    inputBinding:
      position: 11
      prefix: "--sortUsing"
    doc: |
      Indicates which method should be used for sorting. For each row the method is computed.
      For region_length, a dashed line is drawn at the end of the region (reference point TSS
      and center) or the beginning of the region (reference point TES) as appropriate.

  average_type_summary_plot:
    type:
    - "null"
    - type: enum
      name: "average_type_summary_plot"
      symbols: ["mean", "median", "min", "max", "std", "sum"]
    inputBinding:
      position: 12
      prefix: "--averageTypeSummaryPlot"
    doc: |
      Define the type of statistic that should be plotted in the summary image above the heatmap.
      The options are: “mean”, “median”, “min”, “max”, “sum” and “std”.

  what_to_show:
    type:
    - "null"
    - type: enum
      name: "what_to_show"
      symbols:
      - plot, heatmap and colorbar
      - plot and heatmap
      - heatmap only
      - heatmap and colorbar
    inputBinding:
      position: 13
      prefix: "--whatToShow"
    doc: |
      The default is to include a summary or profile plot on top of the heatmap and a heatmap colorbar.
      Other options are: “plot and heatmap”, “heatmap only”, “heatmap and colorbar”, and the default “plot,
      heatmap and colorbar”.

  x_axis_label:
    type: string?
    inputBinding:
      position: 14
      prefix: "--xAxisLabel"
    doc: |
      Description for the x-axis label

  start_label:
    type: string?
    inputBinding:
      position: 15
      prefix: "--startLabel"
    doc: |
      [only for scale-regions mode] Label shown in the plot for the start of the region.
      Default is TSS (transcription start site), but could be changed to anything, e.g.
      “peak start”. Same for the –endLabel option.

  end_label:
    type: string?
    inputBinding:
      position: 16
      prefix: "--endLabel"
    doc: |
      [only for scale-regions mode] Label shown in the plot for the region end.
      Default is TES (transcription end site).

  ref_point_label:
    type: string?
    inputBinding:
      position: 17
      prefix: "--refPointLabel"
    doc: |
      [only for reference-point mode] Label shown in the plot for the reference-point.
      Default is the same as the reference point selected (e.g. TSS), but could be anything,
      e.g. “peak start”.  

  label_rotation_angle:
    type: int?
    inputBinding:
      position: 18
      prefix: "--labelRotation"
    doc: |
      Rotation of the X-axis labels in degrees.
      The default is 0, positive values denote a counter-clockwise rotation.

  regions_label:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      position: 19
      prefix: "--regionsLabel"
    doc: |
      Labels for the regions plotted in the heatmap. If more than one region is being plotted, a list of
      labels separated by spaces is required. If a label itself contains a space, then quotes are needed.
      For example, –regionsLabel label_1, “label 2”.

  samples_label:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      position: 20
      prefix: "--samplesLabel"
    doc: |
      Labels for the samples plotted. The default is to use the file name of the sample. The sample labels
      should be separated by spaces and quoted if a label itselfcontains a space
      E.g. –samplesLabel label-1 “label 2”

  plot_title:
    type: string?
    inputBinding:
      position: 21
      prefix: "--plotTitle"
    doc: |
      Title of the plot, to be printed on top of the generated image.
      Leave blank for no title.

  y_axisLabel:
    type: string?
    inputBinding:
      position: 22
      prefix: "--yAxisLabel"
    doc: |
      Y-axis label for the top panel.

  y_min:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      position: 23
      prefix: "--yMin"
    doc: |
      Minimum value for the Y-axis. Multiple values, separated by spaces can be set for each profile.
      If the number of yMin values is smaller thanthe number of plots, the values are recycled.

  y_max:
    type:
    - "null"
    - int
    - int[]
    inputBinding:
      position: 24
      prefix: "--yMax"
    doc: |
      Maximum value for the Y-axis. Multiple values, separated by spaces can be set for each profile.
      If the number of yMin values is smaller thanthe number of plots, the values are recycled.

  legend_location:
    type:
    - "null"
    - type: enum
      name: "legend_location"
      symbols:
      - best
      - upper-right
      - upper-left
      - upper-center
      - lower-left
      - lower-right
      - lower-center
      - center
      - center-left
      - center-right
      - none
    inputBinding:
      position: 25
      prefix: "--legendLocation"
    doc: |
      Location for the legend in the summary plot. Note that “none” does not work for the profiler.

  per_group:
    type: boolean?
    inputBinding:
      position: 26
      prefix: "--perGroup"
    doc: |
      The default is to plot all groups of regions by sample. Using this option instead plots all
      samples by group of regions. Note that this is only useful if you have multiple groups of
      regions by sample rather than group.

  plot_file_format:
    type:
    - "null"
    - type: enum
      name: "plot_file_format"
      symbols: ["png", "pdf", "svg", "eps", "plotly"]
    inputBinding:
      position: 27
      prefix: "--plotFileFormat"
    doc: |
      Image format type. If given, this option overrides the image format based on the plotFile ending.
      The available options are: “png”, “eps”, “pdf”, “plotly” and “svg”


outputs:

  heatmap_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
    doc: "Heatmap file"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["plotHeatmap", "--verbose"]


stdout: plot_heatmap_stdout.log
stderr: plot_heatmap_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/deeptools-metadata.yaml

label: "plotHeatmap - tool creates a heatmap for scores associated with genomic regions"
s:name: "plotHeatmap - tool creates a heatmap for scores associated with genomic regions"
s:alternateName: "plotHeatmap - tool creates a heatmap for scores associated with genomic regions"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/deeptools-plotheatmap.cwl
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
  This tool creates a heatmap for scores associated with genomic regions.
  The program requires a matrix file generated by the tool computeMatrix.


s:about: |
  usage: plotHeatmap [--matrixFile MATRIXFILE] --outFileName OUTFILENAME
                    [--outFileSortedRegions FILE] [--outFileNameMatrix FILE]
                    [--interpolationMethod STR] [--dpi DPI] [--kmeans KMEANS]
                    [--hclust HCLUST] [--silhouette] [--help] [--version]
                    [--plotType {lines,fill,se,std}]
                    [--sortRegions {descend,ascend,no,keep}]
                    [--sortUsing {mean,median,max,min,sum,region_length}]
                    [--sortUsingSamples SORTUSINGSAMPLES [SORTUSINGSAMPLES ...]]
                    [--linesAtTickMarks]
                    [--clusterUsingSamples CLUSTERUSINGSAMPLES [CLUSTERUSINGSAMPLES ...]]
                    [--averageTypeSummaryPlot {mean,median,min,max,std,sum}]
                    [--missingDataColor MISSINGDATACOLOR]
                    [--colorMap COLORMAP [COLORMAP ...]] [--alpha ALPHA]
                    [--colorList COLORLIST [COLORLIST ...]]
                    [--colorNumber COLORNUMBER] [--zMin ZMIN [ZMIN ...]]
                    [--zMax ZMAX [ZMAX ...]] [--heatmapHeight HEATMAPHEIGHT]
                    [--heatmapWidth HEATMAPWIDTH]
                    [--whatToShow {plot, heatmap and colorbar,plot and heatmap,heatmap only,heatmap and colorbar}]
                    [--boxAroundHeatmaps BOXAROUNDHEATMAPS]
                    [--xAxisLabel XAXISLABEL] [--startLabel STARTLABEL]
                    [--endLabel ENDLABEL] [--refPointLabel REFPOINTLABEL]
                    [--labelRotation LABEL_ROTATION]
                    [--regionsLabel REGIONSLABEL [REGIONSLABEL ...]]
                    [--samplesLabel SAMPLESLABEL [SAMPLESLABEL ...]]
                    [--plotTitle PLOTTITLE] [--yAxisLabel YAXISLABEL]
                    [--yMin YMIN [YMIN ...]] [--yMax YMAX [YMAX ...]]
                    [--legendLocation {best,upper-right,upper-left,upper-center,lower-left,lower-right,lower-center,center,center-left,center-right,none}]
                    [--perGroup] [--plotFileFormat] [--verbose]

  This tool creates a heatmap for scores associated with genomic regions. The
  program requires a matrix file generated by the tool ``computeMatrix``.

  Required arguments:
    --matrixFile MATRIXFILE, -m MATRIXFILE
                          Matrix file from the computeMatrix tool. (default:
                          None)
    --outFileName OUTFILENAME, -out OUTFILENAME, -o OUTFILENAME
                          File name to save the image to. The file ending will
                          be used to determine the image format. The available
                          options are: "png", "eps", "pdf" and "svg", e.g.,
                          MyHeatmap.png. (default: None)

  Output options:
    --outFileSortedRegions FILE
                          File name into which the regions are saved after
                          skipping zeros or min/max threshold values. The order
                          of the regions in the file follows the sorting order
                          selected. This is useful, for example, to generate
                          other heatmaps while keeping the sorting of the first
                          heatmap. Example: Heatmap1sortedRegions.bed (default:
                          None)
    --outFileNameMatrix FILE
                          If this option is given, then the matrix of values
                          underlying the heatmap will be saved using this name,
                          e.g. MyMatrix.tab. (default: None)
    --interpolationMethod STR
                          If the heatmap image contains a large number of
                          columns is usually better to use an interpolation
                          method to produce better results (see https://matplotl
                          ib.org/examples/images_contours_and_fields/interpolati
                          on_methods.html). Be default, plotHeatmap uses the
                          method `nearest` if the number of columns is 1000 or
                          less. Otherwise it uses the bilinear method. This
                          default behaviour can be changed by using any of the
                          following options: "nearest", "bilinear", "bicubic",
                          "gaussian" (default: auto)
    --dpi DPI             Set the DPI to save the figure. (default: 200)

  Clustering arguments:
    --kmeans KMEANS       Number of clusters to compute. When this option is
                          set, the matrix is split into clusters using the
                          k-means algorithm. Only works for data that is not
                          grouped, otherwise only the first group will be
                          clustered. If more specific clustering methods are
                          required, then save the underlying matrix and run the
                          clustering using other software. The plotting of the
                          clustering may fail with an error if a cluster has
                          very few members compared to the total number or
                          regions. (default: None)
    --hclust HCLUST       Number of clusters to compute. When this option is
                          set, then the matrix is split into clusters using the
                          hierarchical clustering algorithm, using "ward
                          linkage". Only works for data that is not grouped,
                          otherwise only the first group will be clustered.
                          --hclust could be very slow if you have >1000 regions.
                          In those cases, you might prefer --kmeans or if more
                          clustering methods are required you can save the
                          underlying matrix and run the clustering using other
                          software. The plotting of the clustering may fail with
                          an error if a cluster has very few members compared to
                          the total number of regions. (default: None)
    --silhouette          Compute the silhouette score for regions. This is only
                          applicable if clustering has been performed. The
                          silhouette score is a measure of how similar a region
                          is to other regions in the same cluster as opposed to
                          those in other clusters. It will be reported in the
                          final column of the BED file with regions. The
                          silhouette evaluation can be very slow when you have
                          morethan 100 000 regions. (default: False)

  Optional arguments:
    --help, -h            show this help message and exit
    --version             show program's version number and exit
    --plotType {lines,fill,se,std}
                          "lines" will plot the profile line based on the
                          average type selected. "fill" fills the region between
                          zero and the profile curve. The fill in color is semi
                          transparent to distinguish different profiles. "se"
                          and "std" color the region between the profile and the
                          standard error or standard deviation of the data.
                          (default: lines)
    --sortRegions {descend,ascend,no,keep}
                          Whether the heatmap should present the regions sorted.
                          The default is to sort in descending order based on
                          the mean value per region. Note that "keep" and "no"
                          are the same thing. (default: descend)
    --sortUsing {mean,median,max,min,sum,region_length}
                          Indicate which method should be used for sorting. For
                          each row the method is computed. For region_length, a
                          dashed line is drawn at the end of the region
                          (reference point TSS and center) or the beginning of
                          the region (reference point TES) as appropriate.
                          (default: mean)
    --sortUsingSamples SORTUSINGSAMPLES [SORTUSINGSAMPLES ...]
                          List of sample numbers (order as in matrix), that are
                          used for sorting by --sortUsing, no value uses all
                          samples, example: --sortUsingSamples 1 3 (default:
                          None)
    --linesAtTickMarks    Draw dashed lines from all tick marks through the
                          heatmap. This is then similar to the dashed line draw
                          at region bounds when using a reference point and
                          --sortUsing region_length (default: False)
    --clusterUsingSamples CLUSTERUSINGSAMPLES [CLUSTERUSINGSAMPLES ...]
                          List of sample numbers (order as in matrix), that are
                          used for clustering by --kmeans or --hclust if not
                          given, all samples are taken into account for
                          clustering. Example: --ClusterUsingSamples 1 3
                          (default: None)
    --averageTypeSummaryPlot {mean,median,min,max,std,sum}
                          Define the type of statistic that should be plotted in
                          the summary image above the heatmap. The options are:
                          "mean", "median", "min", "max", "sum" and "std".
                          (default: mean)
    --missingDataColor MISSINGDATACOLOR
                          If --missingDataAsZero was not set, such cases will be
                          colored in black by default. Using this parameter, a
                          different color can be set. A value between 0 and 1
                          will be used for a gray scale (black is 0). For a list
                          of possible color names see: http://packages.python.or
                          g/ete2/reference/reference_svgcolors.html. Other
                          colors can be specified using the #rrggbb notation.
                          (default: black)
    --colorMap COLORMAP [COLORMAP ...]
                          Color map to use for the heatmap. If more than one
                          heatmap is being plotted the color of each heatmap can
                          be enter individually (e.g. `--colorMap Reds Blues`).
                          Color maps are recycled if the number of color maps is
                          smaller than the number of heatmaps being plotted.
                          Available values can be seen here:
                          http://matplotlib.org/users/colormaps.html The
                          available options are: 'Accent', 'Blues', 'BrBG',
                          'BuGn', 'BuPu', 'CMRmap', 'Dark2', 'GnBu', 'Greens',
                          'Greys', 'OrRd', 'Oranges', 'PRGn', 'Paired',
                          'Pastel1', 'Pastel2', 'PiYG', 'PuBu', 'PuBuGn',
                          'PuOr', 'PuRd', 'Purples', 'RdBu', 'RdGy', 'RdPu',
                          'RdYlBu', 'RdYlGn', 'Reds', 'Set1', 'Set2', 'Set3',
                          'Spectral', 'Wistia', 'YlGn', 'YlGnBu', 'YlOrBr',
                          'YlOrRd', 'afmhot', 'autumn', 'binary', 'bone', 'brg',
                          'bwr', 'cividis', 'cool', 'coolwarm', 'copper',
                          'cubehelix', 'flag', 'gist_earth', 'gist_gray',
                          'gist_heat', 'gist_ncar', 'gist_rainbow',
                          'gist_stern', 'gist_yarg', 'gnuplot', 'gnuplot2',
                          'gray', 'hot', 'hsv', 'icefire', 'inferno', 'jet',
                          'magma', 'mako', 'nipy_spectral', 'ocean', 'pink',
                          'plasma', 'prism', 'rainbow', 'rocket', 'seismic',
                          'spring', 'summer', 'tab10', 'tab20', 'tab20b',
                          'tab20c', 'terrain', 'twilight', 'twilight_shifted',
                          'viridis', 'vlag', 'winter' (default: ['RdYlBu'])
    --alpha ALPHA         The alpha channel (transparency) to use for the
                          heatmaps. The default is 1.0 and values must be
                          between 0 and 1. (default: 1.0)
    --colorList COLORLIST [COLORLIST ...]
                          List of colors to use to create a colormap. For
                          example, if `--colorList black,yellow,blue` is set
                          (colors separated by comas) then a color map that
                          starts with black, continues to yellow and finishes in
                          blue is created. If this option is selected, it
                          overrides the --colorMap chosen. The list of valid
                          color names can be seen here:
                          http://matplotlib.org/examples/color/named_colors.html
                          Hex colors are valid (e.g #34a2b1). If individual
                          colors for different heatmaps need to be specified
                          they need to be separated by space as for example:
                          `--colorList "white,#cccccc" "white,darkred"` As for
                          --colorMap, the color lists are recycled if their
                          number is smaller thatn the number ofplotted heatmaps.
                          The number of transitions is defined by the
                          --colorNumber option. (default: None)
    --colorNumber COLORNUMBER
                          N.B., --colorList is required for an effect. This
                          controls the number of transitions from one color to
                          the other. If --colorNumber is the number of colors in
                          --colorList then there will be no transitions between
                          the colors. (default: 256)
    --zMin ZMIN [ZMIN ...], -min ZMIN [ZMIN ...]
                          Minimum value for the heatmap intensities. Multiple
                          values, separated by spaces can be set for each
                          heatmap. If the number of zMin values is smaller
                          thanthe number of heatmaps the values are recycled.
                          (default: None)
    --zMax ZMAX [ZMAX ...], -max ZMAX [ZMAX ...]
                          Maximum value for the heatmap intensities. Multiple
                          values, separated by spaces can be set for each
                          heatmap. If the number of zMax values is smaller
                          thanthe number of heatmaps the values are recycled.
                          (default: None)
    --heatmapHeight HEATMAPHEIGHT
                          Plot height in cm. The default for the heatmap height
                          is 28. The minimum value is 3 and the maximum is 100.
                          (default: 28)
    --heatmapWidth HEATMAPWIDTH
                          Plot width in cm. The default value is 4 The minimum
                          value is 1 and the maximum is 100. (default: 4)
    --whatToShow {plot, heatmap and colorbar,plot and heatmap,heatmap only,heatmap and colorbar}
                          The default is to include a summary or profile plot on
                          top of the heatmap and a heatmap colorbar. Other
                          options are: "plot and heatmap", "heatmap only",
                          "heatmap and colorbar", and the default "plot, heatmap
                          and colorbar". (default: plot, heatmap and colorbar)
    --boxAroundHeatmaps BOXAROUNDHEATMAPS
                          By default black boxes are plot around heatmaps. This
                          can be turned off by setting --boxAroundHeatmaps no
                          (default: yes)
    --xAxisLabel XAXISLABEL, -x XAXISLABEL
                          Description for the x-axis label. (default: gene
                          distance (bp))
    --startLabel STARTLABEL
                          [only for scale-regions mode] Label shown in the plot
                          for the start of the region. Default is TSS
                          (transcription start site), but could be changed to
                          anything, e.g. "peak start". Same for the --endLabel
                          option. See below. (default: TSS)
    --endLabel ENDLABEL   [only for scale-regions mode] Label shown in the plot
                          for the region end. Default is TES (transcription end
                          site). (default: TES)
    --refPointLabel REFPOINTLABEL
                          [only for reference-point mode] Label shown in the
                          plot for the reference-point. Default is the same as
                          the reference point selected (e.g. TSS), but could be
                          anything, e.g. "peak start". (default: None)
    --labelRotation LABEL_ROTATION
                          Rotation of the X-axis labels in degrees. The default
                          is 0, positive values denote a counter-clockwise
                          rotation. (default: 0.0)
    --regionsLabel REGIONSLABEL [REGIONSLABEL ...], -z REGIONSLABEL [REGIONSLABEL ...]
                          Labels for the regions plotted in the heatmap. If more
                          than one region is being plotted, a list of labels
                          separated by spaces is required. If a label itself
                          contains a space, then quotes are needed. For example,
                          --regionsLabel label_1, "label 2". (default: None)
    --samplesLabel SAMPLESLABEL [SAMPLESLABEL ...]
                          Labels for the samples plotted. The default is to use
                          the file name of the sample. The sample labels should
                          be separated by spaces and quoted if a label
                          itselfcontains a space E.g. --samplesLabel label-1
                          "label 2" (default: None)
    --plotTitle PLOTTITLE, -T PLOTTITLE
                          Title of the plot, to be printed on top of the
                          generated image. Leave blank for no title. (default: )
    --yAxisLabel YAXISLABEL, -y YAXISLABEL
                          Y-axis label for the top panel. (default: )
    --yMin YMIN [YMIN ...]
                          Minimum value for the Y-axis. Multiple values,
                          separated by spaces can be set for each profile. If
                          the number of yMin values is smaller thanthe number of
                          plots, the values are recycled. (default: None)
    --yMax YMAX [YMAX ...]
                          Maximum value for the Y-axis. Multiple values,
                          separated by spaces can be set for each profile. If
                          the number of yMin values is smaller thanthe number of
                          plots, the values are recycled. (default: None)
    --legendLocation {best,upper-right,upper-left,upper-center,lower-left,lower-right,lower-center,center,center-left,center-right,none}
                          Location for the legend in the summary plot. Note that
                          "none" does not work for the profiler. (default: best)
    --perGroup            The default is to plot all groups of regions by
                          sample. Using this option instead plots all samples by
                          group of regions. Note that this is only useful if you
                          have multiple groups of regions. by sample rather than
                          group. (default: False)
    --plotFileFormat      Image format type. If given, this option overrides the
                          image format based on the plotFile ending. The
                          available options are: "png", "eps", "pdf", "plotly"
                          and "svg" (default: None)
    --verbose             If set, warning messages and additional information
                          are given. (default: False)

  An example usage is: plotHeatmap -m <matrix file>