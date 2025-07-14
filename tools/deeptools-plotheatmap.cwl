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
      prefix: --matrixFile
    doc: |
      Matrix file from the computeMatrix tool
  output_filename:
    type: string
    inputBinding:
      position: 6
      prefix: --outFileName
    doc: |
      File name to save the image to. The file ending will be used to determine the image format.
      The available options are: “png”, “eps”, “pdf” and “svg”, e.g., MyHeatmap.png
  interpolation_method:
    type:
    - 'null'
    - type: enum
      name: interpolation_method
      symbols:
      - auto
      - nearest
      - bilinear
      - bicubic
      - gaussian
    inputBinding:
      position: 7
      prefix: --interpolationMethod
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
      prefix: --dpi
    doc: Set the DPI to save the figure.
  plot_type:
    type:
    - 'null'
    - type: enum
      name: plot_type
      symbols:
      - lines
      - fill
      - se
      - std
    inputBinding:
      position: 9
      prefix: --plotType
    doc: |
      “lines” will plot the profile line based on the average type selected.
      “fill” fills the region between zero and the profile curve. The fill in
      color is semi transparent to distinguish different profiles.
      “se” and “std” color the region between the profile and the standard error
      or standard deviation of the data.
  sort_regions:
    type:
    - 'null'
    - type: enum
      name: sort_regions
      symbols:
      - descend
      - ascend
      - 'no'
      - keep
    inputBinding:
      position: 10
      prefix: --sortRegions
    doc: |
      Whether the heatmap should present the regions sorted. The default is to sort in
      descending order based on the mean value per region. Note that “keep” and “no” are
      the same thing.
  sort_using:
    type:
    - 'null'
    - type: enum
      name: sort_using
      symbols:
      - mean
      - median
      - max
      - min
      - sum
      - region_length
    inputBinding:
      position: 11
      prefix: --sortUsing
    doc: |
      Indicates which method should be used for sorting. For each row the method is computed.
      For region_length, a dashed line is drawn at the end of the region (reference point TSS
      and center) or the beginning of the region (reference point TES) as appropriate.
  average_type_summary_plot:
    type:
    - 'null'
    - type: enum
      name: average_type_summary_plot
      symbols:
      - mean
      - median
      - min
      - max
      - std
      - sum
    inputBinding:
      position: 12
      prefix: --averageTypeSummaryPlot
    doc: |
      Define the type of statistic that should be plotted in the summary image above the heatmap.
      The options are: “mean”, “median”, “min”, “max”, “sum” and “std”.
  what_to_show:
    type:
    - 'null'
    - type: enum
      name: what_to_show
      symbols:
      - plot, heatmap and colorbar
      - plot and heatmap
      - heatmap only
      - heatmap and colorbar
    inputBinding:
      position: 13
      prefix: --whatToShow
    doc: |
      The default is to include a summary or profile plot on top of the heatmap and a heatmap colorbar.
      Other options are: “plot and heatmap”, “heatmap only”, “heatmap and colorbar”, and the default “plot,
      heatmap and colorbar”.
  x_axis_label:
    type: string?
    inputBinding:
      position: 14
      prefix: --xAxisLabel
    doc: |
      Description for the x-axis label
  start_label:
    type: string?
    inputBinding:
      position: 15
      prefix: --startLabel
    doc: |
      [only for scale-regions mode] Label shown in the plot for the start of the region.
      Default is TSS (transcription start site), but could be changed to anything, e.g.
      “peak start”. Same for the –endLabel option.
  end_label:
    type: string?
    inputBinding:
      position: 16
      prefix: --endLabel
    doc: |
      [only for scale-regions mode] Label shown in the plot for the region end.
      Default is TES (transcription end site).
  ref_point_label:
    type: string?
    inputBinding:
      position: 17
      prefix: --refPointLabel
    doc: "[only for reference-point mode] Label shown in the plot for the reference-point.\nDefault is the same as the reference point selected (e.g. TSS), but could be anything,\ne.g. “peak start”.  \n"
  label_rotation_angle:
    type: int?
    inputBinding:
      position: 18
      prefix: --labelRotation
    doc: |
      Rotation of the X-axis labels in degrees.
      The default is 0, positive values denote a counter-clockwise rotation.
  regions_label:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      position: 19
      prefix: --regionsLabel
    doc: |
      Labels for the regions plotted in the heatmap. If more than one region is being plotted, a list of
      labels separated by spaces is required. If a label itself contains a space, then quotes are needed.
      For example, –regionsLabel label_1, “label 2”.
  samples_label:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      position: 20
      prefix: --samplesLabel
    doc: |
      Labels for the samples plotted. The default is to use the file name of the sample. The sample labels
      should be separated by spaces and quoted if a label itselfcontains a space
      E.g. –samplesLabel label-1 “label 2”
  plot_title:
    type: string?
    inputBinding:
      position: 21
      prefix: --plotTitle
    doc: |
      Title of the plot, to be printed on top of the generated image.
      Leave blank for no title.
  y_axisLabel:
    type: string?
    inputBinding:
      position: 22
      prefix: --yAxisLabel
    doc: |
      Y-axis label for the top panel.
  y_min:
    type:
    - 'null'
    - int
    - int[]
    inputBinding:
      position: 23
      prefix: --yMin
    doc: |
      Minimum value for the Y-axis. Multiple values, separated by spaces can be set for each profile.
      If the number of yMin values is smaller thanthe number of plots, the values are recycled.
  y_max:
    type:
    - 'null'
    - int
    - int[]
    inputBinding:
      position: 24
      prefix: --yMax
    doc: |
      Maximum value for the Y-axis. Multiple values, separated by spaces can be set for each profile.
      If the number of yMin values is smaller thanthe number of plots, the values are recycled.
  legend_location:
    type:
    - 'null'
    - type: enum
      name: legend_location
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
      prefix: --legendLocation
    doc: |
      Location for the legend in the summary plot. Note that “none” does not work for the profiler.
  per_group:
    type: boolean?
    inputBinding:
      position: 26
      prefix: --perGroup
    doc: |
      The default is to plot all groups of regions by sample. Using this option instead plots all
      samples by group of regions. Note that this is only useful if you have multiple groups of
      regions by sample rather than group.
  plot_file_format:
    type:
    - 'null'
    - type: enum
      name: plot_file_format
      symbols:
      - png
      - pdf
      - svg
      - eps
      - plotly
    inputBinding:
      position: 27
      prefix: --plotFileFormat
    doc: |
      Image format type. If given, this option overrides the image format based on the plotFile ending.
      The available options are: “png”, “eps”, “pdf”, “plotly” and “svg”
outputs:
  heatmap_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
    doc: Heatmap file
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- plotHeatmap
- --verbose
stdout: plot_heatmap_stdout.log
stderr: plot_heatmap_stderr.log
label: plotHeatmap - tool creates a heatmap for scores associated with genomic regions
doc: |
  This tool creates a heatmap for scores associated with genomic regions.
  The program requires a matrix file generated by the tool computeMatrix.
