cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/deeptools:v0.0.1
inputs:
  score_files:
    type:
    - File
    - File[]
    inputBinding:
      position: 5
      prefix: --scoreFileName
    doc: |
      BigWig file(s) containing the scores to be plotted
  regions_files:
    type:
    - File
    - File[]
    inputBinding:
      position: 6
      prefix: --regionsFileName
    doc: |
      File name or names, in BED format, containing the regions to plot
  reference_point:
    type:
    - 'null'
    - type: enum
      name: reference
      symbols:
      - TSS
      - TES
      - center
    inputBinding:
      position: 7
      prefix: --referencePoint
    doc: |
      The reference point for the plotting could be either the region start (TSS),
      the region end (TES) or the center of the region. Note that regardless of what
      you specify, plotHeatmap/plotProfile will default to using “TSS” as the label.
      Default: TSS
  before_region_start_length:
    type: int?
    inputBinding:
      position: 8
      prefix: --beforeRegionStartLength
    doc: |
      Distance upstream of the reference-point selected.
      Default: 500
  after_region_start_length:
    type: int?
    inputBinding:
      position: 9
      prefix: --afterRegionStartLength
    doc: |
      Distance downstream of the reference-point selected.
      Default: 1500
  nan_after_end:
    type: boolean?
    inputBinding:
      position: 10
      prefix: --nanAfterEnd
    doc: |
      If set, any values after the region end are discarded. This is useful to visualize
      the region end when not using the scale-regions mode and when the reference-point
      is set to the TSS.
  bin_size:
    type: int?
    inputBinding:
      position: 11
      prefix: --binSize
    doc: |
      Length, in bases, of the non-overlapping bins for averaging the score over
      the regions length.
      Default: 10
  sort_regions:
    type:
    - 'null'
    - type: enum
      name: sort
      symbols:
      - descend
      - ascend
      - 'no'
      - keep
    inputBinding:
      position: 12
      prefix: --sortRegions
    doc: |
      Whether the output file should present the regions sorted. The default is to
      not sort the regions. Note that this is only useful if you plan to plot the
      results yourself and not, for example, with plotHeatmap, which will override this.
      Note also that unsorted output will be in whatever order the regions happen to
      be processed in and not match the order in the input files. If you require the
      output order to match that of the input regions, then either specify “keep” or
      use computeMatrixOperations to resort the results file.
      Default: keep
  sort_using:
    type:
    - 'null'
    - type: enum
      name: sort_type
      symbols:
      - mean
      - median
      - max
      - min
      - sum
      - region_length
    inputBinding:
      position: 13
      prefix: --sortUsing
    doc: |
      Indicate which method should be used for sorting. The value is computed for
      each row. Note that the region_length option will lead to a dotted line
      within the heatmap that indicates the end of the regions.
      Default: mean
  average_type_bins:
    type:
    - 'null'
    - type: enum
      name: average
      symbols:
      - mean
      - median
      - min
      - max
      - std
      - sum
    inputBinding:
      position: 14
      prefix: --averageTypeBins
    doc: |
      Define the type of statistic that should be used over the bin size range.
      The options are: “mean”, “median”, “min”, “max”, “sum” and “std”.
      Default: mean
  missing_data_as_zero:
    type: boolean?
    inputBinding:
      position: 15
      prefix: --missingDataAsZero
    doc: |
      If set, missing data (NAs) will be treated as zeros. The default is to ignore such cases,
      which will be depicted as black areas in a heatmap. (see the –missingDataColor argument
      of the plotHeatmap command for additional options)
  skip_zeros:
    type: boolean?
    inputBinding:
      position: 16
      prefix: --skipZeros
    doc: |
      Whether regions with only scores of zero should be included or not.
      Default is to include them
  min_threshold:
    type: float?
    inputBinding:
      position: 17
      prefix: --minThreshold
    doc: |
      Numeric value. Any region containing a value that is less than or equal to this will be skipped.
      This is useful to skip, for example, genes where the read count is zero for any of the bins.
      This could be the result of unmappable areas and can bias the overall results.
      Default: None
  max_threshold:
    type: float?
    inputBinding:
      position: 18
      prefix: --maxThreshold
    doc: |
      Numeric value. Any region containing a value greater than or equal to this will be skipped.
      The maxThreshold is useful to skip those few regions with very high read counts (e.g. micro satellites)
      that may bias the average values.
      Default: None
  samples_label:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      position: 19
      prefix: --samplesLabel
    doc: |
      Labels for the samples. This will then be passed to plotHeatmap and plotProfile.
      The default is to use the file name of the sample. The sample labels should be
      separated by spaces and quoted if a label itselfcontains a space
      E.g. –samplesLabel label-1 “label 2”
  blacklisted_regions:
    type: File?
    inputBinding:
      position: 20
      prefix: --blackListFileName
    doc: |
      A BED file containing regions that should be excluded from all analyses. Currently
      this works by rejecting genomic chunks that happen to overlap an entry. Consequently,
      for BAM files, if a read partially overlaps a blacklisted region or a fragment spans
      over it, then the read/fragment might still be considered
  output_filename:
    type: string
    inputBinding:
      position: 21
      prefix: --outFileName
    doc: |
      File name to save the gzipped matrix file needed by the “plotHeatmap” and “plotProfile” tools
  threads:
    type: int?
    inputBinding:
      position: 22
      prefix: --numberOfProcessors
    doc: |
      Number of processors to use
outputs:
  scores_matrix:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)
    doc: |
      Scores per genome regions matrix,
      File that can be used with plotHeatmap and plotProfiles
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- computeMatrix
- reference-point
- --verbose
stdout: compute_matrix_stdout.log
stderr: compute_matrix_stderr.log
label: computeMatrix - prepares an intermediate file that can be used with plotHeatmap and plotProfiles
doc: |
  Tool calculates scores per genome regions and prepares an intermediate file that can be used
  with plotHeatmap and plotProfiles. Typically, the genome regions are genes, but any other
  regions defined in a BED file can be used. computeMatrix accepts multiple score files
  (bigWig format) and multiple regions files (BED format). This tool can also be used to filter
  and sort regions according to their score.
