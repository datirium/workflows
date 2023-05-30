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
      prefix: "--scoreFileName"
    doc: |
      BigWig file(s) containing the scores to be plotted

  regions_files:
    type:
    - File
    - File[]
    inputBinding:
      position: 6
      prefix: "--regionsFileName"
    doc: |
      File name or names, in BED format, containing the regions to plot

  reference_point:
    type:
    - "null"
    - type: enum
      name: "reference"
      symbols: ["TSS", "TES", "center"]
    inputBinding:
      position: 7
      prefix: "--referencePoint"
    doc: |
      The reference point for the plotting could be either the region start (TSS),
      the region end (TES) or the center of the region. Note that regardless of what
      you specify, plotHeatmap/plotProfile will default to using “TSS” as the label.
      Default: TSS

  before_region_start_length:
    type: int?
    inputBinding:
      position: 8
      prefix: "--beforeRegionStartLength"
    doc: |
      Distance upstream of the reference-point selected.
      Default: 500

  after_region_start_length:
    type: int?
    inputBinding:
      position: 9
      prefix: "--afterRegionStartLength"
    doc: |
      Distance downstream of the reference-point selected.
      Default: 1500

  nan_after_end:
    type: boolean?
    inputBinding:
      position: 10
      prefix: "--nanAfterEnd"
    doc: |
      If set, any values after the region end are discarded. This is useful to visualize
      the region end when not using the scale-regions mode and when the reference-point
      is set to the TSS.

  bin_size:
    type: int?
    inputBinding:
      position: 11
      prefix: "--binSize"
    doc: |
      Length, in bases, of the non-overlapping bins for averaging the score over
      the regions length.
      Default: 10

  sort_regions:
    type:
    - "null"
    - type: enum
      name: "sort"
      symbols: ["descend", "ascend", "no", "keep"]
    inputBinding:
      position: 12
      prefix: "--sortRegions"
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
    - "null"
    - type: enum
      name: "sort_type"
      symbols: ["mean", "median", "max", "min", "sum", "region_length"]
    inputBinding:
      position: 13
      prefix: "--sortUsing"
    doc: |
      Indicate which method should be used for sorting. The value is computed for
      each row. Note that the region_length option will lead to a dotted line
      within the heatmap that indicates the end of the regions.
      Default: mean

  average_type_bins:
    type:
    - "null"
    - type: enum
      name: "average"
      symbols: ["mean", "median", "min", "max", "std", "sum"]
    inputBinding:
      position: 14
      prefix: "--averageTypeBins"
    doc: |
      Define the type of statistic that should be used over the bin size range.
      The options are: “mean”, “median”, “min”, “max”, “sum” and “std”.
      Default: mean

  missing_data_as_zero:
    type: boolean?
    inputBinding:
      position: 15
      prefix: "--missingDataAsZero"
    doc: |
      If set, missing data (NAs) will be treated as zeros. The default is to ignore such cases,
      which will be depicted as black areas in a heatmap. (see the –missingDataColor argument
      of the plotHeatmap command for additional options)

  skip_zeros:
    type: boolean?
    inputBinding:
      position: 16
      prefix: "--skipZeros"
    doc: |
      Whether regions with only scores of zero should be included or not.
      Default is to include them

  min_threshold:
    type: float?
    inputBinding:
      position: 17
      prefix: "--minThreshold"
    doc: |
      Numeric value. Any region containing a value that is less than or equal to this will be skipped.
      This is useful to skip, for example, genes where the read count is zero for any of the bins.
      This could be the result of unmappable areas and can bias the overall results.
      Default: None

  max_threshold:
    type: float?
    inputBinding:
      position: 18
      prefix: "--maxThreshold"
    doc: |
      Numeric value. Any region containing a value greater than or equal to this will be skipped.
      The maxThreshold is useful to skip those few regions with very high read counts (e.g. micro satellites)
      that may bias the average values.
      Default: None

  samples_label:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      position: 19
      prefix: "--samplesLabel"
    doc: |
      Labels for the samples. This will then be passed to plotHeatmap and plotProfile.
      The default is to use the file name of the sample. The sample labels should be
      separated by spaces and quoted if a label itselfcontains a space
      E.g. –samplesLabel label-1 “label 2”

  blacklisted_regions:
    type: File?
    inputBinding:
      position: 20
      prefix: "--blackListFileName"
    doc: |
      A BED file containing regions that should be excluded from all analyses. Currently
      this works by rejecting genomic chunks that happen to overlap an entry. Consequently,
      for BAM files, if a read partially overlaps a blacklisted region or a fragment spans
      over it, then the read/fragment might still be considered

  output_filename:
    type: string
    inputBinding:
      position: 21
      prefix: "--outFileName"
    doc: |
      File name to save the gzipped matrix file needed by the “plotHeatmap” and “plotProfile” tools

  threads:
    type: int?
    inputBinding:
      position: 22
      prefix: "--numberOfProcessors"
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


baseCommand: ["computeMatrix", "reference-point", "--verbose"]


stdout: compute_matrix_stdout.log
stderr: compute_matrix_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/deeptools-metadata.yaml

label: "computeMatrix - prepares an intermediate file that can be used with plotHeatmap and plotProfiles"
s:name: "computeMatrix - prepares an intermediate file that can be used with plotHeatmap and plotProfiles"
s:alternateName: "computeMatrix - prepares an intermediate file that can be used with plotHeatmap and plotProfiles"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/deeptools-computematrix-referencepoint.cwl
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
  Tool calculates scores per genome regions and prepares an intermediate file that can be used
  with plotHeatmap and plotProfiles. Typically, the genome regions are genes, but any other
  regions defined in a BED file can be used. computeMatrix accepts multiple score files
  (bigWig format) and multiple regions files (BED format). This tool can also be used to filter
  and sort regions according to their score.


s:about: |
  usage: An example usage is:
    computeMatrix reference-point -S <biwig file(s)> -R <bed file> -a 3000 -b 3000

  optional arguments:
    -h, --help            show this help message and exit

  Required arguments:
    --regionsFileName File [File ...], -R File [File ...]
                          File name or names, in BED or GTF format, containing
                          the regions to plot. If multiple bed files are given,
                          each one is considered a group that can be plotted
                          separately. Also, adding a "#" symbol in the bed file
                          causes all the regions until the previous "#" to be
                          considered one group. (default: None)
    --scoreFileName File [File ...], -S File [File ...]
                          bigWig file(s) containing the scores to be plotted.
                          Multiple files should be separated by spaced. BigWig
                          files can be obtained by using the bamCoverage or
                          bamCompare tools. More information about the bigWig
                          file format can be found at
                          http://genome.ucsc.edu/goldenPath/help/bigWig.html
                          (default: None)

  Output options:
    --outFileName OUTFILENAME, -out OUTFILENAME, -o OUTFILENAME
                          File name to save the gzipped matrix file needed by
                          the "plotHeatmap" and "plotProfile" tools. (default:
                          None)
    --outFileNameMatrix FILE
                          If this option is given, then the matrix of values
                          underlying the heatmap will be saved using the
                          indicated name, e.g. IndividualValues.tab.This matrix
                          can easily be loaded into R or other programs.
                          (default: None)
    --outFileSortedRegions BED file
                          File name in which the regions are saved after skiping
                          zeros or min/max threshold values. The order of the
                          regions in the file follows the sorting order
                          selected. This is useful, for example, to generate
                          other heatmaps keeping the sorting of the first
                          heatmap. Example: Heatmap1sortedRegions.bed (default:
                          None)

  Optional arguments:
    --version             show program's version number and exit
    --referencePoint {TSS,TES,center}
                          The reference point for the plotting could be either
                          the region start (TSS), the region end (TES) or the
                          center of the region. Note that regardless of what you
                          specify, plotHeatmap/plotProfile will default to using
                          "TSS" as the label. (Default: TSS)
    --beforeRegionStartLength INT bp, -b INT bp, --upstream INT bp
                          Distance upstream of the reference-point selected.
                          (Default: 500)
    --afterRegionStartLength INT bp, -a INT bp, --downstream INT bp
                          Distance downstream of the reference-point selected.
                          (Default: 1500)
    --nanAfterEnd         If set, any values after the region end are discarded.
                          This is useful to visualize the region end when not
                          using the scale-regions mode and when the reference-
                          point is set to the TSS. (default: False)
    --binSize BINSIZE, -bs BINSIZE
                          Length, in bases, of the non-overlapping bins for
                          averaging the score over the regions length. (Default:
                          10)
    --sortRegions {descend,ascend,no,keep}
                          Whether the output file should present the regions
                          sorted. The default is to not sort the regions. Note
                          that this is only useful if you plan to plot the
                          results yourself and not, for example, with
                          plotHeatmap, which will override this. Note also that
                          unsorted output will be in whatever order the regions
                          happen to be processed in and not match the order in
                          the input files. If you require the output order to
                          match that of the input regions, then either specify
                          "keep" or use computeMatrixOperations to resort the
                          results file. (Default: keep)
    --sortUsing {mean,median,max,min,sum,region_length}
                          Indicate which method should be used for sorting. The
                          value is computed for each row.Note that the
                          region_length option will lead to a dotted line within
                          the heatmap that indicates the end of the regions.
                          (Default: mean)
    --sortUsingSamples SORTUSINGSAMPLES [SORTUSINGSAMPLES ...]
                          List of sample numbers (order as in matrix), that are
                          used for sorting by --sortUsing, no value uses all
                          samples, example: --sortUsingSamples 1 3 (default:
                          None)
    --averageTypeBins {mean,median,min,max,std,sum}
                          Define the type of statistic that should be used over
                          the bin size range. The options are: "mean", "median",
                          "min", "max", "sum" and "std". The default is "mean".
                          (Default: mean)
    --missingDataAsZero   If set, missing data (NAs) will be treated as zeros.
                          The default is to ignore such cases, which will be
                          depicted as black areas in a heatmap. (see the
                          --missingDataColor argument of the plotHeatmap command
                          for additional options). (default: False)
    --skipZeros           Whether regions with only scores of zero should be
                          included or not. Default is to include them. (default:
                          False)
    --minThreshold MINTHRESHOLD
                          Numeric value. Any region containing a value that is
                          less than or equal to this will be skipped. This is
                          useful to skip, for example, genes where the read
                          count is zero for any of the bins. This could be the
                          result of unmappable areas and can bias the overall
                          results. (Default: None)
    --maxThreshold MAXTHRESHOLD
                          Numeric value. Any region containing a value greater
                          than or equal to this will be skipped. The
                          maxThreshold is useful to skip those few regions with
                          very high read counts (e.g. micro satellites) that may
                          bias the average values. (Default: None)
    --blackListFileName BED file, -bl BED file
                          A BED file containing regions that should be excluded
                          from all analyses. Currently this works by rejecting
                          genomic chunks that happen to overlap an entry.
                          Consequently, for BAM files, if a read partially
                          overlaps a blacklisted region or a fragment spans over
                          it, then the read/fragment might still be considered.
                          (default: None)
    --samplesLabel SAMPLESLABEL [SAMPLESLABEL ...]
                          Labels for the samples. This will then be passed to
                          plotHeatmap and plotProfile. The default is to use the
                          file name of the sample. The sample labels should be
                          separated by spaces and quoted if a label
                          itselfcontains a space E.g. --samplesLabel label-1
                          "label 2" (default: None)
    --smartLabels         Instead of manually specifying labels for the input
                          bigWig and BED/GTF files, this causes deepTools to use
                          the file name after removing the path and extension.
                          (default: False)
    --quiet, -q           Set to remove any warning or processing messages.
                          (default: False)
    --verbose             Being VERY verbose in the status messages. --quiet
                          will disable this. (default: False)
    --scale SCALE         If set, all values are multiplied by this number.
                          (Default: 1)
    --numberOfProcessors INT, -p INT
                          Number of processors to use. Type "max/2" to use half
                          the maximum number of processors or "max" to use all
                          available processors. (Default: 1)

  GTF/BED12 options:
    --metagene            When either a BED12 or GTF file are used to provide
                          regions, perform the computation on the merged exons,
                          rather than using the genomic interval defined by the
                          5-prime and 3-prime most transcript bound (i.e.,
                          columns 2 and 3 of a BED file). If a BED3 or BED6 file
                          is used as input, then columns 2 and 3 are used as an
                          exon. (Default: False)
    --transcriptID TRANSCRIPTID
                          When a GTF file is used to provide regions, only
                          entries with this value as their feature (column 3)
                          will be processed as transcripts. (Default:
                          transcript)
    --exonID EXONID       When a GTF file is used to provide regions, only
                          entries with this value as their feature (column 3)
                          will be processed as exons. CDS would be another
                          common value for this. (Default: exon)
    --transcript_id_designator TRANSCRIPT_ID_DESIGNATOR
                          Each region has an ID (e.g., ACTB) assigned to it,
                          which for BED files is either column 4 (if it exists)
                          or the interval bounds. For GTF files this is instead
                          stored in the last column as a key:value pair (e.g.,
                          as 'transcript_id "ACTB"', for a key of transcript_id
                          and a value of ACTB). In some cases it can be
                          convenient to use a different identifier. To do so,
                          set this to the desired key. (Default: transcript_id)

  deepBlue arguments:
    Options used only for remote bedgraph/wig files hosted on deepBlue

    --deepBlueURL DEEPBLUEURL
                          For remote files bedgraph/wiggle files hosted on
                          deepBlue, this specifies the server URL. The default
                          is "http://deepblue.mpi-inf.mpg.de/xmlrpc", which
                          should not be changed without good reason. (default:
                          http://deepblue.mpi-inf.mpg.de/xmlrpc)
    --userKey USERKEY     For remote files bedgraph/wiggle files hosted on
                          deepBlue, this specifies the user key to use for
                          access. The default is "anonymous_key", which suffices
                          for public datasets. If you need access to a
                          restricted access/private dataset, then request a key
                          from deepBlue and specify it here. (default:
                          anonymous_key)
    --deepBlueTempDir DEEPBLUETEMPDIR
                          If specified, temporary files from preloading datasets
                          from deepBlue will be written here (note, this
                          directory must exist). If not specified, where ever
                          temporary files would normally be written on your
                          system is used. (default: None)
    --deepBlueKeepTemp    If specified, temporary bigWig files from preloading
                          deepBlue datasets are not deleted. A message will be
                          printed noting where these files are and what sample
                          they correspond to. These can then be used if you wish
                          to analyse the same sample with the same regions
                          again. (default: False)