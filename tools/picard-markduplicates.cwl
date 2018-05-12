cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          return inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".bam";
        };
  - var default_output_metrics_filename = function() {
          return inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+".metrics";
        };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/picard:v2.9.0
  dockerFile: >
    $import: ./dockerfiles/picard-Dockerfile


inputs:
  java_opts:
    type:
      - "null"
      - string
    inputBinding:
      position: 1
      shellQuote: false
    doc: |
      JVM arguments should be a quoted, space separated list (e.g. "-Xms128m -Xmx512m")

  max_seq_for_disk_read_ends_map:
    type:
      - "null"
      - int
    doc: |
      This option is obsolete. ReadEnds will always be spilled to disk.
      Default value: 50000.
    inputBinding:
      position: 5
      prefix: 'MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP='
      separate: false

  max_file_handles_for_read_ends_map:
    type:
      - "null"
      - int
    doc: |
      Maximum number of file handles to keep open when spilling read ends to disk.
      Set this number a little lower than the per-process maximum number of file that may be open.
      This number can be found by executing the 'ulimit -n' command on a Unix system.
      Default value: 8000.
    inputBinding:
      position: 6
      prefix: 'MAX_FILE_HANDLES_FOR_READ_ENDS_MAP='
      separate: false

  sorting_coll_size_ratio:
    type:
      - "null"
      - double
    doc: |
      This number, plus the maximum RAM available to the JVM, determine the memory footprint used by some
      of the sorting collections. If you are running out of memory, try reducing this number.
      Default value: 0.25.
    inputBinding:
      position: 7
      prefix: 'SORTING_COLLECTION_SIZE_RATIO='
      separate: false

  barcode_tag:
    type:
      - "null"
      - string
    doc: |
      Barcode SAM tag (ex. BC for 10X Genomics)
      Default value: null.
    inputBinding:
      position: 8
      prefix: 'BARCODE_TAG='
      separate: false

  read_one_barcode_tag:
    type:
      - "null"
      - string
    doc: |
      Read one barcode SAM tag (ex. BX for 10X Genomics)
      Default value: null.
    inputBinding:
      position: 9
      prefix: 'READ_ONE_BARCODE_TAG='
      separate: false

  read_two_barcode_tag:
    type:
      - "null"
      - string
    doc: |
      Read two barcode SAM tag (ex. BX for 10X Genomics).
      Default value: null.
    inputBinding:
      position: 10
      prefix: 'READ_TWO_BARCODE_TAG='
      separate: false

  remove_seq_dup:
    type:
      - "null"
      - boolean
    doc: |
      If true remove 'optical' duplicates and other duplicates that appear to have arisen from the sequencing process
      instead of the library preparation process, even if REMOVE_DUPLICATES is false. If REMOVE_DUPLICATES is true,
      all duplicates are removed and this option is ignored.
      Default value: false.
    inputBinding:
      position: 11
      prefix: 'REMOVE_SEQUENCING_DUPLICATES='
      separate: false

  tagging_policy:
    type:
      - "null"
      - type: enum
        name: "tagging_policy"
        symbols: ['DontTag','OpticalOnly','All']
    doc: |
      Determines how duplicate types are recorded in the DT optional attribute.
      Possible values: {DontTag, OpticalOnly, All}
      Default value: DontTag.
    inputBinding:
      position: 12
      prefix: 'TAGGING_POLICY='
      separate: false

  input_file:
    type:
      - File
    doc: |
      One or more input SAM or BAM files to analyze.
      Must be coordinate sorted.
    inputBinding:
      position: 13
      prefix: 'INPUT='
      separate: false

  output_filename:
    type:
      - string
    doc: |
      The output file to write marked records to Required.
    inputBinding:
      position: 14
      prefix: 'OUTPUT='
      separate: false
      valueFrom: |
        ${
            if (self == null){
              return default_output_filename();
            } else {
              return self;
            }
        }
    default: null

  output_metrics_filename:
    type:
      - string
    doc: |
      File to write duplication metrics to Required.
    inputBinding:
      position: 15
      prefix: 'METRICS_FILE='
      separate: false
      valueFrom: |
        ${
            if (self == null){
              return default_output_metrics_filename();
            } else {
              return self;
            }
        }
    default: null


  remove_dup:
    type:
      - "null"
      - boolean
    doc: |
      If true do not write duplicates to the output file instead of writing them with appropriate flags set.
      Default value: false.
    inputBinding:
      position: 16
      prefix: 'REMOVE_DUPLICATES='
      separate: false

  assume_sort_order:
    type:
      - "null"
      - type: enum
        name: "assume_sort_order"
        symbols: ['unsorted','queryname','coordinate', 'duplicate']
    doc: |
      If not null, assume that the input file has this order even if the header says otherwise.
      Possible values: {unsorted, queryname, coordinate, duplicate}
      Default value: null.
    inputBinding:
      position: 17
      prefix: 'ASSUME_SORT_ORDER='
      separate: false

  dup_scoring_strategy:
    type:
      - "null"
      - type: enum
        name: "dup_scoring_strategy"
        symbols: ['SUM_OF_BASE_QUALITIES','TOTAL_MAPPED_REFERENCE_LENGTH','RANDOM']
    doc: |
      The scoring strategy for choosing the non-duplicate among candidates.
      Possible values: {SUM_OF_BASE_QUALITIES, TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM}
      Default value: SUM_OF_BASE_QUALITIES.
    inputBinding:
      position: 18
      prefix: 'DUPLICATE_SCORING_STRATEGY='
      separate: false

  program_record_id:
    type:
      - "null"
      - string
    doc: |
      The program record ID for the @PG record(s) created by this program.
      This string may have a suffix appended to avoid collision with other program record IDs.
      Default value: MarkDuplicates.
    inputBinding:
      position: 19
      prefix: 'PROGRAM_RECORD_ID='
      separate: false

  program_group_version:
    type:
      - "null"
      - string
    doc: |
      Value of VN tag of PG record to be created.
      If not specified, the version will be detected automatically.
      Default value: null.
    inputBinding:
      position: 20
      prefix: 'PROGRAM_GROUP_VERSION='
      separate: false

  program_group_cmd_line:
    type:
      - "null"
      - string
    doc: |
      Value of CL tag of PG record to be created.
      If not supplied the command line will be detected automatically.
      Default value: null.
    inputBinding:
      position: 21
      prefix: 'PROGRAM_GROUP_COMMAND_LINE='
      separate: false

  program_group_name:
    type:
      - "null"
      - string
    doc: |
      Value of PN tag of PG record to be created.
      Default value: MarkDuplicates.
    inputBinding:
      position: 22
      prefix: 'PROGRAM_GROUP_NAME='
      separate: false

  comment:
    type:
      - "null"
      - string
    doc: |
      Comment(s) to include in the output file's header.
      Default value: null.
    inputBinding:
      position: 23
      prefix: 'COMMENT='
      separate: false

  read_name_regex:
    type:
      - "null"
      - string
    doc: |
      Regular expression that can be used to parse read names in the incoming SAM file.
      Read names are parsed to extract three variables: tile/region, x coordinate and y coordinate.
      These values are used to estimate the rate of optical duplication in order to give a more accurate
      estimated library size.
      The regular expression should contain three capture groups for the three variables, in order.
      It must match the entire read name. Note that if the default regex is specified, a regex match is not
      actually done, but instead the read name is split on colon character. For 5 element names, the 3rd,
      4th and 5th elements are assumed to be tile, x and y values. For 7 element names (CASAVA 1.8),
      the 5th, 6th, and 7th elements are assumed to be tile, x and y values.
      Default value: .
    inputBinding:
      position: 24
      prefix: 'READ_NAME_REGEX='
      separate: false

  optical_dup_pixel_dist:
    type:
      - "null"
      - int
    doc: |
      The maximum offset between two duplicate clusters in order to consider them optical duplicates.
      The default is appropriate for unpatterned versions of the Illumina platform.
      For the patterned flowcell models, 2500 is moreappropriate.
      For other platforms and models, users should experiment to find what works best.
      Default value: 100.
    inputBinding:
      position: 25
      prefix: 'OPTICAL_DUPLICATE_PIXEL_DISTANCE='
      separate: false

outputs:
  output_file:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == null){
              return default_output_filename();
            } else {
              return inputs.output_filename;
            }
        }

  metrics_file:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.output_metrics_filename == null){
              return default_output_metrics_filename();
            } else {
              return inputs.output_metrics_filename;
            }
        }

baseCommand: java

arguments:
  - valueFrom: |
      ${
        return "-jar";
      }
    position: 2
  - valueFrom: |
      ${
        return "/usr/local/bin/picard.jar";
      }
    position: 3
  - valueFrom: |
      ${
        return "MarkDuplicates";
      }
    position: 4


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/picard-metadata.yaml

s:downloadUrl: https://github.com/SciDAP/workflows/blob/master/tools/picard-markduplicates.cwl
s:codeRepository: https://github.com/SciDAP/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:author:
  class: s:Person
  s:name: Michael Kotliar
  s:email: mailto:michael.kotliar@cchmc.org
  s:sameAs:
  - id: http://orcid.org/0000-0002-6486-3898
  s:worksFor:
  - class: s:Organization
    s:name: Cincinnati Children's Hospital Medical Center
    s:location: 3333 Burnet Ave, Cincinnati, OH 45229-3026
    s:department:
    - class: s:Organization
      s:name: Barski Lab

doc: |
  USAGE: MarkDuplicates [options]

  Documentation: http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

  Identifies duplicate reads.
  This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating
  from a single fragment of DNA.  Duplicates can arise during sample preparation e.g. library construction using PCR.
  See also EstimateLibraryComplexity (https://broadinstitute.github.io/picard/command-line-overview.html#EstimateLibraryComplexity)
  for additional notes on PCR duplication artifacts.  Duplicate reads can also result from a single amplification cluster,
  incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument.
  These duplication artifacts are referred to as optical duplicates.

  The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a
  SAM/BAM file.  An BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes.
  After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that
  ranks reads by the sums of their base-quality scores (default method).
  The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for
  each read.  Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024.
  If you are not familiar with this type of annotation, please see the following blog post
  (https://www.broadinstitute.org/gatk/blog?id=7019) for additional information.

  Although the bitwise flag annotation indicates whether a read was marked as a duplicate, it does not identify the type
  of duplicate.  To do this, a new tag called the duplicate type (DT) tag was recently added as an optional output in
  the 'optional field' section of a SAM/BAM file.  Invoking the TAGGING_POLICY option, you can instruct the program to
  mark all the duplicates (All), only the optical duplicates (OpticalOnly), or no duplicates (DontTag).
  The records within the output of a SAM/BAM file will have values for the 'DT' tag (depending on the invoked
  TAGGING_POLICY), as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ).
  This tool uses the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the primary methods to
  identify and differentiate duplicate types.  Set READ_NAME_REGEX to null to skip optical duplicate detection, e.g.
  for RNA-seq or other data where duplicate sets are extremely large and estimating library complexity is not an aim.
  Note that without optical duplicate counts, library size estimation will be inaccurate.

  MarkDuplicates also produces a metrics file indicating the numbers of duplicates for both single- and paired-end reads.
  The program can take either coordinate-sorted or query-sorted inputs, however the behavior is slightly different.
  When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not
  marked as duplicates.  However, when the input is query-sorted (actually query-grouped), then unmapped mates and
  secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.
  If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.

  Usage example:

  java -jar picard.jar MarkDuplicates \
        I=input.bam \
        O=marked_duplicates.bam \
        M=marked_dup_metrics.txt
  Please see MarkDuplicates (http://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics) for detailed explanations of the output metrics.

  Version: 2.8.3-SNAPSHOT


  Options:

  --help
  -h                            Displays options specific to this tool.

  --stdhelp
  -H                            Displays options specific to this tool AND options common to all Picard command line
                                tools.

  --version                     Displays program version.

  MAX_SEQUENCES_FOR_DISK_READ_ENDS_MAP=Integer
  MAX_SEQS=Integer              This option is obsolete. ReadEnds will always be spilled to disk.  Default value: 50000.
                                This option can be set to 'null' to clear the default value.

  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=Integer
  MAX_FILE_HANDLES=Integer      Maximum number of file handles to keep open when spilling read ends to disk. Set this
                                number a little lower than the per-process maximum number of file that may be open. This
                                number can be found by executing the 'ulimit -n' command on a Unix system.  Default
                                value: 8000. This option can be set to 'null' to clear the default value.

  SORTING_COLLECTION_SIZE_RATIO=Double
                                This number, plus the maximum RAM available to the JVM, determine the memory footprint
                                used by some of the sorting collections.  If you are running out of memory, try reducing
                                this number.  Default value: 0.25. This option can be set to 'null' to clear the default
                                value.

  BARCODE_TAG=String            Barcode SAM tag (ex. BC for 10X Genomics)  Default value: null.

  READ_ONE_BARCODE_TAG=String   Read one barcode SAM tag (ex. BX for 10X Genomics)  Default value: null.

  READ_TWO_BARCODE_TAG=String   Read two barcode SAM tag (ex. BX for 10X Genomics)  Default value: null.

  REMOVE_SEQUENCING_DUPLICATES=Boolean
                                If true remove 'optical' duplicates and other duplicates that appear to have arisen from
                                the sequencing process instead of the library preparation process, even if
                                REMOVE_DUPLICATES is false. If REMOVE_DUPLICATES is true, all duplicates are removed and
                                this option is ignored.  Default value: false. This option can be set to 'null' to clear
                                the default value. Possible values: {true, false}

  TAGGING_POLICY=DuplicateTaggingPolicy
                                Determines how duplicate types are recorded in the DT optional attribute.  Default value:
                                DontTag. This option can be set to 'null' to clear the default value. Possible values:
                                {DontTag, OpticalOnly, All}

  INPUT=String
  I=String                      One or more input SAM or BAM files to analyze. Must be coordinate sorted.  Default value:
                                null. This option may be specified 0 or more times.

  OUTPUT=File
  O=File                        The output file to write marked records to  Required.

  METRICS_FILE=File
  M=File                        File to write duplication metrics to  Required.

  REMOVE_DUPLICATES=Boolean     If true do not write duplicates to the output file instead of writing them with
                                appropriate flags set.  Default value: false. This option can be set to 'null' to clear
                                the default value. Possible values: {true, false}

  ASSUME_SORTED=Boolean
  AS=Boolean                    If true, assume that the input file is coordinate sorted even if the header says
                                otherwise. Deprecated, used ASSUME_SORT_ORDER=coordinate instead.  Default value: false.
                                This option can be set to 'null' to clear the default value. Possible values: {true,
                                false}  Cannot be used in conjuction with option(s) ASSUME_SORT_ORDER (ASO)

  ASSUME_SORT_ORDER=SortOrder
  ASO=SortOrder                 If not null, assume that the input file has this order even if the header says otherwise.
                                Default value: null. Possible values: {unsorted, queryname, coordinate, duplicate}
                                Cannot be used in conjuction with option(s) ASSUME_SORTED (AS)

  DUPLICATE_SCORING_STRATEGY=ScoringStrategy
  DS=ScoringStrategy            The scoring strategy for choosing the non-duplicate among candidates.  Default value:
                                SUM_OF_BASE_QUALITIES. This option can be set to 'null' to clear the default value.
                                Possible values: {SUM_OF_BASE_QUALITIES, TOTAL_MAPPED_REFERENCE_LENGTH, RANDOM}

  PROGRAM_RECORD_ID=String
  PG=String                     The program record ID for the @PG record(s) created by this program. Set to null to
                                disable PG record creation.  This string may have a suffix appended to avoid collision
                                with other program record IDs.  Default value: MarkDuplicates. This option can be set to
                                'null' to clear the default value.

  PROGRAM_GROUP_VERSION=String
  PG_VERSION=String             Value of VN tag of PG record to be created. If not specified, the version will be
                                detected automatically.  Default value: null.

  PROGRAM_GROUP_COMMAND_LINE=String
  PG_COMMAND=String             Value of CL tag of PG record to be created. If not supplied the command line will be
                                detected automatically.  Default value: null.

  PROGRAM_GROUP_NAME=String
  PG_NAME=String                Value of PN tag of PG record to be created.  Default value: MarkDuplicates. This option
                                can be set to 'null' to clear the default value.

  COMMENT=String
  CO=String                     Comment(s) to include in the output file's header.  Default value: null. This option may
                                be specified 0 or more times.

  READ_NAME_REGEX=String        Regular expression that can be used to parse read names in the incoming SAM file. Read
                                names are parsed to extract three variables: tile/region, x coordinate and y coordinate.
                                These values are used to estimate the rate of optical duplication in order to give a more
                                accurate estimated library size. Set this option to null to disable optical duplicate
                                detection, e.g. for RNA-seq or other data where duplicate sets are extremely large and
                                estimating library complexity is not an aim. Note that without optical duplicate counts,
                                library size estimation will be inaccurate. The regular expression should contain three
                                capture groups for the three variables, in order. It must match the entire read name.
                                Note that if the default regex is specified, a regex match is not actually done, but
                                instead the read name  is split on colon character. For 5 element names, the 3rd, 4th and
                                5th elements are assumed to be tile, x and y values. For 7 element names (CASAVA 1.8),
                                the 5th, 6th, and 7th elements are assumed to be tile, x and y values.  Default value:
                                <optimized capture of last three ':' separated fields as numeric values>. This option can
                                be set to 'null' to clear the default value.

  OPTICAL_DUPLICATE_PIXEL_DISTANCE=Integer
                                The maximum offset between two duplicate clusters in order to consider them optical
                                duplicates. The default is appropriate for unpatterned versions of the Illumina platform.
                                For the patterned flowcell models, 2500 is moreappropriate. For other platforms and
                                models, users should experiment to find what works best.  Default value: 100. This option
                                can be set to 'null' to clear the default value.
