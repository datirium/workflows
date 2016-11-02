#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: scidap/sratoolkit:v2.8.0
  dockerFile: >
    $import: sratoolkit-Dockerfile

inputs:

  inputFiles:
    type: File[]
    inputBinding:
      position: 60
    doc: |
      Input file[s]

  split-spot:
    type: boolean?
    inputBinding:
      position: 2
      prefix: "--split-spot"
    doc: |
      Split spots into individual reads

  minSpotId:
    type: int?
    inputBinding:
      position: 3
      prefix: "--minSpotId"
    doc: |
      Minimum spot id

  maxSpotId:
    type: int?
    inputBinding:
      position: 4
      prefix: "--maxSpotId"
    doc: |
      Maximum spot id

  spot-groups:
    type:
      - "null"
      - type: array
        items: string
    inputBinding:
      position: 5
      itemSeparator: ","
      prefix: "--spot-groups"
    doc: |
      Filter by SPOT_GROUP (member): name[,...]

  clip:
    type: boolean?
    inputBinding:
      position: 6
      prefix: "--clip"
    doc: |
      Apply left and right clips

  minReadLen:
    type: int?
    inputBinding:
      position: 7
      prefix: "--minReadLen"
    doc: |
      Filter by sequence length >= <len>

  read-filter:
    type:
      - "null"
      - type: enum
        name: "filter"
        symbols: ["pass","reject","criteria","redacted"]
    inputBinding:
      position: 8
      prefix: "--read-filter"
    doc: |
      Split into files by READ_FILTER value optionally filter by value: pass|reject|criteria|redacted


  qual-filter:
    type: boolean?
    inputBinding:
      position: 9
      prefix: "--qual-filter"
    doc: |
      Filter used in early 1000 Genomes data: no sequences starting or ending with >= 10N

  qual-filter-1:
    type: boolean?
    inputBinding:
      position: 10
      prefix: "--qual-filter-1"
    doc: |
      Filter used in current 1000 Genomes data

  aligned:
    type: boolean?
    inputBinding:
      position: 11
      prefix: "--aligned"
    doc: |
      Dump only aligned sequences

  unaligned:
    type: boolean?
    inputBinding:
      position: 12
      prefix: "--unaligned"
    doc: |
      Dump only unaligned sequences

  aligned-region:
    type: string?
    inputBinding:
      position: 13
      prefix: "--aligned-region"
    doc: |
      Filter by position on genome.
      Name can either be accession.version
      (ex:NC_000001.10) or file specific name
      (ex:"chr1" or "1"). "from" and "to" are 1-based coordinates

  matepair-distance:
    type:
      - "null"
      - type: enum
        name: "distance"
        symbols: ["from-to","unknown"]
    inputBinding:
      position: 14
      prefix: "--matepair-distance"
    doc: |
      Filter by distance beiween matepairs.
      Use "unknown" to find matepairs split
      between the references. Use from-to to limit
      matepair distance on the same reference

  skip-technical:
    type: boolean?
    inputBinding:
      position: 15
      prefix: "--skip-technical"
    doc: |
      Dump only biological reads

#  outdir:
#    type: Directory?
#    inputBinding:
#      position: 16
#      prefix: "--outdir"
#    doc: |
#       Output directory, default is working directory '.'

#  stdout:
#    type: boolean?
#    inputBinding:
#      position: 17
#      prefix: "--stdout"
#    doc: |
#      Output to stdout, all split data become
#      joined into single stream

  gzip:
    type: boolean?
    inputBinding:
      position: 18
      prefix: "--gzip"
    doc: |
      Compress output using gzip

  bzip2:
    type: boolean?
    inputBinding:
      position: 19
      prefix: "--bzip2"
    doc: |
      Compress output using bzip2

  split-files:
    type: boolean?
    inputBinding:
      position: 20
      prefix: "--split-files"
    doc: |
      Dump each read into separate file.
      Files will receive suffix corresponding to read number

  split-3:
    type: boolean?
    inputBinding:
      position: 21
      prefix: "--split-3"
    doc: |
      Legacy 3-file splitting for mate-pairs:
      First biological reads satisfying dumping
      conditions are placed in files *_1.fastq and
      *_2.fastq If only one biological read is
      present it is placed in *.fastq Biological
      reads and above are ignored.

  spot-group:
    type: boolean?
    inputBinding:
      position: 22
      prefix: "--spot-group"
    doc: |
      Split into files by SPOT_GROUP (member name)

  group-in-dirs:
    type: boolean?
    inputBinding:
      position: 23
      prefix: "--group-in-dirs"
    doc: |
      Split into subdirectories instead of files

  keep-empty-files:
    type: boolean?
    inputBinding:
      position: 24
      prefix: "--keep-empty-files"
    doc: |
      Do not delete empty files

  dumpcs:
    type: string?
    inputBinding:
      position: 25
      prefix: "--dumpcs"
    doc: |
      Formats sequence using color space (default
      for SOLiD),"cskey" may be specified for
      translation

  dumpbase:
    type: boolean?
    inputBinding:
      position: 26
      prefix: "--dumpbase"
    doc: |
      Formats sequence using base space (default
      for other than SOLiD).

  offset:
    type: int?
    inputBinding:
      position: 27
      prefix: "--offset"
    doc: |
      Offset to use for quality conversion, default is 33

  fasta:
    type: int?
    inputBinding:
      position: 28
      prefix: "--fasta"
    doc: |
      FASTA only, no qualities, optional line
      wrap width (set to zero for no wrapping)

  suppress-qual-for-cskey:
    type: boolean?
    inputBinding:
      position: 29
      prefix: "--suppress-qual-for-cskey"
    doc: |
      supress quality-value for cskey

  origfmt:
    type: boolean?
    inputBinding:
      position: 30
      prefix: "--origfmt"
    doc: |
      Defline contains only original sequence name

  readids:
    type: boolean?
    inputBinding:
      position: 31
      prefix: "--readids"
    doc: |
      Append read id after spot id as 'accession.spot.readid' on defline

  helicos:
    type: boolean?
    inputBinding:
      position: 32
      prefix: "--helicos"
    doc: |
      Helicos style defline

  defline-seq:
    type: string?
    inputBinding:
      position: 33
      prefix: "--defline-seq"
    doc: |
      Defline format specification for sequence.

  defline-qual:
    type: string?
    inputBinding:
      position: 34
      prefix: "--defline-qual"
    doc: |
      Defline format specification for quality.
      <fmt> is string of characters and/or
      variables. The variables can be one of: $ac
      - accession, $si spot id, $sn spot
      name, $sg spot group (barcode), $sl spot
      length in bases, $ri read number, $rn
      read name, $rl read length in bases. '[]'
      could be used for an optional output: if
      all vars in [] yield empty values whole
      group is not printed. Empty value is empty
      string or for numeric variables. Ex:
      @$sn[_$rn]/$ri '_$rn' is omitted if name
      is empty

  disable-multithreading:
    type: boolean?
    inputBinding:
      position: 35
      prefix: "--disable-multithreading"
    doc: |
      disable multithreading

outputs:
  outputDir:
    type: Directory
    outputBinding:
      glob: "."

doc: |
  Usage
  fastq-dump [options] <path> [<path>...]

baseCommand: [fastq-dump]
