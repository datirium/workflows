cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InitialWorkDirRequirement
  listing: |
    ${
      return [
        {"class": "Directory",
         "basename": "default",
         "listing": [inputs.bam_file],
         "writable": true}
      ]
    }
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_folder = function() { if (inputs.output_folder){ return inputs.output_folder.replace(/\t|\s|\[|\]|\>|\<|,|\./g, "_"); } else { return inputs.bam_file.basename.split('.')[0]; } };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/homer:v0.0.2
inputs:
  bam_file:
    type: File
    doc: Alignment file, BAM
  output_folder:
    type: string?
    doc: Name of the directory to save outputs
  fragment_size:
    type:
    - 'null'
    - int
    - string
    inputBinding:
      position: 5
      prefix: -fragLength
    doc: |
      Set fragment size.
      By default is estimated as if it was single end ChIP-Seq experiment.
      Possible values:
        "#" - int value to be used as fragment size
        "given" - use read lengths
        "pe" - calculate from paired end read coordinates
  total_reads:
    type:
    - 'null'
    - int
    - string
    inputBinding:
      position: 6
      prefix: -totalReads
    doc: |
      Set total reads number for downstream normalization.
      Default: autocalculated, equal to uniquely mapped reads number
      Possible values:
        "#" - int value to be used as total reads number
        "all" - autocalculated, equal to uniquely + multi mapped reads number
  min_length:
    type: int?
    inputBinding:
      position: 7
      prefix: -minlen
    doc: |
      Discard reads smaller then
  max_length:
    type: int?
    inputBinding:
      position: 8
      prefix: -maxlen
    doc: |
      Discard reads bigger then
outputs:
  output_tag_folder:
    type: Directory
    outputBinding:
      glob: $(default_output_folder())
    doc: Tag directory
baseCommand:
- makeTagDirectory
arguments:
- valueFrom: $(default_output_folder())
- valueFrom: $("default/" + inputs.bam_file.basename)
doc: |
  Tool runs makeTagDirectory that basically parses through the alignment file and splits the tags into separate
  files based on the chromosomes.

  Multiple alignment files are not supported. Alignment file's format is restricted to be only BAM.

  Output is placed in a folder with the name derived from the input BAM file's basename.

  Skipped arguments:

    Rely on the default value:
      -format           - format will be autodetected
      -precision        - the default value is used

    Not required general functionality:
      -d
      -single
      -force5th
      -t
      -flip
      -tbp

    Not required GC-bias options:
      -genome
      -checkGC
      -normGC
      -normFixedOligo
      -minNormRatio
      -maxNormRatio
      -iterNorm
      -filterReads

    Not required HiC options:
      -removePEbg
      -restrictionSite
      -removeSpikes
      -bowtiePE
      -directional
label: homer-make-tag-directory
