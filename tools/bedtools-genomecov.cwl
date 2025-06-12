cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { var ext = (inputs.depth == "-bg" || inputs.depth == "-bga")?".bedGraph":".tab"; return inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.') + ext; };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0
inputs:
  input_file:
    type: File
    inputBinding:
      position: 16
      valueFrom: |
        ${
          var prefix = ((/.*\.bam$/i).test(inputs.input_file.path))?'-ibam':'-i';
          return [prefix, inputs.input_file.path];
        }
    doc: |
      The input file can be in BAM format (Note: BAM must be sorted by position) or <bed/gff/vcf>.
      Prefix is selected on the base of input file extension
  chrom_length_file:
    type: File?
    inputBinding:
      position: 17
      prefix: -g
    doc: |
      Input genome file. Needed only when -i flag. The genome file is tab delimited <chromName><TAB><chromSize>
  depth:
    type:
    - 'null'
    - type: enum
      name: depth
      symbols:
      - -bg
      - -bga
      - -d
      - -dz
    inputBinding:
      position: 5
    doc: |
      Report the depth type. By default, bedtools genomecov will compute a histogram of coverage
      for the genome file provided (intputs.chrom_length_file)
  scale:
    type: float?
    inputBinding:
      position: 6
      prefix: -scale
    doc: |
      Scale the coverage by a constant factor.
      Each coverage value is multiplied by this factor before being reported.
      Useful for normalizing coverage by, e.g., reads per million (RPM).
      - Default is 1.0; i.e., unscaled.
      - (FLOAT)
  mapped_reads_number:
    type: int?
    inputBinding:
      position: 7
      prefix: -scale
      valueFrom: |
        ${
          if (inputs.scale){
            return null;
          } else if (inputs.mapped_reads_number) {
            return 1000000/inputs.mapped_reads_number;
          } else {
            return null;
          }
        }
    doc: |
      Optional parameter to calculate scale as 1000000/mapped_reads_number if inputs.scale is not provided
  split:
    type: boolean?
    inputBinding:
      position: 8
      prefix: -split
    doc: |
      treat "split" BAM or BED12 entries as distinct BED intervals.
      when computing coverage.
      For BAM files, this uses the CIGAR "N" and "D" operations
      to infer the blocks for computing coverage.
      For BED12 files, this uses the BlockCount, BlockStarts, and BlockEnds
      fields (i.e., columns 10,11,12).
  strand:
    type: string?
    inputBinding:
      position: 9
      prefix: -strand
    doc: |
      Calculate coverage of intervals from a specific strand.
      With BED files, requires at least 6 columns (strand is column 6).
      - (STRING): can be + or -
  pairchip:
    type: boolean?
    inputBinding:
      position: 10
      prefix: -pc
    doc: |
      pair-end chip seq experiment
  du:
    type: boolean?
    inputBinding:
      position: 11
      prefix: -du
    doc: |
      Change strand af the mate read (so both reads from the same strand) useful for strand specific.
      Works for BAM files only
  fragment_size:
    type: int?
    inputBinding:
      position: 12
      prefix: -fs
    doc: |
      Set fixed fragment size
  max:
    type: int?
    inputBinding:
      position: 13
      prefix: -max
    doc: |
      Combine all positions with a depth >= max into
      a single bin in the histogram. Irrelevant
      for -d and -bedGraph
      - (INTEGER)
  m5:
    type: boolean?
    inputBinding:
      position: 14
      prefix: '-5'
    doc: |
      Calculate coverage of 5-end positions (instead of entire interval)
  m3:
    type: boolean?
    inputBinding:
      position: 15
      prefix: '-3'
    doc: |
      Calculate coverage of 3-end positions (instead of entire interval)
  output_filename:
    type: string?
    doc: |
      Name for generated output file
outputs:
  genome_coverage_file:
    type: stdout
    doc: |
      Generated genome coverage output file
stdout: |
  ${
      return inputs.output_filename ? inputs.output_filename : default_output_filename();
  }
baseCommand:
- bedtools
- genomecov
doc: |
  Tool calculates genome coverage from input bam/bed/gff/vcf using `bedtools genomecov`

  Depending on `input_file` extension additional prefix is used: if `*.bam` use `-ibam`, else use `-i`.

  `scale` and `mapped_reads_number` inputs result in the same parameter `-scale`. If `scale` is not provided, check if
  `mapped_reads_number` is not null and calculate `-scale` as `1000000/mapped_reads_number`. If both inputs are
  null, `bedtools genomecov` will use its default scaling value.

  `default_output_filename` function returns default output filename and is used when `output_filename` is not provided.
  Default output file extention is `.tab`. If bedGraph should be generated (check flags `inputs.depth`), extension is
  updated to `.bedGraph`. Default basename of the output file is generated on the base of `input_file` basename.
label: bedtools-genomecov
