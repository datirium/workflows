cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ResourceRequirement
  ramMin: 7024
  coresMin: 1
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { return inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.') + ".fastxstat" };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/fastx_toolkit:v0.0.14
inputs:
  input_file:
    type: File
    inputBinding:
      position: 10
      prefix: -i
    doc: |
      FASTA/Q input file. If FASTA file is given, only nucleotides distribution is calculated (there's no quality info)
  new_output_format:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 5
      prefix: -N
    doc: "New output format (with more information per nucleotide/cycle).\ncycle (previously called 'column') = cycle number\nmax-count\nFor each nucleotide in the cycle (ALL/A/C/G/T/N):\n    count   = number of bases found in this column.\n    min     = Lowest quality score value found in this column.\n    max     = Highest quality score value found in this column.\n    sum     = Sum of quality score values for this column.\n    mean    = Mean quality score value for this column.\n    Q1\t= 1st quartile quality score.\n    med\t= Median quality score.\n    Q3\t= 3rd quartile quality score.\n    IQR\t= Inter-Quartile range (Q3-Q1).\n    lW\t= 'Left-Whisker' value (for boxplotting).\n    rW\t= 'Right-Whisker' value (for boxplotting).\n"
  output_filename:
    type:
    - 'null'
    - string
    inputBinding:
      position: 11
      prefix: -o
      valueFrom: |
        ${
            if (self == ""){
              return default_output_filename();
            } else {
              return self;
            }
        }
    default: ''
    doc: |
      Output file to store generated statistics. If not provided - return from default_output_filename function
outputs:
  statistics_file:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == ""){
              return default_output_filename();
            } else {
              return inputs.output_filename;
            }
        }
    doc: Generated statistics file
baseCommand:
- fastx_quality_stats
doc: |
  Tool calculates statistics on the base of FASTQ file quality scores.
  If `output_filename` is not provided call function `default_output_filename` to return default output file name
  generated as `input_file` basename + `.fastxstat` extension.
label: fastx-quality-stats
