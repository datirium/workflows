cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ResourceRequirement
    ramMin: 7024
    coresMin: 1
  - class: InlineJavascriptRequirement
    expressionLib:
    - var default_output_filename = function() {
            return inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.') + ".fastxstat"
          };


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
      - "null"
      - boolean
    inputBinding:
      position: 5
      prefix: '-N'
    doc: |
      New output format (with more information per nucleotide/cycle).
      cycle (previously called 'column') = cycle number
      max-count
      For each nucleotide in the cycle (ALL/A/C/G/T/N):
          count   = number of bases found in this column.
          min     = Lowest quality score value found in this column.
          max     = Highest quality score value found in this column.
          sum     = Sum of quality score values for this column.
          mean    = Mean quality score value for this column.
          Q1	= 1st quartile quality score.
          med	= Median quality score.
          Q3	= 3rd quartile quality score.
          IQR	= Inter-Quartile range (Q3-Q1).
          lW	= 'Left-Whisker' value (for boxplotting).
          rW	= 'Right-Whisker' value (for boxplotting).

  output_filename:
    type:
      - "null"
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
    default: ""
    doc: |
      Output file to store generated statistics. If not provided - return from default_output_filename function

outputs:

  error_msg:
    type: File?
    outputBinding:
      glob: "error_msg.txt"

  error_report:
    type: File?
    outputBinding:
      glob: "error_report.txt"

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

baseCommand: [fastx_quality_stats]


label: "fastx-quality-stats"
doc: |
  Tool calculates statistics on the base of FASTQ file quality scores.
  If `output_filename` is not provided call function `default_output_filename` to return default output file name
  generated as `input_file` basename + `.fastxstat` extension.

  usage: fastx_quality_stats [-h] [-N] [-i INFILE] [-o OUTFILE]
  Part of FASTX Toolkit 0.0.14 by A. Gordon (assafgordon@gmail.com)
     [-h] = This helpful help screen.
     [-i INFILE]  = FASTQ input file. default is STDIN.
     [-o OUTFILE] = TEXT output file. default is STDOUT.
     [-N]         = New output format (with more information per nucleotide/cycle).
  The *OLD* output TEXT file will have the following fields (one row per column):
  	column	= column number (1 to 36 for a 36-cycles read solexa file)
  	count   = number of bases found in this column.
  	min     = Lowest quality score value found in this column.
  	max     = Highest quality score value found in this column.
  	sum     = Sum of quality score values for this column.
  	mean    = Mean quality score value for this column.
  	Q1	= 1st quartile quality score.
  	med	= Median quality score.
  	Q3	= 3rd quartile quality score.
  	IQR	= Inter-Quartile range (Q3-Q1).
  	lW	= 'Left-Whisker' value (for boxplotting).
  	rW	= 'Right-Whisker' value (for boxplotting).
  	A_Count	= Count of 'A' nucleotides found in this column.
  	C_Count	= Count of 'C' nucleotides found in this column.
  	G_Count	= Count of 'G' nucleotides found in this column.
  	T_Count	= Count of 'T' nucleotides found in this column.
  	N_Count = Count of 'N' nucleotides found in this column.
  	max-count = max. number of bases (in all cycles)
  The *NEW* output format:
  	cycle (previously called 'column') = cycle number
  	max-count
    For each nucleotide in the cycle (ALL/A/C/G/T/N):
  		count   = number of bases found in this column.
  		min     = Lowest quality score value found in this column.
  		max     = Highest quality score value found in this column.
  		sum     = Sum of quality score values for this column.
  		mean    = Mean quality score value for this column.
  		Q1	= 1st quartile quality score.
  		med	= Median quality score.
  		Q3	= 3rd quartile quality score.
  		IQR	= Inter-Quartile range (Q3-Q1).
  		lW	= 'Left-Whisker' value (for boxplotting).
  		rW	= 'Right-Whisker' value (for boxplotting).
