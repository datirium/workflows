cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { if (inputs.output_filename == ""){ var root = inputs.intervals_file.basename.split('.').slice(0,-1).join('.'); return (root == "")?inputs.intervals_file.basename+".fa":root+".fa"; } else { return inputs.output_filename; } };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0
inputs:
  genome_fasta_file:
    type: File
    secondaryFiles:
    - .fai
    inputBinding:
      position: 5
      prefix: -fi
    doc: Genome file in FASTA format
  intervals_file:
    type: File
    inputBinding:
      position: 6
      prefix: -bed
    doc: Intervals file defined in a BED/GFF/VCF format
  output_filename:
    type: string?
    default: ''
    doc: Output file name
outputs:
  sequences_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: Sequences file
baseCommand:
- bedtools
- getfasta
stdout: $(default_output_filename())
doc: |
  Extracts sequences from a FASTA file for each of the intervals defined in a BED/GFF/VCF file. Only selected parameters are implemented.
label: bedtools-getfasta
