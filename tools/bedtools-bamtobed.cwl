cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { if (inputs.output_filename == ""){ var root = inputs.bam_file.basename.split('.').slice(0,-1).join('.'); return (root == "")?inputs.bam_file.basename+".bed":root+".bed"; } else { return inputs.output_filename; } };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/bedtools2:v2.26.0
inputs:
  bam_file:
    type: File
    inputBinding:
      position: 5
      prefix: -i
    doc: Input BAM file (not necessary sorted or indexed)
  output_filename:
    type: string?
    default: ''
    doc: Output BED filename
outputs:
  bed_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: Sequences file
baseCommand:
- bedtools
- bamtobed
stdout: $(default_output_filename())
doc: |
  Converts BAM to BED. All Options are not implemented.
label: bedtools-bamtobed
