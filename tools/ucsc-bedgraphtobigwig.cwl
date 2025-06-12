cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { var basename = inputs.bedgraph_file.location.split('/').slice(-1)[0]; var root = basename.split('.').slice(0,-1).join('.'); var ext = ".bigWig"; return (root == "")?basename+ext:root+ext; };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/ucscuserapps:v358
inputs:
  bedgraph_file:
    type: File
    inputBinding:
      position: 10
    doc: |
      Four column bedGraph file: <chrom> <start> <end> <value>
  chrom_length_file:
    type: File
    inputBinding:
      position: 11
    doc: |
      Two-column chromosome length file: <chromosome name> <size in bases>
  unc:
    type: boolean?
    inputBinding:
      position: 5
      prefix: -unc
    doc: |
      Disable compression
  items_per_slot:
    type: int?
    inputBinding:
      separate: false
      position: 6
      prefix: -itemsPerSlot=
    doc: |
      Number of data points bundled at lowest level. Default 1024
  block_size:
    type: int?
    inputBinding:
      separate: false
      position: 7
      prefix: -blockSize=
    doc: |
      Number of items to bundle in r-tree.  Default 256
  output_filename:
    type: string?
    inputBinding:
      position: 12
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
      If set, writes the output bigWig file to output_filename,
      otherwise generates filename from default_output_filename()
outputs:
  bigwig_file:
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
baseCommand:
- bedGraphToBigWig
doc: |
  Tool converts bedGraph to bigWig file.

  `default_output_filename` function returns filename for generated bigWig if `output_filename` is not provided.
  Default filename is generated on the base of `bedgraph_file` basename with the updated to `*.bigWig` extension.
label: ucsc-bedgraphtobigwig
