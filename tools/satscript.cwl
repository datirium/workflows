cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_prefix = function(ext) {
        ext = ext || "";
        if (inputs.output_prefix == null){
          let root = inputs.bam_file.basename.split('.').slice(0,-1).join('.');
          return (root == "")?inputs.bam_file.basename+ext:root+ext;
        } else {
          return inputs.output_prefix;
        }
    };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/satscript:v0.0.1


inputs:

  bam_file:
    type: File
    inputBinding:
      position: 5
      prefix: "-b"
    doc: "Path to the BAM file"

  macs_log:
    type: File
    inputBinding:
      position: 6
      prefix: "-m"
    doc: "Path to the MACS2 log file"

  percentage:
    type:
      - "null"
      - float
      - float[]
    inputBinding:
      position: 7
      prefix: "-p"
    doc: "Target percentage"

  output_prefix:
    type:
      - "null"
      - string
    inputBinding:
      position: 8
      prefix: "-o"
      valueFrom: $(get_output_prefix("_default_"))
    default: null
    doc: "Output filename prefix"

  output_suffixes:
    type:
      - "null"
      - string[]
    inputBinding:
      position: 9
      prefix: "-s"
    default: ["reads.png", "islands.png", "surface.png", "saturation.txt"]
    doc: |
      Output suffixes for reads, islands, surface and saturation files.

  res_dpi:
    type:
      - "null"
      - int
    inputBinding:
      position: 10
      prefix: "-r"
    doc: "Output picture file resolution, dpi"

outputs:

  reads_file:
    type: File
    outputBinding:
      glob: $("*"+inputs.output_suffixes[0])

  islands_file:
    type: File
    outputBinding:
      glob: $("*"+inputs.output_suffixes[1])

  surface_file:
    type: File
    outputBinding:
      glob: $("*"+inputs.output_suffixes[2])

  saturation_file:
    type: File
    outputBinding:
      glob: $("*"+inputs.output_suffixes[3])


baseCommand: ["SatScript"]