cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_prefix = function(ext) {
        ext = ext || "";
        if (inputs.output_prefix == null){
          let root = inputs.bambai_pair.basename.split('.').slice(0,-1).join('.');
          return (root == "")?inputs.bambai_pair.basename+ext:root+ext;
        } else {
          return inputs.output_prefix;
        }
    };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/plugin-plot-dna:v0.0.1


inputs:

  islands_file:
    type: File
    inputBinding:
      position: 5
      prefix: "-i"
    doc: "Islands file (Iaintersect TSV or MACS2 XLS output)"

  bambai_pair:
    type: File
    inputBinding:
      position: 6
      prefix: "-b"
    secondaryFiles:
      - .bai
    doc: "Indexed BAM + BAI files"

  output_prefix:
    type: string?
    inputBinding:
      position: 9
      prefix: "-o"
      valueFrom: $(get_output_prefix("_default_"))
    default: null
    doc: "Output file prefix"


outputs:

  png_file:
    type: File[]
    outputBinding:
      glob: "*.png"

  pileup_file:
    type: File
    outputBinding:
      glob: "*pileup.tsv"

  length_file:
    type: File
    outputBinding:
      glob: "*length.tsv"

  reads_file:
    type: File
    outputBinding:
      glob: "*reads.tsv"

baseCommand: ["plot_dna.R"]