cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var ext = function() { return inputs.bai?'.bai':(inputs.csi?'.csi':'.bai'); };
  - var default_output_filename = function() { if (inputs.output_filename){ return inputs.output_filename; } else if (inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(-1)[0] == 'cram'){ return inputs.input_file.location.split('/').slice(-1)[0]+'.crai'; } else { return inputs.input_file.location.split('/').slice(-1)[0] + ext(); } };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4
inputs:
  input_file:
    type: File
    inputBinding:
      position: 8
    doc: |
      Input bam file
  csi_interval:
    type: int?
    inputBinding:
      position: 6
      prefix: -m
    doc: |
      Minimum interval size for CSI indices to 2^INT [14]
  threads:
    type: int?
    inputBinding:
      position: 7
      prefix: -@
    doc: |
      Number of threads
  output_filename:
    type: string?
    doc: |
      Output filename
  csi:
    type: boolean?
    doc: |
      Generate CSI-format index for BAM files
  bai:
    type: boolean?
    doc: |
      Generate BAI-format index for BAM files [default]
outputs:
  index_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: Index file
baseCommand:
- samtools
- index
arguments:
- valueFrom: $(inputs.bai?'-b':(inputs.csi?'-c':'-b'))
  position: 5
- valueFrom: $(default_output_filename())
  position: 9
doc: |
  Tool to index input BAM/CRAM file set as `input_file`.

  `default_output_filename` function returns filename for generated index file. If `output_filename` is set, return
  it unchanged. In opposite case, index file basename is set to be equal to `input_file` basename and extension
  depends on `input_file` extension and `csi`/`bai` flags. If `input_file` is `*.cram` - produce CRAI, if not -
  evaluate inputs `csi` and `bai` and produce index according to following logics (function `ext`):
    `csi` &&  `bai`  => BAI
   !`csi` && !`bai ` => BAI
    `csi` && !`bai ` => CSI
  If `input_file` has `*.cram` extension, flags `csi` and `bai` are ignored. `output_filename` extension doesn't
  influence on output index file type.
label: samtools-index
