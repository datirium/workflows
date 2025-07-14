cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/umi_tools:0.5.5--py36h470a237_0
inputs:
  fastq_file_1:
    type: File
    inputBinding:
      prefix: -I
    doc: Input FASTQ file 1
  fastq_file_2:
    type: File?
    inputBinding:
      prefix: --read2-in=
      separate: false
    doc: Input FASTQ file 2
  extract_method:
    type:
      type: enum
      symbols:
      - string
      - regex
    inputBinding:
      prefix: --extract-method=
      separate: false
    default: regex
    doc: How to extract the umi +/- cell barcodes, choose from string or regex
  bc_pattern_1:
    type: string
    inputBinding:
      prefix: --bc-pattern=
      separate: false
    doc: Barcode pattern 1
  bc_pattern_2:
    type: string?
    inputBinding:
      prefix: --bc-pattern2=
      separate: false
    doc: Barcode pattern 2
  output_filename_1:
    type: string
    inputBinding:
      prefix: -S
    doc: Output filename 1
  output_filename_2:
    type: string?
    inputBinding:
      prefix: --read2-out=
      separate: false
    doc: Output filename 2
outputs:
  umi_fastq_file_1:
    type: File
    outputBinding:
      glob: $(inputs.output_filename_1)
  umi_fastq_file_2:
    type: File?
    outputBinding:
      glob: $(inputs.output_filename_2)
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- umi_tools
- extract
stdout: umi_tools_extract_stdout.log
stderr: umi_tools_extract_stderr.log
doc: |
  Extract UMI barcode from a read and add it to the read name, leaving
  any sample barcode in place. Can deal with paired end reads and UMIs
  split across the paired ends. Can also optionally extract cell
  barcodes and append these to the read name also. See the section below
  for an explanation for how to encode the barcode pattern(s) to
  specficy the position of the UMI +/- cell barcode.
label: umi-tools-extract
