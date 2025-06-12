cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap-merge:v0.0.3
inputs:
  feature_files:
    type:
    - File
    - File[]
    inputBinding:
      prefix: --features
    doc: |
      CSV/TSV feature files to be merged
  feature_aliases:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      prefix: --aliases
    doc: |
      Unique aliases for feature files to be used as names for the --report columns
      in the merged file.
      Default: basenames of files provided in --features without extensions
  mergeby:
    type:
    - 'null'
    - string
    - string[]
    inputBinding:
      prefix: --mergeby
    doc: |
      Column names to merge feature files by.
      Default: RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand
  output_prefix:
    type: string?
    inputBinding:
      prefix: --output
    doc: |
      Output file prefix.
      Default: merged_
outputs:
  merged_file:
    type: File
    outputBinding:
      glob: '*report.tsv'
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- run_merge.R
stdout: merge_stdout.log
stderr: merge_stderr.log
label: Feature merge - merges feature files based on the specified columns
doc: "Tool merges input feature files based on the columns provided in --mergeby\ninput. All input feature CSV/TSV files should have the header (case-sensitive)\nFormat of the input files is identified based on file's extension\n*.csv - CSV\n*.tsv - TSV\nOtherwise used CSV by default\n\nThe output file's rows order corresponds to the rows order of the first\nCSV/TSV feature file. Output is always saved in TSV format.\n\nOutput file includes only rows intersected by column names set in --mergeby.\nOutput file includes only columns set in --mergeby and --report parameters.\nColumn set in the --report parameter is renamed based on the --aliases or \nbasenames of the --features files.\n"
