cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap-merge:v0.0.3
inputs:
  feature_file:
    type: File
    inputBinding:
      prefix: --raw
    doc: |
      Input CSV/TSV file to be filtered by SQL query
  sql_query:
    type: string
    inputBinding:
      prefix: --sql
    doc: |
      SQL query parameters that will be appended to `SELECT *` statement
  columns:
    type: string?
    inputBinding:
      prefix: --columns
    doc: |
      Comma-separated list of column to be print in the output.
      Default: all
  header:
    type: boolean?
    inputBinding:
      prefix: --header
    doc: |
      Print header in the output.
      Default: false
  output_filename:
    type: string?
    inputBinding:
      prefix: --output
    doc: |
      Filename for filtered output TSV file.
      Default: rootname(feature_file)_filtered.tsv
outputs:
  filtered_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename?inputs.output_filename:"*_filtered.tsv")
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- sql_select.R
stdout: sql_select_stdout.log
stderr: sql_select_stderr.log
label: Feature select - filters TSV/CSV files based on the provided SQL query parameters
doc: |
  Tool filters input TSV/CSV file based on the provided SQL query.
  Value set in sql_query parameter will be appended to the
  "SELECT * FROM raw_data WHERE". Column names in sql_query are case
  insensitive. Format of the input files is identified based on file
  extension
  *.csv - CSV
  *.tsv - TSV
  Otherwise used CSV by default. Output is always saved in TSV format.
