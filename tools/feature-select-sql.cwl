cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap-merge:v0.0.2


inputs:

  feature_file:
    type: File
    inputBinding:
      prefix: "--raw"
    doc: |
      Input CSV/TSV file to be filtered by SQL query

  sql_query:
    type: string
    inputBinding:
      prefix: "--sql"
    doc: |
      SQL query parameters that will be appended to `SELECT *` statement

  columns:
    type: string?
    inputBinding:
      prefix: "--columns"
    doc: |
      Comma-separated list of column to be print in the output.
      Default: all

  header:
    type: boolean?
    inputBinding:
      prefix: "--header"
    doc: |
      Print header in the output.
      Default: false

  output_filename:
    type: string?
    inputBinding:
      prefix: "--output"
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


baseCommand: [sql_select.R]


stdout: sql_select_stdout.log
stderr: sql_select_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Feature select - filters TSV/CSV files based on the provided SQL query parameters"
s:name: "Feature select - filters TSV/CSV files based on the provided SQL query parameters"
s:alternateName: "Feature select - filters TSV/CSV files based on the provided SQL query parameters"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/feature-select-sql.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Tool filters input TSV/CSV file based on the provided SQL query.
  Value set in sql_query parameter will be appended to the
  "SELECT * FROM raw_data WHERE". Column names in sql_query are case
  insensitive. Format of the input files is identified based on file
  extension
  *.csv - CSV
  *.tsv - TSV
  Otherwise used CSV by default. Output is always saved in TSV format.


s:about: |
  usage: sql_select.R [-h] -r RAW -s SQL [-o OUTPUT]

  Filter CSV/TSV file based on the provided SQL query

  optional arguments:
    -h, --help            show this help message and exit
    -r RAW, --raw RAW     Input CSV/TSV file
    -s SQL, --sql SQL     SQL query parameters that will be appended to `SELECT
                          *` statement
    -o OUTPUT, --output OUTPUT
                          Filename for filtered output TSV file
