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
      prefix: "--features"
    doc: |
      CSV/TSV feature files to be merged

  feature_aliases:
    type:
      - "null"
      - string
      - string[]
    inputBinding:
      prefix: "--aliases"
    doc: |
      Unique aliases for feature files to be used as names for the --report columns
      in the merged file.
      Default: basenames of files provided in --features without extensions

  mergeby:
    type:
      - "null"
      - string
      - string[]
    inputBinding:
      prefix: "--mergeby"
    doc: |
      Column names to merge feature files by.
      Default: RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output file prefix.
      Default: merged_


outputs:

  merged_file:
    type: File
    outputBinding:
      glob: "*report.tsv"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: [run_merge.R]


stdout: merge_stdout.log
stderr: merge_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Feature merge - merges feature files based on the specified columns"
s:name: "Feature merge - merges feature files based on the specified columns"
s:alternateName: "Feature merge - merges feature files based on the specified columns"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/feature-merge.cwl
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
  Tool merges input feature files based on the columns provided in --mergeby
  input. All input feature CSV/TSV files should have the header (case-sensitive)
  Format of the input files is identified based on file's extension
  *.csv - CSV
  *.tsv - TSV
  Otherwise used CSV by default

  The output file's rows order corresponds to the rows order of the first
  CSV/TSV feature file. Output is always saved in TSV format.

  Output file includes only rows intersected by column names set in --mergeby.
  Output file includes only columns set in --mergeby and --report parameters.
  Column set in the --report parameter is renamed based on the --aliases or 
  basenames of the --features files.

s:about: |
  usage: run_merge.R
        [-h] -f FEATURES [FEATURES ...] [-a [ALIASES [ALIASES ...]]]
        [-m [MERGEBY [MERGEBY ...]]] [-o OUTPUT]

  Merge feature files based on the specified columns

  optional arguments:
    -h, --help            show this help message and exit
    -f FEATURES [FEATURES ...], --features FEATURES [FEATURES ...]
                          CSV/TSV feature files to be merged
    -a [ALIASES [ALIASES ...]], --aliases [ALIASES [ALIASES ...]]
                          Unique aliases for feature files to be used as
                          prefixes for not --mergeby columns in the merged file.
                          Default: basenames of files provided in --features
                          without extensions
    -m [MERGEBY [MERGEBY ...]], --mergeby [MERGEBY [MERGEBY ...]]
                          Column names to merge feature files by. Default:
                          RefseqId, GeneId, Chrom, TxStart, TxEnd, Strand
    -o OUTPUT, --output OUTPUT
                          Output file prefix. Default: merged_
