cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var get_root = function(basename) {
          return basename.split('.').slice(0,1).join('.');
      };


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  fastq_file:
    type:
    - File
    - type: array
      items: File
    label: "FASTQ file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, optionally compressed"


outputs:

  fastqc_report:
    type: File
    outputSource: rename_fastqc_report/target_file
    label: "FastqQC HTML report"
    doc: "FastqQC HTML report"
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"


steps:

  extract_fastq:
    run: ../tools/extract-fastq.cwl
    in:
      compressed_file: fastq_file
    out: [fastq_file]

  run_fastqc:
    run: ../tools/fastqc.cwl
    in:
      reads_file: extract_fastq/fastq_file
    out:
      - html_file

  rename_fastqc_report:
    run: ../tools/rename.cwl
    in:
      source_file: run_fastqc/html_file
      target_filename:
        source: fastq_file
        valueFrom: $(get_root(self.basename)+"_fastqc_report.html")
    out: [target_file]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "FastQC - a quality control tool for high throughput sequence data"
label: "FastQC - a quality control tool for high throughput sequence data"
s:alternateName: "FastQC - a quality control tool for high throughput sequence data"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/fastqc.cwl
s:codeRepository: https://github.com/datirium/workflows
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


# doc:
#   $include: ../descriptions/fastqc.md


doc: |
  FastQC - a quality control tool for high throughput sequence data
  =====================================

  FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming
  from high throughput sequencing pipelines. It provides a modular set of analyses which you can use
  to give a quick impression of whether your data has any problems of which you should be aware before
  doing any further analysis.

  The main functions of FastQC are:

  - Import of data from FastQ files (any variant)
  - Providing a quick overview to tell you in which areas there may be problems
  - Summary graphs and tables to quickly assess your data
  - Export of results to an HTML based permanent report
  - Offline operation to allow automated generation of reports without running the interactive application