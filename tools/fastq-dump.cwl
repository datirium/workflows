cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    http_proxy: $(inputs.http_proxy?inputs.http_proxy:"")
    https_proxy: $(inputs.https_proxy?inputs.https_proxy:"")

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/fastqdwnld:v0.0.1


inputs:

  srr_id:
    type:
    - string
    - type: array
      items: string
    inputBinding:
      position: 60
    doc: |
      SRR identifiers

  split_files:
    type: boolean?
    inputBinding:
      position: 10
      prefix: "--split-files"
    doc: |
      Write reads into separate files. Read 
      number will be suffixed to the file name.  
      NOTE! The `--split-3` option is recommended. 
      In cases where not all spots have the same 
      number of reads, this option will produce 
      files that WILL CAUSE ERRORS in most programs 
      which process split pair fastq files. 

  split_3:
    type: boolean?
    inputBinding:
      position: 11
      prefix: "--split-3"
    doc: |
      3-way splitting for mate-pairs. For each 
      spot, if there are two biological reads 
      satisfying filter conditions, the first is 
      placed in the `*_1.fastq` file, and the 
      second is placed in the `*_2.fastq` file. If 
      there is only one biological read 
      satisfying the filter conditions, it is 
      placed in the `*.fastq` file.All other 
      reads in the spot are ignored.

  http_proxy:
    type: string?
    doc: |
      Optional HTTP proxy settings

  https_proxy:
    type: string?
    doc: |
      Optional HTTPS proxy settings

outputs:

  fastq_files:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*.gz"

  metadata_json:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*.json"

  metadata_xml:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*.xml"

  report_md:
    type: File
    outputBinding:
      glob: "report.md"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["sra_download.sh"]

stdout: sra_download_stdout.log
stderr: sra_download_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Fastq-Dump on Steroids"
s:name: "Fastq-Dump on Steroids"
s:alternateName: "Downloads FASTQ files from the provided SRR identifier"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/fastq-dump.cwl
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
  Fastq-Dump on Steroids

  Downloads FASTQ files from the provided SRR identifier


s:about: |
  Custom script to first prefetch raw data based on the
  provided SRR identifiers and them merge exported from
  then FASTQ file.