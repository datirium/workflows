cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      COMBINED="$0"
      function extract {
        FILE=$1
        COMBINED=$2
        T=`file -b "${FILE}" | awk '{print $1}'`
        case "${T}" in
          "bzip2"|"gzip"|"Zip")
            7z e -so "${FILE}" >> "${COMBINED}"
            ;;
          "ASCII")
            cat "${FILE}" >> "${COMBINED}" || true
            ;;
          *)
            echo "Error: file type unknown"
            rm -f "${COMBINED}"
            exit 1
        esac
      } 
      for FILE in "$@"; do
          echo "Extracting:" $FILE;
          extract "${FILE}" "${COMBINED}"
      done;
    inputBinding:
      position: 5
    doc: |
      Bash script to extract compressed file(s)

  output_filename:
    type: string
    inputBinding:
      position: 6
    doc: |
      Output filename for extracted and optionally
      merged file(s)

  file_to_extract:
    type:
    - File
    - type: array
      items: File
    inputBinding:
      position: 7
    doc: |
      Compressed file(s) to extract


outputs:

  extracted_file:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)


baseCommand: ["bash", "-c"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "extract-fastq"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/extract-7z.cwl
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
  Tool to decompress input file(s).
  If several files are provided, they will be concatenated in
  the order that corresponds to files in input.
  Bash script's logic:
  - check file type, decompress if needed, otherwise just cat
    the content of the file
  - return 1, if file type is not recognized


s:about: |
  Tool to decompress input file(s)
