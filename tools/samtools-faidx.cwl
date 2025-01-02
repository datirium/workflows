cwlVersion: v1.0
class: CommandLineTool


requirements:
  - class: ResourceRequirement
    ramMin: 7620
    coresMin: 1
  - class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4
- class: InitialWorkDirRequirement
  listing: |
    ${
      return  [
                {
                  "entry": inputs.fasta_file,
                  "entryname": inputs.fasta_file.basename,
                  "writable": true
                }
              ]
    }

inputs:

  fasta_file:
    type: File
    inputBinding:
      position: 5
    doc: "Genome FASTA file"


outputs:

  fai_file:
    type: File
    outputBinding:
      glob: "*.fai"
    doc: "FAI index file"


baseCommand: ["samtools", "faidx"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

s:name: "samtools-faidx"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/samtools-faidx.cwl
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
  Generates FAI index file for input FASTA file
  Output file has the same basename, as input file, but with updated `.fai` extension. `samtools faidx` exports
  output file alognside the input file. To prevent tool from failing, `input_file` should be staged into output
  directory using `"writable": true`. Setting `writable: true` makes cwl-runner to make a copy of input file and
  mount it to docker container with `rw` mode as part of `--workdir` (if set to false, the file staged into output
  directory will be mounted to docker container separately with `ro` mode)

s:about: |
  Generates FAI index file for input FASTA file
