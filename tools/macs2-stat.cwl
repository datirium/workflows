cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.2

inputs:
  fragments_old:
    type: int
  fragments_new:
    type: int
  islands_old:
    type: int
  islands_new:
    type: int
  trigger:
    type: boolean
  output_filename:
    type: string?
    default: fragment_stat.tsv

outputs:
  fragment_stat_file:
    type: stdout

stdout: $(inputs.output_filename)

baseCommand: [echo]
arguments:
  - valueFrom: |
      ${
        if (inputs.trigger){
          return [inputs.fragments_new, inputs.fragments_old, inputs.islands_new];
        } else {
          return [inputs.fragments_old, inputs.fragments_old, inputs.islands_old];
        }
      }
    position: 1

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "macs2-stat"
s:downloadUrl: https://raw.githubusercontent.com/SciDAP/workflows/master/tools/macs2-stat.cwl
s:codeRepository: https://github.com/SciDAP/workflows
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
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  Returns tab-separated file with FRAGMENT, FRAGMENTE, ISLANDS value.
  Values are set depending on inputs.trigger value.
  Logic corresponds to https://github.com/cincinnati-childrens-hospital/
  biowardrobe/blob/d8eca2d85f720df2ea533e3b14c14984513eb8f4/scripts/RunDNA.py#L242

s:about: >
  Forms output file to be easily uploaded to BioWardrobe DB
