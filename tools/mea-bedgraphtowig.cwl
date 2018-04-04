cwlVersion: v1.0
class: CommandLineTool


requirements:
- $import: ./metadata/envvar-global.yml
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_filename = function() {
      if (inputs.output_filename == null){
        let ext = '.wig';
        let root = inputs.bedgraph_file.basename.split('.').slice(0,-1).join('.');
        return (root == "")?inputs.bedgraph_file.basename+ext:root+ext;
      } else {
        return inputs.output_filename;
      }
    };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      awk 'BEGIN {
          print "track type=wiggle_0"
      }
      NF == 4 {
          print "fixedStep chrom="$1" start="$2+1" step=1 span=1"
          for(i = 0; i < $3-$2; i++) {
              print $4
          }
      }' "$0"
    inputBinding:
      position: 1

  bedgraph_file:
    type: File
    inputBinding:
      position: 2

  output_filename:
    type: string?
    default: null


outputs:

  wig_file:
    type: File
    outputBinding:
      glob: $(get_output_filename())


stdout: $(get_output_filename())


baseCommand: [bash, '-c']


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "mea-bedgraphtowig"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/mea-bedgraphtowig.cwl
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
        s:email: mailto:michael.kotliar@cchmc.org
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898

doc: |
  AWK script to convert input bedGraph file into Wig with fixedStep format

  If input `output_filename` is not set, call `get_output_filename` function to return default output filename
  based on `bedgraph_file` basename with `wig` extension.

s:about: |
  AWK script is taken from
  https://github.com/julienrichardalbert/MEA/blob/e3de228734bafd957cc2072dd8a6a0e84d554724/src/scripts/createTracks.sh#L92