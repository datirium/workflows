cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var ext = function() {
      return inputs.bai?'.bai':(inputs.csi?'.csi':'.bai');
    };

  - var default_output_filename = function() {
      if (inputs.output_filename){
        return inputs.output_filename;
      } else if (inputs.input_file.location.split('/').slice(-1)[0].split('.').slice(-1)[0] == 'cram'){
        return inputs.input_file.location.split('/').slice(-1)[0]+'.crai';
      } else {
        return inputs.input_file.location.split('/').slice(-1)[0] + ext();
      }
    };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4


inputs:
  input_file:
    type: File
    inputBinding:
      position: 8
    doc: |
      Input bam file

  csi_interval:
    type: int?
    inputBinding:
      position: 6
      prefix: -m
    doc: |
      Minimum interval size for CSI indices to 2^INT [14]

  threads:
    type: int?
    inputBinding:
      position: 7
      prefix: -@
    doc: |
      Number of threads

  output_filename:
    type: string?
    doc: |
      Output filename

  csi:
    type: boolean?
    doc: |
      Generate CSI-format index for BAM files

  bai:
    type: boolean?
    doc: |
      Generate BAI-format index for BAM files [default]

outputs:
  index_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: Index file

baseCommand: [samtools, index]
arguments:
- valueFrom: $(inputs.bai?'-b':(inputs.csi?'-c':'-b'))
  position: 5
- valueFrom: $(default_output_filename())
  position: 9


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

s:name: "samtools-index"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/samtools-index.cwl
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
  Tool to index input BAM/CRAM file set as `input_file`.

  `default_output_filename` function returns filename for generated index file. If `output_filename` is set, return
  it unchanged. In opposite case, index file basename is set to be equal to `input_file` basename and extension
  depends on `input_file` extension and `csi`/`bai` flags. If `input_file` is `*.cram` - produce CRAI, if not -
  evaluate inputs `csi` and `bai` and produce index according to following logics (function `ext`):
    `csi` &&  `bai`  => BAI
   !`csi` && !`bai ` => BAI
    `csi` && !`bai ` => CSI
  If `input_file` has `*.cram` extension, flags `csi` and `bai` are ignored. `output_filename` extension doesn't
  influence on output index file type.

s:about: |
  Usage: samtools index [-bc] [-m INT] <in.bam> [out.index]
  Options:
    -b       Generate BAI-format index for BAM files [default]
    -c       Generate CSI-format index for BAM files
    -m INT   Set minimum interval size for CSI indices to 2^INT [14]
    -@ INT   Sets the number of threads [none]