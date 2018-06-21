cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function(ext) {
        let root = inputs.input_filename.basename.split('.').slice(0,-1).join('.');
        return (root == "")?inputs.input_filename.basename+ext:root+ext;
    };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/iaintersect:v0.0.2

inputs:

  input_filename:
    type:
      - File
    inputBinding:
      position: 1
      prefix: --in=
      separate: false
    doc: |
      Input filename with MACS2 peak calling results, tsv

  annotation_filename:
    type:
      - File
    inputBinding:
      position: 2
      prefix: --a=
      separate: false
    doc: |
      Annotation file, tsv

  output_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 3
      prefix: --out=
      separate: false
      valueFrom: |
        ${
            if (self == null){
              return default_output_filename('_iaintersect.tsv');
            } else {
              return self;
            }
        }
    default: null
    doc: |
      Base output file name, tsv

  log_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 4
      prefix: --log=
      separate: false
      valueFrom: |
        ${
            if (self == null){
              return default_output_filename('_iaintersect.log');
            } else {
              return self;
            }
        }
    default: null
    doc: |
      Log filename

  promoter_bp:
    type:
      - "null"
      - int
    inputBinding:
      position: 5
      prefix: --promoter=
      separate: false
    doc: |
      Promoter region around TSS, base pairs

  upstream_bp:
    type:
      - "null"
      - int
    inputBinding:
      position: 6
      prefix: --upstream=
      separate: false
    doc: |
      Upstream region before promoter, base pairs

  ignore_chr:
    type:
      - "null"
      - string
    inputBinding:
      position: 7
      prefix: --sam_ignorechr=
      separate: false
    doc: |
      The chromosome to be ignored, string

outputs:
  log_file:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.log_filename == null){
            return default_output_filename('_iaintersect.log');
          } else {
            return inputs.log_filename;
          }
        }

  result_file:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.output_filename == null){
            return default_output_filename('_iaintersect.tsv');
          } else {
            return inputs.output_filename;
          }
        }


baseCommand: [iaintersect]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/iaintersect-metadata.yaml

s:name: "iaintersect"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/iaintersect.cwl
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
  Tool assigns each peak obtained from MACS2 to a gene and region (upstream, promoter, exon, intron, intergenic)

  `default_output_filename` function returns output filename with sufix set as `ext` argument. Function is called when
  either `output_filename` or `log_filename` inputs are not provided.

s:about: |
  Usage:
    iaintersect [options] --in=pathToFile --a=pathtoFile --out=pathToFile
    --a                        	Tab-separated annotation file
    --in                       	Input filename with MACS2 peak calling results, xls
    --log                      	Log filename (default is ./logfile_def.log)
    --out                      	Base output file name
    --promoter                 	Promoter region around TSS in bp
    --sam_ignorechr            	The chromosome to be ignored
    --upstream                 	Upstream region before promoter in bp
