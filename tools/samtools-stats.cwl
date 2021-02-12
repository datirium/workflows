cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
          if (inputs.output_filename == ""){
            var root = inputs.bambai_pair.basename.split('.').slice(0,-1).join('.');
            return (root == "")?inputs.bambai_pair.basename+".log":root+".log";
          } else {
            return inputs.output_filename;
          }
        };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/samtools:v1.4


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      samtools stats $0 > all_statistics.log
      cat all_statistics.log | grep SN | cut -f 2- > sn_section.log
      echo -e "insert size\tpairs total\tinward oriented pairs\toutward oriented pairs\tother pairs" > is_section.tsv
      cat all_statistics.log | grep ^IS | cut -f 2- | head -n 1501 >> is_section.tsv
      cat all_statistics.log
    inputBinding:
      position: 4
    doc: "samtools stats with optional filtering criteria passes as $1"

  bambai_pair:
    type: File
    inputBinding:
      position: 5
    secondaryFiles:
    - .bai
    doc: "Coordinate sorted BAM alignment and index BAI files"

  output_filename:
    type: string?
    default: ""
    doc: "Output file name"


outputs:

  log_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "BAM file statistics"

  ext_is_section:
    type: File
    outputBinding:
      glob: "is_section.tsv"
    doc: "BAM file statistics (IS section)"

  raw_total_sequences:
    type: int
    outputBinding:
      loadContents: true
      glob: sn_section.log
      outputEval: |
        ${
          var s = self[0].contents.substring(self[0].contents.indexOf("raw total sequences"));
          return parseInt(s.substring(s.indexOf("raw total sequences")+21, s.indexOf("\n")));
        }
    doc: "Raw Total Sequences"

  reads_mapped:
    type: int
    outputBinding:
      loadContents: true
      glob: sn_section.log
      outputEval: |
        ${
          var s = self[0].contents.substring(self[0].contents.indexOf("reads mapped"));
          return parseInt(s.substring(s.indexOf("reads mapped")+14, s.indexOf("\n")));
        }
    doc: "Reads Mapped"

  average_length:
    type: int
    outputBinding:
      loadContents: true
      glob: sn_section.log
      outputEval: |
        ${
          var s = self[0].contents.substring(self[0].contents.indexOf("average length"));
          return parseInt(s.substring(s.indexOf("average length")+16, s.indexOf("\n")));
        }
    doc: "Reads Average Length"


baseCommand: ["bash", "-c"]
stdout: $(default_output_filename())


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/samtools-metadata.yaml

s:name: "samtools-stats"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/samtools-stats.cwl
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
  Generates statistics for the input BAM file.

s:about: |
  Generates statistics for the input BAM file.
