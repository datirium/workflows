cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_prefix = function() {
        if (inputs.output_prefix) {
          return inputs.output_prefix;
        }
        var root = inputs.metrics_summary_report.basename.split('.').slice(0,-1).join('.');
        var suffix = "_stats";
        return (root == "")?inputs.metrics_summary_report.basename+suffix:root+suffix;
    };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scstats:v0.0.2


inputs:

  metrics_summary_report:
    type: File
    inputBinding:
      position: 6
      prefix: "--metrics"

  output_prefix:
    type: string?
    inputBinding:
      position: 7
      prefix: "--output"
      valueFrom: $(get_output_prefix())
    default: ""


outputs:

  collected_statistics_yaml:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".yaml")

  collected_statistics_tsv:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".tsv")

  collected_statistics_md:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".md")


baseCommand: ["cell_ranger_arc_count_stats.py"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


s:name: "Cell Ranger ARC Count Statistics"
label: "Cell Ranger ARC Count Statistics"
s:alternateName: "Collects statistics from Cell Ranger ARC Count experiment"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/collect-stats-sc-arc-count.cwl
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
  Cell Ranger ARC Count Statistics
  ================================
  
  Collects statistics from Cell Ranger ARC Count experiment


s:about: |
  Collects statistics from Cell Ranger ARC Count experiment
