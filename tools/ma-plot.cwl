cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/visualization:v0.0.6


inputs:

  diff_expr_file:
    type: File
    inputBinding:
      position: 5
    doc: |
      TSV file holding data for the plot

  x_axis_column:
    type: string
    inputBinding:
      position: 6
    doc: |
      Name of column in file for the plots x-axis (ex: "baseMean")

  y_axis_column:
    type: string
    inputBinding:
      position: 7
    doc: |
      Name of column in file for the plots y-axis (ex: "log2FoldChange")

  label_column:
    type: string
    inputBinding:
      position: 8
    doc: |
      Name of column in file for each data points 'name' (ex: "GeneId")


outputs:

  html_data:
    type: Directory
    outputBinding:
      glob: "./volcano_plot/MD-MA_plot"
    doc: |
      Directory html data for MA-plot

  html_file:
    type: File
    outputBinding:
      glob: "./volcano_plot/MD-MA_plot/html_data/index.html"
    doc: |
      HTML index file for MA-plot


baseCommand: ["ma_plot.sh"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "MA-plot"
s:name: "MA-plot"
s:alternateName: "Builds ma-plot from the DESeq output"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/ma-plot.cwl
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
  MA-plot

  Builds ma-plot from the DESeq output
