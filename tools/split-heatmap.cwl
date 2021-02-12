cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
    - entryname: split.py
      entry: |
          #!/usr/bin/env python
          import os
          import sys
          import pandas as pd
          heatmap_file = sys.argv[1]
          output_files = [os.path.splitext(os.path.basename(file))[0]+".json" for file in sys.argv[2:]]
          print (output_files)
          data = pd.read_table(heatmap_file, index_col=0)
          size = int(len(data.columns)/len(output_files))
          data = [data[c] for c in [data.columns[i*size:(i+1)*size] for i in range(len(output_files))]]
          for idx,d in enumerate(data):
              d.columns = data[0].columns
              d.to_json(output_files[idx], orient="split")


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/python-pandas:v0.0.1


inputs:

  heatmap_file:
    type: File
    inputBinding:
      position: 5
    doc: "Homer generated heatmap file"

  output_filename:
    type:
      - string
      - File
      - type: array
        items:
          - string
          - File
    inputBinding:
      position: 6
    doc: "Output filenames"


outputs:

  refactored_heatmap_file:
    type:
      - File[]
    outputBinding:
      glob: "*.json"
    doc: "Refactored heatmap file"


baseCommand: ["python", "split.py"]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "split-heatmap"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/split-heatmap.cwl
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
  Splits Homer generated heatmap file into equal parts (by columns).
  Number of parts is defined by the length of the inputs.output_filename (if not array, then length is 1).
  inputs.output_filename are always updated to get the basename without extension + ".json"

s:about: |
  Splits Homer generated heatmap file into equal parts (by columns).
  Number of parts is defined by the length of the inputs.output_filename (if it's array)
  inputs.output_filename are always updated to get only the basename without extension + ".json"