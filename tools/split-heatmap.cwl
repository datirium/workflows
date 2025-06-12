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
    doc: Homer generated heatmap file
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
    doc: Output filenames
outputs:
  refactored_heatmap_file:
    type:
    - File[]
    outputBinding:
      glob: '*.json'
    doc: Refactored heatmap file
baseCommand:
- python
- split.py
doc: |
  Splits Homer generated heatmap file into equal parts (by columns).
  Number of parts is defined by the length of the inputs.output_filename (if not array, then length is 1).
  inputs.output_filename are always updated to get the basename without extension + ".json"
label: split-heatmap
