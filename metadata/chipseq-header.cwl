cwlVersion: v1.0
class: Workflow

inputs:

  cells:
    type: string
    label: "Cells"
    sd:preview:
      position: 1

  conditions:
    type: string
    label: "Conditions"
    sd:preview:
      position: 3

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 2

  catalog:
    type: string?
    label: "Catalog #"

outputs: []
steps: []

