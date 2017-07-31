cwlVersion: v1.0
class: Workflow

inputs:
  cells:
    type: string
    label: "Cells"
  conditions:
    type: string
    label: "Conditions"
  alias:
    type: string
    label: "Experiment short name"
  catalog:
    type: string?
    label: "Catalog #"
  description:
    type: string?
    'sd:type': 'text'
    label: "Catalog #"


outputs: []
steps: []

