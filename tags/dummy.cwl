cwlVersion: v1.0
class: Workflow

inputs:
  string_input:
    type: string
    default: "foobar"

outputs: []

steps:
  dummy_step:
    run:
      class: ExpressionTool
      requirements:
        - class: InlineJavascriptRequirement
      inputs: []
      outputs: []
      expression: ${return null}
    in: []
    out: []

