cwlVersion: v1.0
class: Workflow

inputs:
  string_input:
    type: string
    default: "test text"

outputs: []

steps:
  echo_step:
    run:
      class: CommandLineTool
      inputs:
        string_input:
          type: string
          inputBinding:
            position: 1
      outputs: []
      baseCommand: [echo]
    in:
      string_input: string_input
    out: []