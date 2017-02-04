cwlVersion: v1.0
class: ExpressionTool

requirements:
  - class: InlineJavascriptRequirement

id: "fastqc_results_trigger"
inputs:
  summary:
    type: File
    inputBinding:
      loadContents: true
outputs:
  trigger: boolean
expression: |
  ${
    var criteria_array = inputs.summary.contents.match(/.*Per base sequence quality.*|.*Per sequence quality scores.*|.*Overrepresented sequences.*/g);
    if (criteria_array.length > 0){
      if (criteria_array.toString().match(/FAIL/g) != null){
        return { "trigger": true };
      } else {
        return { "trigger": false };
      }
    } else {
      return { "trigger": true };
    }
  }