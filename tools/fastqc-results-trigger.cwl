cwlVersion: v1.0
class: ExpressionTool

requirements:
  - class: InlineJavascriptRequirement

id: "fastqc_results_trigger"
inputs:
  summary_file:
    type: File
    inputBinding:
      loadContents: true
  criteria:
    type: string?
    default: ".*Per base sequence quality.*|.*Per sequence quality scores.*|.*Overrepresented sequences.*|.*Adapter Content.*"
outputs:
  trigger: boolean
expression: |
  ${
    var re = new RegExp(inputs.criteria, "g");
    var criteria_array = inputs.summary_file.contents.match(re);
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