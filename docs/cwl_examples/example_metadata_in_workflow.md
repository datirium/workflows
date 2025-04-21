# Example of 'metadata' template for user interface in SciDAP: 

> chipseq-header.cwl: file establishing core/common inputs for chipseq related workflows
```yaml
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

```


> include a reference to the metadata file as ```sd:metadata``` in the workflow to include those inputs
```yaml
'sd:metadata':
    - "../metadata/chipseq-header.cwl"
```
