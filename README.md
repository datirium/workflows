[![Build Status](https://travis-ci.org/datirium/workflows.svg?branch=master)](https://travis-ci.org/datirium/workflows)
# Workflows
CWL workflows for [BioWardrobe](https://biowardrobe.com/) project. This package includes 
 all the original BioWardrobe's pipelines, simplifies import 
for biowardrobe-airflow-analysis. It also includes extra pipelines for new project SciDAP.

## Install 

```sh
pip3 install biowardrobe-cwl-workflows
```
or from sources
```sh
git clone https://github.com/datirium/workflows
pip3 install .
```

## Usage

```python
from biowardrobe_cwl_workflows import available
_path_to_workflow=available(workflow='chipseq-se.cwl')
```

## Extra features in CWL

### Metadata & Upstreams
```yaml
'sd:metadata':
  - "../metadata/chipseq-header.cwl"

'sd:upstream':
  genome_indices: "genome-indices.cwl"
  control_file: "chipseq-se.cwl"

```

To enable users to extend user interface (dynamic form) with extra input fields not required by a workflow ```'sd:metadata'``` field were introduced. It defines a list of workflows where inputs field just used for constructing and storing the input object. 


To simplify selection of already analyzed common data ```’sd:upstream’``` were introduced. It defines the upstream workflows of the process. The process is ready to run when upstream output data is available. 


##### Example of extra input fields for user interface:

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
  description:
    type: string?
    'sd:type': 'text'
    label: "Description"

outputs: []
steps: []

```

### CWL VisualPlugins for output data

```yaml
outputs:
  ...
  fastx_statistics:
    type: File
    label: "FASTQ statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats/statistics_file
    'sd:visualPlugins':
    - line:
      Title: 'Base frequency plot'
      xAxisTitle: 'Nucleotide position'
      yAxisTitle: 'Frequency'
      colors: ["#b3de69", "#99c0db", "#fb8072", "#fdc381", "#888888"]
      data: [$12, $13, $14, $15, $16]
  ...
```