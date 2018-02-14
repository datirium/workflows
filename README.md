# workflows
CWL workflows for [BioWardrobe](https://biowardrobe.com/) project

## Install 

```sh
git clone https://github.com/datirium/workflows
pip3 install .
```

## Usage

```python
from biowardrobe_cwl_workflows import available
_path_to_workflow=available(workflow='chipseq-se.cwl')
```
