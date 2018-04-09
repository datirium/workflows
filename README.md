[![Build Status](https://travis-ci.org/Barski-lab/workflows.svg?branch=master)](https://travis-ci.org/Barski-lab/workflows)
# workflows
CWL workflows for [BioWardrobe](https://biowardrobe.com/) project

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
