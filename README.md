[![Build Status](https://travis-ci.org/datirium/workflows.svg?branch=master)](https://travis-ci.org/datirium/workflows)
# Datirium supported Bioinformatic Workflows
ChIP-Seq, ATAC-Seq, CLIP-Seq, RNA-Seq CWL workflows for use in [Scientific Data Analysis Platform (SciDAP)](https://scidap.com) 
or in [BioWardrobe](https://biowardrobe.com/) project or standalone with 
[cwltool](https://github.com/common-workflow-language/cwltool). 
 All the original BioWardrobe's [pipelines](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0720-3)
   has been rewritten in CWL and new has been added. The package simplifies import 
for biowardrobe-airflow-analysis. The repository pulls automatically into SciDAP platform. 

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

## SciDAP supported features in CWL

### Metadata & Upstreams
```yaml
'sd:metadata':
  - "../metadata/chipseq-header.cwl"

'sd:upstream':
  genome_indices: "genome-indices.cwl"
  control_file: "chipseq-se.cwl"

```

To enable users to extend user interface (dynamic form) with extra input fields not required by a workflow ```'sd:metadata'``` field were introduced. 
It defines a list of workflows where inputs fields just used for constructing and storing the input object. 


To chain workflows in order ```’sd:upstream’``` (defines a list of upstream workflows) were introduced it simplifies 
selection of already analyzed upstream data. 


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

outputs: []
steps: []

```

### VisualPlugins for outputs

Usually workflows output results (especially files) are provided as download links in web interfaces. With SciDAP visualization plugins data can be presented as a plot, as a genome browser or as a table. Keyword `'sd:visualPlugins'` enables SciDAP visualization plugins. `line`, `pie`, `chart`, `igvbrowser` and `syncfusiongrid` types are already available in the platform.


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
  diff_expr_file:
    type: File
    label: "DESeq resutls, TSV"
    format: "http://edamontology.org/format_3475"
    doc: "DESeq generated list of differentially expressed items grouped by isoforms, genes or common TSS"
    outputSource: deseq/diff_expr_file
    'sd:visualPlugins':
    - syncfusiongrid:
        Title: 'Combined DESeq results'
  ...
  bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated BigWig file"
    outputSource: bam_to_bigwig/bigwig_file
    'sd:visualPlugins':
    - igvbrowser:
        id: 'igvbrowser'
        type: 'wig'
        name: "BigWig Track"
        height: 120
  ...  
```

