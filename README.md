[![Build Status](https://travis-ci.com/datirium/workflows.svg?branch=master)](https://travis-ci.com/github/datirium/workflows)


## Bioinformatics Workflows by Datirium LLC


ChIP-Seq, ATAC-Seq, CLIP-Seq, RNA-Seq CWL workflows for use in [Scientific Data Analysis Platform (SciDAP)](https://scidap.com)
or in [BioWardrobe](https://biowardrobe.com/) project or standalone with [cwltool](https://github.com/common-workflow-language/cwltool).

All the original BioWardrobe's [pipelines](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0720-3)
    has been rewritten in CWL and new workflows has been added.  The repository pulls automatically into SciDAP platform.

## Augmented CWL standard for SciDAP
There are 4 additional references that can be given to a workflow for added compatability within SciDAP.
1. [Metadata](#metadata)
2. [Upstreams](#upstreams)
3. [Visual Plugins](#visualplugins-for-an-output-type-file) 
4. [Service Tags](#service-tags-for-workflows)

### Metadata

To extend user interface (dynamic form) with extra input fields not required by a workflow, the ```'sd:metadata'``` field was introduced.
It defines a list of workflow templates where the ```inputs``` object is used for constructing and storing extra fields with an original workflow.

#### Example of 'metadata' template for user interface: 

> chipseq-header.cwl
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
and include file as ```sd:metadata```
```yaml
'sd:metadata':
    - "../metadata/chipseq-header.cwl"
```

### Upstreams

To extend the SciDAP UI to allow for already analysed data the be selectable as inputs, we organize a graph of separate workflows. To link workflows we use ```’sd:upstream’```, which defines a list of upstream workflows who's outputs are accible by this workflow.



```yaml
...
'sd:upstream':
  rnaseq_sample_untreated:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"
  rnaseq_sample_treated:
    - "rnaseq-se.cwl"
    - "rnaseq-pe.cwl"
    - "rnaseq-se-dutp.cwl"
    - "rnaseq-pe-dutp.cwl"
    - "rnaseq-se-dutp-mitochondrial.cwl"
    - "rnaseq-pe-dutp-mitochondrial.cwl"
    - "trim-rnaseq-pe.cwl"
    - "trim-rnaseq-se.cwl"
    - "trim-rnaseq-pe-dutp.cwl"
    - "trim-rnaseq-se-dutp.cwl"

inputs:

  untreated_files:
    type: File[]
    format:
     - "http://edamontology.org/format_3752"
     - "http://edamontology.org/format_3475"
    label: "Untreated input CSV/TSV files"
    doc: "Untreated input CSV/TSV files"
    'sd:upstreamSource': "rnaseq_sample_untreated/rpkm_common_tss"
    'sd:localLabel': true
...
```

### VisualPlugins for an output type file

Usually, workflows' output results (especially files) are provided as download links in web interfaces. With SciDAP visualization plugins, data can be presented as a plot, as a genome browser, as a table, or (in the case of html outputs) to be opened in a new tab. The keyword `'sd:visualPlugins'` enables SciDAP visualization plugins. `line`, `pie`, `chart`, `igvbrowser`, `syncfusiongrid`, and `linkList` types are already available in the platform.


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
        label: "DESeq results, TSV"
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

### Service Tags for workflows
The ```'sd:serviceTag'```keyword enables new workflows to be added for the creation of:
- samples: uses keyword ```'sample'```
- analyses: uses keyword ```'analysis'```
- genelist: uses keywork ```'genelist'```