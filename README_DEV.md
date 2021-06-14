[![Build Status](https://travis-ci.org/datirium/workflows.svg?branch=master)](https://travis-ci.org/datirium/workflows)

## Augmented CWL standard with

### Metadata

To extend user interface (dynamic form) with extra input fields not required by a workflow ```'sd:metadata'``` field were introduced.
It defines a list of workflow templates where ```inputs``` object is used for constructing and storing extra fields with an original workflow.

### Example of 'metadata' template for user interface: 

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

To allow selection of already analysed data as input for a workflow we organize a graph of separate workflows. To link workflows we use ```’sd:upstream’``` which defines a list of upstream workflows.

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