# Example of VisualPlugins on workflow outputs

The keyword `'sd:visualPlugins'` enables SciDAP visualization plugins. 

The `line`, `pie`, `chart`, `igvbrowser`, `syncfusiongrid`, and `linkList` types are already available in the platform.

> Examples:
- [line](#line-chart-example)
- [box plot](#box-plot-example)
- [pie chart](#pie-chart-example)
- [scatter plot](#scatter-plot-example)
- [3D scatter plot](#scatter-plot-3d-example)
- [comparetable?](#comparetable-example)
- [grid view](#syncfusiongrid-example)
- [IGV browser](#igvbrowser-example)
- [markdown](#markdown-view-example)
- [image](#image-example)
- [tableView](#table-view)
- [external/internal link](#link-list-example)
- [MA plot](#ma-plot-example)


## Line chart example

//TODO: insert image of example line chart view

```yaml
...

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
            Title: 'Base frequency plot' # title of plot
            xAxisTitle: 'Nucleotide position' # label of axis
            yAxisTitle: 'Frequency'
            colors: ["#b3de69", "#99c0db", "#fb8072", "#fdc381", "#888888"] # default colors used
            data: [$12, $13, $14, $15, $16] #which columns to include (first column is 1)
```

---

## Box plot example

//TODO insert image

```yaml
...

outputs:
    ...
    fastx_statistics_after:
    type: File
    label: "FASTQ statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats_after/statistics_file
    'sd:visualPlugins':
    - boxplot:
        tab: 'QC Plots'
        Title: 'After Clipper Quality Control'
        xAxisTitle: 'Nucleotide position'
        yAxisTitle: 'Quality score'
        colors: ["#b3de69", "#888888", "#fb8072", "#fdc381", "#99c0db"] # index of color correlates to index in file
        data: [$11, $7, $8, $9, $12] #index of columns to use (starts at 1)
```

---

## Pie chart example

Pie charts are generally included on a sample preview (sample-card when for listing samples in a project)

For adding a plugin to the card, the ```'sd:visualPlugin'``` key is prepended with a ```'sd:preview'``` keyword

//TODO insert image

```yaml
...

outputs:
    ...

  get_formatted_stats:
    type: File?
    label: "Bowtie, STAR and GEEP mapping stats"
    format: "http://edamontology.org/format_2330"
    doc: "Processed and combined Bowtie & STAR aligner and GEEP logs"
    outputSource: get_stat/collected_statistics_tsv
    'sd:preview':
      'sd:visualPlugins':
      - pie:
          colors: ['#b3de69', '#99c0db', '#fdc381', '#fb8072', '#778899']
          data: [$2, $3, $4, $5, $6]
```

---

## Scatter plot example

//TODO insert image

```yaml
...

outputs:
    ...
  atdp_result:
    type: File
    label: "ATDP results"
    format: "http://edamontology.org/format_3475"
    doc: "Average Tag Density generated results"
    outputSource: average_tag_density/result_file
    'sd:visualPlugins':
    - scatter:
        tab: 'QC Plots'
        Title: 'Average Tag Density'
        xAxisTitle: 'Distance From TSS (bases)'
        yAxisTitle: 'Average Tag Density (per bp)'
        colors: ["#b3de69"]
        height: 500 # width of plots is flexed based on amount. height can be set (plots can be opened in full screen)
        data: [$1, $2] # index of columns from file to use (starts at 1)
        comparable: "atdp"
```

---

## Scatter plot 3D example

//TODO insert image

```yaml
...

outputs:
    ...
  pca_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "PCA analysis scores results"
    doc: "PCA analysis scores results exported as TSV"
    outputSource: pca/pca_scores_file
    'sd:visualPlugins':
    - scatter3d:
        tab: '3D Plots'
        Title: 'PCA'
        xAxisTitle: 'PCA1'
        yAxisTitle: 'PCA2'
        zAxisTitle: 'PCA3'
        colors: ["#b3de69", "#888888", "#fb8072"]
        height: 600
        data: [$1, $2, $3, $4]
```

---

## Comparetable example

Inserts columns from file into synfusion grid

> (more verbose version of this plugin available with [syncfusion grids](#syncfusiongrid-example))

//TODO insert image

```yaml
...

outputs:
    ...

```

---

## Syncfusiongrid example

Displays file as excel sheet. Includes filtering options for better comparisons.

//TODO insert image

```yaml
...

outputs:
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
            tab: 'custom tab name'
    ...
```

---

## igvBrowser example

Displays bam/bigwig for exploratory visualization.
//TODO insert image

```yaml
...

outputs:
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

---

## Markdown view example

//TODO insert image

```yaml
...

outputs:
    ...
    samplesheet_md:
        type: File
        label: "Samplesheet with condition labels"
        doc: "Samplesheet with condition labels"
        outputSource: run_rnbeads_diff/samplesheet_overview
        'sd:visualPlugins':
        - markdownView:
            tab: 'Overview'
```

---

## Image example

//TODO insert image

```yaml
...

outputs:
    ...
    mbias_plot_png:
        type: File[]
        label: "Methylation bias plot (PNG)"
        doc: "QC data showing methylation bias across read lengths"
        format: "http://edamontology.org/format_3603"
        outputSource: bismark_extract_methylation/mbias_plot_png
        'sd:visualPlugins':
        - image:
            tab: 'Plots'
            Caption: 'Methylation bias plot'
```

---

## Table View

The tableView plugin allows a SINLGE tsv (and yaml) file to be used for creating a QC-table in sample-comparison.

```yaml
outputs: 
    ...
    get_stat_formatted_log:
        type: File?
        label: "Bowtie & Samtools Rmdup combined formatted log"
        format: "http://edamontology.org/format_3475"
        doc: "Processed and combined Bowtie aligner and Samtools rmdup formatted log"
        outputSource: get_stat/collected_statistics_tsv
        'sd:visualPlugins':
        - tableView:
            vertical: true
            tab: 'Overview'
```

There are 2 important things to note.

1. the **tableView** plugin also requires a yaml file (both generated from any **get_statistic_...** tools)
2. the yaml file should be saved as an output. It's name is irrelevant, but the output name from the get_statistic tool used should end in "yaml" 

---

## Link List example

Allows generated html files to be served by satellites and opened within or outside of the platform.

//TODO insert image

```yaml
...

outputs:
    ...
    volcano_plot_html_file:
        type: File
        outputSource: make_volcano_plot/html_file
        label: "Volcano Plot"
        doc: |
            HTML index file with volcano plot data.
        'sd:visualPlugins':
        - linkList:
            tab: 'Overview'
            target: "_blank" # target: "_this" should cause page to be opened within a component on the platform
```

---

## MA plot example

//TODO insert image

```yaml
...

outputs:
    ...
    diff_expr_features:
        type: File
        outputSource: deseq_multi_factor/diff_expr_features
        label: "TSV file with not filtered differentially expressed features"
        doc: |
            TSV file with not filtered differentially expressed features
        'sd:visualPlugins':
        - MAplot:
            tab: "MA plot" # can be own tab, or in "Plots"
            Title: "Differentially expressed features"
            data: ["$1", "$7", "$8"] # col for:   feature, baseMean, log2FC
            xAxisTitle: "Mean Expression (log[x])"
            yAxisTitle: "log 2 Fold Change"
            height: 500,
            colors: ["#b3de69"]
    ...
```
