<!-- [![Build Status](https://travis-ci.com/datirium/workflows.svg?branch=master)](https://travis-ci.com/github/datirium/workflows) -->

# Bioinformatics Workflows by Datirium LLC

- [Repository folder structure](#repo-structure)
<!-- - [Contributing](#contributing) -->
- [CWL for SciDAP](#augmented-cwl-standard-for-scidap)



The workflows in this repository are capable of running analyses on ([NGS sequencing](https://www.illumina.com/science/technology/next-generation-sequencing.html)) data; ChIP-Seq, ATAC-Seq, CLIP-Seq, and RNA-Seq

Workflows are written in [CWL](https://www.illumina.com/science/technology/next-generation-sequencing.html).

These workflows are compatable within;

- [Scientific Data Analysis Platform (SciDAP)](https://scidap.com)
<!-- - [BioWardrobe](https://biowardrobe.com/) project -->
- or standalone with [cwltool](https://github.com/common-workflow-language/cwltool)

<!-- All the original [BioWardrobe pipelines](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0720-3) have been rewritten in CWL, and new workflows has been added.   -->
    
Workflows available on SciDAP are automatically updated and version controlled based on PR's into this repository.

---

## Repo structure
> in order of importance

```
├── workflows
├── tools
├── dockerfiles
│   ├── scripts
├── docs
│   ├── contributing
│   ├── tutorials
│   ├── cwl_examples
├── tests
│   ├── data
```
<!-- ├── .github -->

- The **workflows** directory is the home for all of the CWL files that makeup individual workflows. Workflows will use **tools** for running individual steps.

- The **tools** directory is the home for all of the CWL files that handle general steps. Tools will utilize docker images for running analyses in containers (to maintain run conditions).

- The **dockerfiles** directory is home to the files used to create images that are used by tools.
  - The **scripts** directory is home to scripts that are included in built docker images. These scripts are often either in R, or bash

<!-- - is home to 3 folders: -->
  <!-- - **contributing**: docs related to contributing to the open-source workflows, or documentation, of this repo. -->
- The **docs** directory:
  - [**cwl_examples**](./docs/cwl_examples/): docs with examples of specific aspects of contributing or development
  - [**tutorials**](./docs/tutorials/github_flow_for_workflows.md): docs with tutorials on different aspects of using/developing/testing workflows

- The **tests** directory is home to json files (jobs) for each workflow/tool. These job files are what are used in order to test individual workflows. Using these tests requires data, which the child folder **data** contains through a repo-reference to [Barski-Labs workflow_tests](https://github.com/Barski-lab/workflows_test/tree/305adf1fd61a08dff5e3f348296a1b246fc8683a)




<!-- - The **.github** directory contains yaml files defining actions that github will take depending on what how the remote repository is updated. -->

<!-- ---

## Contributing

See the [contributing guide](./docs/contributing.md) for detailed instructions on how to get started with our project. -->

---

## Augmented CWL standard for SciDAP
There are 4 additional references/tags that can be included in different parts of a workflow for added compatability within SciDAP.

<!-- 1. [Metadata](#metadata): For establishing inputs on the add_sample form that are shared among many workflows. -->
1. [Upstreams](#upstreams): For designating what workflows generate outputs that some workflow can use as inputs.
2. [Visual Plugins](#visualplugins-for-an-output-type-file): For added visualizations of output data within the SciDAP platform
3. [Service Tags](#service-tags-for-workflows): For differentiating what kind of samples this workflow creates
4. [sql query for input](#sql-for-input): How to allow user to dynamically create sql query based off of options (saved as string for cwl input)
---
### Metadata

To extend user interface (dynamic form) with extra input fields not required by a workflow ```'sd:metadata'``` field were introduced.
It defines a list of workflow templates where ```inputs``` object is used for constructing and storing extra fields with an original workflow.

[Example of 'metadata' template for user interface](./docs/cwl_examples/example_metadata_in_workflow.md)

---
### Upstreams

To allow selection of already analysed data as input for a workflow, we organize a graph of separate workflows. To link workflows we use ```’sd:upstream’```, which defines a list of upstream workflows that this workflow can use for input data.

[example of workflow with upstreams](./docs/cwl_examples/example_upstream_workflow_for_input.md)

---
### VisualPlugins for an output type file

Usually, workflows' output results (especially files) are provided as download links on the SciDAP platform. 

With SciDAP's visualization plugins, output data can be presented as a;
- plot
- a genome (igv) browser
- a table
- or (in the case of html outputs), can be opened in a new tab. 

The keyword `'sd:visualPlugins'` enables SciDAP visualization plugins. 

The `line`, `pie`, `chart`, `igvbrowser`, `syncfusiongrid`, and `linkList` types are already available in the platform.

[Example of visual plugins used on workflow outputs](./docs/cwl_examples/example_visualPlugins_in_workflows.md)

--- 
### Service Tags for workflows
The ```'sd:serviceTag'```keyword enables new workflows to be added for the creation of:
- samples: uses keyword ```'sample'```
- analyses: uses keyword ```'analysis'```
- genelist: uses keywork ```'genelist'```

The service tag on a workflow will determine how samples are listed when viewing a project on the SciDAP platform.

Workflows without a service tag (or with one not recognized) will create samples in a tab called "not in use"

### SQL for Input

```yml
inputs: 
#...
  sql_query:
    type: string
    label: "Filtering parameters"
    doc: "Filtering parameters (WHERE parameters for SQL query)"
    'sd:filtering':
      params:
        columns: ["Refseq_id", "Gene_id", "txStart", "txEnd", "Strand", "Region", "Chr", "Start", "End", "Conc", "Conc1", "Conc2", "Fold", "p-value", "FDR", "Called1", "Called2"]
        types:   ["string", "string", "number", "number", "string", "string", "string", "number", "number", "number", "number", "number", "number", "number", "number","number", "number"]
       
```

will create an sql query based on the values given for any grouping and selection of the columns