cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/cellbrowser:v0.0.2


requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entryname: cellbrowser.conf
    entry: |
      name = "cellbrowser"
      shortLabel="cellbrowser"
      priority = 1
      geneIdType="auto"
      exprMatrix="exprMatrix.tsv.gz"
      meta="meta.csv"
      coords=[
          {
              "file": "tsne.coords.csv",
              "shortLabel": "CellRanger t-SNE"
          },
          {
              "file": "umap.coords.csv",
              "shortLabel": "CellRanger UMAP"
          }
      ]
      markers=[
        {
          "file":"markers.tsv",
          "shortLabel":"Cluster-specific genes"
        }
      ]
      enumFields = ["Barcode"]
      clusterField="Cluster"
      labelField="Cluster"
  - entryname: desc.conf
    entry: |
      title = "CellBrowser"
      abstract = ""
      methods = ""
      biorxiv_url = ""
      custom = {}


inputs:

  bash_script:
    type: string?
    default: |
      #!/bin/bash
      echo "Prepare input data"
      mkdir -p ./cellbrowser_input/analysis ./cellbrowser_input/filtered_feature_bc_matrix
      cp -r $0/* ./cellbrowser_input/analysis/
      cp -r $1/* ./cellbrowser_input/filtered_feature_bc_matrix/
      echo "Run cbImportCellranger"
      cbImportCellranger -i cellbrowser_input -o cellbrowser_output --name cellbrowser
      cd ./cellbrowser_output
      echo "Copy UMAP coordinates files"
      cp ../cellbrowser_input/analysis/umap/2_components/projection.csv umap.coords.csv
      echo "Replace configuration files"
      rm -f cellbrowser.conf desc.conf
      cp ../cellbrowser.conf .
      cp ../desc.conf .
      echo "Run cbBuild"
      cbBuild -o html_data
    inputBinding:
      position: 5
    doc: |
      Bash script to run cbImportCellranger and cbBuild commands

  secondary_analysis_report_folder:
    type: Directory
    inputBinding:
      position: 6
    doc: |
      Folder with secondary analysis results including dimensionality reduction,
      cell clustering, and differential expression produced by Cellranger Count
      or Cellranger Aggr

  filtered_feature_bc_matrix_folder:
    type: Directory
    inputBinding:
      position: 7
    doc: |
      Folder with filtered feature-barcode matrices containing only cellular
      barcodes in MEX format produced by Cellranger Count or Cellranger Aggr


outputs:

  html_data:
    type: Directory
    outputBinding:
      glob: "cellbrowser_output/html_data"

  index_html_file:
    type: File
    outputBinding:
      glob: "cellbrowser_output/html_data/index.html"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]


stdout: cbbuild_stdout.log
stderr: cbbuild_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "cellbrowser-build-cellranger"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/cellbrowser-build-cellranger.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Cincinnati Children's Hospital Medical Center"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: "45229"
    s:streetAddress: "3333 Burnet Ave"
    s:telephone: "+1(513)636-4200"
  s:logo: "https://www.cincinnatichildrens.org/-/media/cincinnati%20childrens/global%20shared/childrens-logo-new.png"
  s:department:
  - class: s:Organization
    s:legalName: "Allergy and Immunology"
    s:department:
    - class: s:Organization
      s:legalName: "Barski Research Lab"
      s:member:
      - class: s:Person
        s:name: Michael Kotliar
        s:email: mailto:misha.kotliar@gmail.com
        s:sameAs:
        - id: http://orcid.org/0000-0002-6486-3898


doc: |
  Converts Cellranger outputs into the data structure supported by UCSC CellBrowser


s:about: |
  Usage: cbImportCellranger [options] -i cellRangerDir -o outputDir - convert the cellranger output to cellbrowser format and create a cellranger.conf file
  Options:
    -h, --help            show this help message and exit
    -d, --debug           show debug messages
    -i INDIR, --inDir=INDIR
                          input folder with the cellranger analysis output. This
                          is the directory with the two directories 'analysis'
                          and 'filtered_gene_bc_matrices'
    -o OUTDIR, --outDir=OUTDIR
                          output directory
    -n DATASETNAME, --name=DATASETNAME
                          name of the dataset. No spaces or special characters.
    -m, --noMat           do not export the matrix again, saves some time if you
                          changed something small since the last run


  Usage: cbBuild [options] -i cellbrowser.conf -o outputDir - add a dataset to the single cell viewer directory
      If you have previously built into the same output directory with the same dataset and the
      expression matrix has not changed its filesize, this will be detected and the expression
      matrix will not be copied again. This means that an update of a few meta data attributes
      is quite quick.
  Options:
    -h, --help            show this help message and exit
    --init                copy sample cellbrowser.conf and desc.conf to current
                          directory
    -d, --debug           show debug messages
    -i INCONF, --inConf=INCONF
                          a cellbrowser.conf file that specifies labels and all
                          input files, default is ./cellbrowser.conf, can be
                          specified multiple times
    -o OUTDIR, --outDir=OUTDIR
                          output directory, default can be set through the env.
                          variable CBOUT or ~/.cellbrowser.conf, current value:
                          none
    -p PORT, --port=PORT  if build is successful, start an http server on this
                          port and serve the result via http://localhost:port
    -r, --recursive       run in all subdirectories of the current directory.
                          Useful when rebuilding a full hierarchy.
    --redo=REDO           do not use cached old data. Can be: 'meta' or 'matrix'
                          (matrix includes meta). 