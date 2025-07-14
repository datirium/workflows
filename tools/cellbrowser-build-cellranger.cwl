cwlVersion: v1.0
class: CommandLineTool
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.41
requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entryname: cellbrowser.conf
    entry: |
      name = "RNA"
      shortLabel = "RNA"
      priority = 1
      geneIdType = "auto"
      exprMatrix = "exprMatrix.tsv.gz"
      meta = "meta.csv"
      coords = [
          {
              "file": "tsne.coords.csv",
              "shortLabel": "CellRanger t-SNE"
          },
          {
              "file": "umap.coords.csv",
              "shortLabel": "CellRanger UMAP"
          }
      ]
      markers = [
          {
            "file": "markers.tsv",
            "shortLabel": "Cluster-specific genes"
          }
      ]
      geneLabel = "Feature"
      radius = 3
      alpha = 0.5
      clusterField = "Cluster"
      labelField = "Cluster"
  - entryname: desc.conf
    entry: |
      title = "RNA"
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
      echo "Removing gene_expression_ part from all of the folder names in analysis"
      du -a ./cellbrowser_input/analysis | cut -f 2 | grep gene_expression | xargs -I{} bash -c 'mv "$1" "${1//gene_expression_/}"' -- {}
      echo "Run cbImportCellranger"
      cbImportCellranger -i cellbrowser_input -o cellbrowser_output --name cellbrowser
      cd ./cellbrowser_output
      echo "Copy UMAP coordinates files"
      cp ../cellbrowser_input/analysis/umap/2_components/projection.csv umap.coords.csv
      echo "Replace configuration files"
      rm -f cellbrowser.conf desc.conf
      cp ../cellbrowser.conf .
      cp ../desc.conf .
      if [[ -n $2 ]]; then
        echo "Aggregation metadata file was provided. Adding initial cell identity classes"
        cat $2 | grep -v "sample_id" | awk '{print NR","$0}' > aggregation_metadata.csv
        cat meta.csv | grep -v "Barcode" > meta_headerless.csv
        echo "Barcode,Cluster,Dataset" > meta.csv
        awk -F, 'NR==FNR {identity[$1]=$2; next} {split($1,barcode,"-"); print $0","identity[barcode[2]]}' aggregation_metadata.csv meta_headerless.csv >> meta.csv
        rm -f aggregation_metadata.csv meta_headerless.csv
      fi
      echo "Run cbBuild"
      cbBuild -o html_data
    inputBinding:
      position: 5
    doc: |
      Bash script to run cbImportCellranger
      and cbBuild commands.
  secondary_analysis_report_folder:
    type: Directory
    inputBinding:
      position: 6
    doc: |
      Folder with secondary analysis results
      including dimensionality reduction, cell
      clustering, and differential expression
      produced by Cellranger Count or Cellranger
      Aggr.
  filtered_feature_bc_matrix_folder:
    type: Directory
    inputBinding:
      position: 7
    doc: |
      Folder with filtered feature-barcode
      matrices containing only cellular
      barcodes in MEX format produced by
      Cellranger Count or Cellranger Aggr.
  aggregation_metadata:
    type: File?
    inputBinding:
      position: 8
    doc: |
      Cellranger aggregation CSV file. If
      provided, the Dataset metadata column
      will be added to the meta.csv.
outputs:
  html_data:
    type: Directory
    outputBinding:
      glob: cellbrowser_output/html_data
  index_html_file:
    type: File
    outputBinding:
      glob: cellbrowser_output/html_data/index.html
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- bash
- -c
stdout: cbbuild_stdout.log
stderr: cbbuild_stderr.log
label: Cell Ranger Count/Aggregate to UCSC Cell Browser
doc: |
  Cell Ranger Count/Aggregate to UCSC Cell Browser

  Exports clustering results from Cell Ranger Count Gene Expression
  and Cell Ranger Aggregate experiments into compatible with UCSC
  Cell Browser format.
