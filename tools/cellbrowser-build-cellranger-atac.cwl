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
      name = "ATAC"
      shortLabel = "ATAC"
      priority = 1
      geneIdType = "auto"
      exprMatrix = "exprMatrix.tsv.gz"
      meta = "meta.csv"
      coords = [
          {
              "file": "tsne.coords.csv",
              "shortLabel": "t-SNE"
          },
          {
              "file": "umap.coords.csv",
              "shortLabel": "UMAP"
          },
          {
              "file": "lsa.coords.csv",
              "shortLabel": "LSA"
          }
      ]
      markers=[
          {
              "file": "markers.tsv",
              "shortLabel": "Cluster-specific peaks"
          }
      ]
      geneLabel = "Feature"
      radius = 3
      alpha = 0.5
      clusterField = "Cluster"
      labelField = "Cluster"
      atacSearch = "genome.current"
  - entryname: desc.conf
    entry: |
      title = "ATAC"
      abstract = ""
      methods = ""
      biorxiv_url = ""
      custom = {}
inputs:
  bash_script:
    type: string?
    default: |
      #!/bin/bash
      echo "Preparing ATAC search file"
      sc_cb_utils_atac_search.R --annotations $2
      echo "Prepare input data"
      mkdir -p ./cellbrowser_input/analysis/clustering/graphclust \
               ./cellbrowser_input/analysis/diffexp/graphclust \
               ./cellbrowser_input/filtered_feature_bc_matrix
      cp -r $0/clustering/graphclust/clusters.csv ./cellbrowser_input/analysis/clustering/graphclust/clusters.csv
      cp -r $0/enrichment/graphclust/differential_expression.csv ./cellbrowser_input/analysis/diffexp/graphclust/differential_expression.csv
      cp -r $0/tsne ./cellbrowser_input/analysis/
      cp -r $0/umap ./cellbrowser_input/analysis/
      cp -r $0/lsa ./cellbrowser_input/analysis/
      cp -r $1/* ./cellbrowser_input/filtered_feature_bc_matrix/
      cd ./cellbrowser_input/filtered_feature_bc_matrix/
      gzip barcodes.tsv
      gzip matrix.mtx
      cat peaks.bed | awk '{print $1":"$2"-"$3"\t"$1":"$2"-"$3"\tPeaks\t"$0}' > features.tsv
      gzip features.tsv
      rm -f peaks.bed
      cd -
      echo "Run cbImportCellranger"
      cbImportCellranger -i cellbrowser_input -o cellbrowser_output --name cellbrowser
      cd ./cellbrowser_output
      echo "Copying coordinates files"
      cp ../cellbrowser_input/analysis/tsne/*/projection.csv tsne.coords.csv
      cp ../cellbrowser_input/analysis/umap/*/projection.csv umap.coords.csv
      cp ../cellbrowser_input/analysis/lsa/*/projection.csv lsa.coords.csv
      echo "Replace configuration files"
      rm -f cellbrowser.conf desc.conf
      cp ../cellbrowser.conf .
      cp ../desc.conf .
      if [[ -n $3 ]]; then
        echo "Aggregation metadata file was provided. Adding initial cell identity classes"
        cat $3 | grep -v "library_id" | awk '{print NR","$0}' > aggregation_metadata.csv
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
      Folder with secondary
      analysis results.
  filtered_feature_bc_matrix_folder:
    type: Directory
    inputBinding:
      position: 7
    doc: |
      Folder with filtered peak-barcode matrices
      containing only cellular barcodes
      in MEX format.
  annotation_gtf_file:
    type: File
    inputBinding:
      position: 8
    doc: |
      GTF annotation file.
  aggregation_metadata:
    type: File?
    inputBinding:
      position: 9
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
label: Cell Ranger ATAC Count/Aggregate to UCSC Cell Browser
doc: |
  Cell Ranger ATAC Count/Aggregate to UCSC Cell Browser

  Exports clustering results from Cell Ranger ATAC Count
  or Cell Ranger ATAC Aggregate experiments into compatible
  with UCSC Cell Browser format
