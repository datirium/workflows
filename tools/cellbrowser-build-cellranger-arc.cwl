cwlVersion: v1.0
class: CommandLineTool
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.41
requirements:
- class: EnvVarRequirement
  envDef:
    CBDATAROOT: $(runtime.outdir)
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entryname: cellbrowser_rna.conf
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
              "shortLabel": "t-SNE"
          },
          {
              "file": "umap.coords.csv",
              "shortLabel": "UMAP"
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
      dataRoot = "../"
  - entryname: cellbrowser_atac.conf
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
      markers = [
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
      dataRoot = "../"
      atacSearch = "genome.current"
  - entryname: cellbrowser.conf
    entry: |
      shortLabel = "Multiple datasets"
  - entryname: desc_rna.conf
    entry: |
      title = "RNA"
      abstract = ""
      methods = ""
      biorxiv_url = ""
      custom = {}
  - entryname: desc_atac.conf
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
      echo "Splitting combined feature-barcode matrix into RNA and ATAC matrices"
      sc_cb_utils_split_mex.R --mex $1 --output temp_sc
      echo "Preparing ATAC search file"
      sc_cb_utils_atac_search.R --annotations $2
      echo "Preparing ATAC data"
      mkdir -p ./atac_input/analysis/clustering/graphclust \
              ./atac_input/analysis/diffexp/graphclust \
              ./atac_input/analysis/tsne/2_components \
              ./atac_input/analysis/umap/2_components \
              ./atac_input/analysis/lsa/2_components
      cp -r $0/clustering/atac/graphclust/clusters.csv ./atac_input/analysis/clustering/graphclust/clusters.csv
      cp -r $0/clustering/atac/graphclust/differential_accessibility.csv ./atac_input/analysis/diffexp/graphclust/differential_expression.csv
      cp -r $0/dimensionality_reduction/atac/tsne_projection.csv ./atac_input/analysis/tsne/2_components/projection.csv
      cp -r $0/dimensionality_reduction/atac/umap_projection.csv ./atac_input/analysis/umap/2_components/projection.csv
      cp -r $0/dimensionality_reduction/atac/lsa_projection.csv ./atac_input/analysis/lsa/2_components/projection.csv
      mkdir -p ./atac_input/filtered_feature_bc_matrix
      cp -r temp_sc_atac/* ./atac_input/filtered_feature_bc_matrix/
      echo "Importing ATAC data"
      cbImportCellranger -i atac_input -o atac --name atac
      cd ./atac
      echo "Copying coordinates files"
      cp ../atac_input/analysis/tsne/2_components/projection.csv tsne.coords.csv
      cp ../atac_input/analysis/umap/2_components/projection.csv umap.coords.csv
      cp ../atac_input/analysis/lsa/2_components/projection.csv lsa.coords.csv
      echo "Replacing configuration files"
      rm -f cellbrowser.conf desc.conf
      cp ../cellbrowser_atac.conf cellbrowser.conf
      cp ../desc_atac.conf desc.conf
      if [[ -n $3 ]]; then
          echo "Aggregation metadata file was provided. Adding initial cell identity classes"
          cat $3 | grep -v "library_id" | awk '{print NR","$0}' > aggregation_metadata.csv
          cat meta.csv | grep -v "Barcode" > meta_headerless.csv
          echo "Barcode,Cluster,Dataset" > meta.csv
          awk -F, 'NR==FNR {identity[$1]=$2; next} {split($1,barcode,"-"); print $0","identity[barcode[2]]}' aggregation_metadata.csv meta_headerless.csv >> meta.csv
          rm -f aggregation_metadata.csv meta_headerless.csv
      fi
      cbBuild -o ../html_data
      cd ..
      echo "Preparing RNA data"
      mkdir -p ./rna_input/analysis/clustering/graphclust \
              ./rna_input/analysis/diffexp/graphclust \
              ./rna_input/analysis/tsne/2_components \
              ./rna_input/analysis/umap/2_components
      cp -r $0/clustering/gex/graphclust/clusters.csv ./rna_input/analysis/clustering/graphclust/clusters.csv
      cp -r $0/clustering/gex/graphclust/differential_expression.csv ./rna_input/analysis/diffexp/graphclust/differential_expression.csv
      cp -r $0/dimensionality_reduction/gex/tsne_projection.csv ./rna_input/analysis/tsne/2_components/projection.csv
      cp -r $0/dimensionality_reduction/gex/umap_projection.csv ./rna_input/analysis/umap/2_components/projection.csv
      mkdir -p ./rna_input/filtered_feature_bc_matrix
      cp -r temp_sc_rna/* ./rna_input/filtered_feature_bc_matrix/
      echo "Importing RNA data"
      cbImportCellranger -i rna_input -o rna --name rna
      cd ./rna
      echo "Copying coordinates files"
      cp ../rna_input/analysis/tsne/2_components/projection.csv tsne.coords.csv
      cp ../rna_input/analysis/umap/2_components/projection.csv umap.coords.csv
      echo "Replacing configuration files"
      rm -f cellbrowser.conf desc.conf
      cp ../cellbrowser_rna.conf cellbrowser.conf
      cp ../desc_rna.conf desc.conf
      if [[ -n $3 ]]; then
          echo "Aggregation metadata file was provided. Adding initial cell identity classes"
          cat $3 | grep -v "library_id" | awk '{print NR","$0}' > aggregation_metadata.csv
          cat meta.csv | grep -v "Barcode" > meta_headerless.csv
          echo "Barcode,Cluster,Dataset" > meta.csv
          awk -F, 'NR==FNR {identity[$1]=$2; next} {split($1,barcode,"-"); print $0","identity[barcode[2]]}' aggregation_metadata.csv meta_headerless.csv >> meta.csv
          rm -f aggregation_metadata.csv meta_headerless.csv
      fi
      cbBuild -o ../html_data
      cd ..
      echo "Cleaning up temporary files"
      rm -rf rna_input atac_input atac rna temp_sc_rna temp_sc_atac
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
      produced by Cellranger ARC Count or
      Cellranger ARC Aggr.
  filtered_feature_bc_matrix_folder:
    type: Directory
    inputBinding:
      position: 7
    doc: |
      Folder with filtered feature-barcode
      matrices containing only cellular
      barcodes in MEX format produced by
      Cellranger ARC Count or Cellranger
      ARC Aggr.
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
      glob: html_data
  index_html_file:
    type: File
    outputBinding:
      glob: html_data/index.html
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- bash
- -c
stdout: cbbuild_stdout.log
stderr: cbbuild_stderr.log
label: Cell Ranger ARC Count/Aggregate to UCSC Cell Browser
doc: |
  Cell Ranger ARC Count/Aggregate to UCSC Cell Browser

  Exports clustering results from Cell Ranger ARC Count
  Chromatin Accessibility and Gene Expression or Cell
  Ranger ARC Aggregate experiments into compatible with
  UCSC Cell Browser format
