cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/cellbrowser:v0.0.2


requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing:
  - entryname: cellbrowser_gex.conf
    entry: |
      name = "GEX"
      shortLabel="GEX"
      priority = 1
      geneIdType="auto"
      geneLabel="Gene"
      exprMatrix="exprMatrix.tsv.gz"
      meta="meta.csv"
      coords=[
          {
              "file": "tsne.coords.csv",
              "shortLabel": "t-SNE"
          },
          {
              "file": "umap.coords.csv",
              "shortLabel": "UMAP"
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
      dataRoot="../"
  - entryname: cellbrowser_atac.conf
    entry: |
      name = "ATAC"
      shortLabel="ATAC"
      priority = 1
      geneIdType="auto"
      geneLabel="Peak"
      exprMatrix="exprMatrix.tsv.gz"
      meta="meta.csv"
      coords=[
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
          "file":"markers.tsv",
          "shortLabel":"Cluster-specific peaks"
      }
      ]
      enumFields = ["Barcode"]
      clusterField="Cluster"
      labelField="Cluster"
      dataRoot="../"
  - entryname: cellbrowser.conf
    entry: |
      shortLabel="Multiome"
  - entryname: desc_gex.conf
    entry: |
      title = "GEX"
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
      cp -r $1/* ./atac_input/filtered_feature_bc_matrix/
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
      if [[ -n $2 ]]; then
          echo "Aggregation metadata file was provided. Adding initial cell identity classes"
          cat $2 | grep -v "library_id" | awk '{print NR","$0}' > aggregation_metadata.csv
          cat meta.csv | grep -v "Barcode" > meta_headerless.csv
          echo "Barcode,Cluster,Identity" > meta.csv
          awk -F, 'NR==FNR {identity[$1]=$2; next} {split($1,barcode,"-"); print $0","identity[barcode[2]]}' aggregation_metadata.csv meta_headerless.csv >> meta.csv
          rm -f aggregation_metadata.csv meta_headerless.csv
      fi
      cd ..
      echo "Preparing GEX data"
      mkdir -p ./gex_input/analysis/clustering/graphclust \
              ./gex_input/analysis/diffexp/graphclust \
              ./gex_input/analysis/tsne/2_components \
              ./gex_input/analysis/umap/2_components
      cp -r $0/clustering/gex/graphclust/clusters.csv ./gex_input/analysis/clustering/graphclust/clusters.csv
      cp -r $0/clustering/gex/graphclust/differential_expression.csv ./gex_input/analysis/diffexp/graphclust/differential_expression.csv
      cp -r $0/dimensionality_reduction/gex/tsne_projection.csv ./gex_input/analysis/tsne/2_components/projection.csv
      cp -r $0/dimensionality_reduction/gex/umap_projection.csv ./gex_input/analysis/umap/2_components/projection.csv
      mkdir -p ./gex_input/filtered_feature_bc_matrix
      cp -r $1/* ./gex_input/filtered_feature_bc_matrix/
      echo "Importing GEX data"
      cbImportCellranger -i gex_input -o gex --name gex
      cd ./gex
      echo "Copying coordinates files"
      cp ../gex_input/analysis/tsne/2_components/projection.csv tsne.coords.csv
      cp ../gex_input/analysis/umap/2_components/projection.csv umap.coords.csv
      echo "Replacing configuration files"
      rm -f cellbrowser.conf desc.conf
      cp ../cellbrowser_gex.conf cellbrowser.conf
      cp ../desc_gex.conf desc.conf
      if [[ -n $2 ]]; then
          echo "Aggregation metadata file was provided. Adding initial cell identity classes"
          cat $2 | grep -v "library_id" | awk '{print NR","$0}' > aggregation_metadata.csv
          cat meta.csv | grep -v "Barcode" > meta_headerless.csv
          echo "Barcode,Cluster,Identity" > meta.csv
          awk -F, 'NR==FNR {identity[$1]=$2; next} {split($1,barcode,"-"); print $0","identity[barcode[2]]}' aggregation_metadata.csv meta_headerless.csv >> meta.csv
          rm -f aggregation_metadata.csv meta_headerless.csv
      fi
      cd ..
      echo "Building"
      cbBuild -r -o html_data
      echo "Cleaning up temporary files"
      rm -rf gex_input atac_input atac gex
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
      cell clustering, and differential expression produced by Cellranger ARC
      Count or Cellranger ARC Aggr

  filtered_feature_bc_matrix_folder:
    type: Directory
    inputBinding:
      position: 7
    doc: |
      Folder with filtered feature-barcode matrices containing only cellular
      barcodes in MEX format produced by Cellranger ARC Count or Cellranger ARC Aggr

  aggregation_metadata:
    type: File?
    inputBinding:
      position: 8
    doc: |
      Cellranger aggregation CSV file. If provided, the Identity metadata
      column will be added to the meta.csv


outputs:

  html_data:
    type: Directory
    outputBinding:
      glob: "html_data"

  index_html_file:
    type: File
    outputBinding:
      glob: "html_data/index.html"

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

s:name: "cellbrowser-build-cellranger-arc"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/cellbrowser-build-cellranger-arc.cwl
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
  Converts Cellranger ARC outputs into the data structure supported by UCSC CellBrowser


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