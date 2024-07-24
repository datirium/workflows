cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/sc-tools:v0.0.39


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

label: "Cell Ranger ATAC Count/Aggregate to UCSC Cell Browser"
s:name: "Cell Ranger ATAC Count/Aggregate to UCSC Cell Browser"
s:alternateName: |
  Exports clustering results from Cell Ranger ATAC Count or Cell Ranger ATAC Aggregate
  experiments into compatible with UCSC Cell Browser format

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/cellbrowser-build-cellranger-atac.cwl
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
  Cell Ranger ATAC Count/Aggregate to UCSC Cell Browser

  Exports clustering results from Cell Ranger ATAC Count
  or Cell Ranger ATAC Aggregate experiments into compatible
  with UCSC Cell Browser format


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
      Gene symbol/annotation files are downloaded to ~/cellbrowserData when
      needed. Config defaults can be specified in ~/.cellbrowser. See
      documentation at https://cellbrowser.readthedocs.io/
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
                          variable CBOUT or ~/.cellbrowser, current value: none
    -p PORT, --port=PORT  if build is successful, start an http server on this
                          port and serve the result via http://localhost:port
    -r, --recursive       run in all subdirectories of the current directory.
                          Useful when rebuilding a full hierarchy. Cannot be
                          used with -p.
    --depth=DEPTH         when using -r: only go this many directories deep
    --redo=REDO           do not use cached old data. Can be: 'meta' or 'matrix'
                          (matrix includes meta).
    --force               ignore errors that usually stop the build and go ahead
                          anyways.