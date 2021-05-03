cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      return [
        {
          "entry": inputs.query_feature_bc_matrices_h5,
          "entryname": inputs.query_feature_bc_matrices_h5.basename,
          "writable": true
        },
        {
          "entry": inputs.reference_marker_heatmap_file,
          "entryname": inputs.reference_marker_heatmap_file.basename,
          "writable": true
        },
        {
          "entry": inputs.reference_annotation_metadata_file,
          "entryname": inputs.reference_annotation_metadata_file.basename,
          "writable": true
        },
        {
          "entry": inputs.reference_expression_matrix_file,
          "entryname": inputs.reference_expression_matrix_file.basename,
          "writable": true
        }
      ]
    }


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/altanalyze:v0.0.5


inputs:

  bash_script:
    type: string?
    default: |
      #!/bin/bash
      
      # Copy altanalyze to the current working directory
      # which is mount with -rw- permissions. Otherwise 
      # we can't download anything because by default
      # container is run by cwltool with --read-only
      
      cp -r /opt/altanalyze .

      GENOME_DATA=$0

      GENOME_DATA_BASENAME=`basename ${GENOME_DATA}`
      ARR=(${GENOME_DATA_BASENAME//__/ })
      ENSEMBL_VERSION=${ARR[0]}
      SPECIES=${ARR[1]}

      echo "Load ${SPECIES} from ${ENSEMBL_VERSION}"
      ln -s ${GENOME_DATA} ./altanalyze/AltDatabase/${ENSEMBL_VERSION}

      python ./altanalyze/AltAnalyze.py --species ${SPECIES} \
      --platform RNASeq --cellHarmony yes ${@:1}
    inputBinding:
      position: 5
    doc: |
      Bash script to run AltAnalyze CellHarmony with provided parameters

  genome_data:
    type: Directory
    inputBinding:
      position: 6
    doc: |
      Ensembl database from the altanalyze-prepare-genome.cwl pipeline.
      --species parameter and Ensembl version will be resolved based on
      this folder basename

  query_feature_bc_matrices_h5:
    type: File
    inputBinding:
      position: 7
      prefix: "--input"
    doc: |
      Query dataset. Feature-barcode matrices in HDF5 format. Output from Cell Ranger

  reference_marker_heatmap_file:
    type: File
    inputBinding:
      position: 8
      prefix: "--reference"
    doc: |
      Reference marker gene heatmap file produced by altanalyze-icgs.cwl pipeline

  reference_annotation_metadata_file:
    type: File
    inputBinding:
      position: 9
      prefix: "--labels"
    doc: |
      Reference annotation metadata file produced by altanalyze-icgs.cwl pipeline

  reference_expression_matrix_file:
    type: File
    inputBinding:
      position: 10
      prefix: "--referenceFull"
    doc: |
      Reference expression matrix file produced by altanalyze-icgs.cwl pipeline

  align_by:
    type:
    - "null"
    - type: enum
      symbols:
      - "centroid"
      - "cell"
    default: "centroid"
    inputBinding:
      prefix: "--referenceType"
      position: 11
    doc: |
      Aligning to cluster centroid or cell

  correlation_threshold:
    type: float?
    default: 0.4
    inputBinding:
      prefix: "--correlationCutoff"
      position: 12
    doc: |
      Pearson correlation threshold

  perform_diff_expression:
    type: boolean?
    default: true
    inputBinding:
      prefix: "--performDiffExp"
      position: 13
      valueFrom: $(self?"yes":"no")
    doc: |
      Perform differential expression analysis

  diff_expr_fold_change_threshold:
    type: float?
    default: 1.5
    inputBinding:
      prefix: "--fold"
      position: 14
    doc: |
      Differential expression fold-change threshold.
      Applied if running differential expression

  diff_expr_p_value_threshold:
    type: float?
    default: 0.05
    inputBinding:
      prefix: "--pval"
      position: 15
    doc: |
      Cutoff value for P-value of P-adjusted-value.
      Applied if running differential expression

  use_adjusted_pvalue:
    type: boolean?
    default: true
    inputBinding:
      prefix: "--adjp"
      position: 16
      valueFrom: $(self?"yes":"no")
    doc: |
      Use adjusted P-value for differentially expressed genes threshold.
      Applied if running differential expression


outputs:

  cellharmony_data:
    type: Directory
    outputBinding: 
      glob: "cellHarmony"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]

stdout: altanalyze_cellharmony_stdout.log
stderr: altanalyze_cellharmony_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


s:name: "altanalyze-cellharmony"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/altanalyze-cellharmony.cwl
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
    - class: s:Organization
      s:legalName: "Salomonis Research Lab"
      s:member:
      - class: s:Person
        s:name: Stuart Hay
        s:email: mailto:haysb91@gmail.com


doc: |
  Runs AltAnalyze cellHarmony for 10X Genomics data

s:about: |
  CellHarmony is a cell-matching algorithm designed to identify a cell's most similar
  analogue in a distinct single-cell RNA-Seq (scRNA-Seq) dataset and find differentially
  expressed genes in each cell population.
