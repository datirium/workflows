cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() { if (inputs.output_filename == ""){ var root = inputs.segment_file.basename.split('.').slice(0,-1).join('.'); return (root == "")?inputs.segment_file.basename+".tsv":root+".tsv"; } else { return inputs.output_filename; } };
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/gat:v0.0.1
inputs:
  segment_file:
    type: File
    inputBinding:
      position: 5
      prefix: -s
    doc: |
      BED file (strictly 3 columns) with sets of intervals whose association is tested with annotation_file
  annotation_file:
    type: File
    inputBinding:
      position: 6
      prefix: -a
    doc: |
      BED file (strictly 3 columns) with sets of intervals that are used for testing association of segment_file
  workspace_file:
    type: File
    inputBinding:
      position: 7
      prefix: -w
    doc: |
      BED file (strictly 3 columns) with genomic regions accessible for simulation
  output_filename:
    type: string?
    inputBinding:
      position: 8
      valueFrom: $(default_output_filename())
      prefix: -S
    default: ''
    doc: |
      Output report file name
  iterations:
    type: int?
    inputBinding:
      position: 9
      prefix: -n
    doc: |
      Number of iterations
  counter:
    type:
    - 'null'
    - type: enum
      name: counter
      symbols:
      - nucleotide-overlap
      - nucleotide-density
      - segment-overlap
      - annotation-overlap
      - segment-midoverlap
      - annotation-midoverlap
    inputBinding:
      position: 10
      prefix: -c
    doc: |
      Set the measure of association that is tested.
      Default: nucleotide-overlap
  threads:
    type: int?
    inputBinding:
      position: 11
      prefix: -t
    doc: |
      Threads number
  seed:
    type: int?
    inputBinding:
      position: 12
      prefix: --random-seed=
      separate: false
    doc: |
      Random seed to initialize random number generator with
outputs:
  report_file:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: |
      Report file
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- gat-run.py
- --ignore-segment-tracks
stdout: gat_stdout.log
stderr: gat_stderr.log
doc: |
  A common question in genomic analysis is whether two sets of genomic intervals overlap significantly.
  This question arises, for example, in the interpretation of ChIP-Seq or RNA-Seq data. The Genomic
  Association Tester (GAT) is a tool for computing the significance of overlap between multiple sets of
  genomic intervals. GAT estimates significance based on simulation.

  Gat implemements a sampling algorithm. Given a chromosome (workspace) and segments of interest, for
  example from a ChIP-Seq experiment, gat creates randomized version of the segments of interest falling
  into the workspace. These sampled segments are then compared to existing genomic annotations.

  Note:
  --ignore-segment-tracks parameter is hardcoded
label: gat-run
