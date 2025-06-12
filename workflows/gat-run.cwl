cwlVersion: v1.0
class: Workflow
requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: InlineJavascriptRequirement
- class: MultipleInputFeatureRequirement
inputs:
  alias:
    type: string
    label: Experiment short name/Alias
    sd:preview:
      position: 1
  segment_file:
    type: File
    format: http://edamontology.org/format_3003
    label: 'Segments of interest (headerless TSV/CSV file with at least three columns: chr start end)'
    doc: |
      Headerless TSV/CSV file (at least 3 columns: chr start end) with sets of intervals whose association
      will be tested with annotation_file
  annotation_file:
    type: File
    format: http://edamontology.org/format_3003
    label: 'Annotations (headerless TSV/CSV file with at least three columns: chr start end)'
    doc: |
      Headerless TSV/CSV file (at least 3 columns: chr start end) with sets of intervals that are used
      for testing association of segment_file
  workspace_file:
    type: File
    format: http://edamontology.org/format_3003
    label: 'Workspace (headerless TSV/CSV file with at least three columns: chr start end)'
    doc: |
      Headerless TSV/CSV file (at least 3 columns: chr start end) with genomic regions
      accessible for simulation
  permutations_count:
    type: int?
    default: 10000
    label: Number of random permutations
    doc: Number of random permutations
    sd:layout:
      advanced: true
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
    default: nucleotide-overlap
    label: Set the measure of association to be tested
    doc: |
      Set the measure of association to be tested
    sd:layout:
      advanced: true
  seed:
    type: int?
    default: 12345
    label: Number of random seed
    doc: Number of random seed
    sd:layout:
      advanced: true
  threads:
    type: int?
    default: 4
    label: Number of threads
    doc: Number of threads for those steps that support multithreading
    sd:layout:
      advanced: true
outputs:
  gat_report:
    type: File
    format: http://edamontology.org/format_3475
    label: Interval overlap statistics report
    doc: Interval overlap statistics report
    outputSource: prepare_report_file/output_file
    sd:visualPlugins:
    - syncfusiongrid:
        tab: Interval overlap statistics
        Title: Interval overlap statistics
  gat_stdout_log:
    type: File
    format: http://edamontology.org/format_2330
    label: GAT stdout log
    doc: GAT stdout log
    outputSource: run_gat/stdout_log
  gat_stderr_log:
    type: File
    format: http://edamontology.org/format_2330
    label: GAT stderr log
    doc: GAT stderr log
    outputSource: run_gat/stderr_log
steps:
  prepare_segment_file:
    run: ../tools/custom-bash.cwl
    in:
      input_file: segment_file
      script:
        default: |
          cat "$0" | tr -d '\r' | tr "," "\t" | cut -f 1-3 | awk NF | sort -u -k1,1 -k2,2n -k3,3n > `basename $0`
    out:
    - output_file
  prepare_annotation_file:
    run: ../tools/custom-bash.cwl
    in:
      input_file: annotation_file
      script:
        default: |
          cat "$0" | tr -d '\r' | tr "," "\t" | cut -f 1-3 | awk NF | sort -u -k1,1 -k2,2n -k3,3n > `basename $0`
    out:
    - output_file
  prepare_workspace_file:
    run: ../tools/custom-bash.cwl
    in:
      input_file: workspace_file
      script:
        default: |
          cat "$0" | tr -d '\r' | tr "," "\t" | cut -f 1-3 | awk NF | sort -u -k1,1 -k2,2n -k3,3n > `basename $0`
    out:
    - output_file
  run_gat:
    run: ../tools/gat-run.cwl
    in:
      segment_file: prepare_segment_file/output_file
      annotation_file: prepare_annotation_file/output_file
      workspace_file: prepare_workspace_file/output_file
      iterations: permutations_count
      counter: counter
      seed: seed
      threads: threads
    out:
    - report_file
    - stdout_log
    - stderr_log
  prepare_report_file:
    run: ../tools/custom-bash.cwl
    in:
      input_file: run_gat/report_file
      script:
        default: |
          cat "$0" | cut -f 3-  > `basename $0`
    out:
    - output_file
label: GAT - Genomic Association Tester
doc: |-
  GAT: Genomic Association Tester
  ==============================================

  A common question in genomic analysis is whether two sets of genomic intervals overlap significantly.
  This question arises, for example, in the interpretation of ChIP-Seq or RNA-Seq data. The Genomic
  Association Tester (GAT) is a tool for computing the significance of overlap between multiple sets of
  genomic intervals. GAT estimates significance based on simulation.

  Gat implemements a sampling algorithm. Given a chromosome (workspace) and segments of interest, for
  example from a ChIP-Seq experiment, gat creates randomized version of the segments of interest falling
  into the workspace. These sampled segments are then compared to existing genomic annotations.

  The sampling method is conceptually simple. Randomized samples of the segments of interest are created
  in a two-step procedure. Firstly, a segment size is selected from to same size distribution as the
  original segments of interest. Secondly, a random position is assigned to the segment. The sampling stops
  when exactly the same number of nucleotides have been sampled. To improve the speed of sampling, segment
  overlap is not resolved until the very end of the sampling procedure. Conflicts are then resolved by
  randomly removing and re-sampling segments until a covering set has been achieved. Because the size of
  randomized segments is derived from the observed segment size distribution of the segments of interest,
  the actual segment sizes in the sampled segments are usually not exactly identical to the ones in the
  segments of interest. This is in contrast to a sampling method that permutes segment positions within
  the workspace.
sd:version: 100
