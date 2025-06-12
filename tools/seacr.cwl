cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
hints:
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/seacr:1.3--hdfd78af_2
inputs:
  script:
    type: string?
    default: |
      #!/bin/bash
      printf "$(date)\nLog file for seacr.cwl tool:\n\n"
      printf "INPUTS:\n"
      echo "\$0 - $0"
      echo "\$1 - $1"
      echo "\$2 - $2"
      echo "\$3 - $3"
      echo "\$4 - $4"
      # run seacr
      SEACR_1.3.sh $0 $1 $2 $3 $4
      # error handling logic
      # make empty output file if not created by seacr (i.e. no peaks called)
      if [[ $? == 1 ]]; then
        echo "exit code for SEACR is 1, check stderr log"
        touch $4.$3.bed
      elif [[ $? == 2 ]]; then
        echo "exit code for SEACR is 2, check stderr log"
        touch $4.$3.bed
      fi
    inputBinding:
      position: 1
  treatment_bedgraph:
    type: File
    inputBinding:
      position: 2
    doc: |
      A sorted bed file from 'bedtools-fragmentcounts.cwl',
      ready for peak calling by SEACR.
  numeric_threshold:
    type: float
    inputBinding:
      position: 3
    doc: |
      A numeric threshold n between 0 and 1 returns the top n
      fraction of peaks based on total signal within peaks.
      Default is "0.01".
  norm_control_to_treatment:
    type: string
    inputBinding:
      position: 4
    doc: "Two options, “norm” denotes normalization of control to\ntreatment data, “non” skips this behavior. \"norm\" is \necommended unless experimental and control data are\nalready rigorously normalized to each other (e.g. via\nspike-in). Default is \"non\". \n"
  peakcalling_mode:
    type: string
    inputBinding:
      position: 5
    doc: |
      Two options, “relaxed” uses a total signal threshold
      between the knee and peak of the total signal curve, and
      corresponds to the “relaxed” mode described in the text,
      whereas “stringent” uses the peak of the curve, and
      corresponds to “stringent” mode. Default is "relaxed".
  output_prefix:
    type: string
    inputBinding:
      position: 6
    doc: |
      Basename of input file that SEACR will use to name the
      output tsv file: <output_prefix>.<peakcalling_mode>.bed
outputs:
  peak_tsv_file:
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + '.' + inputs.peakcalling_mode + '.bed')
    doc: |
      SEACR peak calls in bed formatted file.
  log_file_stdout:
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + '.' + inputs.peakcalling_mode + '.seacr.stdout')
    doc: |
      log for stdout for seacr call
  log_file_stderr:
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + '.' + inputs.peakcalling_mode + '.seacr.stderr')
    doc: |
      log for stderr for seacr call
baseCommand:
- bash
- -c
stdout: $(inputs.output_prefix + '.' + inputs.peakcalling_mode + '.seacr.stdout')
stderr: $(inputs.output_prefix + '.' + inputs.peakcalling_mode + '.seacr.stderr')
doc: "Tool runs peak calling using SEACR_1.3.sh command. When running in stringent mode\n(default), if zero peaks are called from your input data there will be no output\nbed file, and the tool script will fail.\n\nSEACR: Sparse Enrichment Analysis for CUT&RUN and CUT&TAG sequencing library data.\n        Input bedGraph may be from processed single or paired-end library data\n\n  Usage: bash SEACR_1.3.sh <experimental bedgraph>.bg [<control bedgraph>.bg | <FDRthreshold>] [norm | non] [relaxed | stringent] output prefix\n\n  Description of input fields:\n\n      Field 1:    Target data bedgraph file in UCSC bedgraph format\n                  (https://genome.ucsc.edu/goldenpath/help/bedgraph.html)\n                  that omits regions containing 0 signal.\n      Field 2:    Control (IgG) data bedgraph file to generate an empirical\n                  threshold for peak calling. Alternatively, a numeric\n                  threshold n between 0 and 1 returns the top n fraction of\n                  peaks based on total signal within peaks.\n      Field 3:    “norm” denotes normalization of control to target data,\n                  “non” skips this behavior. norm is recommended unless\n                  experimental and control data are already rigorously\n                  normalized to each other (e.g. via spike-in).\n      Field 4:    “relaxed” uses a total signal threshold between the knee\n                  and peak of the total signal curve, and corresponds to the\n                  “relaxed” mode described in the text, whereas “stringent”\n                  uses the peak of the curve, and corresponds to “stringent”\n                  mode.\n      Field 5:    Output prefix\n\n  Output file:\n  <output prefix>.<mode>.bed (Bed file of enriched regions)\n\n  Output data structure: \n  <chr>\t<start>\t<end>\t<AUC>\t<max signal>\t<max signal region>\n\n  Description of output fields:\n  Field 1: Chromosome\n  Field 2: Start coordinate\n  Field 3: End coordinate\n  Field 4: Total signal contained within denoted coordinates\n  Field 5: Maximum bedgraph signal attained at any base pair within denoted coordinates\n  Field 6: Region representing the farthest upstream and farthest downstream bases within\n            the denoted coordinates that are represented by the maximum bedgraph signal\n\n  Examples:\n  bash SEACR_1.3.sh target.bedgraph IgG.bedgraph norm relaxed output\n  Calls enriched regions in target data using normalized IgG control track with relaxed threshold\n\n  bash SEACR_1.3.sh target.bedgraph IgG.bedgraph non relaxed output\n  Calls enriched regions in target data using non-normalized IgG control track with relaxed threshold\n  bash SEACR_1.3.sh target.bedgraph 0.01 non relaxed output\n  Calls enriched regions in target data by selecting the top 1% of regions by area under the curve (AUC)\n\n\nReferences:\n    Link to SEACR publication: https://doi.org/10.1186/s13072-019-0287-4\n    Link to SEACR github repo: https://github.com/FredHutch/SEACR\n        \"To address the problem of oversensitivity in CUT&RUN peak calling, we\n        developed Sparse Enrichment Analysis for CUT&RUN (SEACR), a peak calling\n        algorithm that enforces precision from sparse data by quantifying the\n        global distribution of background signal and using it to set a stringent\n        empirical threshold for peak identity.\"\n\n  Scaling factor based on the following publication:\n      Meers MP, Tenenbaum D, Henikoff S. (2019). Peak calling by Sparse Enrichment\n      Analysis for CUT&RUN chromatin profiling. Epigenetics and Chromatin 12(1):42.\n      https://doi.org/10.1186/s13072-019-0287-4\n  And protocol, step 16:\n      https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=16"
label: seacr
