cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/seacr:1.3--hdfd78af_2


inputs:

  treatment_bedgraph:
    type: File?
    inputBinding:
      position: 1
    doc: |
      A bedgraph formatted file from paired-end sequencing as
      input, which can be generated from read pair BED files
      (i.e. BED coordinates eflecting the 5' and 3' termini
      of each read pair) using bedtools genomecov with the
      "-bg" flag.

  numeric_threshold:
    type: float
    inputBinding:
      position: 2
    doc: |
      A numeric threshold n between 0 and 1 returns the top n
      fraction of peaks based on total signal within peaks.
      Default is "0.01".

  norm_control_to_treatment:
    type: string
    inputBinding:
      position: 3
    doc: |
      Two options, “norm” denotes normalization of control to
      treatment data, “non” skips this behavior. "norm" is 
      ecommended unless experimental and control data are
      already rigorously normalized to each other (e.g. via
      spike-in). Default is "non". 

  peakcalling_mode:
    type: string
    inputBinding:
      position: 4
    doc: |
      Two options, “relaxed” uses a total signal threshold
      between the knee and peak of the total signal curve, and
      corresponds to the “relaxed” mode described in the text,
      whereas “stringent” uses the peak of the curve, and
      corresponds to “stringent” mode. Default is "stringent".

  output_prefix:
    type: string
    inputBinding:
      position: 5
    doc: |
      Basename of input file that SEACR will use to name the
      output tsv file: <output_prefix>.<peakcalling_mode>.bed

outputs:

  peak_tsv_file:
    type: File
    outputBinding:
      glob: $(inputs.output_prefix + '.stringent.bed')
    doc: |
      SEACR peak calls in bed formatted file.

baseCommand: [SEACR_1.3.sh]


doc: |
  Tool runs peak calling using SEACR_1.3.sh command.

  SEACR: Sparse Enrichment Analysis for CUT&RUN

    Usage: bash SEACR_1.3.sh <experimental bedgraph>.bg [<control bedgraph>.bg | <FDR threshold>] [norm | non] [relaxed | stringent] output prefix

    Description of input fields:

        Field 1:    Target data bedgraph file in UCSC bedgraph format
                    (https://genome.ucsc.edu/goldenpath/help/bedgraph.html)
                    that omits regions containing 0 signal.
        Field 2:    Control (IgG) data bedgraph file to generate an empirical
                    threshold for peak calling. Alternatively, a numeric
                    threshold n between 0 and 1 returns the top n fraction of
                    peaks based on total signal within peaks.
        Field 3:    “norm” denotes normalization of control to target data,
                    “non” skips this behavior. norm is recommended unless
                    experimental and control data are already rigorously
                    normalized to each other (e.g. via spike-in).
        Field 4:    “relaxed” uses a total signal threshold between the knee
                    and peak of the total signal curve, and corresponds to the
                    “relaxed” mode described in the text, whereas “stringent”
                    uses the peak of the curve, and corresponds to “stringent”
                    mode.
        Field 5:    Output prefix

    Output file:
    <output prefix>.stringent.bed (Bed file of enriched regions)

    Output data structure: 
    <chr>	<start>	<end>	<AUC>	<max signal>	<max signal region>

    Description of output fields:
    Field 1: Chromosome
    Field 2: Start coordinate
    Field 3: End coordinate
    Field 4: Total signal contained within denoted coordinates
    Field 5: Maximum bedgraph signal attained at any base pair within denoted coordinates
    Field 6: Region representing the farthest upstream and farthest downstream bases within the denoted coordinates that are represented by the maximum bedgraph signal

    Examples:
    bash SEACR_1.3.sh target.bedgraph IgG.bedgraph norm stringent output
    Calls enriched regions in target data using normalized IgG control track with stringent threshold

    bash SEACR_1.3.sh target.bedgraph IgG.bedgraph non relaxed output
    Calls enriched regions in target data using non-normalized IgG control track with relaxed threshold
    bash SEACR_1.3.sh target.bedgraph 0.01 non stringent output
    Calls enriched regions in target data by selecting the top 1% of regions by area under the curve (AUC)


  References:
      Link to SEACR publication: https://doi.org/10.1186/s13072-019-0287-4
      Link to SEACR github repo: https://github.com/FredHutch/SEACR
          "To address the problem of oversensitivity in CUT&RUN peak calling, we
          developed Sparse Enrichment Analysis for CUT&RUN (SEACR), a peak calling
          algorithm that enforces precision from sparse data by quantifying the
          global distribution of background signal and using it to set a stringent
          empirical threshold for peak identity."

    Scaling factor based on the following publication:
        Meers MP, Tenenbaum D, Henikoff S. (2019). Peak calling by Sparse Enrichment
        Analysis for CUT&RUN chromatin profiling. Epigenetics and Chromatin 12(1):42.
        https://doi.org/10.1186/s13072-019-0287-4
    And protocol, step 16:
        https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-e6nvw93x7gmk/v1?step=16