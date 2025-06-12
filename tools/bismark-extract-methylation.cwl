cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: DockerRequirement
  dockerPull: biowardrobe2/bismark:v0.0.2
inputs:
  genome_folder:
    type: Directory
    label: Genome folder
    doc: |
      Genome folder with FASTA (fa, fasta) files.
      Bismark generated indices folder can be used also
    inputBinding:
      position: 4
      prefix: --genome_folder
  bam_file:
    type: File
    label: BAM alignment file
    doc: |
      Bismark generated BAM alignment file
    inputBinding:
      position: 5
  processes:
    type: int?
    label: Number of Bismark instances to run
    doc: |
      Set the number of parallel Bismark instances to run concurrently.
      Each Bismark instance simultainously runs the methylation extractor,
      samtools stream and GZIP streams
    inputBinding:
      position: 1
      prefix: --multicore
  ignore_r1_bp:
    type: int?
    label: Number of bp to ignore from the 5' end of read 1
    doc: |
      Ignore the first <int> bp from the 5' end of Read 1 (or single-end alignment
      files) when processing the methylation call string. This can remove e.g. a
      restriction enzyme site at the start of each read or any other source of
      bias (such as PBAT-Seq data).
    inputBinding:
      position: 2
      prefix: --ignore
  ignore_r2_bp:
    type: int?
    label: Number of bp to ignore from the 5' end of read 2
    doc: |
      Ignore the first <int> bp from the 5' end of Read 2 of paired-end sequencing
      results only. Since the first couple of bases in Read 2 of BS-Seq experiments
      show a severe bias towards non-methylation as a result of end-repairing
      sonicated fragments with unmethylated cytosines (see M-bias plot), it is
      recommended that the first couple of bp of Read 2 are removed before
      starting downstream analysis. Please see the section on M-bias plots in the
      Bismark User Guide for more details.
    inputBinding:
      position: 3
      prefix: --ignore_r2
outputs:
  chg_context_file:
    type: File
    label: CHG methylation call
    doc: CHG methylation call
    outputBinding:
      glob: CHG_context*
  chh_context_file:
    type: File
    label: CHH methylation call
    doc: CHH methylation call
    outputBinding:
      glob: CHH_context*
  cpg_context_file:
    type: File
    label: CpG methylation call
    doc: CpG methylation call
    outputBinding:
      glob: CpG_context*
  mbias_plot:
    type: File
    label: Methylation bias plot
    doc: QC data showing methylation bias across read lengths
    outputBinding:
      glob: '*M-bias.txt'
  mbias_plot_png:
    type: File[]
    label: Methylation bias plot (PNG)
    doc: QC data showing methylation bias across read lengths
    outputBinding:
      glob: '*.png'
  bedgraph_coverage_file:
    type: File
    label: Methylation statuses bedGraph coverage file
    doc: Coverage text file summarising cytosine methylation values in bedGraph format (tab-delimited; 0-based start coords, 1-based end coords)
    outputBinding:
      glob: '*bedGraph.gz'
  bismark_coverage_file:
    type: File
    label: Methylation statuses Bismark coverage file
    doc: Coverage text file summarising cytosine methylation values in Bismark format (tab-delimited, 1-based genomic coords)
    outputBinding:
      glob: '*bismark.cov.gz'
  genome_wide_methylation_report:
    type: File
    label: Genome-wide cytosine methylation report
    doc: Genome-wide methylation report for all cytosines in the genome
    outputBinding:
      glob: '*CpG_report.txt'
  splitting_report:
    type: File
    label: Methylation extraction log
    doc: Log file giving summary statistics about methylation extraction
    outputBinding:
      glob: '*splitting_report.txt'
baseCommand:
- bismark_methylation_extractor
- --comprehensive
- --bedgraph
- --cytosine_report
doc: |
  bismark_methylation_extractor script operates on Bismark result files and extracts the methylation call
  for every single C analysed. The position of every single C will be written out to a new output file,
  depending on its context (CpG, CHG or CHH), whereby methylated Cs will be labelled as forward reads (+),
  non-methylated Cs as reverse reads (-).

  Options used:
  --comprehensive
    If strand-specific methylation is not of interest, all available methylation information can be pooled
    into a single context-dependent file (information from any of the four strands will be pooled). This
    will default to three output files (CpG-context, CHG-context and CHH-context).
  --bedgraph
    The Bismark methylation extractor can optionally also output a file in bedGraph format which uses 0-based
    genomic start and 1- based end coordinates. It will be sorted by chromosomal coordinates and looks like this:
    <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>
  --genome_folder
    Bismark methylation extractor can also output a genome-wide cytosine methylation report. It is also sorted by
    chromosomal coordinates but also contains the sequence context and is in the following format:
    <chromosome> <position> <strand> <count methylated> <count unmethylated> <C-context> <trinucleotide context>
    The main difference to the bedGraph or coverage output is that every cytosine on both the top and bottom strands
    will be considered irrespective of whether they were actually covered by any reads in the experiment or not.
    For this to work one has to also specify the genome that was used for the Bismark alignments using the
    option --genome_folder <path>. As for the bedGraph mode, this will only consider cytosines in CpG context.
  --cytosine_report
    After the conversion to bedGraph has completed, the option '--cytosine_report' produces a
    genome-wide methylation report for all cytosines in the genome. By default, the output uses 1-based
    chromosome coordinates (zero-based start coords are optional) and reports CpG context only (all
    cytosine context is optional). The output considers all Cs on both forward and reverse strands and
    reports their position, strand, trinucleotide content and methylation state (counts are 0 if not
    covered). The cytosine report conversion step is performed by the external module
    'coverage2cytosine'; this script needs to reside in the same folder as the bismark_methylation_extractor
    itself.
label: bismark-extract-methylation
