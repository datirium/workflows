cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
      return inputs.output_filename == "" ? inputs.bed_file.nameroot+".peaks.bed":inputs.output_filename;
    };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/clip-toolkit:v0.0.1


inputs:

  big:
    type: boolean?
    inputBinding:
      prefix: "-big"
    doc: "Big input file"

  separate_strands:
    type: boolean?
    inputBinding:
      prefix: "-ss"
    doc: "Separate the two strands"

  valley_seeking:
    type: boolean?
    inputBinding:
      prefix: "--valley-seeking"
    doc: "Find candidate peaks by valley seeking"

  valley_depth:
    type: float?
    inputBinding:
      prefix: "--valley-depth"
    doc: "Depth of valley if valley seeking (0.9)"

  cluster_boundaries:
    type: string?
    inputBinding:
      prefix: "--out-boundary"
    doc: "Output cluster boundaries"

  half_peak_height_boundaries:
    type: string?
    inputBinding:
      prefix: "--out-half-PH"
    doc: "Output half peak height boundaries"

  min_peak_height:
    type: int?
    inputBinding:
      prefix: "--minPH"
    doc: "Min peak height (2)"

  max_peak_height:
    type: int?
    inputBinding:
      prefix: "--maxPH"
    doc: "Max peak height to calculate p-value(-1, no limit if < 0)"

  skip_out_of_range_peaks:
    type: boolean?
    inputBinding:
      prefix: "--skip-out-of-range-peaks"
    doc: "Skip peaks with PH > maxPH"

  peak_gap:
    type: int?
    inputBinding:
      prefix: "--gap"
    doc: "Merge cluster peaks closer than the gap (-1, no merge if < 0)"

  peak_id_prefix:
    type: string?
    inputBinding:
      prefix: "--prefix"
    doc: "Prefix of peak id (Peak)"

  default_gene_bed:
    type:
    - "null"
    - type: enum
      symbols: ["mm10","hg19"]
    inputBinding:
      prefix: "--dbkey"
    doc: "Species to retrieve the default gene bed file (mm10|hg19)"

  gene_bed_file:
    type: File?
    inputBinding:
      prefix: "--gene"
    doc: "Custom gene bed file for scan statistics (will override --dbkey)"

  use_expression:
    type: boolean?
    inputBinding:
      prefix: "--use-expr"
    doc: "Use expression levels given in the score column in the custom gene bed file for normalization"

  p_value_threshold:
    type: float?
    inputBinding:
      prefix: "-p"
    doc: "threshold of p-value to call peak (0.01)"

  multi_test:
    type: boolean?
    inputBinding:
      prefix: "--multi-test"
    doc: "Do Bonferroni multiple test correction"

  bed_file:
    type: File
    inputBinding:
      position: 30
    doc: "IInput BED file of unique CLIP tags"

  output_filename:
    type: string?
    inputBinding:
      position: 31
      valueFrom: $(default_output_filename())
    default: ""
    doc: "Output BED file name"


outputs:

  peaks_bed:
    type: File
    outputBinding:
      glob: $(default_output_filename())
    doc: "BED file of called peaks"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: [tag2peak.pl]
stdout: tag2peak_stdout.log
stderr: tag2peak_stderr.log


doc: |
    detecting peaks from CLIP data
    Usage: tag2peak.pl [options] <tag.bed> <peak.bed>
     <tag.bed> : BED file of unique CLIP tags, input
     <peak.bed>: BED file of called peaks, output
    Options:
     -big                   : big input file
     -ss                    : separate the two strands
     --valley-seeking       : find candidate peaks by valley seeking
     --valley-depth [float] : depth of valley if valley seeking (0.9)
     --out-boundary [string]: output cluster boundaries
     --out-half-PH  [string]: output half peak height boundaries
     --dbkey        [string]: species to retrieve the default gene bed file (mm10|hg19)
     --gene         [string]: custom gene bed file for scan statistics (will override --dbkey)
     --use-expr             : use expression levels given in the score column in the custom gene bed file for normalization
     -p             [float] : threshold of p-value to call peak (0.01)
     --multi-test           : do Bonferroni multiple test correction
     -minPH         [int]   : min peak height (2)
     -maxPH         [int]   : max peak height to calculate p-value(-1, no limit if < 0)
     --skip-out-of-range-peaks: skip peaks with PH > maxPH
     -gap           [int]   : merge cluster peaks closer than the gap (-1, no merge if < 0)
     --prefix       [string]: prefix of peak id (Peak)
     -c             [dir]   : cache dir
     --keep-cache           : keep cache when the job done
     -v                     : verbose
