cwlVersion: v1.0
class: CommandLineTool
hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/manorm:v0.0.4
inputs:
  peak_file_first:
    type: File
    inputBinding:
      prefix: --p1
    doc: Peaks file of sample 1
  peak_file_second:
    type: File
    inputBinding:
      prefix: --p2
    doc: Peaks file of sample 2
  peak_format:
    type:
    - 'null'
    - type: enum
      name: peak_format
      symbols:
      - bed
      - bed3-summit
      - macs
      - macs2
      - narrowpeak
      - broadpeak
    inputBinding:
      prefix: --pf
    doc: |
      "The format of peak files.
       Default BED"
  read_file_first:
    type: File
    inputBinding:
      prefix: --r1
    doc: Reads file of sample 1
  read_file_second:
    type: File
    inputBinding:
      prefix: --r2
    doc: Reads file of sample 2
  read_format:
    type:
    - 'null'
    - type: enum
      name: read_format
      symbols:
      - bed
      - bedpe
      - sam
      - bam
    inputBinding:
      prefix: --rf
    doc: |
      "The format of read files.
       Default BED"
  shift_size_first:
    type: int?
    inputBinding:
      prefix: --s1
    doc: |
      "Reads shift size of sample 1. This value is used to shift reads towards 3' direction
       to determine the precise binding site. Set as half of the fragment length.
       Default 100"
  shift_size_second:
    type: int?
    inputBinding:
      prefix: --s2
    doc: |
      "Reads shift size of sample 2. This value is used to shift reads towards 5' direction
       to determine the precise binding site. Set as half of the fragment length.
       Default 100"
  sample_name_first:
    type: string?
    inputBinding:
      prefix: --n1
    doc: "\"Name of sample 1, which is used in output files. If not specified, \n the name of the peak file will be used as the sample name\"\n"
  sample_name_second:
    type: string?
    inputBinding:
      prefix: --n2
    doc: |
      "Name of sample 2, which is used in output files. If not specified,
       the name of the peak file will be used as the sample name"
  simulations_number:
    type: int?
    inputBinding:
      prefix: --n-random
    doc: |
      "Number of random simulations to test the enrichment of peak
       overlap between two samples. Set to 0 to disable the testing.
       Default: 10"
  m_value_cutoff:
    type: float?
    inputBinding:
      prefix: --m-cutoff
    doc: |
      "Absolute M-value (log2-ratio) cutoff to define biased (differential binding) peaks.
       Default: 1.0"
  p_value_cutoff:
    type: float?
    inputBinding:
      prefix: --p-cutoff
    doc: |
      "P-value cutoff to define biased peaks.
       Default: 0.01"
  paired_end:
    type: boolean?
    inputBinding:
      prefix: --pe
    doc: |
      "The middle point of each read pair is used to represent the genomic locus
      of underlying DNA fragment. --s1 and --s2 are ignored with this option on"
  window_size:
    type: int?
    inputBinding:
      prefix: -w
    doc: |
      "Window size to count reads and calculate read densities. 2000 is recommended for
      sharp histone marks like H3K4me3 and H3K27ac, and 1000 for TFs or DNase-seq.
      Default: 2000"
  summit_distance:
    type: int?
    inputBinding:
      prefix: --summit-dis
    doc: |
      "Overlapping common peaks with summit-to-summit distance beyond this are excluded in model fitting.
      This option is used to exclude common peaks that only overlap on the edge of each other.
      Default: -w/--window-size/4"
outputs:
  ma_values_file:
    type: File
    outputBinding:
      glob: '*_all_MAvalues.xls'
    doc: |
      "File contains the M-A values and normalized read density of each
       peak, common peaks from two samples are merged together.
       Coordinates in .xls file is under 1-based coordinate-system"
  above_m_cutoff_peak_file:
    type: File
    outputBinding:
      glob: output_filters/*_M_above_*_biased_peaks.bed
    doc: Above M-value cutoff peak file
  below_m_cutoff_peak_file:
    type: File
    outputBinding:
      glob: output_filters/*_M_below_*_biased_peaks.bed
    doc: Below M-value cutoff peak file
  unbiased_peak_file:
    type: File
    outputBinding:
      glob: output_filters/*_unbiased_peaks.bed
    doc: Unbiased peak file
  m_values_wig_file:
    type: File
    outputBinding:
      glob: output_tracks/*_M_values.wig
    doc: Genome track file for M-values
  a_values_wig_file:
    type: File
    outputBinding:
      glob: output_tracks/*_A_values.wig
    doc: Genome track file for A-values
  p_values_wig_file:
    type: File
    outputBinding:
      glob: output_tracks/*_P_values.wig
    doc: Genome track file for P-values
  ma_before_normalization_plot:
    type: File
    outputBinding:
      glob: output_figures/*_MA_plot_before_normalization.*
    doc: MA-values before normalization plot
  ma_after_normalization_plot:
    type: File
    outputBinding:
      glob: output_figures/*_MA_plot_after_normalization.*
    doc: MA-values after normalization plot
  ma_with_P_value_plot:
    type: File
    outputBinding:
      glob: output_figures/*_MA_plot_with_P_value.*
    doc: MA-values with P-values plot
  read_density_on_common_peaks_plot:
    type: File
    outputBinding:
      glob: output_figures/*_read_density_on_common_peaks.*
    doc: Read density on common peaks plot
  stdout_log:
    type: stdout
  stderr_log:
    type: stderr
baseCommand:
- manorm
stderr: manorm_stderr.log
stdout: manorm_stdout.log
label: MAnorm - quantitative comparison of ChIP-Seq data
doc: |
  --wa argument is skipped
