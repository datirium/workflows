cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/manorm:v0.0.4


inputs:

  peak_file_first:
    type: File
    inputBinding:
      prefix: "--p1"
    doc: "Peaks file of sample 1"

  peak_file_second:
    type: File
    inputBinding:
      prefix: "--p2"
    doc: "Peaks file of sample 2"

  peak_format:
    type:
      - "null"
      - type: enum
        name: "peak_format"
        symbols: ["bed", "bed3-summit", "macs", "macs2", "narrowpeak", "broadpeak"]
    inputBinding:
      prefix: "--pf"
    doc: |
      "The format of peak files.
       Default BED"

  read_file_first:
    type: File
    inputBinding:
      prefix: "--r1"
    doc: "Reads file of sample 1"

  read_file_second:
    type: File
    inputBinding:
      prefix: "--r2"
    doc: "Reads file of sample 2"

  read_format:
    type:
      - "null"
      - type: enum
        name: "read_format"
        symbols: ["bed", "bedpe", "sam", "bam"]
    inputBinding:
      prefix: "--rf"
    doc: |
      "The format of read files.
       Default BED"

  shift_size_first:
    type: int?
    inputBinding:
      prefix: "--s1"
    doc: |
      "Reads shift size of sample 1. This value is used to shift reads towards 3' direction
       to determine the precise binding site. Set as half of the fragment length.
       Default 100"

  shift_size_second:
    type: int?
    inputBinding:
      prefix: "--s2"
    doc: |
      "Reads shift size of sample 2. This value is used to shift reads towards 5' direction
       to determine the precise binding site. Set as half of the fragment length.
       Default 100"

  sample_name_first:
    type: string?
    inputBinding:
      prefix: "--n1"
    doc: |
      "Name of sample 1, which is used in output files. If not specified, 
       the name of the peak file will be used as the sample name"

  sample_name_second:
    type: string?
    inputBinding:
      prefix: "--n2"
    doc: |
      "Name of sample 2, which is used in output files. If not specified,
       the name of the peak file will be used as the sample name"

  simulations_number:
    type: int?
    inputBinding:
      prefix: "--n-random"
    doc: |
      "Number of random simulations to test the enrichment of peak
       overlap between two samples. Set to 0 to disable the testing.
       Default: 10"

  m_value_cutoff:
    type: float?
    inputBinding:
      prefix: "--m-cutoff"
    doc: |
      "Absolute M-value (log2-ratio) cutoff to define biased (differential binding) peaks.
       Default: 1.0"

  p_value_cutoff:
    type: float?
    inputBinding:
      prefix: "--p-cutoff"
    doc: |
      "P-value cutoff to define biased peaks.
       Default: 0.01"

  paired_end:
    type: boolean?
    inputBinding:
      prefix: "--pe"
    doc: |
      "The middle point of each read pair is used to represent the genomic locus
      of underlying DNA fragment. --s1 and --s2 are ignored with this option on"

  window_size:
    type: int?
    inputBinding:
      prefix: "-w"
    doc: |
      "Window size to count reads and calculate read densities. 2000 is recommended for
      sharp histone marks like H3K4me3 and H3K27ac, and 1000 for TFs or DNase-seq.
      Default: 2000"
  
  summit_distance:
    type: int?
    inputBinding:
      prefix: "--summit-dis"
    doc: |
      "Overlapping common peaks with summit-to-summit distance beyond this are excluded in model fitting.
      This option is used to exclude common peaks that only overlap on the edge of each other.
      Default: -w/--window-size/4"


outputs:

  ma_values_file:
    type: File
    outputBinding:
      glob: "*_all_MAvalues.xls"
    doc: |
      "File contains the M-A values and normalized read density of each
       peak, common peaks from two samples are merged together.
       Coordinates in .xls file is under 1-based coordinate-system"

  above_m_cutoff_peak_file:
    type: File
    outputBinding:
      glob: "output_filters/*_M_above_*_biased_peaks.bed"
    doc: "Above M-value cutoff peak file"

  below_m_cutoff_peak_file:
    type: File
    outputBinding:
      glob: "output_filters/*_M_below_*_biased_peaks.bed"
    doc: "Below M-value cutoff peak file"

  unbiased_peak_file:
    type: File
    outputBinding:
      glob: "output_filters/*_unbiased_peaks.bed"
    doc: "Unbiased peak file"

  m_values_wig_file:
    type: File
    outputBinding:
      glob: "output_tracks/*_M_values.wig"
    doc: "Genome track file for M-values"
    
  a_values_wig_file:
    type: File
    outputBinding:
      glob: "output_tracks/*_A_values.wig"
    doc: "Genome track file for A-values"

  p_values_wig_file:
    type: File
    outputBinding:
      glob: "output_tracks/*_P_values.wig"
    doc: "Genome track file for P-values"

  ma_before_normalization_plot:
    type: File
    outputBinding:
      glob: "output_figures/*_MA_plot_before_normalization.*"
    doc: "MA-values before normalization plot"

  ma_after_normalization_plot:
    type: File
    outputBinding:
      glob: "output_figures/*_MA_plot_after_normalization.*"
    doc: "MA-values after normalization plot"

  ma_with_P_value_plot:
    type: File
    outputBinding:
      glob: "output_figures/*_MA_plot_with_P_value.*"
    doc: "MA-values with P-values plot"

  read_density_on_common_peaks_plot:
    type: File
    outputBinding:
      glob: "output_figures/*_read_density_on_common_peaks.*"
    doc: "Read density on common peaks plot"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["manorm"]
stderr: manorm_stderr.log
stdout: manorm_stdout.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/manorm-metadata.yaml

label: "MAnorm - quantitative comparison of ChIP-Seq data"
s:name: "MAnorm - quantitative comparison of ChIP-Seq data"
s:alternateName: "MAnorm is a robust model for quantitative comparison of ChIP-Seq data sets of TFs (transcription factors) or epigenetic modifications"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/manorm.cwl
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
  --wa argument is skipped

s:about: |
  usage: manorm [-h] [-v] --p1 FILE --p2 FILE [--pf FORMAT] --r1 FILE --r2 FILE [--rf FORMAT] [--n1 NAME] [--n2 NAME] [--s1 N] [--s2 N] [--pe] [-w LENGTH] [--summit-dis DISTANCE]
                [--n-random NUM] [-m FLOAT] [-p FLOAT] [-o DIR] [--wa] [--verbose]

  MAnorm -- A robust model for quantitative comparison of ChIP-seq data sets

  Citation:
    Shao Z, Zhang Y, Yuan GC, Orkin SH, Waxman DJ. MAnorm: a robust model 
    for quantitative comparison of ChIP-Seq data sets. Genome biology. 
    2012 Mar 16;13(3):R16. https://doi.org/10.1186/gb-2012-13-3-r16

  optional arguments:
    -h, --help            show this help message and exit
    -v, --version         show program's version number and exit
    --verbose             Enable verbose log messages.

  Input Options:
    --p1 FILE, --peak1 FILE
                          Peak file of sample 1.
    --p2 FILE, --peak2 FILE
                          Peak file of sample 2.
    --pf FORMAT, --peak-format FORMAT
                          Format of the peak files. Support ['bed', 'bed3-summit', 'macs', 'macs2', 'narrowpeak', 'broadpeak']. Default: bed
    --r1 FILE, --read1 FILE
                          Read file of sample 1.
    --r2 FILE, --read2 FILE
                          Read file of sample 2.
    --rf FORMAT, --read-format FORMAT
                          Format of the read files. Support ['bed', 'bedpe', 'sam', 'bam']. Default: bed
    --n1 NAME, --name1 NAME
                          Name of sample 1. If not specified, the peak file name will be used.
    --n2 NAME, --name2 NAME
                          Name of sample 2. If not specified, the peak file name will be used.

  Reads Manipulation:
    --s1 N, --shiftsize1 N
                          Single-end reads shift size for sample 1. Reads are shifted by `N` bp towards 3' direction and the 5' end of each shifted read is used to represent the
                          genomic locus of the DNA fragment. Set to 0.5 * fragment size of the ChIP-seq library. Default: 100
    --s2 N, --shiftsize2 N
                          Single-end reads shift size for sample 2. Default: 100
    --pe, --paired-end    Paired-end mode. The middle point of each read pair is used to represent the genomic locus of the DNA fragment. If specified, `--s1` and `--s2` will be
                          ignored.

  Normalization Model Options:
    -w LENGTH, --window-size LENGTH
                          Window size to count reads and calculate read densities. Set to 2000 is recommended for sharp histone marks like H3K4me3 or H3K27ac and 1000 for TFs or
                          DNase-seq. Default: 2000
    --summit-dis DISTANCE
                          Summit-to-summit distance cutoff for overlapping common peaks. Common peaks with summit distance beyond the cutoff are excluded in model fitting. Default:
                          `window size` / 4
    --n-random NUM        Number of random simulations to test the enrichment of peak overlap between the specified samples. Set to 0 to disable the testing. Default: 10

  Output Options:
    -m FLOAT, --m-cutoff FLOAT
                          Absolute M-value (log2-ratio) cutoff to define the biased (differential binding) peaks. Default: 1.0
    -p FLOAT, --p-cutoff FLOAT
                          P-value cutoff to define the biased peaks. Default: 0.01
    -o DIR, --output-dir DIR
                          Output directory. Default: Current working directory
    --wa, --write-all     Write two extra output files containing the results of the original (unmerged) peaks.
