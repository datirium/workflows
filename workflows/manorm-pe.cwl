cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  first_chipseq_sample:
    - "chipseq-pe.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-pe.cwl"
    - "cutandrun-pe.cwl"
  second_chipseq_sample:
    - "chipseq-pe.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-pe.cwl"
    - "cutandrun-pe.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  peak_file_first:
    type: File
    format: "http://edamontology.org/format_3468"
    label: "ChIP-Seq PE sample 1"
    doc: |
      XLS peak file from sample 1, formatted as MACS2 output
    'sd:upstreamSource': "first_chipseq_sample/macs2_called_peaks"
    'sd:localLabel': true

  peak_file_second:
    type: File
    format: "http://edamontology.org/format_3468"
    label: "ChIP-Seq PE sample 2"
    doc: |
      XLS peak file from sample 2, formatted as MACS2 output
    'sd:upstreamSource': "second_chipseq_sample/macs2_called_peaks"
    'sd:localLabel': true

  broad_peak_file_first:
    type: File?
    format: "http://edamontology.org/format_3614"
    label: "ChIP-Seq PE sample 1"
    doc: |
      Broad peak file from sample 1
    'sd:upstreamSource': "first_chipseq_sample/macs2_broad_peaks"
    'sd:localLabel': true

  broad_peak_file_second:
    type: File?
    format: "http://edamontology.org/format_3614"
    label: "ChIP-Seq PE sample 2"
    doc: |
      Broad peak file from sample 2
    'sd:upstreamSource': "second_chipseq_sample/macs2_broad_peaks"
    'sd:localLabel': true

  bam_file_first:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "BAM file from sample 1"
    doc: |
      BAM alignment file from sample 1
    'sd:upstreamSource': "first_chipseq_sample/bambai_pair"

  bam_file_second:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "BAM file from sample 2"
    doc: |
      BAM alignment file from sample 2
    'sd:upstreamSource': "second_chipseq_sample/bambai_pair"

  annotation_file:
    type: File
    label: "Annotation file"
    format: "http://edamontology.org/format_3475"
    doc: |
      Tab-separated annotation file
    'sd:upstreamSource': "first_chipseq_sample/genome_indices/annotation"

  shift_size_first:
    type: int?
    default: 100
    label: "Reads shift size for sample 1"
    doc: |
      Reads shift size of sample 1. This value is used to shift reads towards 3' direction
      to determine the precise binding site. Set as half of the fragment length.
      Default 100
    'sd:layout':
      advanced: true

  shift_size_second:
    type: int?
    default: 100
    label: "Reads shift size for sample 2"
    doc: |
      Reads shift size of sample 2. This value is used to shift reads towards 5' direction
      to determine the precise binding site. Set as half of the fragment length.
      Default 100
    'sd:layout':
      advanced: true

  m_value_cutoff:
    type: float?
    default: 1
    label: "M-value (log2-ratio) cutoff"
    doc: "Absolute M-value (log2-ratio) cutoff to define biased (differential binding) peaks. Default: 1.0"
    'sd:layout':
      advanced: true

  p_value_cutoff:
    type: float?
    default: 0.01
    label: "P-value cutoff"
    doc: "P-value cutoff to define biased peaks. Default: 0.01"
    'sd:layout':
      advanced: true

  window_size:
    type: int?
    default: 2000
    label: "Window size (2000 is recommended for sharp histone marks like H3K4me3 and H3K27ac)"
    doc: |
      Window size to count reads and calculate read densities. 2000 is recommended for
      sharp histone marks like H3K4me3 and H3K27ac, and 1000 for TFs or DNase-seq.
      Default: 2000
    'sd:layout':
      advanced: true

  promoter_dist:
    type: int?
    default: 1000
    label: "Promoter distance, bp"
    doc: |
      Max distance from gene TSS (in both direction) overlapping
      which the peak will be assigned to the promoter region.
      Default: 1000 bp
    'sd:layout':
      advanced: true

  upstream_dist:
    type: int?
    default: 20000
    label: "Upstream distance, bp"
    doc: | 
      Max distance from the promoter (only in upstream direction) overlapping
      which the peak will be assigned to the upstream region.
      Default: 20,000 bp
    'sd:layout':
      advanced: true


outputs:

  common_peak_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "MAnorm common peak file with assigned genes"
    doc: |
      "File contains nearest gene information, the M-A values and normalized read
       density of each peak, common peaks from two samples are merged together.
       Coordinates in a result file is under 1-based coordinate-system"
    outputSource: restore_columns/output_file
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'Differential Peak Calling'
        Title: 'MAnorm Common Peak Results'

  above_m_cutoff_peak_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Above M-value cutoff peak file"
    doc: "Above M-value cutoff peak file"
    outputSource: manorm/above_m_cutoff_peak_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Above M-value cutoff peaks"
        height: 120

  below_m_cutoff_peak_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Below M-value cutoff peak file"
    doc: "Below M-value cutoff peak file"
    outputSource: manorm/below_m_cutoff_peak_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Below M-value cutoff peaks"
        height: 120

  unbiased_peak_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Unbiased peak file"
    doc: "Unbiased peak file"
    outputSource: manorm/unbiased_peak_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Unbiased peaks"
        height: 120

  m_values_wig_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Genome track file for M-values"
    doc: "Genome track file for M-values"
    outputSource: manorm/m_values_wig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "M-values"
        height: 120

  a_values_wig_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Genome track file for A-values"
    doc: "Genome track file for A-values"
    outputSource: manorm/a_values_wig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "A-values"
        height: 120

  p_values_wig_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Genome track file for P-values"
    doc: "Genome track file for P-values"
    outputSource: manorm/p_values_wig_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "P-values"
        height: 120

  ma_before_normalization_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    label: "MA-values before normalization plot"
    doc: "MA-values before normalization plot"
    outputSource: manorm/ma_before_normalization_plot
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'MA-values before normalization'

  ma_after_normalization_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    label: "MA-values after normalization plot"
    doc: "MA-values after normalization plot"
    outputSource: manorm/ma_after_normalization_plot
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'MA-values after normalization'

  ma_with_P_value_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    label: "MA-values with P-values plot"
    doc: "MA-values with P-values plot"
    outputSource: manorm/ma_with_P_value_plot
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'MA-values with P-values'

  read_density_on_common_peaks_plot:
    type: File
    format: "http://edamontology.org/format_3603"
    label: "Read density on common peaks plot"
    doc: "Read density on common peaks plot"
    outputSource: manorm/read_density_on_common_peaks_plot
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'Read density on common peaks'

  manorm_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "MAnorm stderr log"
    doc: "MAnorm stderr log"
    outputSource: manorm/stderr_log

  manorm_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "MAnorm stdout log"
    doc: "MAnorm stdout log"
    outputSource: manorm/stdout_log

steps:

  manorm:
    run: ../tools/manorm.cwl
    in:
      peak_file_first: 
        source:
        - peak_file_first          # [0]
        - peak_file_second         # [1]
        - broad_peak_file_first    # [2]
        - broad_peak_file_second   # [3]
        valueFrom: |
          ${
            if (self[2] && self[3]){
              return self[2];
            }
            else {
              return self[0];
            }
          }
      peak_file_second:
        source:
        - peak_file_first          # [0]
        - peak_file_second         # [1]
        - broad_peak_file_first    # [2]
        - broad_peak_file_second   # [3]
        valueFrom: |
          ${
            if (self[2] && self[3]){
              return self[3];
            }
            else {
              return self[1];
            }
          }
      peak_format:
        source:
        - peak_file_first          # [0]
        - peak_file_second         # [1]
        - broad_peak_file_first    # [2]
        - broad_peak_file_second   # [3]
        valueFrom: |
          ${
            if (self[2] && self[3]){
              return "broadpeak";
            }
            else {
              return "macs2";
            }
          }
      read_file_first: bam_file_first
      read_file_second: bam_file_second
      read_format:
        default: "bam"
      shift_size_first: shift_size_first
      shift_size_second: shift_size_second
      paired_end:
        default: true
      m_value_cutoff: m_value_cutoff
      p_value_cutoff: p_value_cutoff
      window_size: window_size
      sample_name_first:
        default: "sample_1"
      sample_name_second:
        default: "sample_2"
    out:
      - ma_values_file
      - above_m_cutoff_peak_file
      - below_m_cutoff_peak_file
      - unbiased_peak_file
      - m_values_wig_file
      - a_values_wig_file
      - p_values_wig_file
      - ma_before_normalization_plot
      - ma_after_normalization_plot
      - ma_with_P_value_plot
      - read_density_on_common_peaks_plot
      - stderr_log
      - stdout_log

  filter_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file: manorm/ma_values_file
      script:
        default: >
          cat $0 | grep -v "start" | awk
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t"NR"\t0\t0\t0\t0"}' > `basename $0`
    out: [output_file]

  assign_genes:
      run: ../tools/iaintersect.cwl
      in:
        input_filename: filter_columns/output_file
        annotation_filename: annotation_file
        promoter_bp: promoter_dist
        upstream_bp: upstream_dist
      out: [result_file]

  restore_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file: [assign_genes/result_file, manorm/ma_values_file]
      script:
        default: |
          cat $0 | grep -v "start" | sort -k 11n > sorted_iaintersect_result.tsv
          cat $1 | grep -v "start" > manorm_result.tsv
          echo -e "refseq_id\tgene_id\ttxStart\ttxEnd\tstrand\tchrom\tstart\tend\tlength\tregion\tsummit\tM_value\tA_value\tP_value\tPeak_Group\tnormalized_read_density_in_sample_1\tnormalized_read_density_in_sample_2" > `basename $0`;
          cat sorted_iaintersect_result.tsv | paste - manorm_result.tsv | cut -f 1-9,15,19-25 >> `basename $0`
          rm sorted_iaintersect_result.tsv manorm_result.tsv
    out: [output_file]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "MAnorm PE - quantitative comparison of ChIP-Seq paired-end data"
label: "MAnorm PE - quantitative comparison of ChIP-Seq paired-end data"
s:alternateName: "MAnorm is a robust model for quantitative comparison of ChIP-Seq paired-end data sets of TFs (transcription factors) or epigenetic modifications"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/manorm-pe.cwl
s:codeRepository: https://github.com/datirium/workflows
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
          s:email: mailto:michael.kotliar@cchmc.org
          s:sameAs:
          - id: http://orcid.org/0000-0002-6486-3898


# doc:
#   $include: ../descriptions/manorm-pe.md


doc: |
  What is MAnorm?
  --------------

  MAnorm is a robust model for quantitative comparison of ChIP-Seq data sets of TFs (transcription factors) or epigenetic modifications and you can use it for:

  * Normalization of two ChIP-seq samples
  * Quantitative comparison (differential analysis) of two ChIP-seq samples
  * Evaluating the overlap enrichment of the protein binding sites(peaks)
  * Elucidating underlying mechanisms of cell-type specific gene regulation

  How MAnorm works?
  ----------------

  MAnorm uses common peaks of two samples as a reference to build the rescaling model for normalization, which is based on the empirical assumption that if a chromatin-associated protein has a substantial number of peaks shared in two conditions, the binding at these common regions will tend to be determined by similar mechanisms, and thus should exhibit similar global binding intensities across samples. The observed differences on common peaks are presumed to reflect the scaling relationship of ChIP-Seq signals between two samples, which can be applied to all peaks.

  What do the inputs mean?
  ----------------

  ### General

  **Experiment short name/Alias**

  * short name for you experiment to identify among the others

  **ChIP-Seq PE sample 1**
  * previously analyzed ChIP-Seq paired-end experiment to be used as Sample 1

  **ChIP-Seq PE sample 2**

  * previously analyzed ChIP-Seq paired-end experiment to be used as Sample 2

  **Genome**

  * Reference genome to be used for gene assigning

  ### Advanced

  **Reads shift size for sample 1**

  * This value is used to shift reads towards 3' direction to determine
    the precise binding site. Set as half of the fragment length. Default 100

  **Reads shift size for sample 2**

  * This value is used to shift reads towards 5' direction to determine
    the precise binding site. Set as half of the fragment length. Default 100

  **M-value (log2-ratio) cutoff**

  * Absolute M-value (log2-ratio) cutoff to define biased (differential binding)
    peaks. Default: 1.0

  **P-value cutoff**

  * P-value cutoff to define biased peaks. Default: 0.01

  **Window size**

  * Window size to count reads and calculate read densities. 2000 is recommended for
    sharp histone marks like H3K4me3 and H3K27ac, and 1000 for TFs or DNase-seq.
    Default: 2000