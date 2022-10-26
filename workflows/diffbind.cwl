cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  first_biological_condition:
    - "chipseq-se.cwl"
    - "chipseq-pe.cwl"
    - "cutandrun-pe.cwl" 
    - "trim-chipseq-se.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "trim-atacseq-pe.cwl"
  second_biological_condition:
    - "chipseq-se.cwl"
    - "chipseq-pe.cwl"
    - "cutandrun-pe.cwl" 
    - "trim-chipseq-se.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "trim-atacseq-pe.cwl"
  blocked_condition:
    - "chipseq-se.cwl"
    - "chipseq-pe.cwl"
    - "cutandrun-pe.cwl" 
    - "trim-chipseq-se.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "trim-atacseq-pe.cwl"
  genome_indices:
    - "genome-indices.cwl"

inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  read_files_cond_1:
    type: File[]
    format: "http://edamontology.org/format_2572"
    label: "Biological condition 1 samples. Minimum 2 samples"
    doc: "Read files for condition 1. Minimim 2 files in BAM format"
    'sd:upstreamSource': "first_biological_condition/bambai_pair"
    'sd:localLabel': true

  read_files_cond_2:
    type: File[]
    format: "http://edamontology.org/format_2572"
    label: "Biological condition 2 samples. Minimum 2 samples"
    doc: "Read files for condition 2. Minimim 2 files in BAM format"
    'sd:upstreamSource': "second_biological_condition/bambai_pair"
    'sd:localLabel': true

  peak_files_cond_1:
    type: File[]
    format: "http://edamontology.org/format_3468"
    label: "Biological condition 1 samples. Minimum 2 samples"
    doc: "XLS peak files for condition 1 from MACS2. Minimim 2 files. Order corresponds to read_files_cond_1"
    'sd:upstreamSource': "first_biological_condition/macs2_called_peaks"

  peak_files_cond_2:
    type: File[]
    format: "http://edamontology.org/format_3468"
    label: "Biological condition 2 samples. Minimum 2 samples"
    doc: "XLS peak files for condition 2 from MACS2. Minimim 2 files. Order corresponds to read_files_cond_2"
    'sd:upstreamSource': "second_biological_condition/macs2_called_peaks"

  genome_coverage_files_cond_1:
    type: File[]
    format: "http://edamontology.org/format_3006"
    label: "Genome coverage(s) for biological condition 1"
    doc: "Genome coverage bigWig file(s) for biological condition 1"
    'sd:upstreamSource': "first_biological_condition/bigwig"

  genome_coverage_files_cond_2:
    type: File[]
    format: "http://edamontology.org/format_3006"
    label: "Genome coverage(s) for biological condition 2"
    doc: "Genome coverage bigWig file(s) for biological condition 2"
    'sd:upstreamSource': "second_biological_condition/bigwig"

  narrow_peaks_files_cond_1:
    type:
      - "null"
      - File[]
    format: "http://edamontology.org/format_3613"
    label: "Called peaks for biological condition 1"
    doc: "Narrow peaks file(s) for biological condition 1"
    'sd:upstreamSource': "first_biological_condition/macs2_narrow_peaks"

  narrow_peaks_files_cond_2:
    type:
      - "null"
      - File[]
    format: "http://edamontology.org/format_3613"
    label: "Called peaks for biological condition 2"
    doc: "Narrow peaks file(s) for biological condition 2"
    'sd:upstreamSource': "second_biological_condition/macs2_narrow_peaks"

  broad_peaks_files_cond_1:
    type:
      - "null"
      - File[]
    format: "http://edamontology.org/format_3614"
    label: "Called peaks for biological condition 1"
    doc: "Broad peaks file(s) for biological condition 1"
    'sd:upstreamSource': "first_biological_condition/macs2_broad_peaks"

  broad_peaks_files_cond_2:
    type:
      - "null"
      - File[]
    format: "http://edamontology.org/format_3614"
    label: "Called peaks for biological condition 2"
    doc: "Broad peaks file(s) for biological condition 2"
    'sd:upstreamSource': "second_biological_condition/macs2_broad_peaks"

  name_cond_1:
    type: string?
    default: "condition_1"
    label: "Condition 1 name, single word with letters and numbers only"
    doc: "Condition 1 name, single word with letters and numbers only"
    'sd:layout':
      advanced: true

  name_cond_2:
    type: string?
    default: "condition_2"
    label: "Condition 2 name, single word with letters and numbers only"
    doc: "Condition 2 name, single word with letters and numbers only"
    'sd:layout':
      advanced: true

  sample_names_cond_1:
    type:
      - "null"
      - string[]
    label: "Biological condition 1 sample names"
    doc: "Aliases for biological condition 1 samples to make the legend for generated plots. Order corresponds to the read_files_cond_1"
    'sd:upstreamSource': "first_biological_condition/alias"

  sample_names_cond_2:
    type:
      - "null"
      - string[]
    label: "Biological condition 2 sample names"
    doc: "Aliases for biological condition 2 samples to make the legend for generated plots. Order corresponds to the read_files_cond_2"
    'sd:upstreamSource': "second_biological_condition/alias"

  blocked_attributes:
    type:
      - "null"
      - string[]
    default: null
    label: "Blocking attributes for multi-factor analysis. Minimum 2"
    doc: |
      Blocking attributes for multi-factor analysis. Minimum 2.
      Either names from --name1 or/and --name2 or array of strings that can be parsed by R to bool.
      In the later case the order and size should correspond to [--read1]+[--read2].
      Default: not applied
    'sd:upstreamSource': "blocked_condition/alias"
    'sd:localLabel': true

  blocked_file:
    type: File?
    label: "Blocking attribute headerless TSV/CSV file for multi-factor analysis with columns to set name and group. If this inputs is set, blocking attributes above are ignored"
    format: "http://edamontology.org/format_2330"
    doc: |
      Blocking attribute metadata file for multi-factor analysis. Headerless TSV/CSV file.
      First column - names from --name1 and --name2, second column - group name. --block is ignored

  annotation_file:
    type: File
    label: "Genome annotation"
    format: "http://edamontology.org/format_3475"
    doc: "Genome annotation file in TSV format"
    'sd:upstreamSource': "genome_indices/annotation"

  chrom_length_file:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Chromosome length file"
    doc: "Chromosome length file"
    'sd:upstreamSource': "genome_indices/chrom_length"   

  fragmentsize:
    type: int?
    default: 125
    label: "Reads extension size, bp"
    doc: "Extended each read from its endpoint along the appropriate strand. Default: 125bp"
    'sd:layout':
      advanced: true

  promoter_dist:
    type: int?
    default: 1000
    label: "Promoter distance, bp"
    doc: "Max distance from gene TSS (in both direction) overlapping which the peak will be assigned to the promoter region. Default: 1000 bp"
    'sd:layout':
      advanced: true

  upstream_dist:
    type: int?
    default: 20000
    label: "Upstream distance, bp"
    doc: "Max distance from the promoter (only in upstream direction) overlapping which the peak will be assigned to the upstream region. Default: 20,000 bp"
    'sd:layout':
      advanced: true

  min_overlap:
    type: int?
    default: 2
    label: "Minimum peakset overlap"
    doc: "Min peakset overlap. Only include peaks in at least this many peaksets when generating consensus peakset. Default: 2"
    'sd:layout':
      advanced: true

  min_read_counts:
    type: int?
    default: 0
    label: "Minimum read counts. Exclude intervals where MAX read counts for all samples < specified value"
    doc: "Min read counts. Exclude all merged intervals where the MAX raw read counts among all of the samples is smaller than the specified value. Default: 0"
    'sd:layout':
      advanced: true

  use_common:
    type: boolean?
    default: false
    label: "Use common peaks within each condition. Ignore Minimum peakset overlap"
    doc: "Derive consensus peaks only from the common peaks within each condition. Min peakset overlap is ignored. Default: false"
    'sd:layout':
      advanced: true

  remove_duplicates:
    type: boolean?
    default: false
    label: "Remove duplicated reads"
    doc: "Remove reads that map to exactly the same genomic position. Default: false"
    'sd:layout':
      advanced: true

  cutoff_value:
    type: float?
    default: 0.05
    label: "P-value or FDR cutoff for reported results"
    doc: "P-value or FDR cutoff for reported results"
    'sd:layout':
      advanced: true

  cutoff_param:
    type:
      - "null"
      - type: enum
        name: "cutoff"
        symbols: ["pvalue", "fdr"]
    default: "pvalue"
    label: "Parameter to which cutoff should be applied"
    doc: "Parameter to which cutoff should be applied (fdr or pvalue). Default: fdr"
    'sd:layout':
      advanced: true

  analysis_method:
    type:
      - "null"
      - type: enum
        name: "method"
        symbols: ["deseq2", "edger"]
    default: "edger"
    label: "Analysis method"
    doc: "Method by which to analyze differential binding affinity. Default: deseq2"
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 4
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  genome_coverage_cond_1:
    type: File[]
    format: "http://edamontology.org/format_3006"
    label: "Genome coverage(s) for biological condition 1"
    doc: "Genome coverage bigWig file(s) for biological condition 1"
    outputSource: pipe/coverage_files_cond_1
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "Genome coverage 1"
        height: 120

  genome_coverage_cond_2:
    type: File[]
    format: "http://edamontology.org/format_3006"
    label: "Genome coverage(s) for biological condition 2"
    doc: "Genome coverage bigWig file(s) for biological condition 2"
    outputSource: pipe/coverage_files_cond_2
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'wig'
        name: "Genome coverage 2"
        height: 120

  narrow_peaks_cond_1:
    type:
      - "null"
      - File[]
    format: "http://edamontology.org/format_3613"
    label: "Called peaks for biological condition 1"
    doc: "Narrow peaks file(s) for biological condition 1"
    outputSource: pipe/n_peaks_files_cond_1
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Called peaks 1"
        displayMode: "COLLAPSE"
        height: 40

  narrow_peaks_cond_2:
    type:
      - "null"
      - File[]
    format: "http://edamontology.org/format_3613"
    label: "Called peaks for biological condition 2"
    doc: "Narrow peaks file(s) for biological condition 2"
    outputSource: pipe/n_peaks_files_cond_2
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Called peaks 2"
        displayMode: "COLLAPSE"
        height: 40

  broad_peaks_cond_1:
    type:
      - "null"
      - File[]
    format: "http://edamontology.org/format_3614"
    label: "Called peaks for biological condition 1"
    doc: "Broad peaks file(s) for biological condition 1"
    outputSource: pipe/b_peaks_files_cond_1
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Called peaks 1"
        displayMode: "COLLAPSE"
        height: 40

  broad_peaks_cond_2:
    type:
      - "null"
      - File[]
    format: "http://edamontology.org/format_3614"
    label: "Called peaks for biological condition 2"
    doc: "Broad peaks file(s) for biological condition 2"
    outputSource: pipe/b_peaks_files_cond_2
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        name: "Called peaks 2"
        displayMode: "COLLAPSE"
        height: 40

  diffbind_report_file:
    type: File
    format: "http://edamontology.org/format_3475"
    label: "Differential binding analysis results"
    doc: "Differential binding analysis results exported as TSV"
    outputSource: restore_columns/output_file
    'sd:visualPlugins':
      - syncfusiongrid:
          tab: 'Differential Peak Calling'
          Title: 'Differential Binding Analysis Results'

  diffbind_bed_file:
    type: File
    format: "http://edamontology.org/format_3004"
    label: "Estimated differential peaks"
    doc: "Estimated differential peaks, bigBed"
    outputSource: bed_to_bigbed/bigbed_file
    'sd:visualPlugins':
    - igvbrowser:
        tab: 'IGV Genome Browser'
        id: 'igvbrowser'
        type: 'annotation'
        format: 'bigbed'
        name: "Differential peaks"
        height: 40

  diffbind_peak_correlation_heatmap:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "Peak overlap correlation heatmap"
    doc: "Peak overlap correlation heatmap"
    outputSource: diffbind/peak_overlap_corr_heatmap
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'Peak overlap correlation heatmap'

  diffbind_peak_correlation_heatmap_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "Peak overlap correlation heatmap"
    doc: "Peak overlap correlation heatmap"
    outputSource: diffbind/peak_overlap_corr_heatmap_pdf

  diffbind_counts_correlation_heatmap:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "Raw counts correlation heatmap"
    doc: "Raw counts correlation heatmap"
    outputSource: diffbind/raw_counts_corr_heatmap
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'Raw counts correlation heatmap'

  diffbind_counts_correlation_heatmap_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "Raw counts correlation heatmap"
    doc: "Raw counts correlation heatmap"
    outputSource: diffbind/raw_counts_corr_heatmap_pdf

  diffbind_consensus_peak_venn_diagram:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "Consensus peak Venn Diagram"
    doc: "Consensus peak Venn Diagram"
    outputSource: diffbind/consensus_peak_venn_diagram
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'Consensus peak Venn Diagram'

  diffbind_consensus_peak_venn_diagram_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "Consensus peak Venn Diagram"
    doc: "Consensus peak Venn Diagram"
    outputSource: diffbind/consensus_peak_venn_diagram_pdf

  diffbind_all_peak_overlap_rate_plot:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "All peak overlap rate plot"
    doc: "All peak overlap rate plot"
    outputSource: diffbind/all_peak_overlap_rate_plot
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'All peak overlap rate plot'  

  diffbind_all_peak_overlap_rate_plot_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "All peak overlap rate plot"
    doc: "All peak overlap rate plot"
    outputSource: diffbind/all_peak_overlap_rate_plot_pdf

  diffbind_peak_overlap_rate_plot_cond_1:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "Condition 1 peak overlap rate plot"
    doc: "Condition 1 peak overlap rate plot"
    outputSource: diffbind/peak_overlap_rate_plot_cond_1
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'Condition 1 peak overlap rate plot'

  diffbind_peak_overlap_rate_plot_cond_1_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "Condition 1 peak overlap rate plot"
    doc: "Condition 1 peak overlap rate plot"
    outputSource: diffbind/peak_overlap_rate_plot_cond_1_pdf

  diffbind_peak_overlap_rate_plot_cond_2:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "Condition 2 peak overlap rate plot"
    doc: "Condition 2 peak overlap rate plot"
    outputSource: diffbind/peak_overlap_rate_plot_cond_2
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'Condition 2 peak overlap rate plot'

  diffbind_peak_overlap_rate_plot_cond_2_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "Condition 2 peak overlap rate plot"
    doc: "Condition 2 peak overlap rate plot"
    outputSource: diffbind/peak_overlap_rate_plot_cond_2_pdf

  diffbind_all_data_correlation_heatmap:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "Not filtered normalized counts correlation heatmap"
    doc: "Not filtered normalized counts correlation heatmap"
    outputSource: select_files/selected_all_norm_counts_corr_heatmap
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'Not filtered normalized counts correlation heatmap'

  diffbind_all_data_correlation_heatmap_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "Not filtered normalized counts correlation heatmap"
    doc: "Not filtered normalized counts correlation heatmap"
    outputSource: select_files/selected_all_norm_counts_corr_heatmap_pdf

  diffbind_db_sites_correlation_heatmap:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "Normalized counts correlation heatmap for significantly differentially bound sites"
    doc: "Normalized counts correlation heatmap for significantly differentially bound sites"
    outputSource: select_files/selected_diff_filtered_norm_counts_corr_heatmap
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'Normalized counts correlation heatmap for significantly differentially bound sites'

  diffbind_db_sites_correlation_heatmap_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "Normalized counts correlation heatmap for significantly differentially bound sites"
    doc: "Normalized counts correlation heatmap for significantly differentially bound sites"
    outputSource: select_files/selected_diff_filtered_norm_counts_corr_heatmap_pdf

  diffbind_db_sites_binding_heatmap:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "Binding heatmap for significantly differentially bound sites"
    doc: "Binding heatmap for significantly differentially bound sites"
    outputSource: select_files/selected_binding_heatmap
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'Binding heatmap for significantly differentially bound sites'

  diffbind_db_sites_binding_heatmap_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "Binding heatmap for significantly differentially bound sites"
    doc: "Binding heatmap for significantly differentially bound sites"
    outputSource: select_files/selected_binding_heatmap_pdf

  diffbind_pca_plot:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "PCA plot for significantly differentially bound sites"
    doc: "PCA plot for significantly differentially bound sites"
    outputSource: select_files/selected_diff_filtered_pca_plot
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'PCA plot for significantly differentially bound sites'

  diffbind_pca_plot_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "PCA plot for significantly differentially bound sites"
    doc: "PCA plot for significantly differentially bound sites"
    outputSource: select_files/selected_diff_filtered_pca_plot_pdf

  diffbind_all_pca_plot:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "PCA plot for all bound sites"
    doc: "PCA plot for all bound sites"
    outputSource: select_files/selected_all_pca_plot
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'PCA plot for all bound sites'

  diffbind_all_pca_plot_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "PCA plot for all bound sites"
    doc: "PCA plot for all bound sites"
    outputSource: select_files/selected_all_pca_plot_pdf

  diffbind_ma_plot:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "MA plot for significantly differentially bound sites"
    doc: "MA plot for significantly differentially bound sites"
    outputSource: select_files/selected_ma_plot   
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'MA plot for significantly differentially bound sites'

  diffbind_ma_plot_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "MA plot for significantly differentially bound sites"
    doc: "MA plot for significantly differentially bound sites"
    outputSource: select_files/selected_ma_plot_pdf

  diffbind_volcano_plot:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "Volcano plot for significantly differentially bound sites"
    doc: "Volcano plot for significantly differentially bound sites"
    outputSource: select_files/selected_volcano_plot
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'Volcano plot for for significantly differentially bound sites'

  diffbind_volcano_plot_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "Volcano plot for significantly differentially bound sites"
    doc: "Volcano plot for significantly differentially bound sites"
    outputSource: select_files/selected_volcano_plot_pdf

  diffbind_boxplot_plot:
    type: File?
    format: "http://edamontology.org/format_3603"
    label: "Box plots of read distributions for significantly differentially bound sites"
    doc: "Box plots of read distributions for significantly differentially bound sites"
    outputSource: select_files/selected_boxplot  
    'sd:visualPlugins':
      - image:
          tab: 'Plots'
          Caption: 'Box plots of read distributions for significantly differentially bound sites'

  diffbind_boxplot_plot_pdf:
    type: File?
    format: "http://edamontology.org/format_3508"
    label: "Box plots of read distributions for significantly differentially bound sites"
    doc: "Box plots of read distributions for significantly differentially bound sites"
    outputSource: select_files/selected_boxplot_pdf

  diffbind_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "diffbind stdout log"
    doc: "diffbind stdout log"
    outputSource: diffbind/stdout_log

  diffbind_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "diffbind stderr log"
    doc: "diffbind stderr log"
    outputSource: diffbind/stderr_log


steps:

  pipe:
    run:
      cwlVersion: v1.0
      class: ExpressionTool
      inputs:
        genome_coverage_files_cond_1:
          type: File[]
        genome_coverage_files_cond_2:
          type: File[]
        narrow_peaks_files_cond_1:
          type:
            - "null"
            - File[]
        narrow_peaks_files_cond_2:
          type:
            - "null"
            - File[]
        broad_peaks_files_cond_1:
          type:
            - "null"
            - File[]
        broad_peaks_files_cond_2:
          type:
            - "null"
            - File[]
      outputs:
        coverage_files_cond_1:
          type: File[]
        coverage_files_cond_2:
          type: File[]
        n_peaks_files_cond_1:
          type:
            - "null"
            - File[]
        n_peaks_files_cond_2:
          type:
            - "null"
            - File[]
        b_peaks_files_cond_1:
          type:
            - "null"
            - File[]
        b_peaks_files_cond_2:
          type:
            - "null"
            - File[]
      expression: |
        ${
          return {
            "coverage_files_cond_1": inputs.genome_coverage_files_cond_1,
            "coverage_files_cond_2": inputs.genome_coverage_files_cond_2,
            "n_peaks_files_cond_1": inputs.narrow_peaks_files_cond_1,
            "n_peaks_files_cond_2": inputs.narrow_peaks_files_cond_2,
            "b_peaks_files_cond_1": inputs.broad_peaks_files_cond_1,
            "b_peaks_files_cond_2": inputs.broad_peaks_files_cond_2
          };
        }
    in:
      genome_coverage_files_cond_1: genome_coverage_files_cond_1
      genome_coverage_files_cond_2: genome_coverage_files_cond_2
      narrow_peaks_files_cond_1: narrow_peaks_files_cond_1
      narrow_peaks_files_cond_2: narrow_peaks_files_cond_2
      broad_peaks_files_cond_1: broad_peaks_files_cond_1
      broad_peaks_files_cond_2: broad_peaks_files_cond_2
    out:
      - coverage_files_cond_1
      - coverage_files_cond_2
      - n_peaks_files_cond_1
      - n_peaks_files_cond_2
      - b_peaks_files_cond_1
      - b_peaks_files_cond_2

  diffbind:
    run: ../tools/diffbind.cwl
    in:
      read_files_cond_1: read_files_cond_1
      read_files_cond_2: read_files_cond_2
      peak_files_cond_1: peak_files_cond_1
      peak_files_cond_2: peak_files_cond_2
      name_cond_1: name_cond_1
      name_cond_2: name_cond_2
      sample_names_cond_1: sample_names_cond_1
      sample_names_cond_2: sample_names_cond_2
      cutoff_value: cutoff_value
      cutoff_param: cutoff_param
      fragmentsize: fragmentsize
      remove_duplicates: remove_duplicates
      analysis_method: analysis_method
      blocked_attributes: blocked_attributes
      blocked_file: blocked_file
      min_overlap: min_overlap
      min_read_counts: min_read_counts
      use_common: use_common
      threads: threads
      peakformat:
        default: "macs"
    out:
      - diff_filtered_report_deseq
      - diff_filtered_report_deseq_blocked
      - diff_filtered_report_edger
      - diff_filtered_report_edger_blocked
      - boxplot_deseq
      - boxplot_deseq_blocked
      - boxplot_edger
      - boxplot_edger_blocked
      - boxplot_deseq_pdf
      - boxplot_deseq_blocked_pdf
      - boxplot_edger_pdf
      - boxplot_edger_blocked_pdf
      - volcano_plot_deseq
      - volcano_plot_deseq_blocked
      - volcano_plot_edger
      - volcano_plot_edger_blocked
      - volcano_plot_deseq_pdf
      - volcano_plot_deseq_blocked_pdf
      - volcano_plot_edger_pdf
      - volcano_plot_edger_blocked_pdf
      - ma_plot_deseq
      - ma_plot_deseq_blocked
      - ma_plot_edger
      - ma_plot_edger_blocked
      - ma_plot_deseq_pdf
      - ma_plot_deseq_blocked_pdf
      - ma_plot_edger_pdf
      - ma_plot_edger_blocked_pdf
      - diff_filtered_pca_plot_deseq
      - diff_filtered_pca_plot_deseq_blocked
      - diff_filtered_pca_plot_edger
      - diff_filtered_pca_plot_edger_blocked
      - diff_filtered_pca_plot_deseq_pdf
      - diff_filtered_pca_plot_deseq_blocked_pdf
      - diff_filtered_pca_plot_edger_pdf
      - diff_filtered_pca_plot_edger_blocked_pdf
      - all_pca_plot_deseq
      - all_pca_plot_deseq_blocked
      - all_pca_plot_edger
      - all_pca_plot_edger_blocked
      - all_pca_plot_deseq_pdf
      - all_pca_plot_deseq_blocked_pdf
      - all_pca_plot_edger_pdf
      - all_pca_plot_edger_blocked_pdf
      - binding_heatmap_deseq
      - binding_heatmap_deseq_blocked
      - binding_heatmap_edger
      - binding_heatmap_edger_blocked
      - binding_heatmap_deseq_pdf
      - binding_heatmap_deseq_blocked_pdf
      - binding_heatmap_edger_pdf
      - binding_heatmap_edger_blocked_pdf
      - diff_filtered_norm_counts_corr_heatmap_deseq
      - diff_filtered_norm_counts_corr_heatmap_deseq_blocked
      - diff_filtered_norm_counts_corr_heatmap_edger
      - diff_filtered_norm_counts_corr_heatmap_edger_blocked
      - diff_filtered_norm_counts_corr_heatmap_deseq_pdf
      - diff_filtered_norm_counts_corr_heatmap_deseq_blocked_pdf
      - diff_filtered_norm_counts_corr_heatmap_edger_pdf
      - diff_filtered_norm_counts_corr_heatmap_edger_blocked_pdf
      - all_norm_counts_corr_heatmap_deseq
      - all_norm_counts_corr_heatmap_deseq_blocked
      - all_norm_counts_corr_heatmap_edger
      - all_norm_counts_corr_heatmap_edger_blocked
      - all_norm_counts_corr_heatmap_deseq_pdf
      - all_norm_counts_corr_heatmap_deseq_blocked_pdf
      - all_norm_counts_corr_heatmap_edger_pdf
      - all_norm_counts_corr_heatmap_edger_blocked_pdf
      - consensus_peak_venn_diagram
      - raw_counts_corr_heatmap
      - peak_overlap_corr_heatmap
      - all_peak_overlap_rate_plot
      - peak_overlap_rate_plot_cond_1
      - peak_overlap_rate_plot_cond_2
      - consensus_peak_venn_diagram_pdf
      - raw_counts_corr_heatmap_pdf
      - peak_overlap_corr_heatmap_pdf
      - all_peak_overlap_rate_plot_pdf
      - peak_overlap_rate_plot_cond_1_pdf
      - peak_overlap_rate_plot_cond_2_pdf
      - stdout_log
      - stderr_log

  filter_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file: 
        source: [analysis_method, blocked_attributes, diffbind/diff_filtered_report_deseq, diffbind/diff_filtered_report_deseq_blocked, diffbind/diff_filtered_report_edger, diffbind/diff_filtered_report_edger_blocked]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      script:
        default: >
          cat $0 | grep -v "Start" | awk
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
      input_file:
        source: [analysis_method, blocked_attributes, diffbind/diff_filtered_report_deseq, diffbind/diff_filtered_report_deseq_blocked, diffbind/diff_filtered_report_edger, diffbind/diff_filtered_report_edger_blocked, assign_genes/result_file]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return [self[6], self[2]];
                } else {
                  return [self[6], self[3]];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return [self[6], self[4]];
                } else {
                  return [self[6], self[5]];
                }
              }
          }
      script:
        default: |
          cat $0 | grep -v "start" | sort -k 11n | cut -f 1-5,15 > iaintersect_result.tsv
          cat $1 | grep -v "Start" > diffbind_result.tsv
          HEADER=`head -n 1 $1`;
          echo -e "Refseq_id\tGene_id\ttxStart\ttxEnd\tStrand\tRegion\t${HEADER}" > `basename $0`;
          cat iaintersect_result.tsv | paste - diffbind_result.tsv >> `basename $0`
          rm iaintersect_result.tsv diffbind_result.tsv
    out: [output_file]

  convert_to_bed:
    run: ../tools/custom-bash.cwl
    in:
      input_file: restore_columns/output_file
      script:
        default: |
          cat "$0" | awk -F "\t" 'NR==1 {for (i=1; i<=NF; i++) {ix[$i]=i} } NR>1 {color="255,0,0"; if ($ix["Fold"]<0) color="0,255,0"; print $ix["Chr"]"\t"$ix["Start"]"\t"$ix["End"]"\tPv="$ix["p-value"]";FDR="$ix["FDR"]"\t"1000"\t"$ix["Strand"]"\t"$ix["Start"]"\t"$ix["End"]"\t"color}' > `basename $0`
    out: [output_file]

  sort_bed:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: convert_to_bed/output_file
      key:
        default: ["1,1","2,2n"]
    out: [sorted_file]

  bed_to_bigbed:
    run: ../tools/ucsc-bedtobigbed.cwl
    in:
      input_bed: sort_bed/sorted_file
      bed_type:
        default: "bed4+5"
      chrom_length_file: chrom_length_file
      output_filename:
        source: sort_bed/sorted_file
        valueFrom: $(self.basename.split('.').slice(0,-1).join('.') + ".bigBed")
    out: [bigbed_file]

  select_files:
    in:
      input_boxplot:
        source: [analysis_method, blocked_attributes, diffbind/boxplot_deseq, diffbind/boxplot_deseq_blocked, diffbind/boxplot_edger, diffbind/boxplot_edger_blocked]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_boxplot_pdf:
        source: [analysis_method, blocked_attributes, diffbind/boxplot_deseq_pdf, diffbind/boxplot_deseq_blocked_pdf, diffbind/boxplot_edger_pdf, diffbind/boxplot_edger_blocked_pdf]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_volcano_plot:
        source: [analysis_method, blocked_attributes, diffbind/volcano_plot_deseq, diffbind/volcano_plot_deseq_blocked, diffbind/volcano_plot_edger, diffbind/volcano_plot_edger_blocked]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_volcano_plot_pdf:
        source: [analysis_method, blocked_attributes, diffbind/volcano_plot_deseq_pdf, diffbind/volcano_plot_deseq_blocked_pdf, diffbind/volcano_plot_edger_pdf, diffbind/volcano_plot_edger_blocked_pdf]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_ma_plot:
        source: [analysis_method, blocked_attributes, diffbind/ma_plot_deseq, diffbind/ma_plot_deseq_blocked, diffbind/ma_plot_edger, diffbind/ma_plot_edger_blocked]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_ma_plot_pdf:
        source: [analysis_method, blocked_attributes, diffbind/ma_plot_deseq_pdf, diffbind/ma_plot_deseq_blocked_pdf, diffbind/ma_plot_edger_pdf, diffbind/ma_plot_edger_blocked_pdf]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_diff_filtered_pca_plot:
        source: [analysis_method, blocked_attributes, diffbind/diff_filtered_pca_plot_deseq, diffbind/diff_filtered_pca_plot_deseq_blocked, diffbind/diff_filtered_pca_plot_edger, diffbind/diff_filtered_pca_plot_edger_blocked]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_diff_filtered_pca_plot_pdf:
        source: [analysis_method, blocked_attributes, diffbind/diff_filtered_pca_plot_deseq_pdf, diffbind/diff_filtered_pca_plot_deseq_blocked_pdf, diffbind/diff_filtered_pca_plot_edger_pdf, diffbind/diff_filtered_pca_plot_edger_blocked_pdf]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_all_pca_plot:
        source: [analysis_method, blocked_attributes, diffbind/all_pca_plot_deseq, diffbind/all_pca_plot_deseq_blocked, diffbind/all_pca_plot_edger, diffbind/all_pca_plot_edger_blocked]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_all_pca_plot_pdf:
        source: [analysis_method, blocked_attributes, diffbind/all_pca_plot_deseq_pdf, diffbind/all_pca_plot_deseq_blocked_pdf, diffbind/all_pca_plot_edger_pdf, diffbind/all_pca_plot_edger_blocked_pdf]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_binding_heatmap:
        source: [analysis_method, blocked_attributes, diffbind/binding_heatmap_deseq, diffbind/binding_heatmap_deseq_blocked, diffbind/binding_heatmap_edger, diffbind/binding_heatmap_edger_blocked]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_binding_heatmap_pdf:
        source: [analysis_method, blocked_attributes, diffbind/binding_heatmap_deseq_pdf, diffbind/binding_heatmap_deseq_blocked_pdf, diffbind/binding_heatmap_edger_pdf, diffbind/binding_heatmap_edger_blocked_pdf]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_diff_filtered_norm_counts_corr_heatmap:
        source: [analysis_method, blocked_attributes, diffbind/diff_filtered_norm_counts_corr_heatmap_deseq, diffbind/diff_filtered_norm_counts_corr_heatmap_deseq_blocked, diffbind/diff_filtered_norm_counts_corr_heatmap_edger, diffbind/diff_filtered_norm_counts_corr_heatmap_edger_blocked]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_diff_filtered_norm_counts_corr_heatmap_pdf:
        source: [analysis_method, blocked_attributes, diffbind/diff_filtered_norm_counts_corr_heatmap_deseq_pdf, diffbind/diff_filtered_norm_counts_corr_heatmap_deseq_blocked_pdf, diffbind/diff_filtered_norm_counts_corr_heatmap_edger_pdf, diffbind/diff_filtered_norm_counts_corr_heatmap_edger_blocked_pdf]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_all_norm_counts_corr_heatmap:
        source: [analysis_method, blocked_attributes, diffbind/all_norm_counts_corr_heatmap_deseq, diffbind/all_norm_counts_corr_heatmap_deseq_blocked, diffbind/all_norm_counts_corr_heatmap_edger, diffbind/all_norm_counts_corr_heatmap_edger_blocked]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
      input_all_norm_counts_corr_heatmap_pdf:
        source: [analysis_method, blocked_attributes, diffbind/all_norm_counts_corr_heatmap_deseq_pdf, diffbind/all_norm_counts_corr_heatmap_deseq_blocked_pdf, diffbind/all_norm_counts_corr_heatmap_edger_pdf, diffbind/all_norm_counts_corr_heatmap_edger_blocked_pdf]
        valueFrom: |
          ${
              if (self[0] == "deseq2") {
                if (self[1] == null || self[1].length == 0){
                  return self[2];
                } else {
                  return self[3];
                }
              } else {
                if (self[1] == null || self[1].length == 0){
                  return self[4];
                } else {
                  return self[5];
                }
              }
          }
    out:
      - selected_boxplot
      - selected_volcano_plot
      - selected_ma_plot
      - selected_diff_filtered_pca_plot
      - selected_all_pca_plot
      - selected_binding_heatmap
      - selected_diff_filtered_norm_counts_corr_heatmap
      - selected_all_norm_counts_corr_heatmap
      - selected_boxplot_pdf
      - selected_volcano_plot_pdf
      - selected_ma_plot_pdf
      - selected_diff_filtered_pca_plot_pdf
      - selected_all_pca_plot_pdf
      - selected_binding_heatmap_pdf
      - selected_diff_filtered_norm_counts_corr_heatmap_pdf
      - selected_all_norm_counts_corr_heatmap_pdf
    run:
      cwlVersion: v1.0
      class: ExpressionTool
      requirements:
        - class: InlineJavascriptRequirement
      inputs:
        input_boxplot:
          type: File?
        input_boxplot_pdf:
          type: File?
        input_volcano_plot:
          type: File?
        input_volcano_plot_pdf:
          type: File?
        input_ma_plot:
          type: File?
        input_ma_plot_pdf:
          type: File?
        input_diff_filtered_pca_plot:
          type: File?
        input_diff_filtered_pca_plot_pdf:
          type: File?        
        input_all_pca_plot:
          type: File?
        input_all_pca_plot_pdf:
          type: File?        
        input_binding_heatmap:
          type: File?
        input_binding_heatmap_pdf:
          type: File?          
        input_diff_filtered_norm_counts_corr_heatmap:
          type: File?
        input_diff_filtered_norm_counts_corr_heatmap_pdf:
          type: File?
        input_all_norm_counts_corr_heatmap:
          type: File?
        input_all_norm_counts_corr_heatmap_pdf:
          type: File?
      outputs:
        selected_boxplot:
          type: File?
        selected_volcano_plot:
          type: File?
        selected_ma_plot:
          type: File?
        selected_diff_filtered_pca_plot:
          type: File?
        selected_all_pca_plot:
          type: File?
        selected_binding_heatmap:
          type: File?
        selected_diff_filtered_norm_counts_corr_heatmap:
          type: File?
        selected_all_norm_counts_corr_heatmap:
          type: File?
        selected_boxplot_pdf:
          type: File?
        selected_volcano_plot_pdf:
          type: File?
        selected_ma_plot_pdf:
          type: File?
        selected_diff_filtered_pca_plot_pdf:
          type: File?
        selected_all_pca_plot_pdf:
          type: File?
        selected_binding_heatmap_pdf:
          type: File?
        selected_diff_filtered_norm_counts_corr_heatmap_pdf:
          type: File?
        selected_all_norm_counts_corr_heatmap_pdf:
          type: File?
      expression: |
        ${
          return {
            "selected_boxplot": inputs.input_boxplot,
            "selected_volcano_plot": inputs.input_volcano_plot,
            "selected_ma_plot": inputs.input_ma_plot,
            "selected_diff_filtered_pca_plot": inputs.input_diff_filtered_pca_plot,
            "selected_all_pca_plot": inputs.input_all_pca_plot,
            "selected_binding_heatmap": inputs.input_binding_heatmap,
            "selected_diff_filtered_norm_counts_corr_heatmap": inputs.input_diff_filtered_norm_counts_corr_heatmap,
            "selected_all_norm_counts_corr_heatmap": inputs.input_all_norm_counts_corr_heatmap,
            "selected_boxplot_pdf": inputs.input_boxplot_pdf,
            "selected_volcano_plot_pdf": inputs.input_volcano_plot_pdf,
            "selected_ma_plot_pdf": inputs.input_ma_plot_pdf,
            "selected_diff_filtered_pca_plot_pdf": inputs.input_diff_filtered_pca_plot_pdf,
            "selected_all_pca_plot_pdf": inputs.input_all_pca_plot_pdf,
            "selected_binding_heatmap_pdf": inputs.input_binding_heatmap_pdf,
            "selected_diff_filtered_norm_counts_corr_heatmap_pdf": inputs.input_diff_filtered_norm_counts_corr_heatmap_pdf,
            "selected_all_norm_counts_corr_heatmap_pdf": inputs.input_all_norm_counts_corr_heatmap_pdf
          };
        }

      
$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "DiffBind - Differential Binding Analysis of ChIP-Seq or CUTß&RUN/Tag Peak Data"
s:name: "DiffBind - Differential Binding Analysis of ChIP-Seq or CUT&RUN/Tag Peak Data"
s:alternateName: "Compute differentially bound sites from multiple ChIP-seq or CUT&RUN/Tag experiments using affinity (quantitative) and occupancy data."

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/diffbind.cwl
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
#   $include: ../descriptions/diffbind.md


doc: |
  Differential Binding Analysis of ChIP-Seq or CUT&RUN/Tag Peak Data
  ---------------------------------------------------
  
  DiffBind processes ChIP-Seq or CUT&RUN/Tag data enriched for genomic loci where specific
  protein/DNA binding occurs, including peak sets identified by peak caller tools and aligned
  sequence read datasets. It is designed to work with multiple peak sets simultaneously,
  representing different ChIP or CUT&RUN/Tag experiments (antibodies, transcription factor
  and/or histone marks, experimental conditions, replicates) as well as managing the results
  of multiple peak callers.

  For more information please refer to:
  -------------------------------------
  Ross-Innes CS, Stark R, Teschendorff AE, Holmes KA, Ali HR, Dunning MJ, Brown GD, Gojis O,
  Ellis IO, Green AR, Ali S, Chin S, Palmieri C, Caldas C, Carroll JS (2012). “Differential
  oestrogen receptor binding is associated with clinical outcome in breast cancer.” Nature,
  481, -4.