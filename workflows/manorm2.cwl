cwlVersion: v1.0
class: Workflow


requirements:
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: MultipleInputFeatureRequirement
  - class: InlineJavascriptRequirement
    expressionLib:
    - var split_features = function(line) {
          function get_unique(value, index, self) {
            return self.indexOf(value) === index && value != "";
          }
          let splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
          return (splitted_line && !!splitted_line.length)?splitted_line:null;
      };


"sd:upstream":
  control:
  - "cutandrun-macs2-pe.cwl"
  - "cutandrun-seacr-pe.cwl"
  - "trim-chipseq-se.cwl"
  - "trim-chipseq-pe.cwl"
  - "trim-atacseq-se.cwl"
  - "trim-atacseq-pe.cwl"
  treatment:
  - "cutandrun-macs2-pe.cwl"
  - "cutandrun-seacr-pe.cwl"
  - "trim-chipseq-se.cwl"
  - "trim-chipseq-pe.cwl"
  - "trim-atacseq-se.cwl"
  - "trim-atacseq-pe.cwl"
  genome_indices:
  - "genome-indices.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name"
    sd:preview:
      position: 1

  read_files_cond_1:
    type: File[]
    label: "Control group sample(s)"
    doc: |
      Coordinate sorted and indexed BAM files with
      aligned reads from the samples that belong
      to the control group.
    "sd:upstreamSource": "control/bambai_pair"
    "sd:localLabel": true

  read_files_cond_2:
    type: File[]
    label: "Treatment group sample(s)"
    doc: |
      Coordinate sorted and indexed BAM files with
      aligned reads from the samples that belong
      to the treatment group.
    "sd:upstreamSource": "treatment/bambai_pair"
    "sd:localLabel": true

  narrow_peak_files_cond_1:
    type:
    - "null"
    - File[]
    label: "Control group sample(s)"
    doc: |
      Narrow peak files with the peaks called from
      the samples that belong to the control group.
    "sd:upstreamSource": "control/macs2_narrow_peaks"

  narrow_peak_files_cond_2:
    type:
    - "null"
    - File[]
    label: "Treatment group sample(s)"
    doc: |
      Narrow peak files with the peaks called from
      the samples that belong to the treatment group.
    "sd:upstreamSource": "treatment/macs2_narrow_peaks"

  broad_peak_files_cond_1:
    type:
    - "null"
    - File[]
    label: "Control group sample(s)"
    doc: |
      Broad peak files with the peaks called from
      the samples that belong to the control group.
    "sd:upstreamSource": "control/macs2_broad_peaks"

  broad_peak_files_cond_2:
    type:
    - "null"
    - File[]
    label: "Treatment group sample(s)"
    doc: |
      Broad peak files with the peaks called from
      the samples that belong to the treatment group.
    "sd:upstreamSource": "treatment/macs2_broad_peaks"

  genome_cov_files_cond_1:
    type: File[]
    label: "Control group sample(s)"
    doc: |
      Genome coverage files with the normalized
      number of fragments at each base from the
      samples that belong to the control group.
    "sd:upstreamSource": "control/bigwig"

  genome_cov_files_cond_2:
    type: File[]
    label: "Treatment group sample(s)"
    doc: |
      Genome coverage files with the normalized
      number of fragments at each base from the
      samples that belong to the treatment group.
    "sd:upstreamSource": "treatment/bigwig"

  summit_files_cond_1:
    type:
    - "null"
    - File[]
    label: "Control group sample(s)"
    doc: |
      BED files with the summits of the peaks
      called from the samples that belong to
      the control group. If not provided, the
      peak center is taken as the summit.
    "sd:upstreamSource": "control/macs2_peak_summits"

  summit_files_cond_2:
    type:
    - "null"
    - File[]
    label: "Treatment group sample(s)"
    doc: |
      BED files with the summits of the peaks
      called from the samples that belong to
      the treatment group. If not provided,
      the peak center is taken as the summit.
    "sd:upstreamSource": "treatment/macs2_peak_summits"

  sample_names_cond_1:
    type: string[]
    label: "Control group sample(s)"
    doc: |
      Names of the samples that belong to
      the control group.
    "sd:upstreamSource": "control/alias"

  sample_names_cond_2:
    type: string[]
    label: "Treatment group sample(s)"
    doc: |
      Names of the samples that belong to
      the treatment group.
    "sd:upstreamSource": "treatment/alias"

  condition_1:
    type: string?
    default: "control"
    label: "Control group sample(s)"
    doc: |
      Alternative name for the
      control group.

  condition_2:
    type: string?
    default: "treatment"
    label: "Treatment group sample(s)"
    doc: |
      Alternative name for the
      treatment group.

  paired_end:
    type: boolean?
    default: false
    label: "Consider all reads as paired-end"
    doc: |
      Consider all reads as paired-end. When
      counting reads in the reference genomic
      bins, use the middle point of each DNA
      fragment. When running in the paired-end
      mode the reads positions are not shifted.
      Default: treat all reads as single-end.

  annotation_file:
    type: File
    label: "Genome type"
    doc: |
      Genome annotation file for
      the nearest genes assignment.
    "sd:upstreamSource": "genome_indices/annotation"

  chrom_length_file:
    type: File
    label: "Genome type"
    doc: | 
      Chromosome length file for
      generating bigBed tracks.
    "sd:upstreamSource": "genome_indices/chrom_length"

  maximum_peak_number:
    type: int?
    default: 0
    label: "Maximum number of peaks (optional)"
    doc: |
      The maximum number of the most significant
      peaks to select from each peak file when
      constructing reference genomic bins. The
      top significant peaks are selected based
      on the score column which is calculated
      by MACS2 either as int(-10*log10qvalue) or
      as int(-10*log10qvalue), depending on the
      cutoff used for peak calling.
      Default: 0 - keep all peaks.

  minimum_peak_gap:
    type: int?
    default: 150
    label: "Minimum peak gap"
    doc: |
      Peaks remained after optional filtering by
      the maximum number of the most significant
      peaks will be merged if the distance between
      them is smaller than the provided value.
      Merging is first applied per sample and then
      to all peaks together before splitting them
      into the reference genomic bins of a selected
      size. Default: 150

  normalization_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "pseudo-reference"
      - "baseline"
      - "hierarchical"
    default: "pseudo-reference"
    label: "Normalization method"
    doc: |
      Normalization method applied to the
      raw read counts. pseudo-reference -
      normalize each sample to the pseudo
      reference that includes the average
      intesities from all samples. A
      reference genomic bin is occupied by
      the pseudo reference if it was occupied
      by at least one sample that the reference
      was constructed from. Each sample is
      MA-normalized to the pseudo reference
      using the common genomic bins between
      the reference and a sample. baseline -
      normalize each sample to the one whose
      log2 size factor is closest to 0.
      hierarchical - similar to the baseline
      but first all samples are normalized
      within their groups, and than two groups
      are normalized between each other.
      Default: pseudo-reference

  maximum_padj:
    type: float?
    default: 0.05
    label: "Maximum adjusted P-value"
    doc: |
      Filtering threshold to report only
      differentially bound sites with the
      adjusted P-value less than or equal
      to the provided value. Default: 0.05

  batch_metadata_file:
    type: File?
    label: "Batch metadata file (optional)"
    doc: |
      Optional headerless TSV/CSV file for
      batch effect correction. First column
      should include sample names from the
      both control and treatment groups in
      an arbitrary order. The second column
      should include categories to which
      each sample should be assigned.
      Default: do not apply batch correction.

  blacklist_regions_file:
    type: File?
    label: "Genomic blacklist regions file (optional)"
    doc: |
      BED file with the genomic blacklist
      regions. Any reference genomic bin
      overlapping a blacklist region will
      be removed from the analysis.
      Default: include all reference
      genomic bins.

  remove_duplicated_reads:
    type: boolean?
    default: false
    label: "Remove duplicated reads"
    doc: |
      Remove duplicated reads identified by
      their coordinates. The location of a
      single-end read is determined by its
      strand and 5' coordinate. For a paired-
      end read, the DNA fragment position is
      used instead. Default: include all reads.
    "sd:layout":
      advanced: true

  bin_size:
    type: int?
    default: 2000
    label: "The size of the reference genomic bins"
    doc: |
      The size of non-overlapping reference
      genomic bins used for generating a table
      of read per peak counts. 2000 bp is
      recommended for sharp histone marks like
      H3K4me3 and H3K27ac, and 1000 bp for TFs
      or DNase-seq. Default: 2000
    "sd:layout":
      advanced: true

  fixed_bin_size:
    type: boolean?
    default: false
    label: "Force equal size for all reference genomic bins"
    doc: |
      Force all reference genomic bins be
      exaclty the same size (as selected).
      Default: when a merged peak is split
      into the reference genomic bins, their
      sizes do not exceed the boundaries of
      the original peak.
    "sd:layout":
      advanced: true

  shift_size:
    type: int?
    default: 100
    label: "Shift size for single-end reads"
    doc: |
      Shift the positions of the 5' ends
      of a single-end reads downstream on
      the selected value. Use the resulting
      points for counting reads in the
      reference genomic bins. Set as half
      of the DNA fragment size. Ignored if
      analysis is run in the paired-end mode.
      Default 100
    "sd:layout":
      advanced: true

  exclude_chromosomes:
    type: string?
    default: null
    label: "Chromosomes to be exluded from the analysis"
    doc: |
      A comma- or space-separated list of
      chromosomes to be exluded from the
      analysis. Default: include all
      chromosomes
    "sd:layout":
      advanced: true

  minimum_overlap:
    type: float?
    default: 1
    label: "Minimum peak overlap between the samples"
    doc: |
      Filtering threshold to keep only those
      reference genomic bins that are present
      in at least this many samples within
      each group. If this threshold has a value
      between zero and one, only those peaks
      will be included that are present in at
      least this fraction of samples.
      Default: 1
    "sd:layout":
      advanced: true

  promoter_dist:
    type: int?
    default: 1000
    label: "Promoter distance, bp"
    doc: |
      Max distance from the gene TSS (in
      both direction) overlapping which
      the differentially bound site will
      be assigned to the promoter region.
      Default: 1000 bp
    "sd:layout":
      advanced: true

  upstream_dist:
    type: int?
    default: 20000
    label: "Upstream distance, bp"
    doc: | 
      Max distance from the promoter (only
      in upstream direction) overlapping
      which the differentially bound site
      will be assigned to the upstream region.
      Default: 20,000 bp
    "sd:layout":
      advanced: true

  threads:
    type:
    - "null"
    - type: enum
      symbols:
      - "1"
      - "2"
      - "3"
      - "4"
      - "5"
      - "6"
    default: "4"
    label: "Cores/CPUs"
    doc: |
      Parallelization parameter to define the
      number of cores/CPUs that can be utilized
      simultaneously.
      Default: 4
    "sd:layout":
      advanced: true


outputs:

  peak_profile_bins_xls:
    type: File
    outputSource: manorm/peak_profile_bins_xls
    label: "Read counts and enrichment status in each sample"
    doc: |
      Read counts and enrichment status for
      each deduced reference genomic bin in
      each sample.

  coverage_files_cond_1:
    type: File[]
    label: "Control group sample(s)"
    doc: |
      Genome coverage files with the normalized
      number of fragments at each base from the
      samples that belong to the control group.
    outputSource: pipe/coverage_files_cond_1
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "Coverage (control)"
        height: 120

  coverage_files_cond_2:
    type: File[]
    label: "Treatment group sample(s)"
    doc: |
      Genome coverage files with the normalized
      number of fragments at each base from the
      samples that belong to the treatment group.
    outputSource: pipe/coverage_files_cond_2
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "wig"
        name: "Coverage (treatment)"
        height: 120

  n_peak_files_cond_1:
    type:
    - "null"
    - File[]
    label: "Control group sample(s)"
    doc: |
      Narrow peak files with the peaks called
      from the samples that belong to the
      control group.
    outputSource: pipe/n_peak_files_cond_1
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "Peaks (control)"
        displayMode: "COLLAPSE"
        height: 40

  n_peak_files_cond_2:
    type:
    - "null"
    - File[]
    label: "Treatment group sample(s)"
    doc: |
      Narrow peak files with the peaks called
      from the samples that belong to the
      treatment group.
    outputSource: pipe/n_peak_files_cond_2
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "Peaks (treatment)"
        displayMode: "COLLAPSE"
        height: 40

  b_peak_files_cond_1:
    type:
    - "null"
    - File[]
    label: "Control group sample(s)"
    doc: |
      Broad peak files with the peaks called from
      the samples that belong to the control group.
    outputSource: pipe/b_peak_files_cond_1
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "Peaks (control)"
        displayMode: "COLLAPSE"
        height: 40

  b_peak_files_cond_2:
    type:
    - "null"
    - File[]
    label: "Treatment group sample(s)"
    doc: |
      Broad peak files with the peaks called from
      the samples that belong to the treatment group.
    outputSource: pipe/b_peak_files_cond_2
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        name: "Peaks (treatment)"
        displayMode: "COLLAPSE"
        height: 40

  diff_rgns_bigbed:
    type: File
    label: "Differentially bound sites (bigBed format)"
    doc: |
      Differentially bound sites.
      bigBed format.
    outputSource: bed_to_bigbed/bigbed_file
    "sd:visualPlugins":
    - igvbrowser:
        tab: "Genome Browser"
        id: "igvbrowser"
        type: "annotation"
        format: "bigbed"
        name: "Differentially bound sites"
        height: 40

  smpl_corr_raw_plot_png:
    type: File?
    label: "Read counts correlation between the samples (raw)"
    doc: |
      Read counts correlation between the samples.
      On the basis of the raw read counts within
      the reference genomic bins.
      PNG format.
    outputSource: manorm/smpl_corr_raw_plot_png
    "sd:visualPlugins":
    - image:
        tab: "Exploratory plots"
        Caption: "Read counts correlation between the samples (raw)"

  smpl_corr_crtd_plot_png:
    type: File?
    label: "Read counts correlation between the samples (batch corrected)"
    doc: |
      Read counts correlation between the samples.
      On the basis of the batch corrected raw read
      counts within the reference genomic bins
      PNG format.
    outputSource: manorm/smpl_corr_crtd_plot_png
    "sd:visualPlugins":
    - image:
        tab: "Exploratory plots"
        Caption: "Read counts correlation between the samples (batch corrected)"

  smpl_corr_norm_plot_png:
    type: File?
    label: "Read counts correlation between the samples (normalized)"
    doc: |
      Read counts correlation between the samples.
      On the basis of the optionally batch corrected
      normalized read counts within the reference
      genomic bins.
      PNG format.
    outputSource: manorm/smpl_corr_norm_plot_png
    "sd:visualPlugins":
    - image:
        tab: "Exploratory plots"
        Caption: "Read counts correlation between the samples (normalized)"

  smpl_vrlp_plot_png:
    type: File?
    label: "Peaks overlap between the samples"
    doc: |
      Peaks overlap between the samples. On the
      basis of the occupied by each sample
      reference genomic bins.
      PNG format.
    outputSource: manorm/smpl_vrlp_plot_png
    "sd:visualPlugins":
    - image:
        tab: "Exploratory plots"
        Caption: "Peaks overlap between the samples"

  cnd_vrlp_plot_png:
    type: File?
    label: "Peaks overlap between the biological conditions"
    doc: |
      Peaks overlap between the biological conditions.
      On the basis of the occupied by each biological
      condition reference genomic bins.
      PNG format.
    outputSource: manorm/cnd_vrlp_plot_png
    "sd:visualPlugins":
    - image:
        tab: "Exploratory plots"
        Caption: "Peaks overlap between the biological conditions"

  ma_corr_plot_png:
    type: File?
    label: "Correlation between M and A values"
    doc: |
      Correlation between M and A values across the
      common peak regions of either each pair of
      biological conditions or each pair of samples.
      PNG format.
    outputSource: manorm/ma_corr_plot_png
    "sd:visualPlugins":
    - image:
        tab: "Exploratory plots"
        Caption: "Correlation between M and A values"

  diff_vlcn_plot_png:
    type: File?
    label: "Volcano plot for differentially bound sites"
    doc: |
      Volcano plot for differentially bound sites.
      PNG format.
    outputSource: manorm/diff_vlcn_plot_png
    "sd:visualPlugins":
    - image:
        tab: "Differential plots"
        Caption: "Volcano plot for differentially bound sites"

  diff_ma_plot_png:
    type: File?
    label: "MA-plot for differentially bound sites"
    doc: |
      MA-plot for differentially bound sites.
      PNG format.
    outputSource: manorm/diff_ma_plot_png
    "sd:visualPlugins":
    - image:
        tab: "Differential plots"
        Caption: "MA-plot for differentially bound sites"

  pca_1_2_plot_png:
    type: File?
    label: "Read counts PCA (PC1/PC2)."
    doc: |
      Read counts PCA (PC1/PC2).
      PNG format.
    outputSource: manorm/pca_1_2_plot_png
    "sd:visualPlugins":
    - image:
        tab: "Exploratory plots"
        Caption: "PCA (1,2) of not filtered normalized counts"

  pca_2_3_plot_png:
    type: File?
    label: "Read counts PCA (PC2/PC3)"
    doc: |
      Read counts PCA (PC2/PC3).
      PNG format.
    outputSource: manorm/pca_2_3_plot_png
    "sd:visualPlugins":
    - image:
        tab: "Exploratory plots"
        Caption: "PCA (2,3) of not filtered normalized counts"

  mds_plot_html:
    type: File?
    outputSource: manorm/mds_plot_html
    label: "MDS plot"
    doc: |
      MDS plot of optionally batch corrected
      normalized read counts within the
      reference genomic bins.
      HTML format.
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  diff_rgns_tsv:
    type: File
    label: "Differentially bound sites with assigned nearest genes"
    doc: |
      Differentially bound sites, not filtered
      by adjusted P-value threshold, with the
      assigned nearest genes.
      TSV format.
    outputSource: restore_columns/output_file
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Differentially bound sites"
        Title: "Differentially bound sites"

  diff_rgns_labeled_tsv:
    type: File
    label: "Differentially bound sites with labels"
    doc: |
      Differentially bound sites, not filtered
      by adjusted P-value threshold, with the
      labels.
      TSV format.
    outputSource: add_label_column/output_file

  volcano_plot_html_file:
    type: File
    label: "Volcano Plot"
    doc: |
      Volcano Plot html index.
    outputSource: make_volcano_plot/html_file
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  volcano_plot_html_data:
    type: Directory
    label: "Volcano Plot (data)"
    doc: |
      Volcano Plot html data.
    outputSource: make_volcano_plot/html_data

  ma_plot_html_file:
    type: File
    label: "MA-plot"
    doc: |
      MA-plot html index.
    outputSource: make_ma_plot/html_file
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  ma_plot_html_data:
    type: Directory
    label: "MA-plot (data)"
    doc: |
      MA-plot html data.
    outputSource: make_ma_plot/html_data

  heatmap_html:
    type: File
    label: "Heatmap of normalized read counts within the reference genomic bins"
    doc: |
      Morpheus heatmap html index.
    outputSource: morpheus_heatmap/heatmap_html
    "sd:visualPlugins":
    - linkList:
        tab: "Overview"
        target: "_blank"

  read_cnts_gct:
    type: File
    label: "Normalized read counts within the reference genomic bins"
    doc: |
      Optionally batch corrected normalized
      read counts within the reference
      genomic bins.
      GCT format
    outputSource: extend_gct/extended_gct

  experiment_info:
    type: File
    label: "Samples order for IGV"
    doc: |
      Markdown file to explain the sample
      order for IGV.
    outputSource: create_metadata/output_file
    "sd:visualPlugins":
    - markdownView:
        tab: "Overview"

  pdf_plots:
    type: File
    outputSource: compress_pdf_plots/compressed_folder
    label: "Compressed folder with all PDF plots"
    doc: |
      Compressed folder with all PDF plots.

  manorm_stdout_log:
    type: File
    label: "MAnorm output log"
    doc: |
      Stdout log from the manorm step.
    outputSource: manorm/stdout_log

  manorm_stderr_log:
    type: File
    label: "MAnorm error log"
    doc: |
      Stderr log from the manorm step.
    outputSource: manorm/stderr_log

  morpheus_stdout_log:
    type: File
    label: "Morpheus output log"
    doc: |
      Stdout log from the morpheus_heatmap step.
    outputSource: morpheus_heatmap/stdout_log

  morpheus_stderr_log:
    type: File
    label: "Morpheus error log"
    doc: |
      Stderr log from the morpheus_heatmap step.
    outputSource: morpheus_heatmap/stderr_log


steps:

  pipe:
    run:
      cwlVersion: v1.0
      class: ExpressionTool
      inputs:
        genome_cov_files_cond_1:
          type: File[]
        genome_cov_files_cond_2:
          type: File[]
        narrow_peak_files_cond_1:
          type:
          - "null"
          - File[]
        narrow_peak_files_cond_2:
          type:
          - "null"
          - File[]
        broad_peak_files_cond_1:
          type:
          - "null"
          - File[]
        broad_peak_files_cond_2:
          type:
          - "null"
          - File[]
      outputs:
        coverage_files_cond_1:
          type: File[]
        coverage_files_cond_2:
          type: File[]
        n_peak_files_cond_1:
          type:
          - "null"
          - File[]
        n_peak_files_cond_2:
          type:
          - "null"
          - File[]
        b_peak_files_cond_1:
          type:
          - "null"
          - File[]
        b_peak_files_cond_2:
          type:
          - "null"
          - File[]
      expression: |
        ${
          var results = {};
          var output_names = [
            "coverage_files_cond_1",
            "coverage_files_cond_2",
            "n_peak_files_cond_1",
            "n_peak_files_cond_2",
            "b_peak_files_cond_1",
            "b_peak_files_cond_2"
          ];
          var sources = [
            inputs.genome_cov_files_cond_1,
            inputs.genome_cov_files_cond_2,
            inputs.narrow_peak_files_cond_1,
            inputs.narrow_peak_files_cond_2,
            inputs.broad_peak_files_cond_1,
            inputs.broad_peak_files_cond_2
          ];
          for (var i = 0; i < sources.length; i++){
            var current_source = sources[i];
            var current_output_name = output_names[i];
            results[current_output_name] = null;

            if (current_source != null && current_source.length > 0){
              for (var j = 0; j < current_source.length; j++){
                    var new_item = current_source[j];
                    new_item["basename"] = "u" + "_" + i + "_" + j+ "_" + new_item.basename;
                    if (results[current_output_name] == null){
                      results[current_output_name] = [new_item];
                    } else {
                      results[current_output_name].push(new_item);
                    }
              }
            }
          }
          return results;
        }
    in:
      genome_cov_files_cond_1: genome_cov_files_cond_1
      genome_cov_files_cond_2: genome_cov_files_cond_2
      narrow_peak_files_cond_1: narrow_peak_files_cond_1
      narrow_peak_files_cond_2: narrow_peak_files_cond_2
      broad_peak_files_cond_1: broad_peak_files_cond_1
      broad_peak_files_cond_2: broad_peak_files_cond_2
    out:
    - coverage_files_cond_1
    - coverage_files_cond_2
    - n_peak_files_cond_1
    - n_peak_files_cond_2
    - b_peak_files_cond_1
    - b_peak_files_cond_2

  manorm:
    run: ../tools/manorm2.cwl
    in:
      read_files_cond_1: read_files_cond_1
      read_files_cond_2: read_files_cond_2
      peak_files_cond_1:
        source:
        - narrow_peak_files_cond_1         # [0]
        - narrow_peak_files_cond_2         # [1]
        - broad_peak_files_cond_1          # [2]
        - broad_peak_files_cond_2          # [3]
        valueFrom: |
          ${
            if (self[2] && self[3]){
              return self[2];
            }
            else {
              return self[0];
            }
          }
      peak_files_cond_2:
        source:
        - narrow_peak_files_cond_1         # [0]
        - narrow_peak_files_cond_2         # [1]
        - broad_peak_files_cond_1          # [2]
        - broad_peak_files_cond_2          # [3]
        valueFrom: |
          ${
            if (self[2] && self[3]){
              return self[3];
            }
            else {
              return self[1];
            }
          }
      sample_names_cond_1: sample_names_cond_1
      sample_names_cond_2: sample_names_cond_2
      summit_files_cond_1: summit_files_cond_1
      summit_files_cond_2: summit_files_cond_2
      condition_1: condition_1
      condition_2: condition_2
      minimum_overlap: minimum_overlap
      maximum_padj: maximum_padj
      batch_metadata_file: batch_metadata_file
      maximum_peak_number:
        source: maximum_peak_number
        valueFrom: |
          ${
            if (self == 0){
              return null;
            }
            else {
              return self;
            }
          }
      minimum_peak_gap: minimum_peak_gap
      bin_size: bin_size
      fixed_bin_size: fixed_bin_size
      blacklist_regions_file: blacklist_regions_file
      remove_duplicated_reads: remove_duplicated_reads
      shift_size: shift_size
      paired_end: paired_end
      exclude_chromosomes:
        source: exclude_chromosomes
        valueFrom: $(split_features(self))
      normalization_method: normalization_method
      export_pdf_plots:
        default: true
      parallel_memory_limit:
        default: 32
      threads:
        source: threads
        valueFrom: $(parseInt(self))
    out:
    - peak_profile_bins_xls
    - diff_rgns_tsv
    - smpl_corr_raw_plot_png
    - smpl_corr_crtd_plot_png
    - smpl_corr_norm_plot_png
    - smpl_vrlp_plot_png
    - cnd_vrlp_plot_png
    - ma_corr_plot_png
    - diff_vlcn_plot_png
    - diff_ma_plot_png
    - pca_1_2_plot_png
    - pca_2_3_plot_png
    - mds_plot_html
    - read_cnts_gct
    - all_plots_pdf
    - stderr_log
    - stdout_log

  folder_pdf_plots:
    run: ../tools/files-to-folder.cwl
    in:
      input_files:
        source:
        - manorm/all_plots_pdf
        valueFrom: $(self.flat().filter(n => n))
      folder_basename:
        default: "pdf_plots"
    out:
    - folder

  compress_pdf_plots:
    run: ../tools/tar-compress.cwl
    in:
      folder_to_compress: folder_pdf_plots/folder
    out:
    - compressed_folder

  filter_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file: manorm/diff_rgns_tsv
      script:
        default: >
          cat $0 | grep -v "Start" | awk
          'BEGIN {print "chr\tstart\tend\tlength\tabs_summit\tpileup\t-log10(pvalue)\tfold_enrichment\t-log10(qvalue)\tname"}
          {print $1"\t"$2"\t"$3"\t"$3-$2+1"\t0\t"NR"\t0\t0\t0\t0"}' > `basename $0`
    out:
    - output_file

  assign_genes:
    run: ../tools/iaintersect.cwl
    in:
      input_filename: filter_columns/output_file
      annotation_filename: annotation_file
      promoter_bp: promoter_dist
      upstream_bp: upstream_dist
    out:
    - result_file

  restore_columns:
    run: ../tools/custom-bash.cwl
    in:
      input_file:
      - assign_genes/result_file
      - manorm/diff_rgns_tsv
      script:
        default: |
          cat $0 | grep -v "start" | sort -k 11n | cut -f 1-5,15 > iaintersect_result.tsv
          cat $1 | grep -v "start" > manorm_result.tsv
          HEADER=`head -n 1 $1`;
          echo -e "refseq_id\tgene_id\ttxStart\ttxEnd\tstrand\tregion\t${HEADER}" > `basename $0`;
          cat iaintersect_result.tsv | paste - manorm_result.tsv >> `basename $0`
          rm iaintersect_result.tsv manorm_result.tsv
    out:
    - output_file

  convert_to_bed:
    run: ../tools/custom-bash.cwl
    in:
      input_file: restore_columns/output_file
      param:
        source: maximum_padj
        valueFrom: $(self + "")                  # to convert it to string
      script:
        default: |
          cat "$0" | awk -F "\t" -v maximum_padj="$1" 'NR==1 {for (i=1; i<=NF; i++) {ix[$i]=i} } NR>1 && $ix["padj"]<=maximum_padj {color="255,0,0"; if ($ix["log2FoldChange"]<0) color="0,255,0"; print $ix["chr"]"\t"$ix["start"]"\t"$ix["end"]"\tpvalue="$ix["pvalue"]+0.0";padj="$ix["padj"]+0.0";log2FC="$ix["log2FoldChange"]"\t"1000"\t"$ix["strand"]"\t"$ix["start"]"\t"$ix["end"]"\t"color}' > `basename $0`
    out:
    - output_file

  sort_bed:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: convert_to_bed/output_file
      key:
        default: ["1,1","2,2n"]
    out:
    - sorted_file

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
    out: 
    - bigbed_file

  add_label_column:
    run: ../tools/custom-bash.cwl
    in:
      input_file: manorm/diff_rgns_tsv
      script:
        default: |
          HEADER=`head -n 1 $0`;
          echo -e "label\t${HEADER}" > diff_rgns_labeled.tsv;
          cat "$0" | grep -v "start" | awk -F "\t" '{print $1":"$2"-"$3"\t"$0}' >> diff_rgns_labeled.tsv
    out:
    - output_file

  make_volcano_plot:
    run: ../tools/volcano-plot.cwl
    in:
      diff_expr_file: add_label_column/output_file
      x_axis_column:
        default: "log2FoldChange"
      y_axis_column:
        default: "padj"
      label_column:
        default: "label"
    out:
    - html_data
    - html_file

  make_ma_plot:
    run: ../tools/ma-plot.cwl
    in:
      diff_expr_file: add_label_column/output_file
      x_axis_column:
        default: "baseMean"
      y_axis_column:
        default: "log2FoldChange"
      label_column:
        default: "label"
    out:
    - html_data
    - html_file

  extend_gct:
    run:
      cwlVersion: v1.0
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/morpheus:v0.0.2
      - class: InitialWorkDirRequirement
        listing:
        - entryname: extend.R
          entry: |
            options(error=function(){traceback(3); quit(save="no", status=1, runLast=FALSE)})
            suppressMessages(library("cmapR"))
            suppressMessages(library("dplyr"))
            suppressMessages(library("tibble"))
            suppressMessages(library("morpheus"))
            suppressMessages(library("argparse"))
            args = commandArgs(trailingOnly=TRUE)
            gct_data <- read.gct(args[1])
            metadata <- read.table(args[2], sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE) %>%
                        mutate(id=paste(chr, paste(start, end, sep="-"), sep=":")) %>%
                        select(id, gene_id, region)
            row_metadata <- gct_data$rowAnnotations %>%
                            rownames_to_column("id") %>%
                            left_join(metadata, by="id") %>%
                            mutate_at("id", as.vector)
            col_metadata <- gct_data$columnAnnotations %>%
                            rownames_to_column("id") %>%
                            mutate_at("id", as.vector)
            gct_data <- new(
                "GCT",
                mat=gct_data$data[row_metadata$id, col_metadata$id],
                rdesc=row_metadata,
                cdesc=col_metadata
            )
            write_gct(ds=gct_data, ofile="extended.gct", appenddim=FALSE)
      inputs:
        input_files:
          type: File[]
          inputBinding:
            position: 5
      outputs:
        extended_gct:
          type: File
          outputBinding:
            glob: "extended.gct"
      baseCommand: ["Rscript", "extend.R"]
    in:
      input_files:
      - manorm/read_cnts_gct
      - restore_columns/output_file
    out:
    - extended_gct

  morpheus_heatmap:
    run: ../tools/morpheus-heatmap.cwl
    in:
     read_counts_gct: extend_gct/extended_gct
    out:
    - heatmap_html
    - stdout_log
    - stderr_log

  create_metadata:
    run: ../tools/custom-bash.cwl
    in:
      input_file:
        source:
        - read_files_cond_1
        - read_files_cond_2
        valueFrom: $(self.flat().filter(n => n))
      param:
        source:
        - sample_names_cond_1
        - sample_names_cond_2
        valueFrom: $(self.flat().filter(n => n))
      script:
        default: |
          #!/bin/bash
          set -- "$0" "$@"
          COUNT=`expr $# / 2`
          echo "| Sample | Index |" > experiment_info.md
          echo "| :-- | --: |" >> experiment_info.md
          j=1
          for i in "${@:$COUNT+1:$#}"; do
            echo "| $i | $j |" >> experiment_info.md
            (( j++ ))
          done;
    out:
    - output_file


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "MAnorm2 for Normalizing and Comparing ChIP-Seq/ATAC-Seq Samples"
s:name: "MAnorm2 for Normalizing and Comparing ChIP-Seq/ATAC-Seq Samples"
s:alternateName: "MAnorm2 for Normalizing and Comparing ChIP-Seq/ATAC-Seq Samples"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/manorm2.cwl
s:codeRepository: https://github.com/Barski-lab/workflows-datirium
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

doc: |
  MAnorm2 for Normalizing and Comparing ChIP-Seq/ATAC-Seq Samples