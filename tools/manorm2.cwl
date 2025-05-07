cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/manorm2:v0.0.2


inputs:

  read_files_cond_1:
    type:
    - File
    - type: array
      items: File
    inputBinding:
      prefix: "--read1"
    doc: |
      Coordinate sorted and indexed BAM files with
      aligned reads from the samples that belong
      to the first biological condition.

  read_files_cond_2:
    type:
    - File
    - type: array
      items: File
    inputBinding:
      prefix: "--read2"
    doc: |
      Coordinate sorted and indexed BAM files with
      aligned reads from the samples that belong
      to the second biological condition.

  peak_files_cond_1:
    type:
    - File
    - type: array
      items: File
    inputBinding:
      prefix: "--peak1"
    doc: |
      Narrow or broad peak files with the peaks
      called from the samples that belong to the
      first biological condition.

  peak_files_cond_2:
    type:
    - File
    - type: array
      items: File
    inputBinding:
      prefix: "--peak2"
    doc: |
      Narrow or broad peak files with the peaks
      called from the samples that belong to the
      second biological condition.

  sample_names_cond_1:
    type:
    - string
    - type: array
      items: string
    inputBinding:
      prefix: "--name1"
    doc: |
      Names of the samples that belong to the first
      biological condition. All values from the --name1
      and --name2 parameters should be unique.

  sample_names_cond_2:
    type:
    - string
    - type: array
      items: string
    inputBinding:
      prefix: "--name2"
    doc: |
      Names of the samples that belong to the second
      biological condition. All values from the --name1
      and --name2 parameters should be unique.

  summit_files_cond_1:
    type:
    - "null"
    - File
    - type: array
      items: File
    inputBinding:
      prefix: "--summit1"
    doc: |
      BED files with the summits of the peaks called
      from the samples that belong to the first
      biological condition. If not provided, the peak
      center is taken as the summit.

  summit_files_cond_2:
    type:
    - "null"
    - File
    - type: array
      items: File
    inputBinding:
      prefix: "--summit2"
    doc: |
      BED files with the summits of the peaks called
      from the samples that belong to the second
      biological condition. If not provided, the peak
      center is taken as the summit.

  condition_1:
    type:
    - "null"
    - string
    inputBinding:
      prefix: "--condition1"
    doc: |
      Name for the first biological condition. The
      direction of comparison is always --condition2
      vs --condition1. Default: control.

  condition_2:
    type:
    - "null"
    - string
    inputBinding:
      prefix: "--condition2"
    doc: |
      Name for the second biological condition. The
      direction of comparison is always --condition2
      vs --condition1. Default: treatment.

  minimum_overlap:
    type:
    - "null"
    - float
    inputBinding:
      prefix: "--minoverlap"
    doc: |
      Filtering threshold to keep only those reference
      genomic bins that are present in at least this
      many samples within the biological condition. If
      this threshold has a value between zero and one,
      only those peaks will be included that are present in
      at least this fraction of datasets. Default: 1

  maximum_padj:
    type:
    - "null"
    - float
    inputBinding:
      prefix: "--padj"
    doc: |
      Filtering threshold to report only differentially
      bound sites with adjusted P-value less than or equal
      to the provided value. Default: 0.05

  batch_metadata_file:
    type:
    - "null"
    - File
    inputBinding:
      prefix: "--batch"
    doc: |
      Optional headerless TSV/CSV file for batch effect
      correction. First column should include values from
      the --name1 and --name2 parameters in any order.
      The second column should include group names to
      which each sample should be assigned. Default: do
      not apply batch correction.

  maximum_peak_number:
    type:
    - "null"
    - int
    inputBinding:
      prefix: "--maxpeaks"
    doc: |
      The maximum number of the most significant peaks
      to select from each peak file when constructing
      reference genomic bins. The top significant peaks
      are selected based on the score column which is
      calculated by MACS2 either as int(-10*log10qvalue)
      or as int(-10*log10qvalue), depending on the cutoff
      used for peak calling. Default: keep all peaks.

  minimum_peak_gap:
    type:
    - "null"
    - int
    inputBinding:
      prefix: "--minpeakgap"
    doc: |
      Peaks remained after optional filtering by --maxpeaks
      parameter will be merged if the distance between them
      is smaller than the provided value. Merging is first
      applied per sample and then to all peaks together
      before splitting them into the reference genomic bins
      of size --binsize. Default: 150

  bin_size:
    type:
    - "null"
    - int
    inputBinding:
      prefix: "--binsize"
    doc: |
      The size of non-overlapping reference genomic bins used
      for generating a table of read per peak counts. 2000 bp
      is recommended for sharp histone marks like H3K4me3 and
      H3K27ac, and 1000 bp for TFs or DNase-seq. Default: 2000

  fixed_bin_size:
    type:
    - "null"
    - boolean
    inputBinding:
      prefix: "--fixbinsize"
    doc: |
      Force all reference genomic bins be exaclty the
      same size as provided in the --binsize parameter.
      Default: when a merged peak is split into the
      reference genomic bins, their sizes do not exceed
      the boundaries of the original peak.

  blacklist_regions_file:
    type:
    - "null"
    - File
    inputBinding:
      prefix: "--blacklist"
    doc: |
      BED file with the genomic blacklist regions. Any
      reference genomic bin overlapping a blacklist region
      will be removed from the analysis. Default: include
      all reference genomic bins.

  remove_duplicated_reads:
    type:
    - "null"
    - boolean
    inputBinding:
      prefix: "--dedup"
    doc: |
      Remove duplicated reads identified by their coordinates.
      The location of a single-end read is determined by its
      strand and 5' coordinate. For a paired-end read, the
      DNA fragment position is used instead. Default: include
      all reads.

  shift_size:
    type:
    - "null"
    - int
    inputBinding:
      prefix: "--shiftsize"
    doc: |
      Shift the positions of the 5' ends of a single-end
      reads downstream on the selected value. Use the resulting
      points for counting reads in the reference genomic bins.
      Set as half of the DNA fragment size. Ignored if --paired
      parameter is provided. Default 100

  paired_end:
    type:
    - "null"
    - boolean
    inputBinding:
      prefix: "--paired"
    doc: |
      Consider all reads as paired-end. When counting reads
      in the reference genomic bins, use the middle point of
      each DNA fragment. --shiftsize parameters is ignored.
      Default: treat all reads as single-end.

  exclude_chromosomes:
    type:
    - "null"
    - string
    - type: array
      items: string
    inputBinding:
      prefix: "--exclude"
    doc: |
      Define the chromosomes to be exluded from the analysis.
      Default: include all chromosomes

  normalization_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "pseudo-reference"
      - "baseline"
      - "hierarchical"
    inputBinding:
      prefix: "--norm"
    doc: |
      Normalization method applied to the raw read counts.
      pseudo-reference - normalize each sample to the pseudo
      reference that includes the average intesities from all
      samples. A reference genomic bin is occupied by the
      pseudo reference if it was occupied by at least one
      sample that the reference was constructed from. Each
      sample is MA-normalized to the pseudo reference using
      the common genomic bins between the reference and a
      sample. baseline - normalize each sample to the one
      whose log2 size factor is closest to 0. hierarchical -
      similar to the baseline but first all samples are
      normalized within the biological conditions, than two
      biological conditions are normalized between each other.
      Default: pseudo-reference

  export_pdf_plots:
    type:
    - "null"
    - boolean
    inputBinding:
      prefix: "--pdf"
    doc: |
      Export plots in PDF. Default: false

  output_prefix:
    type:
    - "null"
    - string
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix. Default: ./manorm

  parallel_memory_limit:
    type:
    - "null"
    - int
    inputBinding:
      prefix: "--memory"
    doc: |
      Maximum memory in GB allowed to be shared between
      the workers when using multiple --cpus.
      Default: 32

  threads:
    type:
    - "null"
    - int
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of cores/cpus to use. Default: 1

  seed:
    type:
    - "null"
    - int
    inputBinding:
      prefix: "--seed"
    doc: |
      Seed number for random values.
      Default: 42


outputs:

  peak_profile_bins_xls:
    type: File
    outputBinding:
      glob: "peak_profile_bins.xls"
    doc: |
      Read counts and enrichment status for
      each deduced reference genomic bin in
      each sample.

  diff_rgns_tsv:
    type: File
    outputBinding:
      glob: "*_diff_rgns.tsv"
    doc: |
      Differentially bound sites, not filtered
      by adjusted P-value threshold.
      TSV format.

  smpl_corr_raw_plot_png:
    type:
    - "null"
    - File
    outputBinding:
      glob: "*_smpl_corr_raw.png"
    doc: |
      Read counts correlation between the samples.
      On the basis of the raw read counts within
      the reference genomic bins.
      PNG format.

  smpl_corr_crtd_plot_png:
    type:
    - "null"
    - File
    outputBinding:
      glob: "*_smpl_corr_crtd.png"
    doc: |
      Read counts correlation between the samples.
      On the basis of the batch corrected raw read
      counts within the reference genomic bins
      PNG format.

  smpl_corr_norm_plot_png:
    type:
    - "null"
    - File
    outputBinding:
      glob: "*_smpl_corr_norm.png"
    doc: |
      Read counts correlation between the samples.
      On the basis of the optionally batch corrected
      normalized read counts within the reference
      genomic bins.
      PNG format.

  smpl_vrlp_plot_png:
    type:
    - "null"
    - File
    outputBinding:
      glob: "*_smpl_vrlp.png"
    doc: |
      Peaks overlap between the samples. On the
      basis of the occupied by each sample
      reference genomic bins.
      PNG format.

  cnd_vrlp_plot_png:
    type:
    - "null"
    - File
    outputBinding:
      glob: "*_cnd_vrlp.png"
    doc: |
      Peaks overlap between the biological conditions.
      On the basis of the occupied by each biological
      condition reference genomic bins.
      PNG format.

  ma_corr_plot_png:
    type:
    - "null"
    - File
    outputBinding:
      glob: "*_ma_corr.png"
    doc: |
      Correlation between M and A values across the
      common peak regions of either each pair of
      biological conditions or each pair of samples.
      PNG format.

  diff_vlcn_plot_png:
    type:
    - "null"
    - File
    outputBinding:
      glob: "*_diff_vlcn.png"
    doc: |
      Volcano plot for differentially bound sites.
      PNG format.

  diff_ma_plot_png:
    type:
    - "null"
    - File
    outputBinding:
      glob: "*_diff_ma.png"
    doc: |
      MA-plot for differentially bound sites.
      PNG format.

  pca_1_2_plot_png:
    type:
    - "null"
    - File
    outputBinding:
      glob: "*_pca_1_2.png"
    doc: |
      Read counts PCA (PC1/PC2)
      PNG format.

  pca_2_3_plot_png:
    type:
    - "null"
    - File
    outputBinding:
      glob: "*_pca_2_3.png"
    doc: |
      Read counts PCA (PC2/PC3)
      PNG format.

  mds_plot_html:
    type:
    - "null"
    - File
    outputBinding:
      glob: "*_mds_plot.html"
    doc: |
      MDS plot of optionally batch corrected
      normalized read counts within the
      reference genomic bins.
      HTML format.

  read_cnts_gct:
    type:
    - "null"
    - File
    outputBinding:
      glob: "*_read_cnts.gct"
    doc: |
      Optionally batch corrected normalized
      read counts within the reference
      genomic bins.
      GCT format

  all_plots_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*.pdf"
    doc: |
      All generated plots.
      PDF format.

  human_log:
    type: File?
    outputBinding:
      glob: "error_report.txt"
    doc: |
      Human readable error log.
      TXT format.

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["run_manorm2.R"]
stderr: manorm_stderr.log
stdout: manorm_stdout.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "MAnorm2 for Normalizing and Comparing ChIP-Seq/ATAC-Seq Samples"
s:name: "MAnorm2 for Normalizing and Comparing ChIP-Seq/ATAC-Seq Samples"
s:alternateName: "MAnorm2 for Normalizing and Comparing ChIP-Seq/ATAC-Seq Samples"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/manorm2.cwl
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
  MAnorm2 for Normalizing and Comparing ChIP-Seq/ATAC-Seq Samples

s:about: |
  usage: run_manorm2.R [-h] --read1 READ1 [READ1 ...] --read2
                            READ2 [READ2 ...] --peak1 PEAK1
                            [PEAK1 ...] --peak2 PEAK2 [PEAK2 ...]
                            [--summit1 [SUMMIT1 ...]]
                            [--summit2 [SUMMIT2 ...]] --name1 NAME1
                            [NAME1 ...] --name2 NAME2 [NAME2 ...]
                            [--condition1 CONDITION1]
                            [--condition2 CONDITION2]
                            [--minoverlap MINOVERLAP] [--padj PADJ]
                            [--batch BATCH] [--maxpeaks MAXPEAKS]
                            [--minpeakgap MINPEAKGAP]
                            [--binsize BINSIZE] [--fixbinsize]
                            [--blacklist BLACKLIST] [--dedup]
                            [--shiftsize SHIFTSIZE] [--paired]
                            [--exclude [EXCLUDE ...]]
                            [--norm {pseudo-reference,baseline,hierarchical}]
                            [--pdf] [--output OUTPUT] [--cpus CPUS]
                            [--memory MEMORY] [--tmpdir TMPDIR]
                            [--seed SEED]

  MAnorm2 for Normalizing and Comparing ChIP-seq Samples

  options:
    -h, --help            show this help message and exit
    --read1 READ1 [READ1 ...]
                          Coordinate sorted and indexed BAM files with aligned
                          reads from the samples that belong to the first
                          biological condition.
    --read2 READ2 [READ2 ...]
                          Coordinate sorted and indexed BAM files with aligned
                          reads from the samples that belong to the second
                          biological condition.
    --peak1 PEAK1 [PEAK1 ...]
                          Narrow or broad peak files with the peaks called from
                          the samples that belong to the first biological
                          condition.
    --peak2 PEAK2 [PEAK2 ...]
                          Narrow or broad peak files with the peaks called from
                          the samples that belong to the second biological
                          condition.
    --summit1 [SUMMIT1 ...]
                          BED files with the summits of the peaks called from
                          the samples that belong to the first biological
                          condition. If not provided, the peak center is taken
                          as the summit.
    --summit2 [SUMMIT2 ...]
                          BED files with the summits of the peaks called from
                          the samples that belong to the second biological
                          condition. If not provided, the peak center is taken
                          as the summit.
    --name1 NAME1 [NAME1 ...]
                          Names of the samples that belong to the first
                          biological condition. All values from the --name1 and
                          --name2 parameters should be unique.
    --name2 NAME2 [NAME2 ...]
                          Names of the samples that belong to the second
                          biological condition. All values from the --name1 and
                          --name2 parameters should be unique.
    --condition1 CONDITION1
                          Name for the first biological condition. The direction
                          of comparison is always --condition2 vs --condition1.
                          Default: control.
    --condition2 CONDITION2
                          Name for the second biological condition. The
                          direction of comparison is always --condition2 vs
                          --condition1. Default: treatment.
    --minoverlap MINOVERLAP
                          Filtering threshold to keep only those reference
                          genomic bins that are present in at least this many
                          samples within the biological condition. If this
                          threshold has a value between zero and one, only those
                          peaks will be included that are present in at least
                          this fraction of datasets. Default: 1
    --padj PADJ           Filtering threshold to report only differentially
                          bound sites with adjusted P-value less than or equal
                          to the provided value. Default: 0.05
    --batch BATCH         Optional headerless TSV/CSV file for batch effect
                          correction. First column should include values from
                          the --name1 and --name2 parameters in any order. The
                          second column should include group names to which each
                          sample should be assigned. Default: do not apply batch
                          correction.
    --maxpeaks MAXPEAKS   The maximum number of the most significant peaks to
                          select from each peak file when constructing reference
                          genomic bins. The top significant peaks are selected
                          based on the score column which is calculated by MACS2
                          either as int(-10*log10qvalue) or as
                          int(-10*log10qvalue), depending on the cutoff used for
                          peak calling. Default: keep all peaks.
    --minpeakgap MINPEAKGAP
                          Peaks remained after optional filtering by --maxpeaks
                          parameter will be merged if the distance between them
                          is smaller than the provided value. Merging is first
                          applied per sample and then to all peaks together
                          before splitting them into the reference genomic bins
                          of size --binsize. Default: 150
    --binsize BINSIZE     The size of non-overlapping reference genomic bins
                          used for generating a table of read per peak counts.
                          2000 bp is recommended for sharp histone marks like
                          H3K4me3 and H3K27ac, and 1000 bp for TFs or DNase-seq.
                          Default: 2000
    --fixbinsize          Force all reference genomic bins be exaclty the same
                          size as provided in the --binsize parameter. Default:
                          when a merged peak is split into the reference genomic
                          bins, their sizes do not exceed the boundaries of the
                          original peak.
    --blacklist BLACKLIST
                          BED file with the genomic blacklist regions. Any
                          reference genomic bin overlapping a blacklist region
                          will be removed from the analysis. Default: include
                          all reference genomic bins.
    --dedup               Remove duplicated reads identified by their
                          coordinates. The location of a single-end read is
                          determined by its strand and 5' coordinate. For a
                          paired-end read, the DNA fragment position is used
                          instead. Default: include all reads.
    --shiftsize SHIFTSIZE
                          Shift the positions of the 5' ends of a single-end
                          reads downstream on the selected value. Use the
                          resulting points for counting reads in the reference
                          genomic bins. Set as half of the DNA fragment size.
                          Ignored if --paired parameter is provided. Default 100
    --paired              Consider all reads as paired-end. When counting reads
                          in the reference genomic bins, use the middle point of
                          each DNA fragment. --shiftsize parameters is ignored.
                          Default: treat all reads as single-end.
    --exclude [EXCLUDE ...]
                          Define the chromosomes to be exluded from the
                          analysis. Default: include all chromosomes
    --norm {pseudo-reference,baseline,hierarchical}
                          Normalization method applied to the raw read counts.
                          pseudo-reference - normalize each sample to the pseudo
                          reference that includes the average intesities from
                          all samples. A reference genomic bin is occupied by
                          the pseudo reference if it was occupied by at least
                          one sample that the reference was constructed from.
                          Each sample is MA-normalized to the pseudo reference
                          using the common genomic bins between the reference
                          and a sample. baseline - normalize each sample to the
                          one whose log2 size factor is closest to 0.
                          hierarchical - similar to the baseline but first all
                          samples are normalized within the biological
                          conditions, than two biological conditions are
                          normalized between each other. Default: pseudo-
                          reference
    --pdf                 Export plots in PDF. Default: false
    --output OUTPUT       Output prefix. Default: ./manorm
    --cpus CPUS           Number of cores/cpus to use. Default: 1
    --memory MEMORY       Maximum memory in GB allowed to be shared between the
                          workers when using multiple --cpus. Default: 32
    --tmpdir TMPDIR       Directory to keep temporary files. Default: either
                          /tmp or defined by environment variables TMPDIR, TMP,
                          TEMP.
    --seed SEED           Seed number for random values. Default: 42
