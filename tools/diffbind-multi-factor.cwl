cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: DockerRequirement
  dockerPull: biowardrobe2/diffbind:v0.0.16
- class: InlineJavascriptRequirement
- class: InitialWorkDirRequirement
  listing: |
    ${
      var listing = [];
      for (var i = 0; i < inputs.alignment_files.length; i++){
        var alignment_file = inputs.alignment_files[i];
        var prefix = "u" + i + "_";
        alignment_file.basename = prefix + alignment_file.basename;
        if (alignment_file.secondaryFiles && alignment_file.secondaryFiles.length > 0){
          for (var j = 0; j < alignment_file.secondaryFiles.length; j++){
            var secondary_file = alignment_file.secondaryFiles[j];
            secondary_file.basename = prefix + secondary_file.basename;
            listing.push(secondary_file);
          }
          delete alignment_file.secondaryFiles;
        }
        listing.push(alignment_file);
      }
      return listing;
    }


inputs:

  alignment_files:
    type: File[]
    secondaryFiles:
    - .bai
    inputBinding:
      prefix: "--alignments"
    doc:
      Sorted and indexed alignment files in bam format

  peak_files:
    type: File[]
    inputBinding:
      prefix: "--peaks"
    doc:
      Peak files in the MACS2 xls format. Number and order of the
      files should correspond to the files provided in --alignments
      parameter.

  dataset_names:
    type: string[]
    inputBinding:
      prefix: "--aliases"
    doc: |
      Unique names for datasets provided in --alignments and --peaks
      parameters, no special characters or spaces are allowed. Number
      and order of the names should correspond to the values provided
      in --alignments and --peaks parameters.

  metadata_file:
    type: File
    inputBinding:
      prefix: "--metadata"
    doc: |
      TSV/CSV metadata file to describe datasets provided in --alignments
      and --peaks parameters. First column should have the name 'sample',
      all other columns names should be selected from the following list:
      Tissue, Factor, Condition, Treatment, Replicate. The values from the
      'sample' column should correspond to the values provided in --aliases
      parameter. For a proper --contrast intepretation, values defined in
      each metadata column should not be used in any of the other columns.
      All metadata columns are treated as factors (no covariates are supported).

  scoreby:
    type:
    - "null"
    - type: enum
      symbols:
      - "pvalue"
      - "qvalue"
    inputBinding:
      prefix: "--scoreby"
    doc: |
      Score metrics to exclude low quality peaks based on the
      threshold provided in the --score parameter.
      Default: pvalue

  score_threshold:
    type: float?
    inputBinding:
      prefix: "--score"
    doc: |
      Filtering threshold to keep only those peaks where the metric selected
      in --scoreby parameter is less than or equal to the provided value.
      Default: 0.05

  rpkm_threshold:
    type: float?
    inputBinding:
      prefix: "--minrpkm"
    doc: |
      Filtering threshold to keep only those peaks where the max RPKM for
      all datasets is bigger than or equal to the provided value.
      Default: 1

  rec_summits:
    type: int?
    inputBinding:
      prefix: "--summits"
    doc: |
      Width in bp to extend peaks around their summits in both directions
      and replace the original ones. Set it to 100 bp for ATAC-Seq and 200
      bp for ChIP-Seq datasets. To skip peaks extension and replacement, set
      it to negative value. Default: 200 bp (results in 401 bp wide peaks)

  overlap_threshold:
    type: float?
    inputBinding:
      prefix: "--minoverlap"
    doc: |
      Filtering threshold to keep only those peaks that are present in at
      least this many datasets when generating consensus set of peaks used
      in differential analysis. If this threshold has a value between zero
      and one, only those peaks will be included that are present in at least
      this proportion of datasets. When combined with --groupby parameter,
      --minoverlap threshold is applied per group, then union of the resulted
      peaks are used in the differential analysis. Default: 2

  groupby:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--groupby"
    doc: |
      Column(s) from the metadata table to define datasets grouping. --minoverlap
      filtering threshold will be applied within each group independently. Union
      of the resulted peaks from each of the groups will be used in the differential
      analysis. Default: apply --minoverlap filtering threshold for all datasets
      jointly

  design_formula:
    type: string
    inputBinding:
      prefix: "--design"
    doc: |
      Design formula comprised of the metadata columns names.
      It should start with ~

  contrast:
    type: string?
    inputBinding:
      prefix: "--contrast"
    doc: |
      Contrast applied to the analysis results when calculating log2 fold changes.
      It should be formatted as a mathematical formula of values present in the
      metadata table. It is a required parameter if --method is set to edger. If not
      provided and --method is set to deseq2, the last term from the design formula
      will be used.

  base_levels:
    type:
    - "null"
    - string
    - string[]
    inputBinding:
      prefix: "--base"
    doc: |
      Base levels for each of the metadata columns. Number and order of the provided
      values should correspond to the metadata columns. Default: define base levels
      alphabetically.

  analysis_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "deseq2"
      - "edger"
    inputBinding:
      prefix: "--method"
    doc: |
      Method used in the differential binding analysis. Should be equal to
      either edger or deseq2.
      Default: deseq2

  normalization_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "auto"
      - "rle"
      - "tmm"
      - "lib"
    inputBinding:
      prefix: "--norm"
    doc: |
      Normalization technique applied to the read counts before running differential
      binding analysis. When set to auto selects rle for deseq2 and tmm for edger.
      Default: auto

  padj_threshold:
    type: float?
    inputBinding:
      prefix: "--padj"
    doc: |
      Filtering threshold to report only differentially bound sites with adjusted
      P-value less than or equal to the provided value.
      Default: 0.05

  cluster_method:
    type:
    - "null"
    - type: enum
      symbols:
      - "row"
      - "column"
      - "both"
    inputBinding:
      prefix: "--cluster"
    doc: |
      Hopach clustering method to be run on normalized read counts.
      Default: do not run clustering

  row_distance:
    type:
    - "null"
    - type: enum
      symbols:
      - "cosangle"
      - "abscosangle"
      - "euclid"
      - "abseuclid"
      - "cor"
      - "abscor"
    inputBinding:
      prefix: "--rowdist"
    doc: |
      Distance metric for HOPACH row clustering. Ignored if --cluster is not
      provided.
      Default: cosangle

  column_distance:
    type:
    - "null"
    - type: enum
      symbols:
      - "cosangle"
      - "abscosangle"
      - "euclid"
      - "abseuclid"
      - "cor"
      - "abscor"
    inputBinding:
      prefix: "--columndist"
    doc: |
      Distance metric for HOPACH column clustering. Ignored if --cluster is not
      provided.
      Default: euclid

  center_row:
    type: boolean?
    inputBinding:
      prefix: "--center"
    doc: |
      Apply mean centering for normalized read counts prior to running
      clustering by row. Ignored when --cluster is not row or both.
      Default: do not centered

  export_pdf_plots:
    type: boolean?
    inputBinding:
      prefix: "--pdf"
    doc: |
      Export plots in PDF.
      Default: false

  output_prefix:
    type: string?
    inputBinding:
      prefix: "--output"
    doc: |
      Output prefix for generated files
      Default: ./diffbind

  threads:
    type: int?
    inputBinding:
      prefix: "--cpus"
    doc: |
      Number of cores/cpus to use.
      Default: 1


outputs:

  pk_vrlp_s_plot_png:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_pk_vrlp_s*.png"
    doc: |
      Peakset overlap rate
      PNG format

  pk_vrlp_s_plot_pdf:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*_pk_vrlp_s*.pdf"
    doc: |
      Peakset overlap rate
      PDF format

  all_pk_scr_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_all_pk_scr_corr.png"
    doc: |
      Datasets correlation (all peaks)
      PNG format

  all_pk_scr_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_all_pk_scr_corr.pdf"
    doc: |
      Datasets correlation (all peaks)
      PDF format

  cns_pk_scr_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_cns_pk_scr_corr.png"
    doc: |
      Datasets correlation (optionally
      recentered consensus peaks)
      PNG format

  cns_pk_scr_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_cns_pk_scr_corr.pdf"
    doc: |
      Datasets correlation (optionally
      recentered consensus peaks)
      PDF format

  rw_rds_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_rw_rds_corr.png"
    doc: |
      Datasets correlation (raw reads)
      PNG format

  rw_rds_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_rw_rds_corr.pdf"
    doc: |
      Datasets correlation (raw reads)
      PDF format

  nr_rds_corr_plot_png:
    type: File?
    outputBinding:
      glob: "*_nr_rds_corr.png"
    doc: |
      Datasets correlation (normalized reads)
      PNG format

  nr_rds_corr_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_nr_rds_corr.pdf"
    doc: |
      Datasets correlation (normalized reads)
      PDF format

  diff_vlcn_plot_png:
    type: File?
    outputBinding:
      glob: "*_diff_vlcn.png"
    doc: |
      Volcano plot for differentially bound sites
      PNG format

  diff_vlcn_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_diff_vlcn.pdf"
    doc: |
      Volcano plot for differentially bound sites
      PDF format

  diff_ma_plot_png:
    type: File?
    outputBinding:
      glob: "*_diff_ma.png"
    doc: |
      MA-plot for differentially bound sites
      PNG format

  diff_ma_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_diff_ma.pdf"
    doc: |
      MA-plot for differentially bound sites
      PDF format

  nr_rds_pca_1_2_plot_png:
    type: File?
    outputBinding:
      glob: "*_nr_rds_pca_1_2.png"
    doc: |
      PCA (1,2) of not filtered normalized counts
      PNG format

  nr_rds_pca_1_2_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_nr_rds_pca_1_2.pdf"
    doc: |
      PCA (1,2) of not filtered normalized counts
      PDF format

  nr_rds_pca_2_3_plot_png:
    type: File?
    outputBinding:
      glob: "*_nr_rds_pca_2_3.png"
    doc: |
      PCA (2,3) of not filtered normalized counts
      PNG format

  nr_rds_pca_2_3_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_nr_rds_pca_2_3.pdf"
    doc: |
      PCA (2,3) of not filtered normalized counts
      PDF format

  nr_rds_mds_html:
    type: File?
    outputBinding:
      glob: "*_nr_rds_mds.html"
    doc: |
      MDS plot of normalized counts.
      HTML format

  pk_prfl_plot_png:
    type: File?
    outputBinding:
      glob: "*_pk_prfl.png"
    doc: |
      Peak profiles
      PNG format

  pk_prfl_plot_pdf:
    type: File?
    outputBinding:
      glob: "*_pk_prfl.pdf"
    doc: |
      Peak profiles
      PDF format

  diff_sts_tsv:
    type: File?
    outputBinding:
      glob: "*_diff_sts.tsv"
    doc: |
      Differentially bound sites. Not filtered.
      TSV format

  nr_rds_gct:
    type: File?
    outputBinding:
      glob: "*_nr_rds.gct"
    doc: |
      Normalized filtered by padj read counts
      GCT format

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["run_diffbind_manual.R"]
stdout: diffbind_manual_stdout.log
stderr: diffbind_manual_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "DiffBind Multi-factor Analysis"
s:name: "DiffBind Multi-factor Analysis"
s:alternateName: "Runs DiffBind multi-factor analysis with manual control over major parameters"


s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/diffbind-multi-factor.cwl
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
  DiffBind Multi-factor Analysis

  Runs DiffBind multi-factor analysis with manual control over major parameters


s:about: |
  usage: run_diffbind_manual.R [-h] --alignments ALIGNMENTS [ALIGNMENTS ...]
                              --peaks PEAKS [PEAKS ...] --aliases ALIASES
                              [ALIASES ...] --metadata METADATA
                              [--scoreby {pvalue,qvalue}] [--score SCORE]
                              [--minrpkm MINRPKM] [--summits SUMMITS]
                              [--minoverlap MINOVERLAP]
                              [--groupby [GROUPBY ...]] --design DESIGN
                              [--contrast CONTRAST] [--base [BASE ...]]
                              [--method {edger,deseq2}]
                              [--norm {auto,rle,tmm,lib}] [--padj PADJ]
                              [--cluster {row,column,both}]
                              [--rowdist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}]
                              [--columndist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}]
                              [--center] [--pdf] [--output OUTPUT]
                              [--cpus CPUS]

  DiffBind Multi-factor Analysis

  options:
    -h, --help            show this help message and exit
    --alignments ALIGNMENTS [ALIGNMENTS ...]
                          Sorted and indexed alignment files in bam format.
    --peaks PEAKS [PEAKS ...]
                          Peak files in the MACS2 xls format. Number and order
                          of the files should correspond to the files provided
                          in --alignments parameter.
    --aliases ALIASES [ALIASES ...]
                          Unique names for datasets provided in --alignments and
                          --peaks parameters, no special characters or spaces
                          are allowed. Number and order of the names should
                          correspond to the values provided in --alignments and
                          --peaks parameters.
    --metadata METADATA   TSV/CSV metadata file to describe datasets provided in
                          --alignments and --peaks parameters. First column
                          should have the name 'sample', all other columns names
                          should be selected from the following list: Tissue,
                          Factor, Condition, Treatment, Replicate. The values
                          from the 'sample' column should correspond to the
                          values provided in --aliases parameter. For a proper
                          --contrast intepretation, values defined in each
                          metadata column should not be used in any of the other
                          columns. All metadata columns are treated as factors
                          (no covariates are supported).
    --scoreby {pvalue,qvalue}
                          Score metrics to exclude low quality peaks based on
                          the threshold provided in the --score parameter.
                          Default: pvalue
    --score SCORE         Filtering threshold to keep only those peaks where the
                          metric selected in --scoreby parameter is less than or
                          equal to the provided value. Default: 0.05
    --minrpkm MINRPKM     Filtering threshold to keep only those peaks where the
                          max RPKM for all datasets is bigger than or equal to
                          the provided value. Default: 1
    --summits SUMMITS     Width in bp to extend peaks around their summits in
                          both directions and replace the original ones. Set it
                          to 100 bp for ATAC-Seq and 200 bp for ChIP-Seq
                          datasets. To skip peaks extension and replacement, set
                          it to negative value. Default: 200 bp (results in 401
                          bp wide peaks)
    --minoverlap MINOVERLAP
                          Filtering threshold to keep only those peaks that are
                          present in at least this many datasets when generating
                          consensus set of peaks used in differential analysis.
                          If this threshold has a value between zero and one,
                          only those peaks will be included that are present in
                          at least this proportion of datasets. When combined
                          with --groupby parameter, --minoverlap threshold is
                          applied per group, then union of the resulted peaks
                          are used in the differential analysis. Default: 2
    --groupby [GROUPBY ...]
                          Column(s) from the metadata table to define datasets
                          grouping. --minoverlap filtering threshold will be
                          applied within each group independently. Union of the
                          resulted peaks from each of the groups will be used in
                          the differential analysis. Default: apply --minoverlap
                          filtering threshold for all datasets jointly
    --design DESIGN       Design formula comprised of the metadata columns
                          names. It should start with ~.
    --contrast CONTRAST   Contrast applied to the analysis results when
                          calculating log2 fold changes. It should be formatted
                          as a mathematical formula of values present in the
                          metadata table. It is a required parameter if --method
                          is set to edger. If not provided and --method is set
                          to deseq2, the last term from the design formula will
                          be used.
    --base [BASE ...]     Base levels for each of the metadata columns. Number
                          and order of the provided values should correspond to
                          the metadata columns. Default: define base levels
                          alphabetically.
    --method {edger,deseq2}
                          Method used in the differential binding analysis.
                          Should be equal to either edger or deseq2. Default:
                          deseq2
    --norm {auto,rle,tmm,lib}
                          Normalization technique applied to the read counts
                          before running differential binding analysis. When set
                          to auto selects rle for deseq2 and tmm for edger.
                          Default: auto
    --padj PADJ           Filtering threshold to report only differentially
                          bound sites with adjusted P-value less than or equal
                          to the provided value. Default: 0.05
    --cluster {row,column,both}
                          Hopach clustering method to be run on normalized read
                          counts. Default: do not run clustering
    --rowdist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}
                          Distance metric for HOPACH row clustering. Ignored if
                          --cluster is not provided. Default: cosangle
    --columndist {cosangle,abscosangle,euclid,abseuclid,cor,abscor}
                          Distance metric for HOPACH column clustering. Ignored
                          if --cluster is not provided. Default: euclid
    --center              Apply mean centering for normalized read counts prior
                          to running clustering by row. Ignored when --cluster
                          is not row or both. Default: do not centered
    --pdf                 Export plots in PDF. Default: false
    --output OUTPUT       Output prefix for generated files
    --cpus CPUS           Number of cores/cpus to use. Default: 1