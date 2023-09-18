cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement


hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-qiime2:stable


inputs:

  metadata_file:
    type: File
    inputBinding:
      prefix: "-r"
    doc: |
      sample metadata file, sample-id [col1] should be identical to sample names being aggregated

  pcoa_label:
    type: string?
    inputBinding:
      prefix: "-c"
    doc: |
      Must be identical to one of the headers of the metadata file. Values under this metadata header must be INT. Required for PCoA analysis.

  sampling_depth:
    type: int?
    inputBinding:
      prefix: "-d"
    doc: |
      Required for differential abundance analyses (along with group and taxonomic level). This step will subsample the counts in each sample without replacement so that each sample in the resulting table has a total count of INT. If the total count for any sample(s) are smaller than this value, those samples will be dropped from further analysis. It's recommend making your choice by reviewing the rarefaction plot. Choose a value that is as high as possible (so you retain more sequences per sample) while excluding as few samples as possible.

  diff_group:
    type: string?
    inputBinding:
      prefix: "-g"
    doc: |
      Required for differential abundance analyses (along with sampling depth and taxonomic level). Group/experimental condition column name from sample metadata file. Must be identical to one of the headers of the sample-metadata file. The corresponding column should only have two groups/conditions.

  taxonomic_level:
    type: string?
    inputBinding:
      prefix: "-l"
    doc: |
      Required for differential abundance analyses (along with sampling depth and group). Collapses the OTU table at the taxonomic level of interest for differential abundance analysis with ANCOM. Default: Genus

  fastq_r1_array:
    type: File[]
    inputBinding:
      prefix: "-a"
      itemSeparator: ","
    doc: |
      Array of forward read data in FASTQ format from SciDAP upstream qiime2-sample-pe workflow.

  fastq_r2_array:
    type: File[]
    inputBinding:
      prefix: "-b"
      itemSeparator: ","
    doc: |
      Array of reverse read data in FASTQ format from SciDAP upstream qiime2-sample-pe workflow.

  trimLeftF:
    type: int
    inputBinding:
      prefix: "-j"
      itemSeparator: ","
    doc: |
      trims the first J bases from the 5' end of each forward sequence

  trimLeftR:
    type: int
    inputBinding:
      prefix: "-k"
      itemSeparator: ","
    doc: |
      trims the first K bases from the 5' end of each reverse sequence

  truncLenF:
    type: int[]
    inputBinding:
      prefix: "-m"
      itemSeparator: ","
    doc: |
      Clips the forward read starting M bases from the 5' end (before trimming). If base quality is OK for entire read, value should be set to the expected number of Illumina cycles for R1.

  truncLenR:
    type: int[]
    inputBinding:
      prefix: "-n"
      itemSeparator: ","
    doc: |
      Clips the reverse read starting N bases from the 5' end (before trimming).  If base quality is OK for entire read, value should be set to the expected number of Illumina cycles for R2.

  threads:
    type: int
    inputBinding:
      prefix: "-t"
    doc: |
      Number of threads for parallel processing.

outputs:

  overview:
    type: File
    outputBinding:
      glob: overview.md
    doc: |
      overview of inputs

  fastq_summary_file:
    type: File
    outputBinding:
      glob: demux.qzv
    doc: |
      summary of input read data

  denoising_stats:
    type: File
    outputBinding:
      glob: stats.qzv
    doc: |
      Detect and correct (where possible) Illumina amplicon sequence data. This quality control process will additionally filter any phiX reads (commonly present in marker gene Illumina sequence data) that are identified in the sequencing data, and will filter chimeric sequences.

  alpha_rarefaction:
    type: File
    outputBinding:
      glob: alpha-rarefaction.qzv
    doc: |
      plot of OTU rarefaction

  taxa_bar_plots:
    type: File
    outputBinding:
      glob: taxa-bar-plots.qzv
    doc: |
      bar plot for exploring the taxonomic composition of the sample

  pcoa_uwunifrac:
    type: File?
    outputBinding:
      glob: core-metrics-results/pcoa-unweighted-unifrac-emperor.qzv
    doc: |
      PCoA using unweighted unifrac method

  pcoa_braycurtis:
    type: File?
    outputBinding:
      glob: core-metrics-results/pcoa-bray-curtis-emperor.qzv
    doc: |
      PCoA using bray curtis method

  heatmap:
    type: File?
    outputBinding:
      glob: heatmap.qzv
    doc: |
      Gneiss differential abundance analysis hierarchical clustering and heatmap

  ancom_family:
    type: File?
    outputBinding:
      glob: ancom-Family.qzv
    doc: |
      ANCOM family taxonomic level differential abundance analysis with volcano plot

  ancom_genus:
    type: File?
    outputBinding:
      glob: ancom-Genus.qzv
    doc: |
      ANCOM genus taxonomic level differential abundance analysis with volcano plot

  ancom_species:
    type: File?
    outputBinding:
      glob: ancom-Species.qzv
    doc: |
      ANCOM species taxonomic level differential abundance analysis with volcano plot

  log_file_stdout:
    type: stdout

  log_file_stderr:
    type: stderr


baseCommand: ["run_qiime2_aggregate.sh"]
stdout: qiime2_aggregate-stdout.log
stderr: qiime2_aggregate-stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "vc-germline"
s:downloadUrl: https://github.com/datirium/workflows/blob/master/tools/qiime2-aggregate.cwl
s:codeRepository: https://github.com/datirium/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:creator:
- class: s:Organization
  s:legalName: "Datirium LLC"
  s:location:
  - class: s:PostalAddress
    s:addressCountry: "USA"
    s:addressLocality: "Cincinnati"
    s:addressRegion: "OH"
    s:postalCode: ""
    s:streetAddress: ""
    s:telephone: ""
  s:logo: "https://avatars.githubusercontent.com/u/33202955?s=200&v=4"
  s:department:
  - class: s:Organization
    s:legalName: "Datirium LLC"
    s:department:
    - class: s:Organization
      s:legalName: "Bioinformatics"
      s:member:
      - class: s:Person
        s:name: Robert Player
        s:email: mailto:support@datirium.com
        s:sameAs:
        - id: https://orcid.org/0000-0001-5872-259X


doc: |
  Shell wrapper for import and quantitation of 2+ paired-end 16S sequencing samples that (each) have already been processed in the SciDAP platform (using the qiime2-sample-pe.cwl workflow).
  Alpha rarefaction and taxonomic classification plots are also output for the aggregated samples.
  Taxonomy classification is performed using a Naive Bayes classifier trained on the Greengenes2 database "gg_2022_10_backbone_full_length.nb.qza".
  Generally, this workflow follows the "moving-pictures" turorial: https://docs.qiime2.org/2023.5/tutorials/moving-pictures/
  If "Metadata header name for PCoA axis label" is provided, principle coordinates analysis (PCoA) will be performed using the unweighted unifrac and bray curtis methods. 3D plots are produced with PCo1, PCo2, and the provided axis label on the x, y, and z axes.
  If the sampling depth and metadata header for differential analysis are provided, differential abundance analysis will be performed using Gneiss and ANCOM methods at the family, genus, and species taxonomic levels. A unsupervised hierarchical clustering heatmap (Gneiss) and volcano plot (ANCOM) are produced at the taxonomic level between the specified group.

      Primary output files:
      - overview.md, list of inputs
      - demux.qzv, summary visualizations of imported data
      - alpha-rarefaction.qzv, plot of OTU rarefaction
      - taxa-bar-plots.qzv, relative frequency of taxomonies barplot
      - table.qza, table containing how many sequences are associated with each sample and with each feature (OTU)
      Optional output files:
      - pcoa-unweighted-unifrac-emperor.qzv, PCoA using unweighted unifrac method
      - pcoa-bray-curtis-emperor.qzv, PCoA using bray curtis method
      - heatmap.qzv, output from gneiss differential abundance analysis using unsupervised correlation-clustering method (this will define the partitions of microbes that commonly co-occur with each other using Ward hierarchical clustering)
      - ancom-\$LEVEL.qzv, output from ANCOM differential abundance analysis at family, genus, and species taxonomic levels (includes volcano plot)

  PARAMS:
      SECTION 1: general
      -h  help         show this message
      -t  INT          number of threads
      -r  FILE         sample metadata file (sample-id [col1] should be identical to sample names being aggregated)
      -a  FILE ARRAY   array (csv) of paths to read1 fastq files
      -b  FILE ARRAY   array (csv) of paths to read2 fastq files
      -j  INT          trims the first J bases from the 5' end of each forward read
      -k  INT          trims the first K bases from the 5' end of each reverse read
      -m  INT          clips the forward read starting M bases from the 5' end (before trimming)
      -n  INT          clips the reverse read starting N bases from the 5' end (before trimming)
      -c  STR          custom axis label
  Must be identical to one of the headers of the sample-metadata file. The corresponding column may only contain INT data.
      -d  INT          rarefaction sampling depth (required for differential abundance execution)
  This step will subsample the counts in each sample without replacement so that each sample in the resulting table has a total count of INT. If the total count for any sample(s) are smaller than this value, those samples will be dropped from the diversity analysis. It's recommend making your choice by reviewing the rarefaction plot. Choose a value that is as high as possible (so you retain more sequences per sample) while excluding as few samples as possible.
      -g  STR          group or experimental condition column name from sample metadata file (required for differential abundance execution)
  Must be identical to one of the headers of the sample-metadata file. The corresponding column should only have two groups/conditions.

  NOTES:
    Example sample metadata file (-r):
    Spacing should be single-tab separated (below things are lined up for clarity)

  sample-id       sample-name     condition       time
  SRR25508255     sample1 c1      1
  SRR25508256     sample2 c1      1
  SRR25508257     sample3 c2      2
  SRR25508258     sample4 c2      2
  
  ____________________________________________________________________________________________________
  References:
      Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumugam M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodríguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen LB, Rivers A, Robeson MS, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, Ul-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Vogtmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, and Caporaso JG. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857. https://doi.org/10.1038/s41587-019-0209-9
