cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: InlineJavascriptRequirement
- class: ShellCommandRequirement
hints:
- class: DockerRequirement
  dockerPull: robertplayer/scidap-qiime2:stable
inputs:
  samplename:
    type: string
    inputBinding:
      prefix: -s
    doc: |
      Name of sample, should be unique within a project/experiment.
  read1file:
    type: File
    inputBinding:
      prefix: -a
    doc: |
      Read1 data in FASTQ format, received after paired end sequencing.
  read2file:
    type: File
    inputBinding:
      prefix: -b
    doc: |
      Read2 data in FASTQ format, received after paired end sequencing.
  trimLeftF:
    type: int
    inputBinding:
      prefix: -j
    doc: |
      trims the first J bases from the 5' end of each forward sequence
  trimLeftR:
    type: int
    inputBinding:
      prefix: -k
    doc: |
      trims the first K bases from the 5' end of each reverse sequence
  truncLenF:
    type: int
    inputBinding:
      prefix: -m
    doc: |
      Clips the forward read starting M bases from the 5' end (before trimming). If base quality is OK for entire read, value should be set to the expected number of Illumina cycles for R1.
  truncLenR:
    type: int
    inputBinding:
      prefix: -n
    doc: |
      Clips the reverse read starting N bases from the 5' end (before trimming).  If base quality is OK for entire read, value should be set to the expected number of Illumina cycles for R2.
  threads:
    type: int
    inputBinding:
      prefix: -t
    doc: |
      Number of threads for parallel processing.
outputs:
  overview:
    type: File
    outputBinding:
      glob: overview.md
    doc: |
      overview of inputs
  fastq_summary:
    type: File?
    outputBinding:
      glob: demux.qzv
    doc: |
      summary of input read data
  alpha_rarefaction:
    type: File?
    outputBinding:
      glob: alpha-rarefaction.qzv
    doc: |
      plot of OTU rarefaction
  taxa_bar_plots:
    type: File?
    outputBinding:
      glob: taxa-bar-plots.qzv
    doc: |
      bar plot for exploring the taxonomic composition of the sample
  log_file_stdout:
    type: stdout
  log_file_stderr:
    type: stderr
baseCommand:
- run_qiime2_sample_pe.sh
stdout: qiime2_sample_pe-stdout.log
stderr: qiime2_sample_pe-stderr.log
doc: |
  Shell wrapper for import and quantitation of one paired-end 16S sequencing sample using DADA2. Alpha rarefaction and taxonomic classification plots are also output.
  Taxonomy classification is performed using a Naive Bayes classifier trained on the Greengenes2 database "gg_2022_10_backbone_full_length.nb.qza".
  Generally, this workflow follows the "moving-pictures" turorial: https://docs.qiime2.org/2023.5/tutorials/moving-pictures/

      Output files:
      - overview.md, list of inputs
      - demux.qzv, summary visualizations of imported data
      - alpha-rarefaction.qzv, plot of OTU rarefaction
      - taxa-bar-plots.qzv, relative frequency of taxomonies barplot


  PARAMS:
      SECTION 1: general
      -h  help   show this message
      -t  INT    number of threads
      -s  STR    sample name
      -a  FILE   path to read1 fastq file
      -b  FILE   path to read2 fastq file
      -j  J      trims the first J bases from the 5' end of each forward read
      -k  K      trims the first K bases from the 5' end of each reverse read
      -m  M      clips the forward read starting M bases from the 5' end (before trimming)
      -n  N      clips the reverse read starting N bases from the 5' end (before trimming)

  ____________________________________________________________________________________________________
  References:
      Bolyen E, Rideout JR, Dillon MR, Bokulich NA, Abnet CC, Al-Ghalith GA, Alexander H, Alm EJ, Arumugam M, Asnicar F, Bai Y, Bisanz JE, Bittinger K, Brejnrod A, Brislawn CJ, Brown CT, Callahan BJ, Caraballo-Rodríguez AM, Chase J, Cope EK, Da Silva R, Diener C, Dorrestein PC, Douglas GM, Durall DM, Duvallet C, Edwardson CF, Ernst M, Estaki M, Fouquier J, Gauglitz JM, Gibbons SM, Gibson DL, Gonzalez A, Gorlick K, Guo J, Hillmann B, Holmes S, Holste H, Huttenhower C, Huttley GA, Janssen S, Jarmusch AK, Jiang L, Kaehler BD, Kang KB, Keefe CR, Keim P, Kelley ST, Knights D, Koester I, Kosciolek T, Kreps J, Langille MGI, Lee J, Ley R, Liu YX, Loftfield E, Lozupone C, Maher M, Marotz C, Martin BD, McDonald D, McIver LJ, Melnik AV, Metcalf JL, Morgan SC, Morton JT, Naimey AT, Navas-Molina JA, Nothias LF, Orchanian SB, Pearson T, Peoples SL, Petras D, Preuss ML, Pruesse E, Rasmussen LB, Rivers A, Robeson MS, Rosenthal P, Segata N, Shaffer M, Shiffer A, Sinha R, Song SJ, Spear JR, Swafford AD, Thompson LR, Torres PJ, Trinh P, Tripathi A, Turnbaugh PJ, Ul-Hasan S, van der Hooft JJJ, Vargas F, Vázquez-Baeza Y, Vogtmann E, von Hippel M, Walters W, Wan Y, Wang M, Warren J, Weber KC, Williamson CHD, Willis AD, Xu ZZ, Zaneveld JR, Zhang Y, Zhu Q, Knight R, and Caporaso JG. 2019. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nature Biotechnology 37: 852–857. https://doi.org/10.1038/s41587-019-0209-9
label: vc-germline
