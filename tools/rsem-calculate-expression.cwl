cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
        if (Array.isArray(inputs.upstream_read_file) && inputs.upstream_read_file.length | 0){
          return inputs.upstream_read_file[0].basename.split('.').slice(0,-1).join('.');
        } else
          if (inputs.upstream_read_file != null){
            return inputs.upstream_read_file.basename.split('.').slice(0,-1).join('.');
          } else
            if (Array.isArray(inputs.downstream_read_file) && inputs.downstream_read_file.length | 0){
              return inputs.downstream_read_file[0].basename.split('.').slice(0,-1).join('.');
            } else
              if (inputs.downstream_read_file != null){
                return inputs.downstream_read_file.basename.split('.').slice(0,-1).join('.');
              } else
                if (inputs.input_aligned != null){
                  return inputs.input_aligned.basename.split('.').slice(0,-1).join('.');
                } else
                  {
                    return null;
                  }
    };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/rsem:v1.3.0


inputs:
  upstream_read_file:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      position: 100
    doc: |
      Comma-separated list of files containing single-end reads
      or upstream reads for paired-end data.
      By default, these files are assumed to be in FASTQ format.
      If the --no-qualities option is specified,
      then FASTA format is expected.

  downstream_read_file:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      position: 101
      itemSeparator: ","
    doc: |
      Comma-separated list of files containing downstream reads
      which are paired with the upstream reads.
      By default, these files are assumed to be in FASTQ format.
      If the --no-qualities option is specified,
      then FASTA format is expected.

  input_aligned:
    type:
      - "null"
      - File
    inputBinding:
      position: 102
      prefix: '--alignments'
    doc: |
      SAM/BAM/CRAM formatted input file.
      RSEM requires all alignments of the same read group together.
      For paired-end reads, RSEM also requires the two mates of any alignment be adjacent.
      In addition, RSEM does not allow the SEQ and QUAL fields to be empty.
      See Description section for how to make input file obey RSEM's requirements.

  indices_folder:
    type: Directory
    doc: |
      Path to the folder where all rsem reference files are saved

  output_filename:
    type:
      - "null"
      - string
    inputBinding:
      position: 104
      valueFrom: |
        ${
            if (self == null){
              return default_output_filename();
            } else {
              return self;
            }
        }
    default: null
    doc: |
      The name of the sample analyzed.
      All output files are prefixed by this name
      (e.g., sample_name.genes.results)

#BASIC OPTIONS

  paired_end:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 99
      prefix: '--paired-end'
    doc: |
      Input reads are paired-end reads
      If ommited and both upstream_read_file and downstream_read_file are set - we'll add this prefix
      automatically from the arguments oprator.
      In a case of using input 'input', user should define himself either it's paired-end or single-end data
      (Default: off)

  no_qualities:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 97
      prefix: '--no-qualities'
    doc: |
      Input reads do not contain quality scores
      (Default: off)

  strandedness:
    type:
      - "null"
      - type: enum
        name: "strandedness"
        symbols: ["none","forward","reverse"]
    inputBinding:
      position: 96
      prefix: '--strandedness'
    doc: |
      This option defines the strandedness of the RNA-Seq reads.
      It recognizes three values: 'none', 'forward', and 'reverse'.
      'none' refers to non-strand-specific protocols.
      'forward' means all (upstream) reads are derived from the forward strand.
      'reverse' means all (upstream) reads are derived from the reverse strand.
      If 'forward'/'reverse' is set, the '--norc'/'--nofw' Bowtie/Bowtie 2 option
      will also be enabled to avoid aligning reads to the opposite strand.
      For Illumina TruSeq Stranded protocols, please use 'reverse'.
      (Default: 'none')

  threads:
    type:
      - "null"
      - int
    inputBinding:
      position: 95
      prefix: '--num-threads'
    doc: |
      Number of threads to use. Both Bowtie/Bowtie2, expression estimation and 'samtools sort' will use this many threads.
      (Default: 1)

  fai_file:
    type:
      - "null"
      - File
    inputBinding:
      position: 94
      prefix: '--fai'
    doc: |
      If the header section of input alignment file does not contain reference sequence information,
      this option should be turned on. <file| is a FAI format file containing each reference sequence's name and length.
      Please refer to the SAM official website for the details of FAI format.
      (Default: off)

  bowtie2:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 93
      prefix: '--bowtie2'
    doc: |
      Use Bowtie 2 instead of Bowtie to align reads. Since currently RSEM does not handle indel,
      local and discordant alignments, the Bowtie2 parameters are set in a way to avoid those alignments.
      In particular, we use options '--sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1' by default.
      The last parameter of '--score-min', '-0.1', is the negative of maximum mismatch rate.
      This rate can be set by option '--bowtie2-mismatch-rate'. If reads are paired-end, we additionally use options
      '--no-mixed' and '--no-discordant'.
      (Default: off)

  star:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 92
      prefix: '--star'
    doc: |
      Use STAR to align reads. Alignment parameters are from ENCODE3's STAR-RSEM pipeline.
      To save computational time and memory resources, STAR's Output BAM file is unsorted.
      It is stored in RSEM's temporary directory with name as 'sample_name.bam'.
      Each STAR job will have its own private copy of the genome in memory. (Default: off)

  append_names:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 91
      prefix: '--append-names'
    doc: |
      If gene_name/transcript_name is available, append it to the end of gene_id/transcript_id (separated by '_')
      in files 'sample_name.isoforms.results' and 'sample_name.genes.results'.
      (Default: off)

  seed:
    type:
      - "null"
      - int
    inputBinding:
      position: 90
      prefix: '--seed'
    doc: |
      Set the seed for the random number generators used in calculating posterior mean estimates and credibility intervals.
      The seed must be a non-negative 32 bit integer.
      (Default: off)

  single_cell_prior:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 89
      prefix: '--single-cell-prior'
    doc: |
      By default, RSEM uses Dirichlet(1) as the prior to calculate posterior mean estimates and credibility intervals.
      However, much less genes are expressed in single cell RNA-Seq data. Thus, if you want to compute posterior
      mean estimates and/or credibility intervals and you have single-cell RNA-Seq data, you are recommended
      to turn on this option. Then RSEM will use Dirichlet(0.1) as the prior which encourage the sparsity
      of the expression levels.
      (Default: off)

  calc_pme:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 88
      prefix: '--calc-pme'
    doc: |
      Run RSEM's collapsed Gibbs sampler to calculate posterior mean estimates.
      (Default: off)

  calc_ci:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 87
      prefix: '--calc-ci'
    doc: |
      Calculate 95% credibility intervals and posterior mean estimates.
      The credibility level can be changed by setting '--ci-credibility-level'.
      (Default: off)

  quiet:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 86
      prefix: '--quiet'
    doc: |
      Suppress the output of logging information.
      (Default: off)

# OUTPUT OPTIONS

  sort_bam_by_read_name:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 85
      prefix: '--sort-bam-by-read-name'
    doc: |
      Sort BAM file aligned under transcript coordidate by read name.
      Setting this option on will produce deterministic maximum likelihood estimations from independent runs.
      Note that sorting will take long time and lots of memory.
      (Default: off)

  no_bam_output:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 84
      prefix: '--no-bam-output'
    doc: |
      Do not output any BAM file.
      (Default: off)

  sampling_for_bam:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 83
      prefix: '--sampling-for-bam'
    doc: |
      When RSEM generates a BAM file, instead of outputting all alignments a read has with their posterior probabilities,
      one alignment is sampled according to the posterior probabilities. The sampling procedure includes the alignment
      to the "noise" transcript, which does not appear in the BAM file. Only the sampled alignment has a weight of 1.
      All other alignments have weight 0. If the "noise" transcript is sampled, all alignments appeared in the BAM file
      should have weight 0.
      (Default: off)

  output_genome_bam:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 82
      prefix: '--output-genome-bam'
    doc: |
      Generate a BAM file, 'sample_name.genome.bam', with alignments mapped to genomic coordinates and annotated with
      their posterior probabilities. In addition, RSEM will call samtools (included in RSEM package) to sort and index
      the bam file. 'sample_name.genome.sorted.bam' and 'sample_name.genome.sorted.bam.bai' will be generated.
      (Default: off)

  sort_bam_by_coordinate:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 81
      prefix: '--sort-bam-by-coordinate'
    doc: |
      Sort RSEM generated transcript and genome BAM files by coordinates and build associated indices.
      (Default: off)

  sort_bam_memory_per_thread:
    type:
      - "null"
      - string
    inputBinding:
      position: 80
      prefix: '--sort-bam-memory-per-thread'
    doc: |
      Set the maximum memory per thread that can be used by 'samtools sort'. <string| represents the memory and accepts
      suffices 'K/M/G'. RSEM will pass <string| to the '-m' option of 'samtools sort'.
      Note that the default used here is different from the default used by samtools.
      (Default: 1G)

# ALIGNER OPTIONS

  seed_length:
    type:
      - "null"
      - int
    inputBinding:
      position: 79
      prefix: '--seed-length'
    doc: |
      Seed length used by the read aligner. Providing the correct value is important for RSEM.
      If RSEM runs Bowtie, it uses this value for Bowtie's seed length parameter.
      Any read with its or at least one of its mates' (for paired-end reads) length less than this value will be ignored.
      If the references are not added poly(A) tails, the minimum allowed value is 5, otherwise, the minimum allowed value is 25.
      Note that this script will only check if the value |= 5 and give a warning message if the value < 25 but |= 5.
      (Default: 25)

  phred33_quals:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 78
      prefix: '--phred33-quals'
    doc: |
      Input quality scores are encoded as Phred+33.
      (Default: on)

  phred64_quals:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 77
      prefix: '--phred64-quals'
    doc: |
      Input quality scores are encoded as Phred+64 (default for GA Pipeline ver. |= 1.3).
      (Default: off)

  solexa_quals:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 76
      prefix: '--solexa-quals'
    doc: |
      Input quality scores are solexa encoded (from GA Pipeline ver. < 1.3).
      (Default: off)

  bowtie_n:
    type:
      - "null"
      - int
    inputBinding:
      position: 74
      prefix: '--bowtie-n'
    doc: |
      (Bowtie parameter) max # of mismatches in the seed.
      (Range: 0-3, Default: 2))

  bowtie_e:
    type:
      - "null"
      - int
    inputBinding:
      position: 73
      prefix: '--bowtie-e'
    doc: |
      (Bowtie parameter) max sum of mismatch quality scores across the alignment.
      (Default: 99999999)

  bowtie_m:
    type:
      - "null"
      - int
    inputBinding:
      position: 72
      prefix: '--bowtie-m'
    doc: |
      (Bowtie parameter) suppress all alignments for a read if | <int| valid alignments exist.
      (Default: 200)

  bowtie_chunkmbs:
    type:
      - "null"
      - int
    inputBinding:
      position: 71
      prefix: '--bowtie-chunkmbs'
    doc: |
      (Bowtie parameter) memory allocated for best first alignment calculation
      (Default: 0 - use Bowtie's default)

  bowtie2_mismatch_rate:
    type:
      - "null"
      - double
    inputBinding:
      position: 69
      prefix: '--bowtie2-mismatch-rate'
    doc: |
      (Bowtie 2 parameter) The maximum mismatch rate allowed.
      (Default: 0.1)

  bowtie2_k:
    type:
      - "null"
      - int
    inputBinding:
      position: 68
      prefix: '--bowtie2-k'
    doc: |
      (Bowtie 2 parameter) Find up to <int| alignments per read.
      (Default: 200)

  bowtie2_sensitivity_level:
    type:
      - "null"
      - type: enum
        name: "bowtie2_sensitivity"
        symbols: ["very_fast","fast","sensitive","very_sensitive"]
    inputBinding:
      position: 67
      prefix: '--bowtie2-sensitivity-level'
    doc: |
      (Bowtie 2 parameter) Set Bowtie 2's preset options in --end-to-end mode.
      This option controls how hard Bowtie 2 tries to find alignments.
      <string| must be one of "very_fast", "fast", "sensitive" and "very_sensitive".
      The four candidates correspond to Bowtie 2's "--very-fast", "--fast", "--sensitive" and "--very-sensitive" options.
      (Default: "sensitive" - use Bowtie 2's default)

  star_gzipped_read_file:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 65
      prefix: '--star-gzipped-read-file'
    doc: |
      (STAR parameter) Input read file(s) is compressed by gzip.
      (Default: off)

  star_bzipped_read_file:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 64
      prefix: '--star-bzipped-read-file'
    doc: |
      (STAR parameter) Input read file(s) is compressed by bzip2.
      (Default: off)

  star_output_genome_bam:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 63
      prefix: '--star-output-genome-bam'
    doc: |
      (STAR parameter) Save the BAM file from STAR alignment under genomic coordinate to 'sample_name.STAR.genome.bam'.
      This file is NOT sorted by genomic coordinate. In this file, according to STAR's manual,
      'paired ends of an alignment are always adjacent, and multiple alignments of a read are adjacent as well'.
      (Default: off)

# ADVANCED OPTIONS

  tag:
    type:
      - "null"
      - string
    inputBinding:
      position: 62
      prefix: '--tag'
    doc: |
      The name of the optional field used in the SAM input for identifying a read with too many valid alignments.
      The field should have the format <tagName|:i:<value|, where a <value| bigger than 0 indicates a read with too many alignments.
      (Default: "")

  fragment_length_min:
    type:
      - "null"
      - int
    inputBinding:
      position: 61
      prefix: '--fragment-length-min'
    doc: |
      Minimum read/insert length allowed. This is also the value for the Bowtie/Bowtie2 -I option.
      (Default: 1)

  fragment_length_max:
    type:
      - "null"
      - int
    inputBinding:
      position: 60
      prefix: '--fragment-length-max'
    doc: |
      Maximum read/insert length allowed. This is also the value for the Bowtie/Bowtie 2 -X option.
      (Default: 1000)

  fragment_length_mean:
    type:
      - "null"
      - double
    inputBinding:
      position: 59
      prefix: '--fragment-length-mean'
    doc: |
      (single-end data only) The mean of the fragment length distribution, which is assumed to be a Gaussian.
      (Default: -1, which disables use of the fragment length distribution)

  fragment_length_sd:
    type:
      - "null"
      - double
    inputBinding:
      position: 58
      prefix: '--fragment-length-sd'
    doc: |
      (single-end data only) The standard deviation of the fragment length distribution, which is assumed to be a Gaussian.
      (Default: 0, which assumes that all fragments are of the same length, given by the rounded value of --fragment-length-mean)

  estimate_rspd:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 57
      prefix: '--estimate-rspd'
    doc: |
      Set this option if you want to estimate the read start position distribution (RSPD) from data. Otherwise, RSEM will use a uniform RSPD.
      (Default: off)

  num_rspd_bins:
    type:
      - "null"
      - int
    inputBinding:
      position: 56
      prefix: '--num-rspd-bins'
    doc: |
      Number of bins in the RSPD. Only relevant when '--estimate-rspd' is specified. Use of the default setting is recommended.
      (Default: 20)

  gibbs_burnin:
    type:
      - "null"
      - int
    inputBinding:
      position: 55
      prefix: '--gibbs-burnin'
    doc: |
      The number of burn-in rounds for RSEM's Gibbs sampler. Each round passes over the entire data set once.
      If RSEM can use multiple threads, multiple Gibbs samplers will start at the same time and all samplers share the same burn-in number.
      (Default: 200)

  gibbs_number_of_samples:
    type:
      - "null"
      - int
    inputBinding:
      position: 54
      prefix: '--gibbs-number-of-samples'
    doc: |
      TThe total number of count vectors RSEM will collect from its Gibbs samplers.
      (Default: 1000)

  gibbs_sampling_gap:
    type:
      - "null"
      - int
    inputBinding:
      position: 53
      prefix: '--gibbs-sampling-gap'
    doc: |
      The number of rounds between two succinct count vectors RSEM collects.
      If the count vector after round N is collected, the count vector after round N + <int| will also be collected.
      (Default: 1)

  ci_credibility_level:
    type:
      - "null"
      - double
    inputBinding:
      position: 52
      prefix: '--ci-credibility-level'
    doc: |
      The credibility level for credibility intervals.
      (Default: 0.95)

  ci_memory:
    type:
      - "null"
      - int
    inputBinding:
      position: 51
      prefix: '--ci-memory'
    doc: |
      Maximum size (in memory, MB) of the auxiliary buffer used for computing credibility intervals (CI).
      (Default: 1024)

  ci_number_of_samples_per_count_vector:
    type:
      - "null"
      - int
    inputBinding:
      position: 50
      prefix: '--ci-number-of-samples-per-count-vector'
    doc: |
      The number of read generating probability vectors sampled per sampled count vector.
      The crebility intervals are calculated by first sampling P(C | D) and then sampling P(Theta | C)
      for each sampled count vector. This option controls how many Theta vectors are sampled per sampled count vector.
      (Default: 50)

  keep_intermediate_files:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 49
      prefix: '--keep-intermediate-files'
    doc: |
      Keep temporary files generated by RSEM. RSEM creates a temporary directory, 'sample_name.temp', into which
      it puts all intermediate output files. If this directory already exists, RSEM overwrites all files generated
      by previous RSEM runs inside of it. By default, after RSEM finishes, the temporary directory is deleted.
      Set this option to prevent the deletion of this directory and the intermediate files inside of it.
      (Default: off)

  temporary_folder:
    type:
      - "null"
      - string
    inputBinding:
      position: 48
      prefix: '--temporary-folder'
    doc: |
      Set where to put the temporary files generated by RSEM.
      If the folder specified does not exist, RSEM will try to create it.
      (Default: sample_name.temp)

  time:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 47
      prefix: '--time'
    doc: |
      Output time consumed by each step of RSEM to 'sample_name.time'.
      (Default: off)

# PRIOR-ENHANCED RSEM OPTIONS

  run_prsem:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 46
      prefix: '--run-pRSEM'
    doc: |
      Running prior-enhanced RSEM (pRSEM). Prior parameters, i.e. isoform's initial pseudo-count for RSEM's Gibbs sampling,
      will be learned from input RNA-seq data and an external data set. When pRSEM needs and only needs ChIP-seq peak
      information to partition isoforms (e.g. in pRSEM's default partition model), either ChIP-seq peak file (with the
      '--chipseq-peak-file' option) or ChIP-seq FASTQ files for target and input and the path for Bowtie executables
      are required (with the '--chipseq-target-read-files <string|', '--chipseq-control-read-files <string|', and
      '--bowtie-path <path| options), otherwise, ChIP-seq FASTQ files for target and control and the path to
      Bowtie executables are required.
      (Default: off)

  chipseq_peak_file:
    type:
      - "null"
      - string
    inputBinding:
      position: 45
      prefix: '--chipseq-peak-file'
    doc: |
      Full path to a ChIP-seq peak file in ENCODE's narrowPeak, i.e. BED6+4, format.
      This file is used when running prior-enhanced RSEM in the default two-partition model.
      It partitions isoforms by whether they have ChIP-seq overlapping with their transcription start site region or not.
      Each partition will have its own prior parameter learned from a training set.
      This file can be either gzipped or ungzipped.
      (Default: "")

  chipseq_target_read_files:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      position: 44
      itemSeparator: ","
      prefix: '--chipseq-target-read-files'
    doc: |
      Comma-separated full path of FASTQ read file(s) for ChIP-seq target.
      This option is used when running prior-enhanced RSEM. It provides information to calculate ChIP-seq peaks and signals.
      The file(s) can be either ungzipped or gzipped with a suffix '.gz' or '.gzip'. The options '--bowtie-path <path|' and
      '--chipseq-control-read-files <string|' must be defined when this option is specified.
      (Default: "")

  chipseq_control_read_files:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      position: 43
      itemSeparator: ","
      prefix: '--chipseq-control-read-files'
    doc: |
      Comma-separated full path of FASTQ read file(s) for ChIP-seq conrol. This option is used when running prior-enhanced RSEM.
      It provides information to call ChIP-seq peaks. The file(s) can be either ungzipped or gzipped with a suffix
      '.gz' or '.gzip'. The options '--bowtie-path <path|' and '--chipseq-target-read-files <string|'
      must be defined when this option is specified.
      (Default: "")

  chipseq_read_files_multi_targets:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      position: 42
      itemSeparator: ","
      prefix: '--chipseq-read-files-multi-targets'
    doc: |
      Comma-separated full path of FASTQ read files for multiple ChIP-seq targets.
      This option is used when running prior-enhanced RSEM, where prior is learned from multiple
      complementary data sets. It provides information to calculate ChIP-seq signals.
      All files can be either ungzipped or gzipped with a suffix '.gz' or '.gzip'.
      When this option is specified, the option '--bowtie-path <path|' must be defined and the option
      '--partition-model <string|' will be set to 'cmb_lgt' automatically.
      (Default: "")

  chipseq_bed_files_multi_targets:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      position: 41
      itemSeparator: ","
      prefix: '--chipseq-bed-files-multi-targets'
    doc: |
      Comma-separated full path of BED files for multiple ChIP-seq targets.
      This option is used when running prior-enhanced RSEM, where prior is learned
      from multiple complementary data sets. It provides information of ChIP-seq signals
      and must have at least the first six BED columns. All files can be either
      ungzipped or gzipped with a suffix '.gz' or '.gzip'. When this option is specified,
      the option '--partition-model <string|' will be set to 'cmb_lgt' automatically.
      (Default: "")

  cap_stacked_chipseq_reads:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 40
      prefix: '--cap-stacked-chipseq-reads'
    doc: |
      Keep a maximum number of ChIP-seq reads that aligned to the same genomic interval.
      This option is used when running prior-enhanced RSEM, where prior is learned from multiple
      complementary data sets. This option is only in use when either '--chipseq-read-files-multi-targets
      <string|' or '--chipseq-bed-files-multi-targets <string|' is specified.
      (Default: off)

  n_max_stacked_chipseq_reads:
    type:
      - "null"
      - int
    inputBinding:
      position: 39
      prefix: '--n-max-stacked-chipseq-reads'
    doc: |
      The maximum number of stacked ChIP-seq reads to keep. This option is used when running prior-enhanced RSEM,
      where prior is learned from multiple complementary data sets. This option is only in use when the option
      '--cap-stacked-chipseq-reads' is set.
      (Default: 5)

  partition_model:
    type:
      - "null"
      - type: enum
        name: "partition_model"
        symbols: ["pk","pk_lgtnopk","lm3","lm4","lm5","lm6","nopk_lm2pk","nopk_lm3pk","nopk_lm4pk","nopk_lm5pk","pk_lm2nopk","pk_lm3nopk","pk_lm4nopk","pk_lm5nopk","cmb_lgt"]
    inputBinding:
      position: 38
      prefix: '--partition-model'
    doc: |
      A keyword to specify the partition model used by prior-enhanced RSEM. It must be one of the following keywords:
      - pk
          Partitioned by whether an isoform has a ChIP-seq peak overlapping with its transcription start site (TSS) region.
          The TSS region is defined as [TSS-500bp, TSS+500bp]. For simplicity, we refer this type of peak as 'TSS peak'
          when explaining other keywords.
      - pk_lgtnopk
          First partitioned by TSS peak. Then, for isoforms in the 'no TSS peak' set, a logistic model is employed to
          further classify them into two partitions.
      - lm3, lm4, lm5, or lm6
          Based on their ChIP-seq signals, isoforms are classified into 3, 4, 5, or 6 partitions by a linear regression model.
      - nopk_lm2pk, nopk_lm3pk, nopk_lm4pk, or nopk_lm5pk
          First partitioned by TSS peak. Then, for isoforms in the 'with TSS peak' set, a linear regression model is
          employed to further classify them into 2, 3, 4, or 5 partitions.
      - pk_lm2nopk, pk_lm3nopk, pk_lm4nopk, or pk_lm5nopk
          First partitioned by TSS peak. Then, for isoforms in the 'no TSS peak' set, a linear regression model is
          employed to further classify them into 2, 3, 4, or 5 partitions.
      - cmb_lgt
          Using a logistic regression to combine TSS signals from multiple complementary data sets and partition training
          set isoform into 'expressed' and 'not expressed'. This partition model is only in use when either '--chipseq-read-files-multi-targets <string|' or '--chipseq-bed-files-multi-targets <string| is specified.

      Parameters for all the above models are learned from a training set. For detailed explainations,
      please see prior-enhanced RSEM's paper.
      (Default: 'pk')


outputs:

  isoform_results_file:
    type:
      - "null"
      - File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == null){
              return default_output_filename() + ".isoforms.results";
            } else {
              return inputs.output_filename + ".isoforms.results";
            }
        }

  gene_results_file:
    type:
      - "null"
      - File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == null){
              return default_output_filename() + ".genes.results";
            } else {
              return inputs.output_filename + ".genes.results";
            }
        }

  alleles_results_file:
    type:
      - "null"
      - File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == null){
              return default_output_filename() + ".alleles.results";
            } else {
              return inputs.output_filename + ".alleles.results";
            }
        }

  genome_bam_file:
    type:
      - "null"
      - File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == null){
              return default_output_filename() + ".genome.bam";
            } else {
              return inputs.output_filename + ".genome.bam";
            }
        }

  transcript_bam_file:
    type:
      - "null"
      - File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == null){
              return default_output_filename() + ".transcript.bam";
            } else {
              return inputs.output_filename + ".transcript.bam";
            }
        }

  transcript_sorted_bam_bai_pair:
    type:
      - "null"
      - File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == null){
              return default_output_filename() + ".transcript.sorted.bam";
            } else {
              return inputs.output_filename + ".transcript.sorted.bam";
            }
        }
    secondaryFiles: ${return self.basename + ".bai"}

  genome_sorted_bam_bai_pair:
    type:
      - "null"
      - File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == null){
              return default_output_filename() + ".genome.sorted.bam";
            } else {
              return inputs.output_filename + ".genome.sorted.bam";
            }
        }
    secondaryFiles: ${return self.basename + ".bai"}

  align_time_file:
    type:
      - "null"
      - File
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == null){
              return default_output_filename() + ".time";
            } else {
              return inputs.output_filename + ".time";
            }
        }

  stat_folder:
    type: Directory
    outputBinding:
      glob: |
        ${
            if (inputs.output_filename == null){
              return default_output_filename() + ".stat";
            } else {
              return inputs.output_filename + ".stat";
            }
        }

  unmapped_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: |
        ${
            let base = inputs.output_filename?inputs.output_filename:default_output_filename();
            return base+".stat/"+base+".cnt";
        }
      outputEval: |
        ${
          return parseInt(self[0].contents.split(/\r?\n/)[0].split(" ")[0]);
        }

  mapped_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: |
        ${
            let base = inputs.output_filename?inputs.output_filename:default_output_filename();
            return base+".stat/"+base+".cnt";
        }
      outputEval: |
        ${
          return parseInt(self[0].contents.split(/\r?\n/)[0].split(" ")[1]);
        }

  multimapped_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: |
        ${
            let base = inputs.output_filename?inputs.output_filename:default_output_filename();
            return base+".stat/"+base+".cnt";
        }
      outputEval: |
        ${
          return parseInt(self[0].contents.split(/\r?\n/)[0].split(" ")[2]);
        }

  total_reads_number:
    type: int
    outputBinding:
      loadContents: true
      glob: |
        ${
            let base = inputs.output_filename?inputs.output_filename:default_output_filename();
            return base+".stat/"+base+".cnt";
        }
      outputEval: |
        ${
          return parseInt(self[0].contents.split(/\r?\n/)[0].split(" ")[3]);
        }


baseCommand: [rsem-calculate-expression]
arguments:
# We should check if user set both upstream_read_file and downstream_read_file
# and didn't set paired_end manually - we will add it
  - valueFrom: |
      ${
        if (inputs.upstream_read_file && inputs.downstream_read_file && !inputs.paired_end){
          return "--paired-end";
        }
        return null;
      }
    position: 99
# We should get the value automatically from indices_folder
  - valueFrom: |
      ${
          for (var i = 0; i < inputs.indices_folder.listing.length; i++) {
              if (inputs.indices_folder.listing[i].basename.split('.').slice(-1)[0] == 'grp'){
                let name = inputs.indices_folder.listing[i].basename.split('.').slice(0,-1).join('.');
                return inputs.indices_folder.listing[i].path.split('/').slice(0,-1).join('/') + '/' + name;
              }
          }
          return null;
      }
    position: 103


$namespaces:
  s: http://schema.org/
$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "rsem-calculate-expression"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/rsem-calculate-expression.cwl
s:codeRepository: https://github.com/Barski-lab/workflows
s:license: http://www.apache.org/licenses/LICENSE-2.0

s:isPartOf:
  class: s:CreativeWork
  s:name: Common Workflow Language
  s:url: http://commonwl.org/

s:mainEntity:
  $import: ./metadata/rsem-metadata.yaml

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
  Tool runs rsem-calculate-expression.

  `reference_name` parameter for RSEM is resolved from `indices_folder` input.
  If `paired_end` input is not set, but both of the `upstream_read_file` and `downstream_read_file` are present,
  set `paired_end` automatically.

  `default_output_filename` function return prefix fot output files generated by RSEM based on `upstream_read_file` or
  `downstream_read_file` basename, if `output_filename` input is not provided

s:about: |
  Usage:
     rsem-calculate-expression [options] upstream_read_file(s) reference_name sample_name
     rsem-calculate-expression [options] --paired-end upstream_read_file(s) downstream_read_file(s) reference_name sample_name
     rsem-calculate-expression [options] --alignments [--paired-end] input reference_name sample_name

