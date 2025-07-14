cwlVersion: v1.0
class: CommandLineTool
requirements:
- class: ResourceRequirement
  ramMin: 7024
  coresMin: 1
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_log_name = function() { var lognames = {}; lognames["pair"] = (inputs.paired && inputs.input_file_pair) ? inputs.input_file_pair.basename+'_trimming_report.txt':null; lognames["single"] = inputs.input_file.basename+'_trimming_report.txt'; return lognames; }
- class: InitialWorkDirRequirement
  listing: |
    ${
      var listing = [inputs.input_file]
      if (inputs.input_file_pair){
        listing.push(inputs.input_file_pair);
      }
      return listing;
    }
hints:
- class: DockerRequirement
  dockerPull: scidap/trimgalore:v0.6.6
inputs:
  bash_script:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = "true" ]
      then
        echo "Run: trimgalore " ${@:1}
        trim_galore "${@:1}"
      else
        echo "Skip run trimgalore"
      fi
    inputBinding:
      position: 1
    doc: |
      Bash function to run trimgalore with all input parameters or skip it if trigger is false
  trigger:
    type: boolean?
    default: true
    inputBinding:
      position: 2
      valueFrom: $(self?"true":"false")
    doc: |
      If true - run trimgalore, if false - return input fastq files, previously staged into output directory.
      Use valueFrom to return string instead of boolean to make sure that value is printed to command line in both cases
  input_file:
    type:
    - File
    inputBinding:
      position: 100
    doc: |
      Input FASTQ file
  input_file_pair:
    type:
    - 'null'
    - File
    inputBinding:
      position: 101
      valueFrom: |
        ${
            if (!inputs.paired){
              return null;
            } else {
              return self;
            }
        }
    doc: |
      Input FASTQ file, if paired-end data is analysed
  quality:
    type:
    - 'null'
    - int
    inputBinding:
      position: 5
      prefix: -q
    doc: |
      Trim low-quality ends from reads in addition to adapter removal.
      Default Phred score: 20.
  phred33:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 6
      prefix: --phred33
    doc: |
      Instructs Cutadapt to use ASCII+33 quality scores as Phred scores
      (Sanger/Illumina 1.9+ encoding) for quality trimming.
      Default: ON.
  phred64:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 7
      prefix: --phred64
    doc: |
      Instructs Cutadapt to use ASCII+64 quality scores as Phred scores
      (Illumina 1.5 encoding) for quality trimming.
  adapter:
    type:
    - 'null'
    - string
    inputBinding:
      position: 10
      prefix: -a
    doc: |
      Adapter sequence to be trimmed. If not specified explicitly, Trim Galore will
      try to auto-detect whether the Illumina universal, Nextera transposase or Illumina
      small RNA adapter sequence was used. Also see '--illumina', '--nextera' and
      '--small_rna'. If no adapter can be detected within the first 1 million sequences
      of the first file specified Trim Galore defaults to '--illumina'.
  adapter_pair:
    type:
    - 'null'
    - string
    inputBinding:
      position: 11
      prefix: -a2
    doc: |
      Optional adapter sequence to be trimmed off read 2 of paired-end files. This
      option requires '--paired' to be specified as well. If the libraries to be trimmed
      are smallRNA then a2 will be set to the Illumina small RNA 5' adapter automatically
      (GATCGTCGGACT).
  illumina:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 12
      prefix: --illumina
    doc: |
      Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter
      'AGATCGGAAGAGC' instead of the default auto-detection of adapter sequence.
  nextera:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 13
      prefix: --nextera
    doc: |
      Adapter sequence to be trimmed is the first 12bp of the Nextera adapter
      'CTGTCTCTTATA' instead of the default auto-detection of adapter sequence.
  small_rna:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 14
      prefix: --small_rna
    doc: |
      Adapter sequence to be trimmed is the first 12bp of the Illumina Small RNA 3' Adapter
      'TGGAATTCTCGG' instead of the default auto-detection of adapter sequence. Selecting
      to trim smallRNA adapters will also lower the --length value to 18bp. If the smallRNA
      libraries are paired-end then a2 will be set to the Illumina small RNA 5' adapter
      automatically (GATCGTCGGACT) unless -a 2 had been defined explicitly.
  max_length:
    type:
    - 'null'
    - int
    inputBinding:
      position: 15
      prefix: --max_length
    doc: |
      Discard reads that are longer than <INT> bp after trimming. This is only advised for
      smallRNA sequencing to remove non-small RNA sequences.
  stringency:
    type:
    - 'null'
    - int
    inputBinding:
      position: 16
      prefix: --stringency
    doc: |
      Overlap with adapter sequence required to trim a sequence. Defaults to a
      very stringent setting of 1, i.e. even a single bp of overlapping sequence
      will be trimmed off from the 3' end of any read.
  error_rate:
    type:
    - 'null'
    - float
    inputBinding:
      position: 17
      prefix: -e
    doc: |
      Maximum allowed error rate (no. of errors divided by the length of the matching region)
      Default: 0.1
  gzip:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 18
      prefix: --gzip
    doc: |
      Compress the output file with GZIP. If the input files are GZIP-compressed
      the output files will automatically be GZIP compressed as well. As of v0.2.8 the
      compression will take place on the fly.
  dont_gzip:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 19
      prefix: --dont_gzip
    doc: |
      Output files won't be compressed with GZIP. This option overrides --gzip.
  length:
    type:
    - 'null'
    - int
    inputBinding:
      position: 20
      prefix: --length
    doc: |
      Discard reads that became shorter than length INT because of either
      quality or adapter trimming. A value of '0' effectively disables
      this behaviour.
      Default: 20 bp.
  max_n:
    type:
    - 'null'
    - int
    inputBinding:
      position: 21
      prefix: --max_n
    doc: |
      The total number of Ns (as integer) a read may contain before it will be removed altogether.
      In a paired-end setting, either read exceeding this limit will result in the entire
      pair being removed from the trimmed output files.
  trim_n:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 22
      prefix: --trim-n
    doc: |
      Removes Ns from either side of the read. This option does currently not work in RRBS mode.
  no_report_file:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 23
      prefix: --no_report_file
    doc: |
      If specified no report file will be generated.
  suppress_warn:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 24
      prefix: --suppress_warn
    doc: |
      If specified any output to STDOUT or STDERR will be suppressed.
  clip_R1:
    type:
    - 'null'
    - int
    inputBinding:
      position: 25
      prefix: --clip_R1
    doc: |
      Instructs Trim Galore to remove <int> bp from the 5' end of read 1 (or single-end
      reads). This may be useful if the qualities were very poor, or if there is some
      sort of unwanted bias at the 5' end.
      Default: OFF.
  clip_R2:
    type:
    - 'null'
    - int
    inputBinding:
      position: 26
      prefix: --clip_R2
    doc: |
      Instructs Trim Galore to remove <int> bp from the 5' end of read 2 (paired-end reads
      only). This may be useful if the qualities were very poor, or if there is some sort
      of unwanted bias at the 5' end. For paired-end BS-Seq, it is recommended to remove
      the first few bp because the end-repair reaction may introduce a bias towards low
      methylation. Please refer to the M-bias plot section in the Bismark User Guide for
      some examples.
      Default: OFF.
  three_prime_clip_R1:
    type:
    - 'null'
    - int
    inputBinding:
      position: 27
      prefix: --three_prime_clip_R1
    doc: |
      Instructs Trim Galore to remove <int> bp from the 3' end of read 1 (or single-end
      reads) AFTER adapter/quality trimming has been performed. This may remove some unwanted
      bias from the 3' end that is not directly related to adapter sequence or basecall quality.
      Default: OFF.
  three_prime_clip_R2:
    type:
    - 'null'
    - int
    inputBinding:
      position: 28
      prefix: --three_prime_clip_R2
    doc: |
      Instructs Trim Galore to remove <int> bp from the 3' end of read 2 AFTER
      adapter/quality trimming has been performed. This may remove some unwanted bias from
      the 3' end that is not directly related to adapter sequence or basecall quality.
      Default: OFF.
  rrbs:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 29
      prefix: --rrbs
    doc: |
      Specifies that the input file was an MspI digested RRBS sample (recognition
      site: CCGG). Single-end or Read 1 sequences (paired-end) which were adapter-trimmed
      will have a further 2 bp removed from their 3' end. Sequences which were merely
      trimmed because of poor quality will not be shortened further. Read 2 of paired-end
      libraries will in addition have the first 2 bp removed from the 5' end (by setting
      '--clip_r2 2'). This is to avoid using artificial methylation calls from the filled-in
      cytosine positions close to the 3' MspI site in sequenced fragments.
      This option is not recommended for users of the NuGEN ovation RRBS System 1-16
      kit (see below).
  non_directional:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 30
      prefix: --non_directional
    doc: |
      Selecting this option for non-directional RRBS libraries will screen
      quality-trimmed sequences for 'CAA' or 'CGA' at the start of the read
      and, if found, removes the first two basepairs. Like with the option
      '--rrbs' this avoids using cytosine positions that were filled-in
      during the end-repair step. '--non_directional' requires '--rrbs' to
      be specified as well. Note that this option does not set '--clip_r2 2' in
      paired-end mode.
  paired:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 32
      prefix: --paired
      valueFrom: |
        ${
            if (!inputs.input_file_pair){
              return null;
            } else {
              return self;
            }
        }
    doc: |
      This option performs length trimming of quality/adapter/RRBS trimmed reads for
      paired-end files. To pass the validation test, both sequences of a sequence pair
      are required to have a certain minimum length which is governed by the option
      --length (see above). If only one read passes this length threshold the
      other read can be rescued (see option --retain_unpaired). Using this option lets
      you discard too short read pairs without disturbing the sequence-by-sequence order
      of FastQ files which is required by many aligners.
  trim1:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 33
      prefix: --trim1
    doc: |
      Trims 1 bp off every read from its 3' end. This may be needed for FastQ files that
      are to be aligned as paired-end data with Bowtie.
  retain_unpaired:
    type:
    - 'null'
    - boolean
    inputBinding:
      position: 34
      prefix: --retain_unpaired
    doc: |
      If only one of the two paired-end reads became too short, the longer
      read will be written to either '.unpaired_1.fq' or '.unpaired_2.fq'
      output files. The length cutoff for unpaired single-end reads is
      governed by the parameters -r1/--length_1 and -r2/--length_2.
      Default: OFF.
  length_1:
    type:
    - 'null'
    - int
    inputBinding:
      position: 35
      prefix: -r1
    doc: |
      Unpaired single-end read length cutoff needed for read 1 to be written to
      '.unpaired_1.fq' output file. These reads may be mapped in single-end mode.
      Default: 35 bp.
  length_2:
    type:
    - 'null'
    - int
    inputBinding:
      position: 36
      prefix: -r2
    doc: |
      Unpaired single-end read length cutoff needed for read 2 to be written to
      '.unpaired_2.fq' output file. These reads may be mapped in single-end mode.
      Default: 35 bp.
outputs:
  trimmed_file:
    type: File
    outputBinding:
      glob: |
        ${
            if (inputs.trigger == false){
              return inputs.input_file.basename;
            } else if (inputs.paired && inputs.input_file_pair){
              return "*_val_1.fq*";
            } else {
              return "*_trimmed.fq*";
            }
        }
  trimmed_file_pair:
    type: File?
    outputBinding:
      glob: |
        ${
            if (inputs.paired && inputs.input_file_pair && inputs.trigger == false){
              return inputs.input_file_pair.basename;
            } else {
              return "*_val_2.fq*";
            }
        }
  unpaired_file_1:
    type: File?
    outputBinding:
      glob: '*_unpaired_1.fq*'
  unpaired_file_2:
    type: File?
    outputBinding:
      glob: '*_unpaired_2.fq*'
  report_file:
    type: File?
    outputBinding:
      glob: $(default_log_name()['single'])
  report_file_pair:
    type: File?
    outputBinding:
      glob: $(default_log_name()['pair'])
baseCommand:
- bash
- -c
doc: |
  Tool runs Trimgalore - the wrapper around Cutadapt and FastQC to consistently apply adapter and quality trimming
  to FastQ files.

  `default_log_name` function returns names for generated log files (for both paired-end and single-end cases).
  `trim_galore` itself doesn't support setting custom names for output files.

  For paired-end data processing both `input_file_pair` and `paired` should be set. If either of them is not set,
  the other one becomes unset automatically.

  If input trigger was set to false, skip running trimaglore and return unchanged input files
label: trimgalore
