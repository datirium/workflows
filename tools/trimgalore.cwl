cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_log_name = function() {
      let lognames = {};
      lognames["pair"] = (inputs.paired && inputs.input_file_pair) ? inputs.input_file_pair.basename+'_trimming_report.txt':null;
      lognames["single"] = inputs.input_file.basename+'_trimming_report.txt';
      return lognames;
    }

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/trimgalore:v0.4.4

inputs:

  input_file:
    type:
      - File
    inputBinding:
      position: 100
    doc: |
      Input FASTQ file

  input_file_pair:
    type:
      - "null"
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
      - "null"
      - int
    inputBinding:
      position: 5
      prefix: '-q'
    doc: |
      Trim low-quality ends from reads in addition to adapter removal.
      Default Phred score: 20.

  phred33:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 6
      prefix: '--phred33'
    doc: |
      Instructs Cutadapt to use ASCII+33 quality scores as Phred scores
      (Sanger/Illumina 1.9+ encoding) for quality trimming.
      Default: ON.

  phred64:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 7
      prefix: '--phred64'
    doc: |
        Instructs Cutadapt to use ASCII+64 quality scores as Phred scores
        (Illumina 1.5 encoding) for quality trimming.

  adapter:
    type:
      - "null"
      - string
    inputBinding:
      position: 10
      prefix: '-a'
    doc: |
      Adapter sequence to be trimmed. If not specified explicitly, Trim Galore will
      try to auto-detect whether the Illumina universal, Nextera transposase or Illumina
      small RNA adapter sequence was used. Also see '--illumina', '--nextera' and
      '--small_rna'. If no adapter can be detected within the first 1 million sequences
      of the first file specified Trim Galore defaults to '--illumina'.

  adapter_pair:
    type:
      - "null"
      - string
    inputBinding:
      position: 11
      prefix: '-a2'
    doc: |
      Optional adapter sequence to be trimmed off read 2 of paired-end files. This
      option requires '--paired' to be specified as well. If the libraries to be trimmed
      are smallRNA then a2 will be set to the Illumina small RNA 5' adapter automatically
      (GATCGTCGGACT).

  illumina:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 12
      prefix: '--illumina'
    doc: |
        Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter
        'AGATCGGAAGAGC' instead of the default auto-detection of adapter sequence.

  nextera:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 13
      prefix: '--nextera'
    doc: |
        Adapter sequence to be trimmed is the first 12bp of the Nextera adapter
        'CTGTCTCTTATA' instead of the default auto-detection of adapter sequence.

  small_rna:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 14
      prefix: '--small_rna'
    doc: |
        Adapter sequence to be trimmed is the first 12bp of the Illumina Small RNA 3' Adapter
        'TGGAATTCTCGG' instead of the default auto-detection of adapter sequence. Selecting
        to trim smallRNA adapters will also lower the --length value to 18bp. If the smallRNA
        libraries are paired-end then a2 will be set to the Illumina small RNA 5' adapter
        automatically (GATCGTCGGACT) unless -a 2 had been defined explicitly.

  max_length:
    type:
      - "null"
      - int
    inputBinding:
      position: 15
      prefix: '--max_length'
    doc: |
      Discard reads that are longer than <INT> bp after trimming. This is only advised for
      smallRNA sequencing to remove non-small RNA sequences.

  stringency:
    type:
      - "null"
      - int
    inputBinding:
      position: 16
      prefix: '--stringency'
    doc: |
      Overlap with adapter sequence required to trim a sequence. Defaults to a
      very stringent setting of 1, i.e. even a single bp of overlapping sequence
      will be trimmed off from the 3' end of any read.

  error_rate:
    type:
      - "null"
      - float
    inputBinding:
      position: 17
      prefix: '-e'
    doc: |
      Maximum allowed error rate (no. of errors divided by the length of the matching region)
      Default: 0.1

  gzip:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 18
      prefix: '--gzip'
    doc: |
      Compress the output file with GZIP. If the input files are GZIP-compressed
      the output files will automatically be GZIP compressed as well. As of v0.2.8 the
      compression will take place on the fly.

  dont_gzip:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 19
      prefix: '--dont_gzip'
    doc: |
      Output files won't be compressed with GZIP. This option overrides --gzip.

  length:
    type:
      - "null"
      - int
    inputBinding:
      position: 20
      prefix: '--length'
    doc: |
      Discard reads that became shorter than length INT because of either
      quality or adapter trimming. A value of '0' effectively disables
      this behaviour.
      Default: 20 bp.

  max_n:
    type:
      - "null"
      - int
    inputBinding:
      position: 21
      prefix: '--max_n'
    doc: |
      The total number of Ns (as integer) a read may contain before it will be removed altogether.
      In a paired-end setting, either read exceeding this limit will result in the entire
      pair being removed from the trimmed output files.

  trim_n:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 22
      prefix: '--trim-n'
    doc: |
      Removes Ns from either side of the read. This option does currently not work in RRBS mode.

  no_report_file:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 23
      prefix: '--no_report_file'
    doc: |
      If specified no report file will be generated.

  suppress_warn:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 24
      prefix: '--suppress_warn'
    doc: |
      If specified any output to STDOUT or STDERR will be suppressed.

  clip_R1:
    type:
      - "null"
      - int
    inputBinding:
      position: 25
      prefix: '--clip_R1'
    doc: |
      Instructs Trim Galore to remove <int> bp from the 5' end of read 1 (or single-end
      reads). This may be useful if the qualities were very poor, or if there is some
      sort of unwanted bias at the 5' end.
      Default: OFF.

  clip_R2:
    type:
      - "null"
      - int
    inputBinding:
      position: 26
      prefix: '--clip_R2'
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
      - "null"
      - int
    inputBinding:
      position: 27
      prefix: '--three_prime_clip_R1'
    doc: |
      Instructs Trim Galore to remove <int> bp from the 3' end of read 1 (or single-end
      reads) AFTER adapter/quality trimming has been performed. This may remove some unwanted
      bias from the 3' end that is not directly related to adapter sequence or basecall quality.
      Default: OFF.

  three_prime_clip_R2:
    type:
      - "null"
      - int
    inputBinding:
      position: 28
      prefix: '--three_prime_clip_R2'
    doc: |
      Instructs Trim Galore to remove <int> bp from the 3' end of read 2 AFTER
      adapter/quality trimming has been performed. This may remove some unwanted bias from
      the 3' end that is not directly related to adapter sequence or basecall quality.
      Default: OFF.

  rrbs:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 29
      prefix: '--rrbs'
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
      - "null"
      - boolean
    inputBinding:
      position: 30
      prefix: '--non_directional'
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
      - "null"
      - boolean
    inputBinding:
      position: 32
      prefix: '--paired'
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
      - "null"
      - boolean
    inputBinding:
      position: 33
      prefix: '--trim1'
    doc: |
      Trims 1 bp off every read from its 3' end. This may be needed for FastQ files that
      are to be aligned as paired-end data with Bowtie.

# NOTE influence on output
  retain_unpaired:
    type:
      - "null"
      - boolean
    inputBinding:
      position: 34
      prefix: '--retain_unpaired'
    doc: |
      If only one of the two paired-end reads became too short, the longer
      read will be written to either '.unpaired_1.fq' or '.unpaired_2.fq'
      output files. The length cutoff for unpaired single-end reads is
      governed by the parameters -r1/--length_1 and -r2/--length_2.
      Default: OFF.

  length_1:
    type:
      - "null"
      - int
    inputBinding:
      position: 35
      prefix: '-r1'
    doc: |
      Unpaired single-end read length cutoff needed for read 1 to be written to
      '.unpaired_1.fq' output file. These reads may be mapped in single-end mode.
      Default: 35 bp.

  length_2:
    type:
      - "null"
      - int
    inputBinding:
      position: 36
      prefix: '-r2'
    doc: |
      Unpaired single-end read length cutoff needed for read 2 to be written to
      '.unpaired_2.fq' output file. These reads may be mapped in single-end mode.
      Default: 35 bp.


outputs:
  trimmed_file:
    type: File
    outputBinding:
      glob: $((inputs.paired && inputs.input_file_pair) ? "*_val_1.fq*":"*_trimmed.fq*")

  trimmed_file_pair:
    type: File?
    outputBinding:
      glob: "*_val_2.fq*"

  unpaired_file_1:
    type: File?
    outputBinding:
      glob: "*_unpaired_1.fq*"

  unpaired_file_2:
    type: File?
    outputBinding:
      glob: "*_unpaired_2.fq*"

  report_file:
    type: File?
    outputBinding:
      glob: $(default_log_name()['single'])

  report_file_pair:
    type: File?
    outputBinding:
      glob: $(default_log_name()['pair'])

baseCommand: [trim_galore]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:mainEntity:
  $import: ./metadata/trimgalore-metadata.yaml

s:name: "trimgalore"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/trimgalore.cwl
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
  Tool runs Trimgalore - the wrapper around Cutadapt and FastQC to consistently apply adapter and quality trimming
  to FastQ files.

  `default_log_name` function returns names for generated log files (for both paired-end and single-end cases).
  `trim_galore` itself doesn't support setting custom names for output files.

  For paired-end data processing both `input_file_pair` and `paired` should be set. If either of them is not set,
  the other one becomes unset automatically.

s:about: |
  trim_galore [options] <filename(s)>

  Currently NOT supported as cwl tool options:

  -h/--help               Print this help message and exits.

  -v/--version            Print the version information and exits.

  --fastqc                Run FastQC in the default mode on the FastQ file once trimming is complete.

  --fastqc_args "<ARGS>"  Passes extra arguments to FastQC. If more than one argument is to be passed
                          to FastQC they must be in the form "arg1 arg2 etc.". An example would be:
                          --fastqc_args "--nogroup --outdir /home/". Passing extra arguments will
                          automatically invoke FastQC, so --fastqc does not have to be specified
                          separately.

  --path_to_cutadapt </path/to/cutadapt>     You may use this option to specify a path to the Cutadapt executable,
                          e.g. /my/home/cutadapt-1.7.1/bin/cutadapt. Else it is assumed that Cutadapt is in
                          the PATH.

  -o/--output_dir <DIR>   If specified all output will be written to this directory instead of the current
                          directory.

  --keep                  Keep the quality trimmed intermediate file. Default: off, which means
                          the temporary file is being deleted after adapter trimming. Only has
                          an effect for RRBS samples since other FastQ files are not trimmed
                          for poor qualities separately.


  General options:

  -q/--quality <INT>      Trim low-quality ends from reads in addition to adapter removal. For
                          RRBS samples, quality trimming will be performed first, and adapter
                          trimming is carried in a second round. Other files are quality and adapter
                          trimmed in a single pass. The algorithm is the same as the one used by BWA
                          (Subtract INT from all qualities; compute partial sums from all indices
                          to the end of the sequence; cut sequence at the index at which the sum is
                          minimal). Default Phred score: 20.

  --phred33               Instructs Cutadapt to use ASCII+33 quality scores as Phred scores
                          (Sanger/Illumina 1.9+ encoding) for quality trimming. Default: ON.

  --phred64               Instructs Cutadapt to use ASCII+64 quality scores as Phred scores
                          (Illumina 1.5 encoding) for quality trimming.

  -a/--adapter <STRING>   Adapter sequence to be trimmed. If not specified explicitly, Trim Galore will
                          try to auto-detect whether the Illumina universal, Nextera transposase or Illumina
                          small RNA adapter sequence was used. Also see '--illumina', '--nextera' and
                          '--small_rna'. If no adapter can be detected within the first 1 million sequences
                          of the first file specified Trim Galore defaults to '--illumina'.

  -a2/--adapter2 <STRING> Optional adapter sequence to be trimmed off read 2 of paired-end files. This
                          option requires '--paired' to be specified as well. If the libraries to be trimmed
                          are smallRNA then a2 will be set to the Illumina small RNA 5' adapter automatically
                          (GATCGTCGGACT).

  --illumina              Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter
                          'AGATCGGAAGAGC' instead of the default auto-detection of adapter sequence.

  --nextera               Adapter sequence to be trimmed is the first 12bp of the Nextera adapter
                          'CTGTCTCTTATA' instead of the default auto-detection of adapter sequence.

  --small_rna             Adapter sequence to be trimmed is the first 12bp of the Illumina Small RNA 3' Adapter
                          'TGGAATTCTCGG' instead of the default auto-detection of adapter sequence. Selecting
                          to trim smallRNA adapters will also lower the --length value to 18bp. If the smallRNA
                          libraries are paired-end then a2 will be set to the Illumina small RNA 5' adapter
                          automatically (GATCGTCGGACT) unless -a 2 had been defined explicitly.

  --max_length <INT>      Discard reads that are longer than <INT> bp after trimming. This is only advised for
                          smallRNA sequencing to remove non-small RNA sequences.


  --stringency <INT>      Overlap with adapter sequence required to trim a sequence. Defaults to a
                          very stringent setting of 1, i.e. even a single bp of overlapping sequence
                          will be trimmed off from the 3' end of any read.

  -e <ERROR RATE>         Maximum allowed error rate (no. of errors divided by the length of the matching
                          region) (default: 0.1)

  --gzip                  Compress the output file with GZIP. If the input files are GZIP-compressed
                          the output files will automatically be GZIP compressed as well. As of v0.2.8 the
                          compression will take place on the fly.

  --dont_gzip             Output files won't be compressed with GZIP. This option overrides --gzip.

  --length <INT>          Discard reads that became shorter than length INT because of either
                          quality or adapter trimming. A value of '0' effectively disables
                          this behaviour. Default: 20 bp.

                          For paired-end files, both reads of a read-pair need to be longer than
                          <INT> bp to be printed out to validated paired-end files (see option --paired).
                          If only one read became too short there is the possibility of keeping such
                          unpaired single-end reads (see --retain_unpaired). Default pair-cutoff: 20 bp.

  --max_n COUNT           The total number of Ns (as integer) a read may contain before it will be removed altogether.
                          In a paired-end setting, either read exceeding this limit will result in the entire
                          pair being removed from the trimmed output files.

  --trim-n                Removes Ns from either side of the read. This option does currently not work in RRBS mode.

  --no_report_file        If specified no report file will be generated.

  --suppress_warn         If specified any output to STDOUT or STDERR will be suppressed.

  --clip_R1 <int>         Instructs Trim Galore to remove <int> bp from the 5' end of read 1 (or single-end
                          reads). This may be useful if the qualities were very poor, or if there is some
                          sort of unwanted bias at the 5' end. Default: OFF.

  --clip_R2 <int>         Instructs Trim Galore to remove <int> bp from the 5' end of read 2 (paired-end reads
                          only). This may be useful if the qualities were very poor, or if there is some sort
                          of unwanted bias at the 5' end. For paired-end BS-Seq, it is recommended to remove
                          the first few bp because the end-repair reaction may introduce a bias towards low
                          methylation. Please refer to the M-bias plot section in the Bismark User Guide for
                          some examples. Default: OFF.

  --three_prime_clip_R1 <int>     Instructs Trim Galore to remove <int> bp from the 3' end of read 1 (or single-end
                          reads) AFTER adapter/quality trimming has been performed. This may remove some unwanted
                          bias from the 3' end that is not directly related to adapter sequence or basecall quality.
                          Default: OFF.

  --three_prime_clip_R2 <int>     Instructs Trim Galore to remove <int> bp from the 3' end of read 2 AFTER
                          adapter/quality trimming has been performed. This may remove some unwanted bias from
                          the 3' end that is not directly related to adapter sequence or basecall quality.
                          Default: OFF.

  --path_to_cutadapt </path/to/cutadapt>     You may use this option to specify a path to the Cutadapt executable,
                          e.g. /my/home/cutadapt-1.7.1/bin/cutadapt. Else it is assumed that Cutadapt is in
                          the PATH.


  RRBS-specific options (MspI digested material):

  --rrbs                  Specifies that the input file was an MspI digested RRBS sample (recognition
                          site: CCGG). Single-end or Read 1 sequences (paired-end) which were adapter-trimmed
                          will have a further 2 bp removed from their 3' end. Sequences which were merely
                          trimmed because of poor quality will not be shortened further. Read 2 of paired-end
                          libraries will in addition have the first 2 bp removed from the 5' end (by setting
                          '--clip_r2 2'). This is to avoid using artificial methylation calls from the filled-in
                          cytosine positions close to the 3' MspI site in sequenced fragments.
                          This option is not recommended for users of the NuGEN ovation RRBS System 1-16
                          kit (see below).

  --non_directional       Selecting this option for non-directional RRBS libraries will screen
                          quality-trimmed sequences for 'CAA' or 'CGA' at the start of the read
                          and, if found, removes the first two basepairs. Like with the option
                          '--rrbs' this avoids using cytosine positions that were filled-in
                          during the end-repair step. '--non_directional' requires '--rrbs' to
                          be specified as well. Note that this option does not set '--clip_r2 2' in
                          paired-end mode.


  Note for RRBS using the NuGEN Ovation RRBS System 1-16 kit:

  Owing to the fact that the NuGEN Ovation kit attaches a varying number of nucleotides (0-3) after each MspI
  site Trim Galore should be run WITHOUT the option --rrbs. This trimming is accomplished in a subsequent
  diversity trimming step afterwards (see their manual).


  Note for RRBS using MseI:

  If your DNA material was digested with MseI (recognition motif: TTAA) instead of MspI it is NOT necessary
  to specify --rrbs or --non_directional since virtually all reads should start with the sequence
  'TAA', and this holds true for both directional and non-directional libraries. As the end-repair of 'TAA'
  restricted sites does not involve any cytosines it does not need to be treated especially. Instead, simply
  run Trim Galore! in the standard (i.e. non-RRBS) mode.


  Paired-end specific options:

  --paired                This option performs length trimming of quality/adapter/RRBS trimmed reads for
                          paired-end files. To pass the validation test, both sequences of a sequence pair
                          are required to have a certain minimum length which is governed by the option
                          --length (see above). If only one read passes this length threshold the
                          other read can be rescued (see option --retain_unpaired). Using this option lets
                          you discard too short read pairs without disturbing the sequence-by-sequence order
                          of FastQ files which is required by many aligners.

                          Trim Galore! expects paired-end files to be supplied in a pairwise fashion, e.g.
                          file1_1.fq file1_2.fq SRR2_1.fq.gz SRR2_2.fq.gz ... .

  -t/--trim1              Trims 1 bp off every read from its 3' end. This may be needed for FastQ files that
                          are to be aligned as paired-end data with Bowtie. This is because Bowtie (1) regards
                          alignments like this:

                            R1 --------------------------->     or this:    ----------------------->  R1
                            R2 <---------------------------                       <-----------------  R2

                          as invalid (whenever a start/end coordinate is contained within the other read).
                          NOTE: If you are planning to use Bowtie2, BWA etc. you don't need to specify this option.

  --retain_unpaired       If only one of the two paired-end reads became too short, the longer
                          read will be written to either '.unpaired_1.fq' or '.unpaired_2.fq'
                          output files. The length cutoff for unpaired single-end reads is
                          governed by the parameters -r1/--length_1 and -r2/--length_2. Default: OFF.

  -r1/--length_1 <INT>    Unpaired single-end read length cutoff needed for read 1 to be written to
                          '.unpaired_1.fq' output file. These reads may be mapped in single-end mode.
                          Default: 35 bp.

  -r2/--length_2 <INT>    Unpaired single-end read length cutoff needed for read 2 to be written to
                          '.unpaired_2.fq' output file. These reads may be mapped in single-end mode.
                          Default: 35 bp.