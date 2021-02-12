cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/umi_tools:0.5.5--py36h470a237_0


inputs:

  fastq_file_1:
    type: File
    inputBinding:
      prefix: "-I"
    doc: "Input FASTQ file 1"

  fastq_file_2:
    type: File?
    inputBinding:
      prefix: "--read2-in="
      separate: false
    doc: "Input FASTQ file 2"

  extract_method:
    type:
      type: enum
      symbols: ["string", "regex"]
    inputBinding:
      prefix: "--extract-method="
      separate: false
    default: regex
    doc: "How to extract the umi +/- cell barcodes, choose from string or regex"

  bc_pattern_1:
    type: string
    inputBinding:
      prefix: "--bc-pattern="
      separate: false
    doc: "Barcode pattern 1"

  bc_pattern_2:
    type: string?
    inputBinding:
      prefix: "--bc-pattern2="
      separate: false
    doc: "Barcode pattern 2"

  output_filename_1:
    type: string
    inputBinding:
      prefix: "-S"
    doc: "Output filename 1"

  output_filename_2:
    type: string?
    inputBinding:
      prefix: "--read2-out="
      separate: false
    doc: "Output filename 2"


outputs:

  umi_fastq_file_1:
    type: File
    outputBinding:
      glob: $(inputs.output_filename_1)

  umi_fastq_file_2:
    type: File?
    outputBinding:
      glob: $(inputs.output_filename_2)

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: [umi_tools, extract]
stdout: umi_tools_extract_stdout.log
stderr: umi_tools_extract_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "umi-tools-extract"
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
  Extract UMI barcode from a read and add it to the read name, leaving
  any sample barcode in place. Can deal with paired end reads and UMIs
  split across the paired ends. Can also optionally extract cell
  barcodes and append these to the read name also. See the section below
  for an explanation for how to encode the barcode pattern(s) to
  specficy the position of the UMI +/- cell barcode.

s:about: |
  Usage:
  extract.py - Extract UMI from fastq
  ====================================================

  :Author: Ian Sudbery, Tom Smith
  :Release: $Id$
  :Date: |today|
  :Tags: Python UMI

  Purpose
  -------

  Extract UMI barcode from a read and add it to the read name, leaving
  any sample barcode in place. Can deal with paired end reads and UMIs
  split across the paired ends. Can also optionally extract cell
  barcodes and append these to the read name also. See the section below
  for an explanation for how to encode the barcode pattern(s) to
  specficy the position of the UMI +/- cell barcode.


  Filtering and correcting cell barcodes
  --------------------------------------

  umi_tools extract can optionally filter cell barcodes
  (--filter-cell-barcode) against a user-supplied whitelist
  (--whitelist). If a whitelist is not available for your data, e.g
  if you have performed droplet-based scRNA-Seq, you can use the
  whitelist tool.

  Cell barcodes which do not match the whitelist (user-generated or
  automatically generated) can also be optionally corrected using the
  --error-correct-cell option.

  The whitelist should be in the following format (tab-separated):

      AAAAAA      AGAAAA
      AAAATC
      AAACAT
      AAACTA      AAACTN,GAACTA
      AAATAC
      AAATCA      GAATCA
      AAATGT      AAAGGT,CAATGT

  Where column 1 is the whitelisted cell barcodes and column 2 is
  the list (comma-separated) of other cell barcodes which should be
  corrected to the barcode in column 1. If the --error-correct-cell
  option is not used, this column will be ignored. Any additional columns
  in the whitelist input, such as the counts columns from the output of
  umi_tools whitelist, will be ignored.


  Usage:
  ------

  For single ended reads:
          umi_tools extract --extract-method=string
          --bc-pattern=[PATTERN] -L extract.log [OPTIONS]

  reads from stdin and outputs to stdout.

  For paired end reads:
          umi_tools extract --extract-method=string
          --bc-pattern=[PATTERN] --bc-pattern2=[PATTERN]
          --read2-in=[FASTQIN] --read2-out=[FASTQOUT] -L extract.log [OPTIONS]

  reads end one from stdin and end two from FASTQIN and outputs end one to stdin
  and end two to FASTQOUT.


  Using regex and filtering against a whitelist of cell barcodes:
          umi_tools extract --extract-method=regex --filter-cell-barcode
          --bc-pattern=[REGEX] --whitlist=[WHITELIST_TSV]
          -L extract.log [OPTIONS]


  Barcode extraction
  ------------------

  There are two methods enabled to extract the umi barcode (+/- cell
  barcode). For both methods, the patterns should be provided using the
  --bc-pattern and --bc-pattern2 options. The method is specified using
  the --extract-method option

  -'string':
        This should be used where the barcodes are always in the same
        place in the read.

        - N = UMI position (required)
        - C = cell barcode position (optional)
        - X = sample position (optional)

        Bases with Ns and Cs will be extracted and added to the read
        name. The corresponding sequence qualities will be removed from
        the read. Bases with an X will be reattached to the read.

        E.g. If the pattern is NNNNCC,
        Then the read:
        @HISEQ:87:00000000 read1
        AAGGTTGCTGATTGGATGGGCTAG
        DA1AEBFGGCG01DFH00B1FF0B
        +
        will become:
        @HISEQ:87:00000000_TT_AAGG read1
        GCTGATTGGATGGGCTAG
        1AFGGCG01DFH00B1FF0B
        +

        where 'TT' is the cell barcode and 'AAGG' is the UMI.

  -'regex'
        This method allows for more flexible barcode extraction and
        should be used where the cell barcodes are variable in
        length. Alternatively, the regex option can also be used to
        filter out reads which do not contain an expected adapter
        sequence.

        The expected groups in the regex are:

        umi_n = UMI positions, where n can be any value (required)
        cell_n = cell barcode positions, where n can be any value (optional)
        discard_n = positions to discard, where n can be any value (optional)

        UMI positions and cell barcode positions will be extracted and
        added to the read name. The corresponding sequence qualities
        will be removed from the read. Discard bases and the
        corresponding quality scores will be removed from the read. All
        bases matched by other groups or components of the regex will be
        reattached to the read sequence

        For example, the following regex can be used to extract reads
        from the Klein et al inDrop data:

        (?P<cell_1>.{8,12})(?P<discard_1>GAGTGATTGCTTGTGACGCCTT)(?P<cell_2>.{8}
  )(?P<umi_1>.{6})T{3}.*

        Where only reads with a 3' T-tail and GAGTGATTGCTTGTGACGCCTT in
        the correct position to yield two cell barcodes of 8-12 and 8bp
        respectively, and a 6bp UMI will be retained.

        You can also specify fuzzy matching to allow errors. For example if
        the discard group above was specified as below this would enable
        matches with up to 2 errors in the discard_1 group.

        (?P<discard_1>GAGTGATTGCTTGTGACGCCTT){s<=2}

        Note that all UMIs must be the same length for downstream
        processing with dedup, group or count commands


  additional whitelist/extract options
  -------------------------

  --3prime
        By default the barcode is assumed to be on the 5' end of the
        read, but use this option to sepecify that it is on the 3' end
        instead. This option only works with --extract-method=string
        since 3' encoding can be specified explicitly with a regex, e.g
        ".*(?P<umi_1>.{5})$"



  Options:
    --version             show program's version number and exit
    -h, --help            show this help message and exit
    -p PATTERN, --bc-pattern=PATTERN
                          Barcode pattern
    --bc-pattern2=PATTERN2
                          Barcode pattern for paired reads
    --3prime              barcode is on 3' end of read.
    --read2-in=READ2_IN   file name for read pairs
    --read2-out=READ2_OUT
                          file to output processed paired read to
    --read2-stdout        Paired reads, send read2 to stdout, discarding read1
    --quality-filter-threshold=QUALITY_FILTER_THRESHOLD
                          Remove reads where any UMI base quality score falls
                          below this threshold
    --quality-filter-mask=QUALITY_FILTER_MASK
                          If a UMI base has a quality below this threshold,
                          replace the base with 'N'
    --quality-encoding=QUALITY_ENCODING
                          Quality score encoding. Choose from 'phred33'[33-77]
                          'phred64' [64-106] or 'solexa' [59-106]
    --extract-method=EXTRACT_METHOD
                          How to extract the umi +/- cell barcodes, Choose from
                          'string' or 'regex'
    --filter-cell-barcode
                          Filter the cell barcodes
    --error-correct-cell  Correct errors in the cell barcode
    --whitelist=WHITELIST
                          A whitelist of accepted cell barcodes
    --blacklist=BLACKLIST
                          A blacklist of accepted cell barcodes
    --reads-subset=READS_SUBSET
                          Only extract from the first N reads. If N is greater
                          than the number of reads, all reads will be used
    --reconcile-pairs     Allow the presences of reads in read2 input that
                          arenot present in read1 input. This allows cell
                          barcodefiltering of read1s without considering read2s

    profiling options:
      --timeit=TIMEIT_FILE
                          store timeing information in file [none].
      --timeit-name=TIMEIT_NAME
                          name in timing file for this class of jobs [all].
      --timeit-header     add header for timing information [none].

    common options:
      -v LOGLEVEL, --verbose=LOGLEVEL
                          loglevel [1]. The higher, the more output.
      -?                  output short help (command line options only.
      --random-seed=RANDOM_SEED
                          random seed to initialize number generator with
                          [none].

    Input/output options:
      -I FILE, --stdin=FILE
                          file to read stdin from [default = stdin].
      -L FILE, --log=FILE
                          file with logging information [default = stdout].
      -E FILE, --error=FILE
                          file with error information [default = stderr].
      -S FILE, --stdout=FILE
                          file where output is to go [default = stdout].
      --log2stderr        send logging information to stderr [default = False].
      --compresslevel=COMPRESSLEVEL
                          Level of Gzip compression to use. Default (6)
                          matchesGNU gzip rather than python gzip default (which
                          is 9)