cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

hints:
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/umi_tools:0.5.5--py36h470a237_0


inputs:


  dedup_log:
    type: string
    default: "dedup.log"
    label: "filename for a file with logging information"
    inputBinding:
      prefix: -L
      position: 7

  error_log:
    type: string
    default: "dedup_error.log"
    label: "filename for a file with errors information"
    inputBinding:
      prefix: -E
      position: 7

  output_file:
    type: string?
    label: "filename for a file where output is to go"
    inputBinding:
      prefix: -S
      position: 8
      valueFrom: |
        ${
            if (inputs.output_file == ""){
              var _f = inputs.input_file.location.split('/').slice(-1)[0].split('.');
              return _f.slice(0,-1).join('.') + '-dedup.' + _f.slice(-1)[0];
            } else {
              return inputs.output_file;
            }
        }
    default: ""

  input_file:
    type: File
    label: "BAM file to read from"
    inputBinding:
      prefix: -I
      position: 8

  input_pair:
    type: boolean?
    label: "BAM is paired end?"
    doc: |
         BAM is paired end - output both read pairs. This will also
         force the use of the template length to determine reads with
         the same mapping coordinates.
    inputBinding:
      prefix: --paired
      position: 3


outputs:
  output:
    type: File
    outputBinding:
      glob: |
       ${
            if (inputs.output_file == ""){
              var _f = inputs.input_file.location.split('/').slice(-1)[0].split('.');
              return _f.slice(0,-1).join('.') + '-dedup.' + _f.slice(-1)[0];
            } else {
              return inputs.output_file;
            }
        }

  log:
    type: File
    outputBinding:
      glob: $(inputs.dedup_log)

  error_log:
    type: File
    outputBinding:
      glob: $(inputs.error_log)


baseCommand: [umi_tools, dedup]


label: "Deduplicate BAM files based on the first mapping co-ordinate and the UMI attached to the read"

doc: |
  dedup.py - Deduplicate reads that are coded with a UMI
  =========================================================

  :Author: Ian Sudbery, Tom Smith
  :Release: $Id$
  :Date: |today|
  :Tags: Python UMI

  Purpose
  -------

  The purpose of this command is to deduplicate BAM files based
  on the first mapping co-ordinate and the UMI attached to the read.

  Selecting the representative read
  ---------------------------------

  The following criteria are applied to select the read that will be retained
  from a group of duplicated reads:

  1. The read with the lowest number of mapping coordinates (see
  --multimapping-detection-method option)
  2. The read with the highest mapping quality

  Otherwise a read is chosen at random.

s:about: |

  Usage:
  dedup.py - Deduplicate reads that are coded with a UMI
  =========================================================

  :Author: Ian Sudbery, Tom Smith
  :Release: $Id$
  :Date: |today|
  :Tags: Python UMI

  Purpose
  -------

  The purpose of this command is to deduplicate BAM files based
  on the first mapping co-ordinate and the UMI attached to the read.

  Selecting the representative read
  ---------------------------------

  The following criteria are applied to select the read that will be retained
  from a group of duplicated reads:

  1. The read with the lowest number of mapping coordinates (see
  --multimapping-detection-method option)
  2. The read with the highest mapping quality

  Otherwise a read is chosen at random.


  dedup-specific options
  ----------------------

  --output-stats (string, filename_prefix)
         Output edit distance statistics and UMI usage statistics
         using this prefix.

         Output files are:

         "[prefix]_stats_per_umi_per_position.tsv"
             Histogram of counts per position per UMI pre- and post-deduplication

         "[prefix]_stats_per_umi_per.tsv"
             Table of stats per umi. Number of times UMI was observed,
             total counts and median counts, pre- and post-deduplication

         "[prefix]_stats_edit_distance.tsv"
             Edit distance between UMIs at each position. Positions with a
             single UMI are reported seperately. Pre- and post-deduplication and
             inluding null expectations from random sampling of UMIs from the
             UMIs observed across all positions.


  It is assumed that the FASTQ files were processed with extract_umi.py
  before mapping and thus the UMI is the last word of the read name. e.g:

  @HISEQ:87:00000000_AATT

  where AATT is the UMI sequeuence.

  If you have used an alternative method which does not separate the
  read id and UMI with a "_", such as bcl2fastq which uses ":", you can
  specify the separator with the option "--umi-separator=<sep>",
  replacing <sep> with e.g ":".

  Alternatively, if your UMIs are encoded in a tag, you can specify this
  by setting the option --extract-umi-method=tag and set the tag name
  with the --umi-tag option. For example, if your UMIs are encoded in
  the 'UM' tag, provide the following options:
  "--extract-umi-method=tag --umi-tag=UM"

  Finally, if you have used umis to extract the UMI +/- cell barcode,
  you can specify --extract-umi-method=umis

  The start position of a read is considered to be the start of its alignment
  minus any soft clipped bases. A read aligned at position 500 with
  cigar 2S98M will be assumed to start at position 498.


  Input/Output Options
  -------------
  -i, --in-sam/-o, --out-sam
        By default, inputs are assumed to be in BAM format and outputs are
  written
        in BAM format. Use these options to specify the use of SAM format for
        inputs or outputs.

  -I    (string, filename) input file name
        The input file must be sorted and indexed.

  -S    (string, filename) output file name

  -L    (string, filename) log file name

  .. note::
     In order to get a valid sam/bam file you need to redirect logging
     information or turn it off logging via -v 0. You can redirect the
     logging to a file with -L <logfile> or use the --log2stderr option
     to send the logging to stderr.

  --paired
         BAM is paired end - output both read pairs. This will also
         force the use of the template length to determine reads with
         the same mapping coordinates.

  --extract-umi-method (choice)
        How are the barcodes encoded in the read?

        Options are:

        - "read_id" (default)
              Barcodes are contained at the end of the read separated as
              specified with --umi-separator option

        - "tag"
              Barcodes contained in a tag(s), see --umi-tag/--cell-tag
              options

        - "umis"
              Barcodes were extracted using umis (https://github.com/vals/umis)

  --umi-separator (string)
        Separator between read id and UMI. See --extract-umi-method above

  --umi-tag (string)
        Tag which contains UMI. See --extract-umi-method above

  --cell-tag (string)
        Tag which contains cell barcode. See --extract-umi-method above


  UMI grouping options
  ---------------------------

  --method (string, choice)
      count can be run with multiple methods to identify group of reads with
      the same (or similar) UMI(s), from which a single read is
      returned. All methods start by identifying the reads with the same
      mapping position.

      The simplest methods, unique and percentile, group reads with
      the exact same UMI. The network-based methods, cluster, adjacency and
      directional, build networks where nodes are UMIs and edges connect UMIs
      with an edit distance <= threshold (usually 1). The groups of reads
      are then defined from the network in a method-specific manner. For all
      the network-based methods, each read group is equivalent to one read
      count for the gene.

        "unique"
            Reads group share the exact same UMI

        "percentile"
            Reads group share the exact same UMI. UMIs with counts < 1% of the
            median counts for UMIs at the same position are ignored.

        "cluster"
            Identify clusters of connected UMIs (based on hamming distance
            threshold). Each network is a read group

        "adjacency"
            Cluster UMIs as above. For each cluster, select the node (UMI)
            with the highest counts. Visit all nodes one edge away. If all
            nodes have been visited, stop. Otherwise, repeat with remaining
            nodes until all nodes have been visted. Each step
            defines a read group.

        "directional" (default)
            Identify clusters of connected UMIs (based on hamming distance
            threshold) and umi A counts >= (2* umi B counts) - 1. Each
            network is a read group.


  --edit-distance-threshold (int)
         For the adjacency and cluster methods the threshold for the
         edit distance to connect two UMIs in the network can be
         increased. The default value of 1 works best unless the UMI is
         very long (>14bp)


  Single-cell RNA-Seq options
  ---------------------------

  --per-gene (string)
        Reads will be grouped together if they have the same gene.  This
        is useful if your library prep generates PCR duplicates with non
        identical alignment positions such as CEL-Seq. Note this option
        is hardcoded to be on with the count command. I.e counting is
        always performed per-gene. Must be combined with either
        --gene-tag or --per-contig option

  --gene-tag (string)
        Deduplicate per gene. The gene information is encoded in the bam
        read tag specified

  --assigned-status-tag (string)
        BAM tag which describes whether a read is assigned to a
        gene. Defaults to the same value as given for --gene-tag

  --skip-tags-regex (string)
        Used in conjunction with the --assigned-status-tag option. Skip any
  reads
        where the tag matches this regex.
        Default matches anything which starts with "__" or "Unassigned":
        ("^[__|Unassigned]")

  --per-contig
        Deduplicate per contig (field 3 in BAM; RNAME).
        All reads with the same contig will be considered to have the
        same alignment position. This is useful if you have aligned to a
        reference transcriptome with one transcript per gene. If you
        have aligned to a transcriptome with more than one transcript
        per gene, you can supply a map between transcripts and gene
        using the --gene-transcript-map option

  --gene-transcript-map (string)
        File mapping genes to transcripts (tab separated), e.g:

        gene1   transcript1
        gene1   transcript2
        gene2   transcript3

  --per-cell (string)
        Reads will only be grouped together if they have the same cell
        barcode. Can be combined with --per-gene.


  Debug options
  -------------

  --subset (float, [0-1])
        Only consider a fraction of the reads, chosen at random. This is useful
        for doing saturation analyses.

  --chrom (string)
        Only consider a single chromosome. This is useful for debugging purposes


  Group/Dedup options
  -------------------

  --no-sort-output
         By default, output is sorted. This involves the
         use of a temporary unsorted file since reads are considered in
         the order of their start position which may not be the same
         as their alignment coordinate due to soft-clipping and reverse
         alignments. The temp file will be saved in $TMPDIR and deleted
         when it has been sorted to the outfile. Use this option to turn
         off sorting.

  --spliced-is-unique
         Causes two reads that start in the same position on the same
         strand and having the same UMI to be considered unique if one is
  spliced
         and the other is not. (Uses the 'N' cigar operation to test for
         splicing)

  --soft-clip-threshold (int)
         Mappers that soft clip will sometimes do so rather than mapping a
         spliced read if there is only a small overhang over the exon
         junction. By setting this option, you can treat reads with at least
         this many bases soft-clipped at the 3' end as spliced.

  --multimapping-detection-method (string, choice)
         If the sam/bam contains tags to identify multimapping reads, you can
         specify for use when selecting the best read at a given loci.
         Supported tags are "NH", "X0" and "XT". If not specified, the read
         with the highest mapping quality will be selected

  --read-length
        Use the read length as a criteria when deduping, for e.g sRNA-Seq

  --buffer-whole-contig
        forces dedup to parse an entire contig before yielding any reads
        for deduplication. This is the only way to absolutely guarantee
        that all reads with the same start position are grouped together
        for deduplication since dedup uses the start position of the
        read, not the alignment coordinate on which the reads are
        sorted. However, by default, dedup reads for another 1000bp
        before outputting read groups which will avoid any reads being
        missed with short read sequencing (<1000bp)


  Options:
    --version             show program's version number and exit
    -h, --help            show this help message and exit

    dedup-specific options:
      --output-stats=STATS
                          Specify location to output stats

    Input options:
      -i, --in-sam        Input file is in sam format [default=False]
      --extract-umi-method=GET_UMI_METHOD
                          how is the read UMI +/ cell barcode encoded?
                          [default=read_id]
      --umi-separator=UMI_SEP
                          separator between read id and UMI
      --umi-tag=UMI_TAG   tag containing umi
      --umi-tag-split=UMI_TAG_SPLIT
                          split UMI in tag and take the first element
      --umi-tag-delimiter=UMI_TAG_DELIM
                          concatenate UMI in tag separated by delimiter
      --cell-tag=CELL_TAG
                          tag containing cell barcode
      --cell-tag-split=CELL_TAG_SPLIT
                          split cell barcode in tag and take the first
                          elementfor e.g 10X GEM tags
      --cell-tag-delimiter=CELL_TAG_DELIM
                          concatenate cell barcode in tag separated by delimiter
      --paired            paired BAM. [default=False]
      --mapping-quality=MAPPING_QUALITY
                          Minimum mapping quality for a read to be retained
                          [default=0]

    UMI grouping options:
      --method=METHOD     method to use for umi grouping [default=directional]
      --edit-distance-threshold=THRESHOLD
                          Edit distance theshold at which to join two UMIs when
                          grouping UMIs. [default=1]

    Single-cell RNA-Seq options:
      --per-gene          Group/Dedup/Count per gene. Must combine with either
                          --gene-tag or --per-contig
      --gene-tag=GENE_TAG
                          gene is defined by this bam tag [default=none]
      --assigned-status-tag=ASSIGNED_TAG
                          Bam tag describing whether read is assigned to a
                          geneBy defualt, this is set as the same tag as --gene-
                          tag
      --skip-tags-regex=SKIP_REGEX
                          Used with --gene-tag. Ignore reads where the gene-tag
                          matches this regex
      --per-contig        count per contig (field 3 in BAM; RNAME), e.g for
                          transcriptome where contig = gene
      --gene-transcript-map=GENE_TRANSCRIPT_MAP
                          file mapping transcripts to genes (tab separated)
      --per-cell          Group/Dedup/Count per cell

    Group/Dedup options:
      -o, --out-sam       Output alignments in sam format [default=False]
      --no-sort-output    Don't Sort the output
      --buffer-whole-contig
                          Read whole contig before outputting bundles:
                          guarantees that no reads are missed, but increases
                          memory usage
      --multimapping-detection-method=DETECTION_METHOD
                          Some aligners identify multimapping using bam tags.
                          Setting this option to NH, X0 or XT will use these
                          tags when selecting the best read amongst reads with
                          the same position and umi [default=none]
      --spliced-is-unique
                          Treat a spliced read as different to an unspliced one
                          [default=False]
      --soft-clip-threshold=SOFT_CLIP_THRESHOLD
                          number of bases clipped from 5' end beforeread is
                          counted as spliced [default=4]
      --read-length       use read length in addition to position and UMIto
                          identify possible duplicates [default=False]

    debug options:
      --ignore-umi        Ignore UMI and dedup only on position
      --chrom=CHROM       Restrict to one chromosome
      --subset=SUBSET     Use only a fraction of reads, specified by subset

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
