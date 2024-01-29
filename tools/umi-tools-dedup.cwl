cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_filename = function() {
        var ext = inputs.bam_file.basename.split('.').slice(-1)[0];
        var root = inputs.bam_file.basename.split('.').slice(0,-1).join('.');
        return inputs.output_filename?inputs.output_filename:root+"_dedup."+ext;
    };

hints:
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/umi_tools:1.0.1--py38h0213d0e_2


inputs:

  bash_script:
    type: string?
    default: |
      #!/bin/bash
      if [ "$0" = "true" ]; then
        echo "Copy $2 to __temp.bam"
        cp $2 __temp.bam
        samtools sort __temp.bam -o __temp_sorted.bam
        samtools index __temp_sorted.bam
        umi_tools dedup --random-seed=12345 -I __temp_sorted.bam "${@:3}"
        rm -f __temp.bam __temp_sorted.bam __temp_sorted.bam.bai
      else
        echo "Skip umi_tools dedup " ${@:1}
      fi
    inputBinding:
      position: 5
    doc: |
      Bash function to run umi_tools dedup with all input parameters or skip it if trigger is false

  trigger:
    type: boolean?
    default: true
    inputBinding:
      position: 6
      valueFrom: $(self?"true":"false")
    doc: |
      If true - run umi_tools dedup, if false - return bam_file input, previously staged into the output directory

  bam_file:
    type: File
    secondaryFiles:
    - .bai
    inputBinding:
      position: 7
      prefix: "-I"
    doc: "Input BAM file"

  paired_end:
    type: boolean?
    inputBinding:
      position: 8
      prefix: "--paired"
    doc: |
      Inputs BAM file is paired end - output both read pairs.
      This will also force the use of the template length to
      determine reads with the same mapping coordinates.

  output_filename:
    type: string?
    inputBinding:
      position: 9
      prefix: "-S"
      valueFrom: $(default_output_filename())
    default: ""
    doc: "Output filename"

  output_stats:
    type: string?
    inputBinding:
      position: 10
      prefix: "--output-stats="
      separate: false
    default: "umi_tools_stats"
    doc: "Specify location to output stats"

  multimapping_detection_method:
    type: string?
    inputBinding:
      position: 11
      prefix: "--multimapping-detection-method="
      separate: false
    doc: |
      Some aligners identify multimapping using bam tags.
      Setting this option to NH, X0 or XT will use these
      tags when selecting the best read amongst reads with
      the same position and umi [default=none]


outputs:

  dedup_bam_file:
    type: File
    outputBinding:
      glob: |
        ${ return inputs.trigger?default_output_filename():inputs.bam_file.basename }
    secondaryFiles: |
      ${
          if (inputs.bam_file.secondaryFiles && inputs.trigger == false){
            return inputs.bam_file.secondaryFiles;
          } else {
            return "null";
          }
        }

  output_stats:
    type:
      - "null"
      - File[]
    outputBinding:
      glob: $(inputs.output_stats + "*")

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: [bash, '-c']
stdout: umi_tools_dedup_stdout.log
stderr: umi_tools_dedup_stderr.log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:name: "umi-tools-dedup"
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
  Deduplicate BAM files based on the first mapping co-ordinate and the UMI attached to the read
  Only -I, --paired and -S parameters are implemented.


s:about: |
  UMI-Tools: Version 1.0.1

  dedup - Deduplicate reads using UMI and mapping coordinates

  Usage: umi_tools dedup [OPTIONS] [--stdin=IN_BAM] [--stdout=OUT_BAM]

        note: If --stdout is ommited, standard out is output. To
              generate a valid BAM file on standard out, please
              redirect log with --log=LOGFILE or --log2stderr

  For full UMI-tools documentation, see https://umi-tools.readthedocs.io/en/latest/

  Options:
    --version             show program's version number and exit

    dedup-specific options:
      --output-stats=STATS
                          Specify location to output stats

    Barcode extraction options:
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

    UMI grouping options:
      --method=METHOD     method to use for umi grouping [default=directional]
      --edit-distance-threshold=THRESHOLD
                          Edit distance theshold at which to join two UMIs when
                          grouping UMIs. [default=1]
      --spliced-is-unique
                          Treat a spliced read as different to an unspliced one
                          [default=False]
      --soft-clip-threshold=SOFT_CLIP_THRESHOLD
                          number of bases clipped from 5' end before read is
                          counted as spliced [default=4]
      --read-length       use read length in addition to position and UMI to
                          identify possible duplicates [default=False]

    single-cell RNA-Seq options:
      --per-gene          Group/Dedup/Count per gene. Must combine with either
                          --gene-tag or --per-contig
      --gene-tag=GENE_TAG
                          Gene is defined by this bam tag [default=none]
      --assigned-status-tag=ASSIGNED_TAG
                          Bam tag describing whether read is assigned to a gene
                          By defualt, this is set as the same tag as --gene-tag
      --skip-tags-regex=SKIP_REGEX
                          Used with --gene-tag. Ignore reads where the gene-tag
                          matches this regex
      --per-contig        group/dedup/count UMIs per contig (field 3 in BAM;
                          RNAME), e.g for transcriptome where contig = gene
      --gene-transcript-map=GENE_TRANSCRIPT_MAP
                          File mapping transcripts to genes (tab separated)
      --per-cell          group/dedup/count per cell

    group/dedup options:
      --buffer-whole-contig
                          Read whole contig before outputting bundles:
                          guarantees that no reads are missed, but increases
                          memory usage
      --multimapping-detection-method=DETECTION_METHOD
                          Some aligners identify multimapping using bam tags.
                          Setting this option to NH, X0 or XT will use these
                          tags when selecting the best read amongst reads with
                          the same position and umi [default=none]

    SAM/BAM options:
      --mapping-quality=MAPPING_QUALITY
                          Minimum mapping quality for a read to be retained
                          [default=0]
      --unmapped-reads=UNMAPPED_READS
                          How to handle unmapped reads. Options are 'discard',
                          'use' or 'correct' [default=discard]
      --chimeric-pairs=CHIMERIC_PAIRS
                          How to handle chimeric read pairs. Options are
                          'discard', 'use' or 'correct' [default=use]
      --unpaired-reads=UNPAIRED_READS
                          How to handle unpaired reads. Options are 'discard',
                          'use' or 'correct' [default=use]
      --ignore-umi        Ignore UMI and dedup only on position
      --chrom=CHROM       Restrict to one chromosome
      --subset=SUBSET     Use only a fraction of reads, specified by subset
      -i, --in-sam        Input file is in sam format [default=False]
      --paired            paired input BAM. [default=False]
      -o, --out-sam       Output alignments in sam format [default=False]
      --no-sort-output    Don't Sort the output

    input/output options:
      -I FILE, --stdin=FILE
                          file to read stdin from [default = stdin].
      -L FILE, --log=FILE
                          file with logging information [default = stdout].
      -E FILE, --error=FILE
                          file with error information [default = stderr].
      -S FILE, --stdout=FILE
                          file where output is to go [default = stdout].
      --temp-dir=FILE     Directory for temporary files. If not set, the bash
                          environmental variable TMPDIR is used[default = None].
      --log2stderr        send logging information to stderr [default = False].
      --compresslevel=COMPRESSLEVEL
                          Level of Gzip compression to use. Default (6)
                          matchesGNU gzip rather than python gzip default (which
                          is 9)

    profiling options:
      --timeit=TIMEIT_FILE
                          store timeing information in file [none].
      --timeit-name=TIMEIT_NAME
                          name in timing file for this class of jobs [all].
      --timeit-header     add header for timing information [none].

    common options:
      -v LOGLEVEL, --verbose=LOGLEVEL
                          loglevel [1]. The higher, the more output.
      -h, --help          output short help (command line options only).
      --help-extended     Output full documentation
      --random-seed=RANDOM_SEED
                          random seed to initialize number generator with
                          [none].
