cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- class: EnvVarRequirement
  envDef:
    http_proxy: $(inputs.http_proxy?inputs.http_proxy:"")
    https_proxy: $(inputs.https_proxy?inputs.https_proxy:"")

hints:
- class: DockerRequirement
  dockerPull: pegi3s/sratoolkit:3.0.1


inputs:

  script:
    type: string?
    default: |
        #!/bin/bash
        set -- "$0" "$@"

        SRA_IDS=()
        PARAMS=()

        for i in "$@"; do
            if [[ "$i" = "--split-files" ]] || [[ "$i" = "--split-3" ]]; then
                echo "Adding param $i"
                PARAMS+=($i)
            else
                echo "Adding SRR $i"
                SRA_IDS+=($i)
            fi
        done;

        echo "### Single files statistics" > single.md
        echo "### Merged files statistics" > merged.md

        for SRA in ${SRA_IDS[@]}; do
            echo "Downloading $SRA with ${PARAMS[@]}"
            fastq-dump --gzip --log-level info ${PARAMS[@]} $SRA
            j=1
            for FASTQ in $SRA*.gz; do
              echo "#### `basename $FASTQ`" >> single.md
              echo "**`zcat $FASTQ | wc -l`** lines, **`stat -c%s $FASTQ`** bytes, top **5** reads" >> single.md
              echo "\`\`\`" >> single.md
              echo "`zcat $FASTQ | head -n 20`" >> single.md
              echo "\`\`\`" >> single.md
              echo "Adding $FASTQ to read_$j.fastq.gz"
              cat $FASTQ >> read_$j.fastq.gz
              rm -f $FASTQ
              (( j++ ))
            done;
        done;
        
        for MERGED in read*.gz; do
            echo "#### `basename $MERGED`" >> merged.md
            echo "**`zcat $MERGED | wc -l`** lines, **`stat -c%s $MERGED`** bytes, top **5** reads" >> merged.md
        done;

        cat merged.md single.md > report.md
        rm -f merged.md single.md

    inputBinding:
      position: 1
    doc: |
      Bash function to run refgene-sort and atdp

  srr_id:
    type:
    - string
    - type: array
      items: string
    inputBinding:
      position: 60
    doc: |
      SRR identifiers

  split_files:
    type: boolean?
    inputBinding:
      position: 10
      prefix: "--split-files"
    doc: |
      Write reads into separate files. Read 
      number will be suffixed to the file name.  
      NOTE! The `--split-3` option is recommended. 
      In cases where not all spots have the same 
      number of reads, this option will produce 
      files that WILL CAUSE ERRORS in most programs 
      which process split pair fastq files. 

  split_3:
    type: boolean?
    inputBinding:
      position: 11
      prefix: "--split-3"
    doc: |
      3-way splitting for mate-pairs. For each 
      spot, if there are two biological reads 
      satisfying filter conditions, the first is 
      placed in the `*_1.fastq` file, and the 
      second is placed in the `*_2.fastq` file. If 
      there is only one biological read 
      satisfying the filter conditions, it is 
      placed in the `*.fastq` file.All other 
      reads in the spot are ignored.

  http_proxy:
    type: string?
    doc: |
      Optional HTTP proxy settings

  https_proxy:
    type: string?
    doc: |
      Optional HTTPS proxy settings

outputs:

  fastq_files:
    type:
    - "null"
    - type: array
      items: File
    outputBinding:
      glob: "*.gz"

  report_md:
    type: File
    outputBinding:
      glob: "report.md"

  stdout_log:
    type: stdout

  stderr_log:
    type: stderr


baseCommand: ["bash", "-c"]

stdout: fastq_dump_stdout.log
stderr: fastq_dump_stderr.log

successCodes: [1, 3]


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Fastq-Dump"
s:name: "Fastq-Dump"
s:alternateName: "Downloads FASTQ files from the provided SRR identifier"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/fastq-dump.cwl
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
  Fastq-Dump

  Downloads FASTQ files from the provided SRR identifier


s:about: |
  Usage:
    fastq-dump [options] <path> [<path>...]
    fastq-dump [options] <accession>

  INPUT
    -A|--accession <accession>       Replaces accession derived from <path> in 
                                    filename(s) and deflines (only for single 
                                    table dump) 
    --table <table-name>             Table name within cSRA object, default is 
                                    "SEQUENCE" 

  PROCESSING

  Read Splitting                     Sequence data may be used in raw form or
                                      split into individual reads
    --split-spot                     Split spots into individual reads 

  Full Spot Filters                  Applied to the full spot independently
                                      of --split-spot
    -N|--minSpotId <rowid>           Minimum spot id 
    -X|--maxSpotId <rowid>           Maximum spot id 
    --spot-groups <[list]>           Filter by SPOT_GROUP (member): name[,...] 
    -W|--clip                        Remove adapter sequences from reads 

  Common Filters                     Applied to spots when --split-spot is not
                                      set, otherwise - to individual reads
    -M|--minReadLen <len>            Filter by sequence length >= <len> 
    -R|--read-filter <[filter]>      Split into files by READ_FILTER value 
                                    optionally filter by value: 
                                    pass|reject|criteria|redacted 
    -E|--qual-filter                 Filter used in early 1000 Genomes data: no 
                                    sequences starting or ending with >= 10N 
    --qual-filter-1                  Filter used in current 1000 Genomes data 

  Filters based on alignments        Filters are active when alignment
                                      data are present
    --aligned                        Dump only aligned sequences 
    --unaligned                      Dump only unaligned sequences 
    --aligned-region <name[:from-to]>  Filter by position on genome. Name can 
                                    either be accession.version (ex: 
                                    NC_000001.10) or file specific name (ex: 
                                    "chr1" or "1"). "from" and "to" are 1-based 
                                    coordinates 
    --matepair-distance <from-to|unknown>  Filter by distance between matepairs. 
                                    Use "unknown" to find matepairs split 
                                    between the references. Use from-to to limit 
                                    matepair distance on the same reference 

  Filters for individual reads       Applied only with --split-spot set
    --skip-technical                 Dump only biological reads 

  OUTPUT
    -O|--outdir <path>               Output directory, default is working 
                                    directory '.' ) 
    -Z|--stdout                      Output to stdout, all split data become 
                                    joined into single stream 
    --gzip                           Compress output using gzip: deprecated, not 
                                    recommended 
    --bzip2                          Compress output using bzip2: deprecated, 
                                    not recommended 

  Multiple File Options              Setting these options will produce more
                                      than 1 file, each of which will be suffixed
                                      according to splitting criteria.
    --split-files                    Write reads into separate files. Read 
                                    number will be suffixed to the file name.  
                                    NOTE! The `--split-3` option is recommended. 
                                    In cases where not all spots have the same 
                                    number of reads, this option will produce 
                                    files that WILL CAUSE ERRORS in most programs 
                                    which process split pair fastq files. 
    --split-3                        3-way splitting for mate-pairs. For each 
                                    spot, if there are two biological reads 
                                    satisfying filter conditions, the first is 
                                    placed in the `*_1.fastq` file, and the 
                                    second is placed in the `*_2.fastq` file. If 
                                    there is only one biological read 
                                    satisfying the filter conditions, it is 
                                    placed in the `*.fastq` file.All other 
                                    reads in the spot are ignored. 
    -G|--spot-group                  Split into files by SPOT_GROUP (member name) 
    -R|--read-filter <[filter]>      Split into files by READ_FILTER value 
                                    optionally filter by value: 
                                    pass|reject|criteria|redacted 
    -T|--group-in-dirs               Split into subdirectories instead of files 
    -K|--keep-empty-files            Do not delete empty files 

  FORMATTING

  Sequence
    -C|--dumpcs <[cskey]>            Formats sequence using color space (default 
                                    for SOLiD),"cskey" may be specified for 
                                    translation 
    -B|--dumpbase                    Formats sequence using base space (default 
                                    for other than SOLiD). 

  Quality
    -Q|--offset <integer>            Offset to use for quality conversion, 
                                    default is 33 
    --fasta <[line width]>           FASTA only, no qualities, optional line 
                                    wrap width (set to zero for no wrapping) 
    --suppress-qual-for-cskey        suppress quality-value for cskey 

  Defline
    -F|--origfmt                     Defline contains only original sequence name 
    -I|--readids                     Append read id after spot id as 
                                    'accession.spot.readid' on defline 
    --helicos                        Helicos style defline 
    --defline-seq <fmt>              Defline format specification for sequence. 
    --defline-qual <fmt>             Defline format specification for quality. 
                                    <fmt> is string of characters and/or 
                                    variables. The variables can be one of: $ac 
                                    - accession, $si spot id, $sn spot 
                                    name, $sg spot group (barcode), $sl spot 
                                    length in bases, $ri read number, $rn 
                                    read name, $rl read length in bases. '[]' 
                                    could be used for an optional output: if 
                                    all vars in [] yield empty values whole 
                                    group is not printed. Empty value is empty 
                                    string or for numeric variables. Ex: 
                                    @$sn[_$rn]/$ri '_$rn' is omitted if name 
                                    is empty
  
  OTHER:
    --ngc <path>                     <path> to ngc file 
    --disable-multithreading         disable multithreading 
    -h|--help                        Output brief explanation of program usage 
    -V|--version                     Display the version of the program 
    -L|--log-level <level>           Logging level as number or enum string One 
                                    of (fatal|sys|int|err|warn|info) or (0-5) 
                                    Current/default is warn 
    -v|--verbose                     Increase the verbosity level of the program 
                                    Use multiple times for more verbosity 
    --ncbi_error_report              Control program execution environment 
                                    report generation (if implemented). One of 
                                    (never|error|always). Default is error 
    --legacy-report                  use legacy style 'Written spots' for tool 