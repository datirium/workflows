#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: envvar-global.yml
- class: ShellCommandRequirement
- class: EnvVarRequirement
  envDef:
    - envName: "BOWTIE_INDEXES"
      envValue: $(inputs.ebwtFile.dirname)
- class: InlineJavascriptRequirement
  expressionLib:
  - var parseEbwtFile = function() {

      var dirname = inputs.ebwtFile.location.split("/").slice(0,-1).join("/");
      var basename = inputs.ebwtFile.location.split("/").slice(-1)[0];
      var extname = basename.split(".").slice(-1)[0];
      var ewbtBasename = "";
      if ( (/.*\.rev\.\d\.ebwt$/i).test(basename) || (/.*\.rev\.\d\.ebwtl$/i).test(basename) ){
        ewbtBasename = basename.split(".").slice(0, -3).join (".");
      } else if ( (/.*\.\d\.ebwt$/i).test(basename) || (/.*\.\d\.ebwtl$/i).test(basename) ){
        ewbtBasename = basename.split(".").slice(0, -2).join (".");
      }

      return [dirname, basename, extname, ewbtBasename]
    };


hints:
- class: DockerRequirement
    #dockerImageId: scidap/bowtie:v1.1.2 #not yet ready
  dockerPull: scidap/bowtie:v1.1.2
  dockerFile: >
    $import: ./dockerfiles/bowtie-Dockerfile

inputs:

  ebwtFile:
    type: File
    format:
      - http://edamontology.org/format_3484 # ebwt
      - http://edamontology.org/format_3491 # ebwtl
    doc: >
      One of the index files.
      OLD:
        The basename of the index to be searched.
        The basename is the name of any of the index files up to but not including the final .1.ebwt / .rev.1.ebwt / etc. bowtie looks for
        the specified index first in the current directory,
        then in the indexes subdirectory under the directory where the bowtie executable is located,
        then looks in the directory specified in the BOWTIE_INDEXES environment variable.
    inputBinding:
      position: 81
      valueFrom: $(self.secondaryFiles[0].dirname + "/" + parseEbwtFile()[3])
    secondaryFiles: |
      ${
        return [{"class": "File", "location": parseEbwtFile()[0] + "/" + parseEbwtFile()[3] + ".1." + parseEbwtFile()[2]},
                {"class": "File", "location": parseEbwtFile()[0] + "/" + parseEbwtFile()[3] + ".2." + parseEbwtFile()[2]},
                {"class": "File", "location": parseEbwtFile()[0] + "/" + parseEbwtFile()[3] + ".3." + parseEbwtFile()[2]},
                {"class": "File", "location": parseEbwtFile()[0] + "/" + parseEbwtFile()[3] + ".4." + parseEbwtFile()[2]},
                {"class": "File", "location": parseEbwtFile()[0] + "/" + parseEbwtFile()[3] + ".rev.1." + parseEbwtFile()[2]},
                {"class": "File", "location": parseEbwtFile()[0] + "/" + parseEbwtFile()[3] + ".rev.2." + parseEbwtFile()[2]},
        ]
      }

  filelist:
    type:
      - "null"
      - File
      - type: array
        items: File
    doc: |
      {-1 <m1> -2 <m2> | --12 <r> | <s>}
      <m1>    Comma-separated list of files containing upstream mates (or the
            sequences themselves, if -c is set) paired with mates in <m2>
      <m2>    Comma-separated list of files containing downstream mates (or the
            sequences themselves if -c is set) paired with mates in <m1>
      <r>     Comma-separated list of files containing Crossbow-style reads.  Can be
            a mixture of paired and unpaired.  Specify "-"for stdin.
      <s>     Comma-separated list of files containing unpaired reads, or the
            sequences themselves, if -c is set.  Specify "-"for stdin.
    inputBinding:
      itemSeparator: ","
      position: 83

  filelist_mates:
    type:
      - "null"
      - File
      - type: array
        items: File
    inputBinding:
      itemSeparator: ","
      position: 85

  filelist_crossbow:
    type: File[]?
    inputBinding:
      itemSeparator: ","
      position: 86
      prefix: "-12"
    doc: |
      Comma-separated list of files containing Crossbow-style reads.
      Can be a mixture of paired and unpaired.  Specify "-"for stdin.


  filename:
    type: string
    inputBinding:
      position: 90

  q:
    type:
      - "null"
      - boolean
    doc: "query input files are FASTQ .fq/.fastq (default)\n"
    inputBinding:
      position: 1
      prefix: '-q'

  f:
    type:
      - "null"
      - boolean
    doc: "query input files are (multi-)FASTA .fa/.mfa\n"
    inputBinding:
      position: 1
      prefix: '-f'

  r:
    type:
      - "null"
      - boolean
    doc: "query input files are raw one-sequence-per-line\n"
    inputBinding:
      position: 1
      prefix: '-r'

  c:
    type:
      - "null"
      - boolean
    doc: "query sequences given on cmd line (as <mates>, <singles>)\n"
    inputBinding:
      position: 1
      prefix: '-c'

  C:
    type:
      - "null"
      - boolean
    doc: "reads and index are in colorspace\n"
    inputBinding:
      position: 1
      prefix: '-C'
  Q:
    type:
      - "null"
      - File
    doc: >
      --quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C
    inputBinding:
      position: 1
      prefix: '-Q'
  Q1:
    type:
      - "null"
      - boolean
    doc: |
      --Q2 <file>   same as -Q, but for mate files 1 and 2 respectively
    inputBinding:
      position: 1
      prefix: '--Q1'
  s:
    type:
      - "null"
      - int
    doc: |
      --skip <int>    skip the first <int> reads/pairs in the input
    inputBinding:
      position: 1
      prefix: '-s'
  u:
    type:
      - "null"
      - int
    doc: |
      --qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)
    inputBinding:
      position: 1
      prefix: '-u'
  5p:
    type: int?
    doc: |
      --trim5 <int>   trim <int> bases from 5' (left) end of reads
    inputBinding:
      position: 1
      prefix: '-5'

  3p:
    type: int?
    doc: |
      --trim3 <int>   trim <int> bases from 3' (right) end of reads
    inputBinding:
      position: 1
      prefix: '-3'
  phred33-quals:
    type:
      - "null"
      - boolean
    doc: "input quals are Phred+33 (default)\n"
    inputBinding:
      position: 1
      prefix: '--phred33-quals'
  phred64-quals:
    type:
      - "null"
      - boolean
    doc: "input quals are Phred+64 (same as --solexa1.3-quals)\n"
    inputBinding:
      position: 1
      prefix: '--phred64-quals'
  solexa-quals:
    type:
      - "null"
      - boolean
    doc: "input quals are from GA Pipeline ver. < 1.3\n"
    inputBinding:
      position: 1
      prefix: '--solexa-quals'
  solexa1.3-quals:
    type:
      - "null"
      - boolean
    doc: "input quals are from GA Pipeline ver. >= 1.3\n"
    inputBinding:
      position: 1
      prefix: '--solexa1.3-quals'
  integer-quals:
    type:
      - "null"
      - boolean
    doc: "qualities are given as space-separated integers (not ASCII)\n"
    inputBinding:
      position: 1
      prefix: '--integer-quals'
  large-index:
    type:
      - "null"
      - boolean
    doc: "force usage of a 'large' index, even if a small one is present\nAlignment:\n"
    inputBinding:
      position: 1
      prefix: '--large-index'
  v:
    type:
      - "null"
      - int
    doc: >
      <int>           report end-to-end hits w/ <=v mismatches; ignore qualities

      or
    inputBinding:
      position: 1
      prefix: '-v'
  n:
    type:
      - "null"
      - int
    doc: |
      --seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)
    inputBinding:
      position: 1
      prefix: '-n'
  e:
    type:
      - "null"
      - int
    doc: >
      --maqerr <int>  max sum of mismatch quals across alignment for -n (def:
      70)
    inputBinding:
      position: 1
      prefix: '-e'
  l:
    type:
      - "null"
      - int
    doc: |
      --seedlen <int> seed length for -n (default: 28)
    inputBinding:
      position: 1
      prefix: '-l'
  nomaqround:
    type:
      - "null"
      - boolean
    doc: "disable Maq-like quality rounding for -n (nearest 10 <= 30)\n"
    inputBinding:
      position: 1
      prefix: '--nomaqround'
  I:
    type:
      - "null"
      - int
    doc: |
      --minins <int>  minimum insert size for paired-end alignment (default: 0)
    inputBinding:
      position: 1
      prefix: '-I'
  X:
    type:
      - "null"
      - int
    doc: >
      --maxins <int>  maximum insert size for paired-end alignment (default:
      250)
    inputBinding:
      position: 1
      prefix: '-X'
  fr:
    type:
      - "null"
      - boolean
    doc: |
      --rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)
    inputBinding:
      position: 1
      prefix: '--fr'
  nofw:
    type:
      - "null"
      - boolean
    doc: |
      --norc      do not align to forward/reverse-complement reference strand
    inputBinding:
      position: 1
      prefix: '--nofw'
  maxbts:
    type:
      - "null"
      - int
    doc: |
      <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)
    inputBinding:
      position: 1
      prefix: '--maxbts'
  pairtries:
    type:
      - "null"
      - int
    doc: |
      <int>  max # attempts to find mate for anchor hit (default: 100)
    inputBinding:
      position: 1
      prefix: '--pairtries'
  y:
    type:
      - "null"
      - boolean
    doc: >
      --tryhard       try hard to find valid alignments, at the expense of speed
    inputBinding:
      position: 1
      prefix: '-y'
  chunkmbs:
    type:
      - "null"
      - int
    doc: |
      <int>   max megabytes of RAM for best-first search frames (def: 64)
      Reporting:
    inputBinding:
      position: 1
      prefix: '--chunkmbs'
  k:
    type:
      - "null"
      - int
    doc: |
      <int>           report up to <int> good alignments per read (default: 1)
    inputBinding:
      position: 1
      prefix: '-k'
  a:
    type:
      - "null"
      - boolean
    doc: |
      --all           report all alignments per read (much slower than low -k)
    inputBinding:
      position: 1
      prefix: '-a'
  m:
    type:
      - "null"
      - int
    doc: |
      <int>           suppress all alignments if > <int> exist (def: no limit)
    inputBinding:
      position: 1
      prefix: '-m'
  M:
    type:
      - "null"
      - int
    doc: >
      <int>           like -m, but reports 1 random hit (MAPQ=0); requires
      --best
    inputBinding:
      position: 1
      prefix: '-M'
  best:
    type:
      - "null"
      - boolean
    doc: "hits guaranteed best stratum; ties broken by quality\n"
    inputBinding:
      position: 1
      prefix: '--best'
  strata:
    type:
      - "null"
      - boolean
    doc: "hits in sub-optimal strata aren't reported (requires --best)\nOutput:\n"
    inputBinding:
      position: 1
      prefix: '--strata'
  t:
    type:
      - "null"
      - boolean
    doc: |
      --time          print wall-clock time taken by search phases
    inputBinding:
      position: 1
      prefix: '-t'
  B:
    type:
      - "null"
      - int
    doc: |
      --offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)
    inputBinding:
      position: 1
      prefix: '-B'
  quiet:
    type:
      - "null"
      - boolean
    doc: "print nothing but the alignments\n"
    inputBinding:
      position: 1
      prefix: '--quiet'
  refout:
    type:
      - "null"
      - boolean
    doc: "write alignments to files refXXXXX.map, 1 map per reference\n"
    inputBinding:
      position: 1
      prefix: '--refout'
  refidx:
    type:
      - "null"
      - boolean
    doc: "refer to ref. seqs by 0-based index rather than name\n"
    inputBinding:
      position: 1
      prefix: '--refidx'
  al:
    type:
      - "null"
      - boolean
    doc: |
      <fname>       write aligned reads/pairs to file(s) <fname>
    inputBinding:
      position: 1
      prefix: '--al'
  un:
    type:
      - "null"
      - boolean
    doc: |
      <fname>       write unaligned reads/pairs to file(s) <fname>
    inputBinding:
      position: 1
      prefix: '--un'
  max:
    type:
      - "null"
      - boolean
    doc: |
      <fname>      write reads/pairs over -m limit to file(s) <fname>
    inputBinding:
      position: 1
      prefix: '--max'
  suppress:
    type:
      - "null"
      - boolean
    doc: |
      <cols>  suppresses given columns (comma-delim'ed) in default output
    inputBinding:
      position: 1
      prefix: '--suppress'
  fullref:
    type:
      - "null"
      - boolean
    doc: "write entire ref name (default: only up to 1st space)\nColorspace:\n"
    inputBinding:
      position: 1
      prefix: '--fullref'
  snpphred:
    type:
      - "null"
      - int
    doc: |
      <int>   Phred penalty for SNP when decoding colorspace (def: 30)
      or
    inputBinding:
      position: 1
      prefix: '--snpphred'
  snpfrac:
    type:
      - "null"
      - boolean
    doc: |
      <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred
    inputBinding:
      position: 1
      prefix: '--snpfrac'
  col-cseq:
    type:
      - "null"
      - boolean
    doc: "print aligned colorspace seqs as colors, not decoded bases\n"
    inputBinding:
      position: 1
      prefix: '--col-cseq'
  col-cqual:
    type:
      - "null"
      - boolean
    doc: "print original colorspace quals, not decoded quals\n"
    inputBinding:
      position: 1
      prefix: '--col-cqual'
  col-keepends:
    type:
      - "null"
      - boolean
    doc: "keep nucleotides at extreme ends of decoded alignment\nSAM:\n"
    inputBinding:
      position: 1
      prefix: '--col-keepends'
  sam:
    type:
      - "null"
      - boolean
    doc: |
      --sam           write hits in SAM format
    inputBinding:
      position: 1
      prefix: '-S'
  mapq:
    type:
      - "null"
      - int
    doc: |
      <int>       default mapping quality (MAPQ) to print for SAM alignments
    inputBinding:
      position: 1
      prefix: '--mapq'
  sam-nohead:
    type:
      - "null"
      - boolean
    doc: "supppress header lines (starting with @) for SAM output\n"
    inputBinding:
      position: 1
      prefix: '--sam-nohead'
  sam-nosq:
    type:
      - "null"
      - boolean
    doc: "supppress @SQ header lines for SAM output\n"
    inputBinding:
      position: 1
      prefix: '--sam-nosq'
  sam-RG:
    type:
      - "null"
      - string
    doc: |
      <text>    add <text> (usually "lab=value") to @RG line of SAM header
      Performance:
    inputBinding:
      position: 1
      prefix: '--sam-RG'
  o:
    type:
      - "null"
      - int
    doc: |
      --offrate <int> override offrate of index; must be >= index's offrate
    inputBinding:
      position: 1
      prefix: '-o'
  p:
    type:
      - "null"
      - int
    doc: |
      --threads <int> number of alignment threads to launch (default: 1)
    inputBinding:
      position: 1
      prefix: '-p'
  mm:
    type:
      - "null"
      - boolean
    doc: "use memory-mapped I/O for index; many 'bowtie's can share\n"
    inputBinding:
      position: 1
      prefix: '--mm'
  shmem:
    type:
      - "null"
      - boolean
    doc: "use shared mem for index; many 'bowtie's can share\nOther:\n"
    inputBinding:
      position: 1
      prefix: '--shmem'
  seed:
    type:
      - "null"
      - int
    doc: |
      <int>       seed for random number generator
    inputBinding:
      position: 1
      prefix: '--seed'
  verbose:
    type:
      - "null"
      - boolean
    doc: "verbose output (for debugging)\n"
    inputBinding:
      position: 1
      prefix: '--verbose'

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.filename)

  output_bowtie_log:
    type: File
    outputBinding:
      glob: $(inputs.filename + ".log")

baseCommand:
  - bowtie

arguments:
  - valueFrom: |
      ${
        if (inputs.filelist && inputs.filelist_mates){
          return "-1";
        }
        return null;
      }
    position: 82
  - valueFrom: |
      ${
        if (inputs.filelist && inputs.filelist_mates){
          return "-2";
        }
        return null;
      }
    position: 84
  - valueFrom: $('2> ' + inputs.filename + '.log')
    position: 100000
    shellQuote: false

doc: |
  bowtie.cwl is developed for CWL consortium

  Usage: 
  bowtie [options]* <ebwt> {-1 <m1> -2 <m2> | --12 <r> | <s>} [<hit>]

    <m1>    Comma-separated list of files containing upstream mates (or the
            sequences themselves, if -c is set) paired with mates in <m2>
    <m2>    Comma-separated list of files containing downstream mates (or the
            sequences themselves if -c is set) paired with mates in <m1>
    <r>     Comma-separated list of files containing Crossbow-style reads.  Can be
            a mixture of paired and unpaired.  Specify "-"for stdin.
    <s>     Comma-separated list of files containing unpaired reads, or the
            sequences themselves, if -c is set.  Specify "-"for stdin.
    <hit>   File to write hits to (default: stdout)
  Input:
    -q                 query input files are FASTQ .fq/.fastq (default)
    -f                 query input files are (multi-)FASTA .fa/.mfa
    -r                 query input files are raw one-sequence-per-line
    -c                 query sequences given on cmd line (as <mates>, <singles>)
    -C                 reads and index are in colorspace
    -Q/--quals <file>  QV file(s) corresponding to CSFASTA inputs; use with -f -C
    --Q1/--Q2 <file>   same as -Q, but for mate files 1 and 2 respectively
    -s/--skip <int>    skip the first <int> reads/pairs in the input
    -u/--qupto <int>   stop after first <int> reads/pairs (excl. skipped reads)
    -5/--trim5 <int>   trim <int> bases from 5' (left) end of reads
    -3/--trim3 <int>   trim <int> bases from 3' (right) end of reads
    --phred33-quals    input quals are Phred+33 (default)
    --phred64-quals    input quals are Phred+64 (same as --solexa1.3-quals)
    --solexa-quals     input quals are from GA Pipeline ver. < 1.3
    --solexa1.3-quals  input quals are from GA Pipeline ver. >= 1.3
    --integer-quals    qualities are given as space-separated integers (not ASCII)
    --large-index      force usage of a 'large' index, even if a small one is present
  Alignment:
    -v <int>           report end-to-end hits w/ <=v mismatches; ignore qualities
      or
    -n/--seedmms <int> max mismatches in seed (can be 0-3, default: -n 2)
    -e/--maqerr <int>  max sum of mismatch quals across alignment for -n (def: 70)
    -l/--seedlen <int> seed length for -n (default: 28)
    --nomaqround       disable Maq-like quality rounding for -n (nearest 10 <= 30)
    -I/--minins <int>  minimum insert size for paired-end alignment (default: 0)
    -X/--maxins <int>  maximum insert size for paired-end alignment (default: 250)
    --fr/--rf/--ff     -1, -2 mates align fw/rev, rev/fw, fw/fw (default: --fr)
    --nofw/--norc      do not align to forward/reverse-complement reference strand
    --maxbts <int>     max # backtracks for -n 2/3 (default: 125, 800 for --best)
    --pairtries <int>  max # attempts to find mate for anchor hit (default: 100)
    -y/--tryhard       try hard to find valid alignments, at the expense of speed
    --chunkmbs <int>   max megabytes of RAM for best-first search frames (def: 64)
  Reporting:
    -k <int>           report up to <int> good alignments per read (default: 1)
    -a/--all           report all alignments per read (much slower than low -k)
    -m <int>           suppress all alignments if > <int> exist (def: no limit)
    -M <int>           like -m, but reports 1 random hit (MAPQ=0); requires --best
    --best             hits guaranteed best stratum; ties broken by quality
    --strata           hits in sub-optimal strata aren't reported (requires --best)
  Output:
    -t/--time          print wall-clock time taken by search phases
    -B/--offbase <int> leftmost ref offset = <int> in bowtie output (default: 0)
    --quiet            print nothing but the alignments
    --refout           write alignments to files refXXXXX.map, 1 map per reference
    --refidx           refer to ref. seqs by 0-based index rather than name
    --al <fname>       write aligned reads/pairs to file(s) <fname>
    --un <fname>       write unaligned reads/pairs to file(s) <fname>
    --max <fname>      write reads/pairs over -m limit to file(s) <fname>
    --suppress <cols>  suppresses given columns (comma-delim'ed) in default output
    --fullref          write entire ref name (default: only up to 1st space)
  Colorspace:
    --snpphred <int>   Phred penalty for SNP when decoding colorspace (def: 30)
       or
    --snpfrac <dec>    approx. fraction of SNP bases (e.g. 0.001); sets --snpphred
    --col-cseq         print aligned colorspace seqs as colors, not decoded bases
    --col-cqual        print original colorspace quals, not decoded quals
    --col-keepends     keep nucleotides at extreme ends of decoded alignment
  SAM:
    -S/--sam           write hits in SAM format
    --mapq <int>       default mapping quality (MAPQ) to print for SAM alignments
    --sam-nohead       supppress header lines (starting with @) for SAM output
    --sam-nosq         supppress @SQ header lines for SAM output
    --sam-RG <text>    add <text> (usually "lab=value") to @RG line of SAM header
  Performance:
    -o/--offrate <int> override offrate of index; must be >= index's offrate
    -p/--threads <int> number of alignment threads to launch (default: 1)
    --mm               use memory-mapped I/O for index; many 'bowtie's can share
    --shmem            use shared mem for index; many 'bowtie's can share
  Other:
    --seed <int>       seed for random number generator
    --verbose          verbose output (for debugging)
    --version          print version information and quit
    -h/--help          print this usage message


$namespaces:
  schema: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

schema:mainEntity:
#  $import: https://scidap.com/description/tools/bowtie.yaml
  class: schema:SoftwareSourceCode
  schema:name: "bowtie"
  schema:about: >
    Bowtie is an ultrafast, memory-efficient short read aligner.
    It aligns short DNA sequences (reads) to the human genome at a rate of over 25 million 35-bp reads per hour.
    Bowtie indexes the genome with a Burrows-Wheeler index to keep its memory footprint small: typically about 2.2 GB for the human genome (2.9 GB for paired-end).
  schema:url: http://bowtie-bio.sourceforge.net
  schema:codeRepository: https://github.com/BenLangmead/bowtie.git

  schema:license:
  - https://opensource.org/licenses/GPL-3.0

  schema:targetProduct:
    class: schema:SoftwareApplication
    schema:softwareVersion: "1.1.2"
    schema:applicationCategory: "commandline tool"

  schema:programmingLanguage: "C++"

  schema:publication:
  - class: schema:ScholarlyArticle
    id: http://dx.doi.org/10.1186/gb-2009-10-3-r25

schema:isPartOf:
  class: schema:CreativeWork
  schema:name: "Common Workflow Language"
  schema:url: http://commonwl.org/

schema:author:
  class: schema:Person
  schema:name: "Andrey Kartashov"
  schema:email: mailto:Andrey.Kartashov@cchmc.org
  schema:sameAs:
  - id: http://orcid.org/0000-0001-9102-5681
  schema:worksFor:
  - class: schema:Organization
    schema:name: "Cincinnati Children's Hospital Medical Center"
    schema:location: "3333 Burnet Ave, Cincinnati, OH 45229-3026"
    schema:department:
    - class: schema:Organization
      schema:name: "Barski Lab"
