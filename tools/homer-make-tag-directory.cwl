cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InitialWorkDirRequirement
  listing: |
    ${
      return [
        {"class": "Directory",
         "basename": "default",
         "listing": [inputs.bam_file],
         "writable": true}
      ]
    }
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_folder = function() {
      if (inputs.output_folder){
        return inputs.output_folder.replace(/\t|\s|\[|\]|\>|\<|,|\./g, "_");
      } else {
        return inputs.bam_file.basename.split('.')[0];
      }
    };

hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/homer:v0.0.2


inputs:

  bam_file:
    type: File
    doc: "Alignment file, BAM"

  output_folder:
    type: string?
    doc: "Name of the directory to save outputs"

  fragment_size:
    type:
      - "null"
      - int
      - string
    inputBinding:
      position: 5
      prefix: "-fragLength"
    doc: |
      Set fragment size.
      By default is estimated as if it was single end ChIP-Seq experiment.
      Possible values:
        "#" - int value to be used as fragment size
        "given" - use read lengths
        "pe" - calculate from paired end read coordinates

  total_reads:
    type:
      - "null"
      - int
      - string
    inputBinding:
      position: 6
      prefix: "-totalReads"
    doc: |
      Set total reads number for downstream normalization.
      Default: autocalculated, equal to uniquely mapped reads number
      Possible values:
        "#" - int value to be used as total reads number
        "all" - autocalculated, equal to uniquely + multi mapped reads number

  min_length:
    type: int?
    inputBinding:
      position: 7
      prefix: "-minlen"
    doc: |
      Discard reads smaller then

  max_length:
    type: int?
    inputBinding:
      position: 8
      prefix: "-maxlen"
    doc: |
      Discard reads bigger then


outputs:

  output_tag_folder:
    type: Directory
    outputBinding:
      glob: $(default_output_folder())
    doc: "Tag directory"



baseCommand: ["makeTagDirectory"]
arguments:
  - valueFrom: $(default_output_folder())
  - valueFrom: $("default/" + inputs.bam_file.basename)


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/homer-metadata.yaml

s:name: "homer-make-tag-directory"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/homer-make-tag-directory.cwl
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
  Tool runs makeTagDirectory that basically parses through the alignment file and splits the tags into separate
  files based on the chromosomes.

  Multiple alignment files are not supported. Alignment file's format is restricted to be only BAM.

  Output is placed in a folder with the name derived from the input BAM file's basename.

  Skipped arguments:

    Rely on the default value:
      -format           - format will be autodetected
      -precision        - the default value is used

    Not required general functionality:
      -d
      -single
      -force5th
      -t
      -flip
      -tbp

    Not required GC-bias options:
      -genome
      -checkGC
      -normGC
      -normFixedOligo
      -minNormRatio
      -maxNormRatio
      -iterNorm
      -filterReads

    Not required HiC options:
      -removePEbg
      -restrictionSite
      -removeSpikes
      -bowtiePE
      -directional

s:about: |
  Usage:
    makeTagDirectory <directory> <alignment file 1> [file 2] ... [options]

        Creates a platform-independent 'tag directory' for later analysis.
        Currently BED, eland, bowtie, and sam files are accepted. The program will try to
        automatically detect the alignment format if not specified.  Program will also
        unzip *.gz, *.bz2, and *.zip files and convert *.bam to sam files on the fly
        Existing tag directories can be added or combined to make a new one using -d/-t
        If more than one format is needed and the program cannot auto-detect it properly,
        make separate tag directories by running the program separately, then combine them.
        To perform QC/manipulations on an existing tag directory, add "-update"

        Options:
            -fragLength <# | given | pe> (Set estimated fragment length or use PE length - given: use read lengths)
                By default treats the sample as a single read ChIP-Seq experiment
            -format <X> where X can be: (with column specifications underneath)
                bed - BED format files:
                    (1:chr,2:start,3:end,4:+/- or read name,5:# tags,6:+/-)
                    -force5th (5th column of BED file contains # of reads mapping to position)
                sam - SAM formatted files (use samTools to covert BAMs into SAM if you have BAM)
                    -unique (keep if there is a single best alignment based on mapq)
                        -mapq <#> (Minimum mapq for -unique, default: 10, set negative to use AS:i:/XS:i:)
                    -keepOne (keep one of the best alignments even if others exist)
                    -keepAll (include all alignments in SAM file)
                    -mis (Maximum allowed mismatches, default: no limit, uses MD:Z: tag)
                    -sspe (strand specific, paired-end reads[flips strand of 2nd read to match])
                    -read1/-read2 (only analyze 1st or 2nd read for PE sequencing)
                bowtie - output from bowtie (run with --best -k 2 options)
                    (1:read name,2:+/-,3:chr,4:position,5:seq,6:quality,7:NA,8:misInfo)
                eland_result - output from basic eland
                    (1:read name,2:seq,3:code,4:#zeroMM,5:#oneMM,6:#twoMM,7:chr,
                                8:position,9:F/R,10-:mismatches
                eland_export - output from illumina pipeline (22 columns total)
                    (1-5:read name info,9:sequence,10:quality,11:chr,13:position,14:strand)
                eland_extended - output from illumina pipeline (4 columns total)
                    (1:read name,2:sequence,3:match stats,4:positions[,])
                mCpGbed - encode style mCpG reporting in extended BED format, no auto-detect
                    (1:chr,2:start,3:end,4:name,5:,6:+/-,7:,8:,9:,10:#C,11:#mC)
                allC - Lister style output files detailing the read information about all cytosines
                    (1:chr,2:pos,3:strand,4:context,#mC,#totalC,#unmC
                bismark - Bismark style output files detailing the read information about all cytosines
                    (1:chr,2:pos,3:strand,4:#mC,5:#unmC,6:context,7:triseq
                    -minCounts <#> (minimum number of reads to report mC/C ratios, default: 10)
                    -mCcontext <CG|CHG|CHH|all> (only use C's in this context, default: CG)
                HiCsummary - minimal paired-end read mapping information
                    (1:readname,2:chr1,3:5'pos1,4:strand1,5:chr2,6:5'pos2,7:strand2)
            -flip (flip strand of each read, i.e. might want to use with some RNA-seq)
            -totalReads <#|all|default> (set the effective total number of reads - all includes multimappers)
            -force5th (5th column of BED file contains # of reads mapping to position)
            -d <tag directory> [tag directory 2] ... (add Tag directory to new tag directory)
            -t <tag file> [tag file 2] ... (add tag file i.e. *.tags.tsv to new tag directory)
            -single (Create a single tags.tsv file for all "chromosomes" - i.e. if >100 chromosomes)
            -update (Use current tag directory for QC/processing, do not parse new alignment files)
            -tbp <#> (Maximum tags per bp, default: no maximum)
            -precision <1|2|3> (number of decimal places to use for tag totals, default: 1)
            -minlen <#> and -maxlen <#> (Filter reads with lengths outside this range)

            GC-bias options:
            -genome <genome version> (To see available genomes, use "-genome list")
                -or- (for custom genomes):
            -genome <path-to-FASTA file or directory of FASTA files>

            -checkGC (check Sequence bias, requires "-genome")
                -freqStart <#> (offset to start calculating frequency, default: -50)
                -freqEnd <#> (distance past fragment length to calculate frequency, default: +50)
                -oligoStart <#> (oligo bias start)
                -oligoEnd <#> (oligo bias end)
            -normGC <target GC profile file> (i.e. tagGCcontent.txt file from control experiment)
                Use "-normGC default" to match the genomic GC distribution
            -normFixedOligo <oligoFreqFile> (normalize 5' end bias, "-normFixedOligo default" ok)
            -minNormRatio <#> (Minimum deflation ratio of tag counts, default: 0.25)
            -maxNormRatio <#> (Maximum inflation ratio of tag counts, default: 2.0)
            -iterNorm <#> (Sets -max/minNormRatio to 1 and 0, iteratively normalizes such that the
                resulting distrubtion is no more than #% different than target, i.e. 0.1,default: off)
            -filterReads <seq> <offset> <keep|remove> (filter reads based on oligo sequence in the genome)

        HiC options
            -removePEbg (remove paired end tags within 1.5x fragment length on same chr)
                -PEbgLength <#> (remove PE  reads facing on another within this distance, default: 1.5x fragLen)
            -restrictionSite <seq> (i.e. AAGCTT for HindIII, assign data < 1.5x fragment length to sites)
                Must specify genome sequence directory too. (-rsmis <#> to specify mismatches, def: 0)
                -both, -one, -onlyOne, -none (Keeps reads near restriction sites, default: keep all)
                -removeSelfLigation (removes reads linking same restriction fragment)
                -removeRestrictionEnds (removes reads starting on a restriction fragment)
                -assignMidPoint (will place reads in the middle of HindIII fragments)
                -restrictionSiteLength <#> (maximum distance from restriction site, default: 1.5x fragLen)
            -removeSpikes <size bp> <#> (remove tags from regions with > than # times
                the average tags per size bp, suggest "-removeSpikes 10000 8")
            -bowtiePE (PE alignments in bowtie alignment, assumes last character of read name is 0 or 1)
                (don't need this for sam/bam files)
            -directional (matrix is not symetric, i.e. grid-seq, 1st/2nd read mean something)