cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: InlineJavascriptRequirement


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/homer:v0.0.2


inputs:

  peak_file:
    type: File
    doc: "Homer generated peak file or BED"

  tag_folders:
    type:
      - Directory
      - Directory[]
    inputBinding:
      position: 7
      prefix: "-d"
    doc: "Tag folders from homer-make-tag-directory tool"

  hist_width:
    type:
      - int
      - string
    inputBinding:
      position: 8
      prefix: "-size"
    doc: |
      Possible values:
        "#" - performs analysis on the "#" bp surrounding the peak centers
        "#,#" - performs analysis from "#" to "#" relative to peak center
        "given" - set size to actual coordinates in peak/BED file

  hist_bin_size:
    type: int
    inputBinding:
      position: 9
      prefix: "-hist"
    doc: |
      Bin size, bp. If hist_width is "given" or skipped, this
      parameter will set the number of bins to divide each region into

  export_heatmap:
    type: boolean?
    inputBinding:
      position: 10
      prefix: "-ghist"
    doc: |
      Generate heatmap. Instead of averaging all of the data
      from each peak, keep data from each peak separate

  norm_fpkm:
    type: boolean?
    inputBinding:
      position: 11
      prefix: "-fpkm"
    doc: |
      Normalize read counts to million reads or fragments per kilobase mapped

  norm_raw:
    type: boolean?
    inputBinding:
      position: 12
      prefix: "-raw"
    doc: |
      Do not adjust the tag counts based on total tags sequenced.
      By default all tag counts will be normalized to norm_tag_count

  norm_tag_count:
    type: int?
    inputBinding:
      position: 13
      prefix: "-norm"
    doc: |
      Normalize tags to this tag count, default=1e7, 0=average tag count in all directories

  norm_fragment_size:
    type: int?
    inputBinding:
      position: 14
      prefix: "-normLength"
    doc: |
      Fragment length to normlize to for experiments with different lens. Default: 100bp

  strand:
    type: string?
    inputBinding:
      position: 15
      prefix: "-strand"
    doc: |
      Count tags on specific strands relative to peak. Default: both
      Possible values: +|-

  threads:
    type: int?
    inputBinding:
      position: 16
      prefix: "-cpu"
    doc: |
      Set the number of threads. Default: 1

  histogram_filename:
    type: string
    doc: "Output histogram's filename"


outputs:

  histogram_file:
    type: stdout
    doc: "Output histogram file"


stdout: ${return inputs.histogram_filename;}


baseCommand: ["annotatePeaks.pl"]
arguments:
  - valueFrom: $(inputs.peak_file)
    position: 5
  - valueFrom: $("none")
    position: 6


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

s:mainEntity:
  $import: ./metadata/homer-metadata.yaml

s:name: "homer-annotate-peaks-hist"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/homer-annotate-peaks-hist.cwl
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
  Tool is used to produce histogram or heatmaps only. Rest of the functionality is not implemented intentionally.
  If TSS analysis needed, input peak_file should be centered on TSS, where the 'center' of the peak in the actual TSS.
  For example:
    1	chr4	978796	978796	-
    2	chr4	1052109	1052109	+
    3	chr4	1105422	1105422	-

  Skipped arguments:

    Related to peaks annotation:
      -organism
      -gtf
      -gff
      -gff3
      -gid
      -ann
      -mask
      -p
      -pdist
      -pcount
      -vcf
      -editDistance
      -individuals
      -gene
      -go
      -genomeOntology
      -ratio
      -rlog
      -vst
      -CpG
      -nfr
      -nfrSize
      -gwasCatalog
      -map
      -noann

    Related to tss/tts/rna modes:
      tss
      tts
      rna
      -list
      -cTSS

    Related to motifs:
      -m
      -mscore
      -nmotifs
      -mdist
      -mfasta
      -fm
      -rmrevopp
      -matrix
      -mbed
      -mlogic
      -norevopp

    Related to peak centering:
      -center
      -mirror
      -multi

    Related to genome comparisons
      -cmpGenome
      -cmpLiftover

    Currently not needed functionality:
      -bedGraph
      -wig
      -nuc
      -di
      -histNorm
      -rm
      -log
      -sqrt
      -len
      -pc
      -noblanks
      -homer1
      -homer2

s:about: |
  Usage: annotatePeaks.pl <peak file | tss> <genome version>  [additional options...]

      Available Genomes (required argument): (name,org,directory,default promoter set)
              -- or --
          Custom: provide the path to genome FASTA files (directory or single file)
          If no genome is available, specify 'none'.
          If using FASTA file or none, may want to specify '-organism <...>'

      User defined annotation files (default is UCSC refGene annotation):
          annotatePeaks.pl accepts GTF (gene transfer formatted) files to annotate positions relative
          to custom annotations, such as those from de novo transcript discovery or Gencode.
          -gtf <gtf format file> (Use -gff and -gff3 if appropriate, but GTF is better)
          -gid (by default the GTF file is processed by transcript_id, use this option for gene_id)
          -ann <custom homer annotation file> (created by assignGenomeAnnotation, see website)

      Peak vs. tss/tts/rna mode (works with custom GTF file):
          If the first argument is "tss" (i.e. annotatePeaks.pl tss hg18 ...) then a TSS centric
          analysis will be carried out.  Tag counts and motifs will be found relative to the TSS.
          (no position file needed) ["tts" now works too - e.g. 3' end of gene]
          ["rna" specifies gene bodies, will automaticall set "-size given"]
          NOTE: The default TSS peak size is 4000 bp, i.e. +/- 2kb (change with -size option)
          -list <gene id list> (subset of genes to perform analysis [unigene, gene id, accession,
               probe, etc.], default = all promoters)
          -cTSS <promoter position file i.e. peak file> (should be centered on TSS)

      Primary Annotation Options:
          -mask (Masked repeats, can also add 'r' to end of genome name)
          -m <motif file 1> [motif file 2] ... (list of motifs to find in peaks)
              -mscore (reports the highest log-odds score within the peak)
              -nmotifs (reports the number of motifs per peak)
              -mdist (reports distance to closest motif)
              -mfasta <filename> (reports sites in a fasta file - for building new motifs)
              -fm <motif file 1> [motif file 2] (list of motifs to filter from above)
              -rmrevopp <#> (only count sites found within <#> on both strands once, i.e. palindromic)
              -matrix <prefix> (outputs a motif co-occurrence files:
                  prefix.count.matrix.txt - number of peaks with motif co-occurrence
                  prefix.ratio.matrix.txt - ratio of observed vs. expected  co-occurrence
                  prefix.logPvalue.matrix.txt - co-occurrence enrichment
                  prefix.stats.txt - table of pair-wise motif co-occurrence statistics
                  additional options:
                  -matrixMinDist <#> (minimum distance between motif pairs - to avoid overlap, default: 4)
                  -matrixMaxDist <#> (maximum distance between motif pairs)
              -mbed <filename> (Output motif positions to a BED file to load at UCSC (or -mpeak))
              -mlogic <filename> (will output stats on common motif orientations)
          -d <tag directory 1> [tag directory 2] ... (list of experiment directories to show
              tag counts for) NOTE: -dfile <file> where file is a list of directories in first column
          -bedGraph <bedGraph file 1> [bedGraph file 2] ... (read coverage counts from bedGraph files)
          -wig <wiggle file 1> [wiggle file 2] ... (read coverage counts from wiggle files)
          -p <peak file> [peak file 2] ... (to find nearest peaks)
              -pdist to report only distance (-pdist2 gives directional distance)
              -pcount to report number of peaks within region
          -vcf <VCF file> (annotate peaks with genetic variation infomation, one col per individual)
              -editDistance (Computes the # bp changes relative to reference)
              -individuals <name1> [name2] ... (restrict analysis to these individuals)
              -editDistance (Computes the # bp changes relative to reference)
              -individuals <name1> [name2] ... (restrict analysis to these individuals)
          -gene <data file> ... (Adds additional data to result based on the closest gene.
              This is useful for adding gene expression data.  The file must have a header,
              and the first column must be a GeneID, Accession number, etc.  If the peak
              cannot be mapped to data in the file then the entry will be left empty.
          -go <output directory> (perform GO analysis using genes near peaks)
          -genomeOntology <output directory> (perform genomeOntology analysis on peaks)
              -gsize <#> (Genome size for genomeOntology analysis, default: 2e9)

      Annotation vs. Histogram mode:
          -hist <bin size in bp> (i.e 1, 2, 5, 10, 20, 50, 100 etc.)
          The -hist option can be used to generate histograms of position dependent features relative
          to the center of peaks.  This is primarily meant to be used with -d and -m options to map
          distribution of motifs and ChIP-Seq tags.  For ChIP-Seq peaks for a Transcription factor
          you might want to use the -center option (below) to center peaks on the known motif
          ** If using "-size given", histogram will be scaled to each region (i.e. 0-100%), with
          the -hist parameter being the number of bins to divide each region into.
              Histogram Mode specific Options:
              -nuc (calculated mononucleotide frequencies at each position,
                  Will report by default if extracting sequence for other purposes like motifs)
              -di (calculated dinucleotide frequencies at each position)
              -histNorm <#> (normalize the total tag count for each region to 1, where <#> is the
                  minimum tag total per region - use to avoid tag spikes from low coverage
              -ghist (outputs profiles for each gene, for peak shape clustering)
              -rm <#> (remove occurrences of same motif that occur within # bp)

      Peak Centering: (other options are ignored)
          -center <motif file> (This will re-center peaks on the specified motif, or remove peak
              if there is no motif in the peak.  ONLY recentering will be performed, and all other
              options will be ignored.  This will output a new peak file that can then be reanalyzed
              to reveal fine-grain structure in peaks (It is advised to use -size < 200) with this
              to keep peaks from moving too far (-mirror flips the position)
          -multi (returns genomic positions of all sites instead of just the closest to center)

      Genome comparisons (need genome & liftOver)
          -cmpGenome <genome1> [genome2] (Genomes to compare for sequence/motifs)
          -cmpLiftover <liftover1> [genome2] (Genomes to compare for sequence/motifs)

      Normalization options:
          -fpkm (normalize read counts to million reads or fragments per kilobase mapped)
          -raw (do not adjust the tag counts based on total tags sequenced, -noadj works too)
          -norm <#> (normalize tags to this tag count, default=1e7, 0=average tag count in all directories)
          -normLength <#> (Fragment length to normlize to for experiments with different lens, def: 100)
          -log (output tag counts as log2(x+1+rand) values - for scatter plots)
          -sqrt (output tag counts as sqrt(x+rand) values - for scatter plots)
          -ratio (process tag values as ratios - i.e. chip-seq, or mCpG/CpG)

      Advanced normalization options: (-rlog and -vst require R and DESeq2 to be installed)
          -rlog (quantile/variance normalization on reported genes using DESeq2 rlog funcition, needs R)
          -vst (quantile/variance normalization on reported genes using DESeq2 vst function, needs R)

      Advanced Options:
          -len <#> / -fragLength <#> (Fragment length, default=auto, might want to set to 1 for 5'RNA)
          -size <#> (Peak size[from center of peak], default=inferred from peak file)
              -size #,# (i.e. -size -10,50 count tags from -10 bp to +50 bp from center)
              -size "given" (count tags etc. using the actual regions - for variable length regions)
          -strand <+|-|both> (Count tags on specific strands relative to peak, default: both)
          -pc <#> (maximum number of tags to count per bp, default=0 [no maximum], -tbp <#> works too)
          -CpG (Calculate CpG/GC content)
          -nfr (report nuclesome free region scores instead of tag counts, also -nfrSize <#>)
          -norevopp (do not search for motifs on the opposite strand [works with -center too])
          -gwasCatalog <gwasCatalog file from UCSC> (list overlapping GWAS risk SNPs)
          -pdist (only report distance to nearest peak using -p, not peak name)
          -map <mapping file> (mapping between peak IDs and promoter IDs, overrides closest assignment)
          -noann, -nogene (skip genome annotation step, skip TSS annotation)
          -homer1/-homer2 (by default, the new version of homer [-homer2] is used for finding motifs)
          -cpu <#> (Number of processors to use when possible - only some parts utilize multiple cores)
          -noblanks (remove peaks/rows with missing data)
