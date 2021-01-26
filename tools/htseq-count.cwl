cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: quay.io/biocontainers/htseq:0.13.5--py38h1773678_0


inputs:

  script:
    type: string?
    default: |
      #!/bin/bash
      shopt -s nocaseglob
      set -- "$0" "$@"
      echo "Run htseq-count with the following parameters"
      for i in "$@";
        do echo $i;
      done;
      htseq-count -f bam -r pos "$@" 1> temp.tsv 2> feature_counts_stderr.log
      cat temp.tsv | grep "^__" > feature_counts_stdout.log
      cat temp.tsv | grep -v "^__" > feature_counts_report.tsv
      rm -f temp.tsv
    inputBinding:
      position: 1
    doc: |
      Script runs htseq-count with the provided parameters splitting
      the report file into summary and the actual counts data.

  alignment_bam_file:
    type: File
    secondaryFiles:
    - .bai
    inputBinding:
      position: 5
    doc: |
      Path to the coordinate sorted indexed BAM file

  annotation_gtf_file:
    type: File
    inputBinding:
      position: 6
    doc: |
      GTF annotation file

  strand_specific:
    type:
    - "null"
    - type: enum
      symbols:
      - "yes"
      - "no"
      - "reverse"
    inputBinding:
      position: 7
      prefix: "-s"
    doc: |
      Whether the data is from a strand-specific assay. For stranded=no, a read is
      considered overlapping with a feature regardless of whether it is mapped to
      the same or the opposite strand as the feature. For stranded=yes and single-end
      reads, the read has to be mapped to the same strand as the feature. For paired-end
      reads, the first read has to be on the same strand and the second read on the
      opposite strand. For stranded=reverse, these rules are reversed.
      Default: "yes"

  feature_type:
    type:
    - "null"
    - type: enum
      symbols:
      - "CDS"
      - "exon"
      - "start_codon"
      - "stop_codon"
      - "transcript"
    inputBinding:
      position: 8
      prefix: "-t"
    doc: |
      Feature type (3rd column in GFF file) to be used,
      all features of other type are ignored.
      Default: exon

  feature_id:
    type:
    - "null"
    - type: enum
      symbols:
      - "gene_id"
      - "transcript_id"
      - "exon_id"
      - "gene_name"
    inputBinding:
      position: 9
      prefix: "-i"
    doc: |
      GTF attribute to be used as feature ID. Several GTF lines with the same feature
      ID will be considered as parts of the same feature. The feature ID is used to
      identity the counts in the output table.
      FYI: no use of having here "exon_number", so it was excluded from enum
      Default: gene_id

  additional_id:
    type:
    - "null"
    - type: enum
      symbols:
      - "gene_id"
      - "transcript_id"
      - "exon_number"
      - "exon_id"
      - "gene_name"
    inputBinding:
      position: 10
      prefix: "--additional-attr"
    doc: |
      Additional feature attributes, which will be printed as an additional
      column after the primary attribute column but before the counts column.
      Default: none

  overlapping_mode:
    type:
    - "null"
    - type: enum
      symbols:
      - "union"
      - "intersection-strict"
      - "intersection-nonempty"
    inputBinding:
      position: 11
      prefix: "-m"
    doc: |
      Mode to handle reads overlapping more than one feature
      Default: union

  nonunique_mode:
    type:
    - "null"
    - type: enum
      symbols:
      - "none"
      - "all"
      - "fraction"
      - "random"
    inputBinding:
      position: 12
      prefix: "--nonunique"
    doc: |
      Mode to handle reads that align to or are assigned to more than
      one feature in the selected overlapping_mode
      Options:
        "none": the read (or read pair) is counted as ambiguous and not
          counted for any features. Also, if the read (or read pair) aligns
          to more than one location in the reference, it is scored as
          alignment_not_unique.
        "all": the read (or read pair) is counted as ambiguous and is also
          counted in all features to which it was assigned. Also, if the read
          (or read pair) aligns to more than one location in the reference, it
          is scored as alignment_not_unique and also separately for each location.
        "fraction": the read (or read pair) is counted as ambiguous and is also
          counted fractionally in all features to which it was assigned. For
          example, if the read overlaps with 3 features, it will be counted 1/3
          to each of them.
        "random": the read (or read pair) is counted as ambiguous and is also
          counted uniformly at random to one of the features to which it
          was assigned.
      Default: "none"

  secondary_alignments_mode:
    type:
    - "null"
    - type: enum
      symbols:
      - "score"
      - "ignore"
    inputBinding:
      position: 13
      prefix: "--secondary-alignments"
    doc: |
      Mode to handle secondary alignments (SAM flag 0x100).
      Can be "score" and "ignore"
      Default: "score"

  supplementary_alignments_mode:
    type:
    - "null"
    - type: enum
      symbols:
      - "score"
      - "ignore"
    inputBinding:
      position: 14
      prefix: "--supplementary-alignments"
    doc: |
      Mode to handle supplementary/chimeric alignments (SAM flag 0x800).
      Can be "score" and "ignore"
      Default: "score"


outputs:
  
  feature_counts_report_file:
    type: File
    outputBinding:
      glob: "feature_counts_report.tsv"

  stdout_log:
    type: File
    outputBinding:
      glob: "feature_counts_stdout.log"

  stderr_log:
    type: File
    outputBinding:
      glob: "feature_counts_stderr.log"


baseCommand: ["bash", "-c"]


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

label: "HTSeq: Analysing high-throughput sequencing data"
s:name: "HTSeq: Analysing high-throughput sequencing data"
s:alternateName: "HTSeq: Analysing high-throughput sequencing data"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/htseq-count.cwl
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
  For convenience to use in the workflow that sort and index BAM files by coordinate
  this tools expects coordinate sorted and indexed BAM file as input. For single-read
  dat it won't influence on anything, for paired-end the more memory will be used to
  keep reads while looking for their proper pairs (see --max-reads-in-buffer parameter).

  Current limitations:
    - only one `--additional-attr` is supported
    - skip `--nprocesses` parameter as it's not helpful when we use only one input BAM file

s:about: |
  usage: htseq-count [options] alignment_file gff_file

  This script takes one or more alignment files in SAM/BAM format and a feature file
  in GFF format and calculates for each feature the number of reads mapping to it. See
  http://htseq.readthedocs.io/en/master/count.html for details.

  options:
    -f <format>, --format=<format>¶
      Format of the input data. Possible values are sam (for text SAM files) and bam
      (for binary BAM files). Default is sam.

      DEPRECATED: Modern versions of samtools/htslibs, which HTSeq uses to access
      SAM/BAM/CRAM files, have automatic file type detection. This flag will be removed
      in future versions of htseq-count.

    -r <order>, --order=<order>
      For paired-end data, the alignment have to be sorted either by read name or by
      alignment position. If your data is not sorted, use the samtools sort function
      of samtools to sort it. Use this option, with name or pos for <order> to indicate
      how the input data has been sorted. The default is name.

      If name is indicated, htseq-count expects all the alignments for the reads of a
      given read pair to appear in adjacent records in the input data. For pos, this
      is not expected; rather, read alignments whose mate alignment have not yet been
      seen are kept in a buffer in memory until the mate is found. While, strictly speaking,
      the latter will also work with unsorted data, sorting ensures that most alignment
      mates appear close to each other in the data and hence the buffer is much less likely
      to overflow.

    --max-reads-in-buffer=<number>
      When <alignment_file> is paired end sorted by position, allow only so many reads to
      stay in memory until the mates are found (raising this number will use more memory).
      Has no effect for single end or paired end sorted by name. (default: 30000000)

    -s <yes/no/reverse>, --stranded=<yes/no/reverse>
      Whether the data is from a strand-specific assay (default: yes)
      For stranded=no, a read is considered overlapping with a feature regardless of whether
      it is mapped to the same or the opposite strand as the feature. For stranded=yes and
      single-end reads, the read has to be mapped to the same strand as the feature. For
      paired-end reads, the first read has to be on the same strand and the second read on
      the opposite strand. For stranded=reverse, these rules are reversed.

    -a <minaqual>, --a=<minaqual>
      Skip all reads with MAPQ alignment quality lower than the given minimum value
      (default: 10). MAPQ is the 5th column of a SAM/BAM file and its usage depends on the
      software used to map the reads.

    -t <feature type>, --type=<feature type>
      Feature type (3rd column in GTF file) to be used, all features of other type are ignored
      (default, suitable for RNA-Seq analysis using an Ensembl GTF file: exon)

    -i <id attribute>, --idattr=<id attribute>
      GTF attribute to be used as feature ID. Several GTF lines with the same feature ID will
      be considered as parts of the same feature. The feature ID is used to identity the counts
      in the output table. The default, suitable for RNA-Seq analysis using an Ensembl GTF file,
      is gene_id.

    --additional-attr=<id attributes>
      Additional feature attributes, which will be printed as an additional column after the
      primary attribute column but before the counts column(s). The default is none, a suitable
      value to get gene names using an Ensembl GTF file is gene_name. To use more than one
      additional attribute, repeat the option in the command line more than once, with a single
      attribute each time, e.g. --additional-attr=gene_name --additional_attr=exon_number.

    -m <mode>, --mode=<mode>
      Mode to handle reads overlapping more than one feature. Possible values for <mode> are union,
      intersection-strict and intersection-nonempty (default: union)

    --nonunique=<nonunique mode>
      Mode to handle reads that align to or are assigned to more than one feature in the overlap
      <mode> of choice (see -m option). <nonunique mode> are none and all (default: none)

    --secondary-alignments=<mode>
      Mode to handle secondary alignments (SAM flag 0x100). <mode> can be score and ignore
      (default: score)

    --supplementary-alignments=<mode>
      Mode to handle supplementary/chimeric alignments (SAM flag 0x800). <mode> can be score
      and ignore (default: score)

    -o <samout>, --samout=<samout>
      Write out all SAM alignment records into SAM files (one per input file needed), annotating
      each line with its feature assignment (as an optional field with tag ‘XF’)

    -n <n>, --nprocesses=<n>
      Number of parallel CPU processes to use (default: 1). This option is useful to process
      several input files at once. Each file will use only 1 CPU. It is possible, of course, to
      split a very large input SAM/BAM files into smaller chunks upstream to make use of this option.

    -p <samout_format>, --samout-format=<samout_format>
      Format to use with the –samout option, can be bam or sam (default: sam).

    -q, --quiet
      Suppress progress report and warnings

    -h, --help
      Show a usage summary and exit

    --version
      Show software version and exit

  Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology
  Laboratory (EMBL). (c) 2010. Released under the terms of the GNU General
  Public License v3. Part of the 'HTSeq' framework, version 0.11.2.