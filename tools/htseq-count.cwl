cwlVersion: v1.0
class: CommandLineTool


hints:
- class: DockerRequirement
  dockerPull: biocontainers/htseq:v0.11.2-1-deb-py3_cv1


inputs:

  alignment_bam_file:
    type: File
    secondaryFiles:
    - .bai
    inputBinding:
      position: 2
    doc: "Path to the coordinate sorted indexed BAM file"

  annotation_gtf_file:
    type: File
    label: ""
    inputBinding:
      position: 3
    doc: "GTF annotation file"


outputs:
  
  gene_expression_report:
    type: stdout


baseCommand: ["htseq-count", "-f", "bam", "-r", "pos", "-s", "no", "-i", "gene_id"]


stdout: gene_expression_report.tsv


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
  Use minimum number of parameters. Harcoded to return gene expression (gene_id)
  iformation from coordinate sorted and indexed BAM file. Not strand specific.


s:about: |
  usage: htseq-count [options] alignment_file gff_file

  This script takes one or more alignment files in SAM/BAM format and a feature
  file in GFF format and calculates for each feature the number of reads mapping
  to it. See http://htseq.readthedocs.io/en/master/count.html for details.

  positional arguments:
    samfilenames          Path to the SAM/BAM files containing the mapped reads.
                          If '-' is selected, read from standard input
    featuresfilename      Path to the file containing the features

  optional arguments:
    -h, --help            show this help message and exit
    -f {sam,bam}, --format {sam,bam}
                          type of <alignment_file> data, either 'sam' or 'bam'
                          (default: sam)
    -r {pos,name}, --order {pos,name}
                          'pos' or 'name'. Sorting order of <alignment_file>
                          (default: name). Paired-end sequencing data must be
                          sorted either by position or by read name, and the
                          sorting order must be specified. Ignored for single-
                          end data.
    --max-reads-in-buffer MAX_BUFFER_SIZE
                          When <alignment_file> is paired end sorted by
                          position, allow only so many reads to stay in memory
                          until the mates are found (raising this number will
                          use more memory). Has no effect for single end or
                          paired end sorted by name
    -s {yes,no,reverse}, --stranded {yes,no,reverse}
                          whether the data is from a strand-specific assay.
                          Specify 'yes', 'no', or 'reverse' (default: yes).
                          'reverse' means 'yes' with reversed strand
                          interpretation
    -a MINAQUAL, --minaqual MINAQUAL
                          skip all reads with alignment quality lower than the
                          given minimum value (default: 10)
    -t FEATURETYPE, --type FEATURETYPE
                          feature type (3rd column in GFF file) to be used, all
                          features of other type are ignored (default, suitable
                          for Ensembl GTF files: exon)
    -i IDATTR, --idattr IDATTR
                          GFF attribute to be used as feature ID (default,
                          suitable for Ensembl GTF files: gene_id)
    --additional-attr ADDITIONAL_ATTR
                          Additional feature attributes (default: none, suitable
                          for Ensembl GTF files: gene_name). Use multiple times
                          for each different attribute
    -m {union,intersection-strict,intersection-nonempty}, --mode {union,intersection-strict,intersection-nonempty}
                          mode to handle reads overlapping more than one feature
                          (choices: union, intersection-strict, intersection-
                          nonempty; default: union)
    --nonunique {none,all}
                          Whether to score reads that are not uniquely aligned
                          or ambiguously assigned to features
    --secondary-alignments {score,ignore}
                          Whether to score secondary alignments (0x100 flag)
    --supplementary-alignments {score,ignore}
                          Whether to score supplementary alignments (0x800 flag)
    -o SAMOUTS, --samout SAMOUTS
                          write out all SAM alignment records into SAM files
                          (one per input file needed), annotating each line with
                          its feature assignment (as an optional field with tag
                          'XF')
    -q, --quiet           suppress progress report

  Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology
  Laboratory (EMBL). (c) 2010. Released under the terms of the GNU General
  Public License v3. Part of the 'HTSeq' framework, version 0.11.2.