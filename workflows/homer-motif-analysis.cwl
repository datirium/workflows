cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  genome_indices:
    - "genome-indices.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  regions_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Regions file. Headerless BED file with minimum [chrom start end] columns"
    doc: "Regions of interest. Formatted as headerless BED file with minimum [chrom start end] columns"

  motifs_db:
    type:
      - "null"
      - type: enum
        symbols: ["vertebrates", "insects", "worms", "plants", "yeast", "all"]
    default: "vertebrates"
    label: "Set motifs DB to check against"
    doc: "Set motifs DB to check against"

  chrom_length_file:
    type: File
    format: "http://edamontology.org/format_2330"
    label: "Chromosome length file"
    doc: "Chromosome length file"
    'sd:upstreamSource': "genome_indices/chrom_length"

  genome_fasta_file:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "Reference genome FASTA file"
    doc: "Reference genome FASTA file. Includes all chromosomes in a single file"
    'sd:upstreamSource': "genome_indices/fasta"

  skip_denovo:
    type: boolean?
    default: True
    label: "Skip de novo motif enrichment"
    doc: "Skip de novo motif enrichment"
    'sd:layout':
      advanced: true

  skip_known:
    type: boolean?
    default: False
    label: "Skip known motif enrichment"
    doc: "Skip known motif enrichment"
    'sd:layout':
      advanced: true

  use_binomial:
    type: boolean?
    default: False
    label: "Use binomial distribution instead of hypergeometric to calculate p-values"
    doc: "Use binomial distribution instead of hypergeometric to calculate p-values"
    'sd:layout':
      advanced: true
  
  threads:
    type: int?
    default: 4
    label: "Threads number"
    doc: "Number of threads for those steps that support multithreading"
    'sd:layout':
      advanced: true


outputs:

  homer_found_motifs:
    type: File
    outputSource: find_motifs/compressed_results_folder
    label: "Compressed file with Homer motifs"
    doc: "Homer motifs"

  homer_stdout_log:
    type: File
    format: "http://edamontology.org/format_2330"
    outputSource: find_motifs/stdout_log
    label: "Homer stdout log"
    doc: "Homer stdout log"

  homer_known_motifs:
    type: File?
    format: "http://edamontology.org/format_2331"
    outputSource: find_motifs/known_motifs
    label: "Known motifs html file"
    doc: "Known motifs html file"

  homer_denovo_motifs:
    type: File?
    format: "http://edamontology.org/format_2331"
    outputSource: find_motifs/denovo_motifs
    label: "de novo motifs html file"
    doc: "de novo motifs html file"

  homer_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    outputSource: find_motifs/stderr_log
    label: "Homer stderr log"
    doc: "Homer stderr log"


steps:

  make_unique:
    run: ../tools/custom-bash.cwl
    in:
      input_file: regions_file
      script:
        default: |
          cut -f 1-3 "$0" | sort -k1,1 -k2,2n | sort -u > `basename $0`
    out:
      - output_file

  bedtools_slop:
    run: ../tools/bedtools-slop.cwl
    in:
      bed_file: make_unique/output_file
      chrom_length_file: chrom_length_file
      bi_direction:
        default: 20000
    out:
      - extended_bed_file

  bedtools_sort:
    run: ../tools/linux-sort.cwl
    in:
      unsorted_file: bedtools_slop/extended_bed_file
      key:
        default: ["1,1","2,2n"]
    out:
      - sorted_file

  bedtools_merge:
    run: ../tools/bedtools-merge.cwl
    in:
      bed_file: bedtools_sort/sorted_file
    out:
      - merged_bed_file

  bedtools_subtract:
    run: ../tools/bedtools-subtract.cwl
    in:
      reduced_bed_file: bedtools_merge/merged_bed_file
      subtracted_bed_file: make_unique/output_file
    out:
      - difference_bed_file

  bedtools_shuffle:
    run: ../tools/bedtools-shuffle.cwl
    in:
      bed_file: make_unique/output_file
      chrom_length_file: chrom_length_file
      incl_bed_file: bedtools_subtract/difference_bed_file
      no_overlapping:
        default: True
      max_tries:
        default: 10000
      seed:
        default: 123456789
    out:
      - shuffled_bed_file

  bedtools_get_fasta_target:
    run: ../tools/bedtools-getfasta.cwl
    in:
      intervals_file: make_unique/output_file
      genome_fasta_file: genome_fasta_file
    out:
      - sequences_file

  bedtools_get_fasta_background:
    run: ../tools/bedtools-getfasta.cwl
    in:
      intervals_file: bedtools_shuffle/shuffled_bed_file
      genome_fasta_file: genome_fasta_file
    out:
      - sequences_file

  find_motifs:
    run: ../tools/homer-find-motifs.cwl
    in:
      target_fasta_file: bedtools_get_fasta_target/sequences_file
      background_fasta_file: bedtools_get_fasta_background/sequences_file
      skip_denovo: skip_denovo
      skip_known: skip_known
      use_binomial: use_binomial
      motifs_db: motifs_db
      threads: threads
    out:
      - compressed_results_folder
      - known_motifs
      - denovo_motifs
      - stdout_log
      - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

label: "Motif Finding with HOMER from FASTA files"
s:name: "Motif Finding with HOMER from FASTA files"
s:alternateName: "Motif Finding with HOMER from FASTA files"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/homer-motif-analysis.cwl
s:codeRepository: https://github.com/datirium/workflows
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
          s:email: mailto:michael.kotliar@cchmc.org
          s:sameAs:
          - id: http://orcid.org/0000-0002-6486-3898


# doc:
#   $include: ../descriptions/homer-motif-analysis.md


doc: |
  Motif Finding with HOMER from FASTA files
  ---------------------------------------------------
  
  HOMER contains a novel motif discovery algorithm that was designed for regulatory element analysis
  in genomics applications (DNA only, no protein). It is a differential motif discovery algorithm,
  which means that it takes two sets of sequences and tries to identify the regulatory elements that
  are specifically enriched in on set relative to the other. It uses ZOOPS scoring (zero or one
  occurrence per sequence) coupled with the hypergeometric enrichment calculations (or binomial) to
  determine motif enrichment. HOMER also tries its best to account for sequenced bias in the dataset.
  It was designed with ChIP-Seq and promoter analysis in mind, but can be applied to pretty much any
  nucleic acids motif finding problem.

  For more information please refer to:
  -------------------------------------
  [Official documentation](http://homer.ucsd.edu/homer/motif/)