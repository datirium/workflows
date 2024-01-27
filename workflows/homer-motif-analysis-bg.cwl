cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  genome_indices:
    - "genome-indices.cwl"
  target_peaklist:
    - "filter-diffbind-for-heatmap.cwl"
    - "filter-peaks-for-heatmap.cwl"
    - "filter-peaks-by-overlap.cwl"
    - "genelists-sets.cwl"
  background_peaklist:
    - "filter-diffbind-for-heatmap.cwl"
    - "filter-peaks-for-heatmap.cwl"
    - "filter-peaks-by-overlap.cwl"
    - "genelists-sets.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  target_regions_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Target regions. Headerless BED file with minimum [chrom start end name dummy strand] columns. Optionally, CSV"
    doc: "Target regions. Headerless BED file with minimum [chrom start end unique_id dummy strand] columns. Optionally, CSV"
    'sd:upstreamSource': "target_peaklist/filtered_file"
    'sd:localLabel': true

  background_regions_file:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Background regions. Headerless BED file with minimum [chrom start end name dummy strand] columns. Optionally, CSV"
    doc: "Background regions. Headerless BED file with minimum [chrom start end unique_id dummy strand] columns. Optionally, CSV"
    'sd:upstreamSource': "background_peaklist/filtered_file"
    'sd:localLabel': true

  genome_fasta_file:
    type: File
    format: "http://edamontology.org/format_1929"
    label: "Reference genome FASTA file"
    doc: "Reference genome FASTA file. Includes all chromosomes in a single file"
    'sd:upstreamSource': "genome_indices/fasta_output"

  motifs_db:
    type:
      - "null"
      - type: enum
        symbols: ["vertebrates", "insects", "worms", "plants", "yeast", "all"]
    default: "vertebrates"
    label: "Set motifs DB to check against"
    doc: "Set motifs DB to check against"

  chopify_background_regions:
    type: boolean?
    default: false
    label: "Chop up large background regions to the avg size of target regions"
    doc: "Chop up large background regions to the avg size of target regions"
    'sd:layout':
      advanced: true

  apply_mask_on_genome:
    type: boolean?
    default: true
    label: "Mask all repeats with N"
    doc: "Use the repeat-masked sequence (all repeats will be masked by N)"
    'sd:layout':
      advanced: true

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

  use_hypergeometric:
    type: boolean?
    default: False
    label: "Use hypergeometric for p-values, instead of default binomial. Usefull if the number of background sequences is smaller than target sequences"
    doc: "Use hypergeometric for p-values, instead of default binomial. Usefull if the number of background sequences is smaller than target sequences"
    'sd:layout':
      advanced: true
  
  search_size:
    type: string?
    default: "200"
    label: "Fragment size to use for motif finding. Can be set as <#> or <#,#>"
    doc: |
      Fragment size to use for motif finding.
      <#> - i.e. -size 300 will get sequences from -150 to +150 relative from center
      <#,#> - i.e. -size -100,50 will get sequences from -100 to +50 relative from center
      given - will use the exact regions you give it.
      Default=200
    'sd:layout':
      advanced: true

  motif_length:
    type: string?
    default: "8,10,12"
    label: "Motif length(s) for de novo motif discovery"
    doc: "<#>[,<#>,<#>...] - motif length. Default=8,10,12"
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

  homer_known_motifs:
    type: File?
    format: "http://edamontology.org/format_2331"
    outputSource: find_motifs/known_motifs
    label: "Known motifs"
    doc: "Known motifs html file"
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  homer_denovo_motifs:
    type: File?
    format: "http://edamontology.org/format_2331"
    outputSource: find_motifs/denovo_motifs
    label: "de novo motifs"
    doc: "de novo motifs html file"
    'sd:visualPlugins':
    - linkList:
        tab: 'Overview'
        target: "_blank"

  stdout_log_file:
    type: File
    format: "http://edamontology.org/format_2330"
    outputSource: find_motifs/stdout_log
    label: "Homer stdout log"
    doc: "Homer stdout log"

  stderr_log_file:
    type: File
    format: "http://edamontology.org/format_2330"
    outputSource: find_motifs/stderr_log
    label: "Homer stderr log"
    doc: "Homer stderr log"


steps:

  make_target_regions_unique:
    run: ../tools/custom-bash.cwl
    in:
      input_file: target_regions_file
      script:
        default: |
          cat "$0" | tr -d '\r' | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n | awk '{print $1"\t"$2"\t"$3"\tp"NR"\t"$5"\t"$6}' > `basename $0`
    out:
      - output_file

  make_background_regions_unique:
    run: ../tools/custom-bash.cwl
    in:
      input_file: background_regions_file
      script:
        default: |
          cat "$0" | tr -d '\r' | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n | awk '{print $1"\t"$2"\t"$3"\tp"NR"\t"$5"\t"$6}' > `basename $0`
    out:
      - output_file

  find_motifs:
    run: ../tools/homer-find-motifs-genome.cwl
    in:
      target_regions_file: make_target_regions_unique/output_file
      background_regions_file: make_background_regions_unique/output_file
      genome_fasta_file: genome_fasta_file
      chopify_background_regions: chopify_background_regions
      search_size: search_size
      motif_length: motif_length
      apply_mask_on_genome: apply_mask_on_genome
      use_hypergeometric: use_hypergeometric
      skip_denovo: skip_denovo
      skip_known: skip_known
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
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

label: "Motif Finding with HOMER with custom background regions"
s:name: "Motif Finding with HOMER with custom background regions"
s:alternateName: "Motif Finding with HOMER with custom background regions"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/homer-motif-analysis-bg.cwl
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
#   $include: ../descriptions/homer-motif-analysis-bg.md


doc: |
  Motif Finding with HOMER with custom background regions
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