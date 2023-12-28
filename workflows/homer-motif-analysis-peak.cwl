cwlVersion: v1.0
class: Workflow


requirements:
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement


'sd:upstream':
  genome_indices:
    - "genome-indices.cwl"
  regions_a:
    - "chipseq-se.cwl"
    - "chipseq-pe.cwl"
    - "trim-chipseq-se.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "trim-atacseq-pe.cwl"
  regions_b:
    - "chipseq-se.cwl"
    - "chipseq-pe.cwl"
    - "trim-chipseq-se.cwl"
    - "trim-chipseq-pe.cwl"
    - "trim-atacseq-se.cwl"
    - "trim-atacseq-pe.cwl"

inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  regions_files_a:
    type: File[]
    format: "http://edamontology.org/format_3613"
    label: "Samples to select target regions from"
    doc: "Narrow peak files to select target regions from"
    'sd:upstreamSource': "regions_a/macs2_narrow_peaks"
    'sd:localLabel': true

  regions_files_b:
    type: File[]
    format: "http://edamontology.org/format_3613"
    label: "Samples to select background regions from"
    doc: "Narrow peak files to select background regions from"
    'sd:upstreamSource': "regions_b/macs2_narrow_peaks"
    'sd:localLabel': true

  diff_regions_file_a:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Target regions ranges. Headerless BED file with minimum [chrom start end name dummy strand] columns. Optionally, CSV"
    doc: "Target regions ranges. Headerless BED file with minimum [chrom start end unique_id dummy strand] columns. Optionally, CSV"

  diff_regions_file_b:
    type: File
    format: "http://edamontology.org/format_3003"
    label: "Background regions ranges. Headerless BED file with minimum [chrom start end name dummy strand] columns. Optionally, CSV"
    doc: "Background regions ranges. Headerless BED file with minimum [chrom start end unique_id dummy strand] columns. Optionally, CSV"

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

  min_signal_regions_a:
    type: string?
    default: "0"
    label: "Min signalValue for peaks selected as target regions"
    doc: "Discard all peaks from narrowPeak file of target regions with signalValue smaller than this threshold"
    'sd:layout':
      advanced: true

  min_signal_regions_b:
    type: string?
    default: "0"
    label: "Min signalValue for peaks selected as background regions"
    doc: "Discard all peaks from narrowPeak file of background regions with signalValue smaller than this threshold"
    'sd:layout':
      advanced: true

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

  homer_stdout_log_file:
    type: File
    format: "http://edamontology.org/format_2330"
    outputSource: find_motifs/stdout_log
    label: "Homer stdout log"
    doc: "Homer stdout log"

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

  homer_stderr_log:
    type: File
    format: "http://edamontology.org/format_2330"
    outputSource: find_motifs/stderr_log
    label: "Homer stderr log"
    doc: "Homer stderr log"


steps:

  dedup_and_sort_diff_regions_a:
    run: ../tools/custom-bash.cwl
    in:
      input_file: diff_regions_file_a
      script:
        default: |
          cat "$0" | tr -d '\r' | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n | awk '{print $1"\t"$2"\t"$3"\tp"NR"\t"$5"\t"$6}' > sorted_unique_diff_regions_from_a.bed
    out:
      - output_file

  dedup_and_sort_diff_regions_b:
    run: ../tools/custom-bash.cwl
    in:
      input_file: diff_regions_file_b
      script:
        default: |
          cat "$0" | tr -d '\r' | tr "," "\t" | awk NF | sort -u -k1,1 -k2,2n -k3,3n | awk '{print $1"\t"$2"\t"$3"\tp"NR"\t"$5"\t"$6}' > sorted_unique_diff_regions_from_b.bed
    out:
      - output_file

  concat_dedup_and_sort_regions_a:
    run: ../tools/custom-bash.cwl
    in:
      input_file: regions_files_a
      param: min_signal_regions_a
      script:
        default: |
          set -- "$0" "$@"
          echo "files: " "${@:1:$#-1}"
          TH="${@: -1}"
          echo "Threashold" "$TH"
          cat "${@:1:$#-1}" | tr -d '\r' | tr "," "\t" | awk NF | awk -v th=$TH '{if ($7 >= th) print $0}' | sort -u -k1,1 -k2,2n -k3,3n | awk '{print $1"\t"$2"\t"$3"\tp"NR"\t"$5"\t"$6}' > concatenated_sorted_unique_regions_from_a.bed
    out:
      - output_file

  concat_dedup_and_sort_regions_b:
    run: ../tools/custom-bash.cwl
    in:
      input_file: regions_files_b
      param: min_signal_regions_b
      script:
        default: |
          set -- "$0" "$@"
          echo "files: " "${@:1:$#-1}"
          TH="${@: -1}"
          echo "Threashold" "$TH"
          cat "${@:1:$#-1}" | tr -d '\r' | tr "," "\t" | awk NF | awk -v th=$TH '{if ($7 >= th) print $0}' | sort -u -k1,1 -k2,2n -k3,3n | awk '{print $1"\t"$2"\t"$3"\tp"NR"\t"$5"\t"$6}' > concatenated_sorted_unique_regions_from_b.bed
    out:
      - output_file

  get_overlapped_with_diff_regions_a:
    run: ../tools/bedtools-intersect.cwl
    in:
      file_a: concat_dedup_and_sort_regions_a/output_file
      file_b: dedup_and_sort_diff_regions_a/output_file
      report_from_a_once:
        default: true
      output_filename:
        default: "concatenated_sorted_unique_regions_from_a_overlapped_with_diff.bed"
    out: [intersected_file]

  get_overlapped_with_diff_regions_b:
    run: ../tools/bedtools-intersect.cwl
    in:
      file_a: concat_dedup_and_sort_regions_b/output_file
      file_b: dedup_and_sort_diff_regions_b/output_file
      report_from_a_once:
        default: true
      output_filename:
        default: "concatenated_sorted_unique_regions_from_b_overlapped_with_diff.bed"
    out: [intersected_file]

  merge_overlapped_with_diff_regions_a:
    run: ../tools/bedtools-merge.cwl
    in:
      bed_file: get_overlapped_with_diff_regions_a/intersected_file
      output_filename:
        default: "merged_overlapped_with_diff_concatenated_sorted_unique_regions_from_a.bed"
    out: [merged_bed_file]
  
  merge_overlapped_with_diff_regions_b:
    run: ../tools/bedtools-merge.cwl
    in:
      bed_file: get_overlapped_with_diff_regions_b/intersected_file
      output_filename:
        default: "merged_overlapped_with_diff_concatenated_sorted_unique_regions_from_b.bed"
    out: [merged_bed_file]

  find_motifs:
    run: ../tools/homer-find-motifs-genome.cwl
    in:
      target_regions_file: merge_overlapped_with_diff_regions_a/merged_bed_file
      background_regions_file: merge_overlapped_with_diff_regions_b/merged_bed_file
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

label: "DEPRECATED - Motif Finding with HOMER with target and background regions from peaks"
s:name: "DEPRECATED - Motif Finding with HOMER with target and background regions from peaks"
s:alternateName: "DEPRECATED - Motif Finding with HOMER with target and background regions from peaks"

s:downloadUrl: https://raw.githubusercontent.com/datirium/workflows/master/workflows/homer-motif-analysis-peak.cwl
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
#   $include: ../descriptions/homer-motif-analysis-peak.md


doc: |
  Motif Finding with HOMER with target and background regions from peaks
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