cwlVersion: v1.1
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var split_by_common_delim = function(line) {
        function get_unique(value, index, self) {
          return self.indexOf(value) === index && value != "";
        }
        let splitted_line = line?line.split(/[\s,]+/).filter(get_unique):null;
        return (splitted_line && !!splitted_line.length)?splitted_line:null;
    };


"sd:upstream":
  sc_rnaseq_sample:
  - "cellranger-multi.cwl"
  - "single-cell-preprocess-cellranger.cwl"


inputs:

  alias:
    type: string
    label: "Analysis name"
    sd:preview:
      position: 1

  possorted_genome_bam_bai:
    type: File
    secondaryFiles:
    - .bai
    label: "Cell Ranger Count RNA or RNA+VDJ Sample"
    doc: |
      Any Cell Ranger Count (RNA) or
      Cell Ranger Count (RNA+VDJ) Sample
    "sd:upstreamSource": "sc_rnaseq_sample/possorted_genome_bam_bai"
    "sd:localLabel": true

  genome_fasta_file:
    type: File
    secondaryFiles:
    - .fai
    label: "Cell Ranger Count RNA or RNA+VDJ Sample"
    "sd:upstreamSource": "sc_rnaseq_sample/genome_indices/genome_indices/fasta_output"

  filtered_feature_bc_matrix_folder:
    type: File
    label: "Cell Ranger Count RNA or RNA+VDJ Sample"
    "sd:upstreamSource": "sc_rnaseq_sample/filtered_feature_bc_matrix_folder"

  clusters_count:
    type: int
    label: "Number of clusters to detect (number of donors merged into one single-cell experiment)"
    doc: |
      Number of clusters to detect (number of donors merged into one single-cell experiment)

  barcodes_data:
    type: File?
    label: "Selected cell barcodes (optional)"
    doc: |
      A TSV/CSV file to optionally prefilter
      the single cell data by including only
      the cells with the selected barcodes.
      The provided file should have one cell
      barcode per line and do not include any
      header information.

  ploidy_count:
    type: int?
    default: 2
    label: "Ploidy, must be 1 or 2"
    doc: |
      Ploidy, must be 1 or 2
    "sd:layout":
      advanced: true

  min_alt:
    type: int?
    default: 10
    label: "Min alt to use locus"
    doc: |
      Min alt to use locus
    "sd:layout":
      advanced: true

  min_ref:
    type: int?
    default: 10
    label: "Min ref to use locus"
    doc: |
      Min ref to use locus
    "sd:layout":
      advanced: true

  max_loci:
    type: int?
    default: 2048
    label: "Max loci per cell, affects speed"
    doc: |
      Max loci per cell, affects speed
    "sd:layout":
      advanced: true

  restarts_count:
    type: int?
    default: 100
    label: "Number of restarts in clustering, when there are > 12 clusters we recommend increasing this to avoid local minima"
    doc: |
      Number of restarts in clustering, when there are > 12
      clusters we recommend increasing this to avoid local
      minima
    "sd:layout":
      advanced: true

  known_genotypes_sample_names:
    type: string?
    label: "Which samples in population VCF from known genotypes option represent the donors in your sample"
    doc: |
      Which samples in population VCF from known genotypes
      option represent the donors in your sample
    "sd:layout":
      advanced: true

  skip_remap:
    type: boolean?
    default: false
    label: "Don't remap with minimap2 (not recommended unless Common variant loci VCF file was provided)"
    doc: |
      Don't remap with minimap2 (not recommended unless in
      conjunction with --common_variants)
    "sd:layout":
      advanced: true

  ignore_data_errors:
    type: boolean?
    label: "Ignore data error assertions"
    doc: |
      Set to True to ignore data error assertions
    "sd:layout":
      advanced: true

  threads:
    type: int?
    default: 2
    label: "Threads number to use"
    doc: |
      Threads number
    "sd:layout":
      advanced: true

  common_variants_vcf_file:
    type: File?
    label: "Common variant loci or known variant loci VCF file"
    doc: |
      Common variant loci or known variant loci VCF file,
      must be made vs the same reference fasta
    "sd:layout":
      advanced: true

  known_genotypes_vcf_file:
    type: File?
    label: "Known variants per clone in population VCF file"
    doc: |
      Known variants per clone in population VCF mode, must be .vcf
    "sd:layout":
      advanced: true


outputs:

  genotype_cluster_tsv_file:
    type: File
    outputSource: rna_souporcell/genotype_cluster_tsv_file
    label: "Cellurar barcodes file clustered by genotype"
    doc: |
      Cellurar barcodes file clustered by genotype
    "sd:visualPlugins":
    - syncfusiongrid:
        tab: "Genotypes"
        Title: "Cells clustered by genotype"

  genotype_cluster_vcf_file:
    type: File
    outputSource: rna_souporcell/genotype_cluster_vcf_file
    label: "VCF file with genotypes for each cluster for each variant call"
    doc: |
      VCF file with genotypes for each cluster for each variant call.
      Refer to http://software.broadinstitute.org/software/igv/viewing_vcf_files
      for track description when displaying in IGV.

  ambient_rna_file:
    type: File
    outputSource: rna_souporcell/ambient_rna_file
    label: "Ambient RNA evaluation text file"
    doc: |
      Ambient RNA evaluation text file

  rna_souporcell_stdout_log:
    type: File
    outputSource: rna_souporcell/stdout_log
    label: stdout log generated by souporcell
    doc: |
      stdout log generated by souporcell

  rna_souporcell_stderr_log:
    type: File
    outputSource: rna_souporcell/stderr_log
    label: stderr log generated by souporcell
    doc: |
      stderr log generated by souporcell


steps:

  get_barcodes_tsv_file:
    run:
      cwlVersion: v1.1
      class: CommandLineTool
      hints:
      - class: DockerRequirement
        dockerPull: biowardrobe2/scidap:v0.0.3
      inputs:
        script:
          type: string?
          default: |
            #!/bin/bash
            tar xzf $0
            mv */barcodes.tsv.gz .
            gunzip barcodes.tsv.gz
            if [ -f "$1" ]; then
                echo "Filter by user provided barcodes"
                comm -12 --check-order <(sort barcodes.tsv) <(sort $1) > cell_barcodes.tsv
            else
                echo "Do not filter by user provided barcodes"
                mv barcodes.tsv cell_barcodes.tsv
            fi
          inputBinding:
            position: 5
        filtered_feature_bc_matrix_folder:
          type: File
          inputBinding:
            position: 6
        barcodes_data:
          type: File?
          inputBinding:
            position: 7
      outputs:
        barcodes_tsv_file:
          type: File
          outputBinding:
            glob: "cell_barcodes.tsv"
      baseCommand: ["bash", "-c"]
    in:
      filtered_feature_bc_matrix_folder: filtered_feature_bc_matrix_folder
      barcodes_data: barcodes_data
    out:
    - barcodes_tsv_file

  rna_souporcell:
    run: ../tools/souporcell.cwl
    in:
      possorted_genome_bam_bai: possorted_genome_bam_bai
      barcodes_tsv_file: get_barcodes_tsv_file/barcodes_tsv_file
      genome_fasta_file: genome_fasta_file
      clusters_count: clusters_count
      ploidy_count: ploidy_count
      min_alt: min_alt
      min_ref: min_ref
      max_loci: max_loci
      restarts_count: restarts_count
      common_variants_vcf_file: common_variants_vcf_file
      known_genotypes_vcf_file: known_genotypes_vcf_file
      known_genotypes_sample_names:
        source: known_genotypes_sample_names
        valueFrom: $(split_by_common_delim(self))
      skip_remap: skip_remap
      ignore_data_errors: ignore_data_errors
      threads: threads
    out:
    - genotype_cluster_tsv_file
    - genotype_cluster_vcf_file
    - ambient_rna_file
    - stdout_log
    - stderr_log


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Souporcell Cluster by Genotype for RNA"
s:name: "Souporcell Cluster by Genotype for RNA"
s:alternateName: "Souporcell Cluster by Genotype for RNA"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/souporcell-rna.cwl
s:codeRepository: https://github.com/Barski-lab/workflows-datirium
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
  Souporcell Cluster by Genotype for RNA

  Souporcell: robust clustering of single-cell data by
  genotype without reference genotypes