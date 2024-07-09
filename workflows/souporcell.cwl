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


'sd:upstream':
  sc_rnaseq_sample:
  - "https://github.com/datirium/workflows/workflows/cellranger-arc-count.cwl"
  - "cellranger-arc-count.cwl"


inputs:

  alias:
    type: string
    label: "Experiment short name/Alias"
    sd:preview:
      position: 1

  rna_possorted_genome_bam_bai:
    type: File
    secondaryFiles:
    - .bai
    label: "scRNA-Seq Cell Ranger ARC Experiment"
    doc: |
      GEX position-sorted reads aligned to the genome and transcriptome
      annotated with barcode information in BAM format
    'sd:upstreamSource': "sc_rnaseq_sample/rna_possorted_genome_bam_bai"
    'sd:localLabel': true

  atac_possorted_genome_bam_bai:
    type: File
    secondaryFiles:
    - .bai
    label: "scRNA-Seq Cell Ranger ARC Experiment"
    doc: |
      ATAC position-sorted reads aligned to the genome annotated with
      barcode information in BAM format
    'sd:upstreamSource': "sc_rnaseq_sample/atac_possorted_genome_bam_bai"
    'sd:localLabel': true

  genome_fasta_file:
    type: File
    secondaryFiles:
    - .fai
    label: "scRNA-Seq Cell Ranger ARC Experiment"
    doc: |
      Reference genome FASTA file that includes all chromosomes + fai index file
    'sd:upstreamSource': "sc_rnaseq_sample/genome_indices/genome_indices/fasta_output"
    'sd:localLabel': true

  filtered_feature_bc_matrix_folder:
    type: File
    label: "scRNA-Seq Cell Ranger ARC Experiment"
    doc: |
      Compressed folder with filtered feature barcode matrix stored as a CSC sparse
      matrix in MEX format. The rows consist of all the gene and peak features
      concatenated together (identical to raw feature barcode matrix) and the columns
      are restricted to those barcodes that are identified as cells.
    'sd:upstreamSource': "sc_rnaseq_sample/filtered_feature_bc_matrix_folder"
    'sd:localLabel': true

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

  regions_bed_file:
    type: File?
    label: "Selected regions (optional)"
    doc: |
      A BED file to optionally prefilter
      reads from the both RNA and ATAC BAM
      files. If minimap2 remapping is not
      skipped, filtering by regions will be
      done after it.

  ploidy_count:
    type: int?
    default: 2
    label: "Ploidy, must be 1 or 2"
    doc: |
      Ploidy, must be 1 or 2
    'sd:layout':
      advanced: true

  min_alt:
    type: int?
    default: 10
    label: "Min alt to use locus"
    doc: |
      Min alt to use locus
    'sd:layout':
      advanced: true

  min_ref:
    type: int?
    default: 10
    label: "Min ref to use locus"
    doc: |
      Min ref to use locus
    'sd:layout':
      advanced: true

  max_loci:
    type: int?
    default: 2048
    label: "Max loci per cell, affects speed"
    doc: |
      Max loci per cell, affects speed
    'sd:layout':
      advanced: true

  restarts_count:
    type: int?
    default: 100
    label: "Number of restarts in clustering, when there are > 12 clusters we recommend increasing this to avoid local minima"
    doc: |
      Number of restarts in clustering, when there are > 12
      clusters we recommend increasing this to avoid local
      minima
    'sd:layout':
      advanced: true

  known_genotypes_sample_names:
    type: string?
    label: "Which samples in population VCF from known genotypes option represent the donors in your sample"
    doc: |
      Which samples in population VCF from known genotypes
      option represent the donors in your sample
    'sd:layout':
      advanced: true

  skip_remap:
    type: boolean?
    default: false
    label: "Don't remap with minimap2 (not recommended unless Common variant loci VCF file was provided)"
    doc: |
      Don't remap with minimap2 (not recommended unless in
      conjunction with --common_variants)
    'sd:layout':
      advanced: true

  ignore_data_errors:
    type: boolean?
    label: "Ignore data error assertions"
    doc: |
      Set to True to ignore data error assertions
    'sd:layout':
      advanced: true

  threads:
    type: int?
    default: 2
    label: "Threads number to use"
    doc: |
      Threads number
    'sd:layout':
      advanced: true

  common_variants_vcf_file:
    type: File?
    label: "Common variant loci or known variant loci VCF file"
    doc: |
      Common variant loci or known variant loci VCF file,
      must be made vs the same reference fasta
    'sd:layout':
      advanced: true

  known_genotypes_vcf_file:
    type: File?
    label: "Known variants per clone in population VCF file"
    doc: |
      Known variants per clone in population VCF mode, must be .vcf
    'sd:layout':
      advanced: true


outputs:

  gex_genotype_cluster_tsv_file:
    type: File
    outputSource: rename_outputs/gex_genotype_cluster_tsv_file
    label: "Cellurar barcodes file clustered by genotype (GEX)"
    doc: |
      Cellurar barcodes file clustered by genotype (GEX)
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'GEX'
        Title: 'Cells clustered by genotype (GEX)'

  gex_genotype_cluster_vcf_file:
    type: File
    outputSource: rename_outputs/gex_genotype_cluster_vcf_file
    label: "VCF file with genotypes for each cluster for each variant call (GEX)"
    doc: |
      VCF file with genotypes for each cluster for each variant call (GEX).
      Refer to http://software.broadinstitute.org/software/igv/viewing_vcf_files
      for track description when displaying in IGV.

  gex_ambient_rna_file:
    type: File
    outputSource: rename_outputs/gex_ambient_rna_file
    label: "Ambient RNA evaluation text file (GEX)"
    doc: |
      Ambient RNA evaluation text file (GEX)

  gex_souporcell_stdout_log:
    type: File
    outputSource: rename_outputs/gex_souporcell_stdout_log
    label: stdout log generated by souporcell (GEX)
    doc: |
      stdout log generated by souporcell (GEX)

  gex_souporcell_stderr_log:
    type: File
    outputSource: rename_outputs/gex_souporcell_stderr_log
    label: stderr log generated by souporcell (GEX)
    doc: |
      stderr log generated by souporcell (GEX)


  atac_genotype_cluster_tsv_file:
    type: File
    outputSource: rename_outputs/atac_genotype_cluster_tsv_file
    label: "Cellurar barcodes file clustered by genotype (ATAC)"
    doc: |
      Cellurar barcodes file clustered by genotype (ATAC)
    'sd:visualPlugins':
    - syncfusiongrid:
        tab: 'ATAC'
        Title: 'Cells clustered by genotype (ATAC)'

  atac_genotype_cluster_vcf_file:
    type: File
    outputSource: rename_outputs/atac_genotype_cluster_vcf_file
    label: "VCF file with genotypes for each cluster for each variant call (ATAC)"
    doc: |
      VCF file with genotypes for each cluster for each variant call (ATAC).
      Refer to http://software.broadinstitute.org/software/igv/viewing_vcf_files
      for track description when displaying in IGV.

  atac_ambient_rna_file:
    type: File
    outputSource: rename_outputs/atac_ambient_rna_file
    label: "Ambient RNA evaluation text file (ATAC)"
    doc: |
      Ambient RNA evaluation text file (ATAC)

  atac_souporcell_stdout_log:
    type: File
    outputSource: rename_outputs/atac_souporcell_stdout_log
    label: stdout log generated by souporcell (ATAC)
    doc: |
      stdout log generated by souporcell (ATAC)

  atac_souporcell_stderr_log:
    type: File
    outputSource: rename_outputs/atac_souporcell_stderr_log
    label: stderr log generated by souporcell (ATAC)
    doc: |
      stderr log generated by souporcell (ATAC)

  upset_plot_png:
    type: File
    outputSource: overlap_clusters/upset_plot_png
    label: RNA and ATAC clusters overlap plot
    doc: |
      RNA and ATAC cluster overlap plot
      PNG format
    'sd:visualPlugins':
    - image:
        tab: 'Plots'
        Caption: 'RNA and ATAC cluster overlap'


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
            mv filtered_feature_bc_matrix/barcodes.tsv.gz .
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

  gex_souporcell:
    run: ../tools/souporcell.cwl
    in:
      possorted_genome_bam_bai: rna_possorted_genome_bam_bai
      barcodes_tsv_file: get_barcodes_tsv_file/barcodes_tsv_file
      genome_fasta_file: genome_fasta_file
      regions_bed_file: regions_bed_file
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

  atac_souporcell:
    run: ../tools/souporcell.cwl
    in:
      possorted_genome_bam_bai: atac_possorted_genome_bam_bai
      barcodes_tsv_file: get_barcodes_tsv_file/barcodes_tsv_file
      genome_fasta_file: genome_fasta_file
      regions_bed_file: regions_bed_file
      clusters_count: clusters_count
      ploidy_count: ploidy_count
      min_alt: min_alt
      min_ref: min_ref
      max_loci: max_loci
      restarts_count: restarts_count
      # common_variants_vcf_file: common_variants_vcf_file
      known_genotypes_vcf_file: known_genotypes_vcf_file
      known_genotypes_sample_names:
        source: known_genotypes_sample_names
        valueFrom: $(split_by_common_delim(self))
      skip_remap: skip_remap
      no_umi:
        default: true
      ignore_data_errors: ignore_data_errors
      threads: threads
    out:
    - genotype_cluster_tsv_file
    - genotype_cluster_vcf_file
    - ambient_rna_file
    - stdout_log
    - stderr_log

  rename_outputs:
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
            set -- "$0" "$@"
            for FILE in "$@"; do
              BASENAME=$(basename "$FILE")
              GEX=gex_${BASENAME}
              ATAC=atac_${BASENAME}
              if [ -f "$GEX" ]; then
                  echo "Copying ${FILE} to ${ATAC}"
                  cp $FILE $ATAC
              else
                  echo "Copying ${FILE} to ${GEX}"
                  cp $FILE $GEX
              fi
            done;
          inputBinding:
            position: 5
        gex_files:
          type: File[]
          inputBinding:
            position: 6
        atac_files:
          type: File[]
          inputBinding:
            position: 7
      outputs:
        gex_genotype_cluster_tsv_file:
          type: File
          outputBinding:
            glob: "gex_clusters.tsv"
        gex_genotype_cluster_vcf_file:
          type: File
          outputBinding:
            glob: "gex_cluster_genotypes.vcf"
        gex_ambient_rna_file:
          type: File
          outputBinding:
            glob: "gex_ambient_rna.txt"
        gex_souporcell_stdout_log:
          type: File
          outputBinding:
            glob: "gex_souporcell_stdout.log"
        gex_souporcell_stderr_log:
          type: File
          outputBinding:
            glob: "gex_souporcell_stderr.log"
        atac_genotype_cluster_tsv_file:
          type: File
          outputBinding:
            glob: "atac_clusters.tsv"
        atac_genotype_cluster_vcf_file:
          type: File
          outputBinding:
            glob: "atac_cluster_genotypes.vcf"
        atac_ambient_rna_file:
          type: File
          outputBinding:
            glob: "atac_ambient_rna.txt"
        atac_souporcell_stdout_log:
          type: File
          outputBinding:
            glob: "atac_souporcell_stdout.log"
        atac_souporcell_stderr_log:
          type: File
          outputBinding:
            glob: "atac_souporcell_stderr.log"
      baseCommand: ["bash", "-c"]
    in:
      gex_files:
      - gex_souporcell/genotype_cluster_tsv_file
      - gex_souporcell/genotype_cluster_vcf_file
      - gex_souporcell/ambient_rna_file
      - gex_souporcell/stdout_log
      - gex_souporcell/stderr_log
      atac_files:
      - atac_souporcell/genotype_cluster_tsv_file
      - atac_souporcell/genotype_cluster_vcf_file
      - atac_souporcell/ambient_rna_file
      - atac_souporcell/stdout_log
      - atac_souporcell/stderr_log
    out:
    - gex_genotype_cluster_tsv_file
    - gex_genotype_cluster_vcf_file
    - gex_ambient_rna_file
    - gex_souporcell_stdout_log
    - gex_souporcell_stderr_log
    - atac_genotype_cluster_tsv_file
    - atac_genotype_cluster_vcf_file
    - atac_ambient_rna_file
    - atac_souporcell_stdout_log
    - atac_souporcell_stderr_log

  overlap_clusters:
    run: ../tools/souporcell-cluster-overlap.cwl
    in:
      cluster_files:
      - rename_outputs/gex_genotype_cluster_tsv_file
      - rename_outputs/atac_genotype_cluster_tsv_file
      cluster_names:
        default:
        - "rna"
        - "atac"
    out:
    - upset_plot_png


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf


label: "Souporcell Cluster by Genotype"
s:name: "Souporcell Cluster by Genotype"
s:alternateName: "Souporcell: robust clustering of single-cell data by genotype without reference genotypes"

s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows-datirium/master/workflows/souporcell.cwl
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
  Souporcell Cluster by Genotype
  ==============================

  Souporcell: robust clustering of single-cell data by
  genotype without reference genotypes