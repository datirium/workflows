cwlVersion: v1.0
class: Workflow


requirements:
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement
- class: MultipleInputFeatureRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var default_output_name = function(named_input, suffix, ext) {
        if (Array.isArray(named_input) && named_input.length > 0){
          return named_input[0].location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+suffix+ext;
        } else {
          return named_input.location.split('/').slice(-1)[0].split('.').slice(0,-1).join('.')+suffix+ext;
        }
    };


inputs:

  fastq_files:
    type:
      - File
      - type: array
        items: File
    label: "Input FASTQ file(s)"
    doc: "Input FASTQ file or array of files"

  insilico_star_indices_folder:
    type: Directory
    label: "STAR indices folder for insilico genome"
    doc: "Path to STAR generated indices folder for insilico genome"

  reference_star_indices_folder:
    type: Directory
    label: "STAR indices folder for reference genome"
    doc: "Path to STAR generated indices folder for reference genome"

  reference_chrom_length_file:
    type: File
    label: "Chromosome length file for reference genome"
    doc: "Chromosome length file for reference genome"

  strain1:
    type: string
    label: "I strain name"
    doc: "First strain name"

  strain2:
    type: string
    label: "II strain name"
    doc: "Second strain name"

  strain1_chain_file:
    type: File
    label: "I strain chain file"
    doc: "Chain file to project strain I to reference genome"

  strain2_chain_file:
    type: File
    label: "II strain chain file"
    doc: "Chain file to project strain II to reference genome"

  threads:
    type: int?
    default: 2
    label: "Number of threads"
    doc: "Number of threads for those steps that support multithreading"


outputs:

  strain1_bambai_pair:
    type: File
    outputSource: strain1_project/projected_file
    label: "I strain output BAM"
    doc: "Coordinate sorted BAM file mapped to the first strain genome, projected to reference genome"

  strain2_bambai_pair:
    type: File
    outputSource: strain2_project/projected_file
    label: "II strain output BAM"
    doc: "Coordinate sorted BAM file mapped to the second strain genome, projected to reference genome"

  reference_bambai_pair:
    type: File
    outputSource: reference_samtools_sort/sorted_file
    label: "Reference output BAM"
    doc: "Coordinate sorted BAM file mapped to reference genome"

  strain1_bigwig:
    type: File
    label: "I strain bigWig file"
    doc: "Generated bigWig file for the first strain, projected to reference genome"
    outputSource: strain1_bam_to_bigwig/bigwig_file

  strain2_bigwig:
    type: File
    label: "II strain bigWig file"
    doc: "Generated bigWig file for the second strain, projected to reference genome"
    outputSource: strain2_bam_to_bigwig/bigwig_file

  reference_bigwig:
    type: File
    label: "Reference bigWig file"
    doc: "Generated BigWig file for the reference genome"
    outputSource: reference_bam_to_bigwig/bigwig_file

  insilico_star_final_log:
    type: File
    outputSource: insilico_star_aligner/log_final
    label: "STAR final log for insilico genome"
    doc: "STAR Log.final.out for insilico genome"

  insilico_star_out_log:
    type: File?
    outputSource: insilico_star_aligner/log_out
    label: "STAR log out for insilico genome"
    doc: "STAR Log.out for insilico genome"

  insilico_star_progress_log:
    type: File?
    outputSource: insilico_star_aligner/log_progress
    label: "STAR progress log for insilico genome"
    doc: "STAR Log.progress.out for insilico genome"

  insilico_star_stdout_log:
    type: File?
    outputSource: insilico_star_aligner/log_std
    label: "STAR stdout log for insilico genome"
    doc: "STAR Log.std.out for insilico genome"

  reference_star_final_log:
    type: File
    outputSource: reference_star_aligner/log_final
    label: "STAR final log for reference genome"
    doc: "STAR Log.final.out for reference genome"

  reference_star_out_log:
    type: File?
    outputSource: reference_star_aligner/log_out
    label: "STAR log out for reference genome"
    doc: "STAR Log.out for reference genome"

  reference_star_progress_log:
    type: File?
    outputSource: reference_star_aligner/log_progress
    label: "STAR progress log for reference genome"
    doc: "STAR Log.progress.out for reference genome"

  reference_star_stdout_log:
    type: File?
    outputSource: reference_star_aligner/log_std
    label: "STAR stdout log for reference genome"
    doc: "STAR Log.std.out for reference genome"

steps:

  insilico_star_aligner:
    run: ../tools/star-alignreads.cwl
    in:
      readFilesIn: fastq_files
      genomeDir: insilico_star_indices_folder
      outFileNamePrefix:
        source: [strain1,strain2]
        valueFrom: $(default_output_name(inputs.readFilesIn, "_"+self.join("_"), "."))
      outFilterMultimapNmax:
        default: 1
      outSAMtype:
        default: ["SAM"]
      threads: threads
    out:
      - aligned_file
      - log_final
      - log_out
      - log_progress
      - log_std
      - uniquely_mapped_reads_number

  strain1_sam_filter:
    run: ../tools/custom-bash.cwl
    in:
      input_file: insilico_star_aligner/aligned_file
      script:
        default: 'cat "$0" | grep "$1" | sed "s/$1/chr/g"  > `basename $0`'
      param:
        source: strain1
        valueFrom: $(self+"_")
    out: [output_file]

  strain2_sam_filter:
    run: ../tools/custom-bash.cwl
    in:
      input_file: insilico_star_aligner/aligned_file
      script:
        default: 'cat "$0" | grep "$1" | sed "s/$1/chr/g"  > `basename $0`'
      param:
        source: strain2
        valueFrom: $(self+"_")
    out: [output_file]

  strain1_samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: strain1_sam_filter/output_file
      sort_output_filename:
        source: [fastq_files, strain1]
        valueFrom: $(default_output_name(self[0], "_"+self[1], ".bam"))
      threads: threads
    out: [bam_bai_pair]

  strain2_samtools_sort_index:
    run: ../tools/samtools-sort-index.cwl
    in:
      sort_input: strain2_sam_filter/output_file
      sort_output_filename:
        source: [fastq_files, strain2]
        valueFrom: $(default_output_name(self[0], "_"+self[1], ".bam"))
      threads: threads
    out: [bam_bai_pair]

  strain1_project:
    run: ../tools/crossmap.cwl
    in:
      input_file_type:
        default: "bam"
      chain_file: strain1_chain_file
      input_file: strain1_samtools_sort_index/bam_bai_pair
    out: [projected_file]

  strain2_project:
    run: ../tools/crossmap.cwl
    in:
      input_file_type:
        default: "bam"
      chain_file: strain2_chain_file
      input_file: strain2_samtools_sort_index/bam_bai_pair
    out: [projected_file]

  strain1_bam_to_bigwig:
    run: bam-bedgraph-bigwig.cwl
    in:
      bam_file: strain1_project/projected_file
      chrom_length_file: reference_chrom_length_file
      mapped_reads_number:
        source: [fastq_files, insilico_star_aligner/uniquely_mapped_reads_number]
        valueFrom: |
          ${
            return (Array.isArray(self[0]) && self[0].length>1)?2*self[1]:self[1];
          }
    out: [bigwig_file]

  strain2_bam_to_bigwig:
    run: bam-bedgraph-bigwig.cwl
    in:
      bam_file: strain2_project/projected_file
      chrom_length_file: reference_chrom_length_file
      mapped_reads_number:
        source: [fastq_files, insilico_star_aligner/uniquely_mapped_reads_number]
        valueFrom: |
          ${
            return (Array.isArray(self[0]) && self[0].length>1)?2*self[1]:self[1];
          }
    out: [bigwig_file]

  reference_star_aligner:
    run: ../tools/star-alignreads.cwl
    in:
      readFilesIn: fastq_files
      genomeDir: reference_star_indices_folder
      outFilterMultimapNmax:
        default: 1
      threads: threads
    out:
      - aligned_file
      - uniquely_mapped_reads_number
      - log_final
      - log_out
      - log_progress
      - log_std

  reference_samtools_sort:
    run: ../tools/samtools-sort.cwl
    in:
      sort_input: reference_star_aligner/aligned_file
      sort_output_filename:
        source: fastq_files
        valueFrom: $(default_output_name(self, "", ".bam"))
      threads: threads
    out: [sorted_file]

  reference_bam_to_bigwig:
    run: bam-bedgraph-bigwig.cwl
    in:
      bam_file: reference_samtools_sort/sorted_file
      chrom_length_file: reference_chrom_length_file
      mapped_reads_number:
        source: [fastq_files, reference_star_aligner/uniquely_mapped_reads_number]
        valueFrom: |
          ${
            return (Array.isArray(self[0]) && self[0].length>1)?2*self[1]:self[1];
          }
    out: [bigwig_file]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "allele-vcf-alignreads-se-pe"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/workflows/allele-vcf-alignreads-se-pe.cwl
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
  Workflow maps FASTQ files from `fastq_files` input into reference genome `reference_star_indices_folder` and
  insilico generated `insilico_star_indices_folder` genome (concatenated genome for both `strain1` and `strain2` strains).
  For both genomes STAR is run with `outFilterMultimapNmax` parameter set to 1 to discard all of the multimapped reads.
  For insilico genome SAM file is generated. Then it's splitted into two SAM files based on strain names and then sorted
  by coordinates into the BAM format. For reference genome output BAM file from STAR slignment is also coordinate sorted.

s:about: |
  Workflow corresponds to MEA alignReads command from
  https://github.com/julienrichardalbert/MEA/blob/e3de228734bafd957cc2072dd8a6a0e84d554724/src/scripts/alignReads.sh
  Samtools quality and flag filtering of generated SAM/BAM files are replaced by outFilterMultimapNmax=1 parameter on
  mapping stage. Flag filtering 1540 should be clarified as long as it's not absolutely the same as STAR's implementation
  of outFilterMultimapNmax=1.