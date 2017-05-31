cwlVersion: v1.0
class: Workflow

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: StepInputExpressionRequirement
  - class: InlineJavascriptRequirement
  - class: MultipleInputFeatureRequirement

inputs:

  fastq_input_file:
    type: File
    label: "FASTQ input file"
    format: "http://edamontology.org/format_1930"
    doc: "Reads data in a FASTQ format, received after single end sequencing"

  bowtie_indices_folder:
    type: Directory
    label: "BOWTIE indices folder"
    doc: "Path to BOWTIE generated indices folder"

  clip_3p_end:
    type: int
    label: "Clip from 3p end"
    doc: "Number of bases to clip from the 3p end"

  clip_5p_end:
    type: int
    label: "Clip from 5p end"
    doc: "Number of bases to clip from the 5p end"

  threads:
    type: int?
    label: "Threads"
    doc: "Number of threads for those steps that support multithreading"

  remove_duplicates:
    type: boolean
    label: "Remove duplicates"
    doc: "Calls samtools rmdup to remove duplicates from sortesd BAM file"

  control_file:
    type: File?
    label: "Control BAM file"
    format: "http://edamontology.org/format_2572"
    doc: "Control BAM file file for MACS2 peak calling"

  exp_fragment_size:
    type: int
    label: "Expected fragment size"
    doc: "Expected fragment size for MACS2"

  force_fragment_size:
    type: boolean
    label: "Force fragment size"
    doc: "Force MACS2 to use exp_fragment_size"

  broad_peak:
    type: boolean
    label: "Callpeak broad"
    doc: "Set to call broad peak for MACS2"

  chrom_length:
    type: File
    label: "Chromosome length file"
    format: "http://edamontology.org/format_2330"
    doc: "Chromosome length file"

  genome_size:
    type: string
    label: "Effective genome size"
    doc: "MACS2 effective genome size: hs, mm, ce, dm or number, for example 2.7e9"

outputs:

  bigwig:
    type: File
    format: "http://edamontology.org/format_3006"
    label: "BigWig file"
    doc: "Generated BigWig file"
    outputSource: bam_to_bigwig/outfile

  fastx_statistics:
    type: File
    label: "FASTQ statistics"
    format: "http://edamontology.org/format_2330"
    doc: "fastx_quality_stats generated FASTQ file quality statistics file"
    outputSource: fastx_quality_stats/statistics

  bowtie_log:
    type: File
    label: "BOWTIE alignment log"
    format: "http://edamontology.org/format_2330"
    doc: "BOWTIE generated alignment log"
    outputSource: bowtie_aligner/output_bowtie_log

  samtools_rmdup_log:
    type: File
    label: "Remove duplicates log"
    format: "http://edamontology.org/format_2330"
    doc: "Samtools rmdup generated log"
    outputSource: samtools_rmdup/rmdup_log

  bambai_pair:
    type: File
    format: "http://edamontology.org/format_2572"
    label: "Coordinate sorted BAM alignment file (+index BAI)"
    doc: "Coordinate sorted BAM file and BAI index file"
    outputSource: samtools_sort_index_after_rmdup/bam_bai_pair

  macs_called_peaks:
    type: File?
    label: "Called peaks"
    format: "http://edamontology.org/format_3468"
    doc: "XLS file to include information about called peaks"
    outputSource: macs2_callpeak_forced/peak_xls_file

  macs_narrow_peaks:
    type: File?
    label: "Narrow peaks"
    format: "http://edamontology.org/format_3613"
    doc: "Contains the peak locations together with peak summit, pvalue and qvalue"
    outputSource: macs2_callpeak_forced/narrow_peak_file

  macs_broad_peaks:
    type: File?
    label: "Broad peaks"
    format: "http://edamontology.org/format_3614"
    doc: "Contains the peak locations together with peak summit, pvalue and qvalue"
    outputSource: macs2_callpeak_forced/broad_peak_file

  macs_peak_summits:
    type: File?
    label: "Peak summits"
    format: "http://edamontology.org/format_3003"
    doc: "Contains the peak summits locations for every peaks"
    outputSource: macs2_callpeak_forced/peak_summits_file

  macs_moder_r:
    type: File?
    label: "MACS2 generated R script"
    format: "http://edamontology.org/format_2330"
    doc: "R script to produce a PDF image about the model based on your data"
    outputSource: macs2_callpeak_forced/moder_r_file

  macs_gapped_peak:
    type: File?
    label: "Gapped peak"
    format: "http://edamontology.org/format_3586"
    doc: "Contains both the broad region and narrow peaks"
    outputSource: macs2_callpeak_forced/gapped_peak_file

  macs_log:
    type: File?
    label: "MACS2 log"
    format: "http://edamontology.org/format_2330"
    doc: "MACS2 output log"
    outputSource: macs2_callpeak_forced/macs_log

  get_stat_log:
    type: File?
    label: "Bowtie & Samtools Rmdup combined log"
    format: "http://edamontology.org/format_2330"
    doc: "Processed and combined Bowtie aligner and Samtools rmdup log"
    outputSource: get_stat/output

steps:

  fastx_quality_stats:
    run: ../../tools/fastx-quality-stats.cwl
    in:
      input_file: fastq_input_file
    out: [statistics]

  bowtie_aligner:
    run: ../../tools/bowtie.cwl
    in:
      filelist: fastq_input_file
      indices_folder: bowtie_indices_folder
      clip_3p_end: clip_3p_end
      clip_5p_end: clip_5p_end
      v:
        default: 3
      m:
        default: 1
      best:
        default: true
      strata:
        default: true
      sam:
        default: true
      threads: threads
      q:
        default: true
    out: [output, output_bowtie_log]

  samtools_sort_index:
    run: ../../tools/samtools-sort-index.cwl
    in:
      sort_input: bowtie_aligner/output
      threads: threads
    out: [bam_bai_pair]

  samtools_rmdup:
    run: ../../tools/samtools-rmdup.cwl
    in:
      trigger: remove_duplicates
      input_file: samtools_sort_index/bam_bai_pair
      single_end:
        default: true
    out: [rmdup_output, rmdup_log]

  samtools_sort_index_after_rmdup:
    run: ../../tools/samtools-sort-index.cwl
    in:
      trigger: remove_duplicates
      sort_input: samtools_rmdup/rmdup_output
      threads: threads
    out: [bam_bai_pair]

  macs2_callpeak:
    run: ../../tools/macs2-callpeak.cwl
    in:
      treatment: samtools_sort_index_after_rmdup/bam_bai_pair
      control: control_file
      nolambda:
        source: control_file
        valueFrom: |
          ${
            return !Boolean(self);
          }
      genome_size: genome_size
      mfold:
        default: "4 40"
      verbose:
        default: 3
      nomodel:
        source: force_fragment_size
        valueFrom: |
          ${
            return Boolean(self);
          }
      extsize:
        source: [force_fragment_size, exp_fragment_size]
        valueFrom: |
          ${
            if (self[0]){
              return self[1];
            } else {
              return null;
            }
          }
      bw:
        source: [force_fragment_size, exp_fragment_size]
        valueFrom: |
          ${
            if (!self[0]){
              return self[1];
            } else {
              return null;
            }
          }
      broad: broad_peak
      call_summits:
        source: broad_peak
        valueFrom: $(!self)
      keep_dup:
        default: auto
      q_value:
        default: 0.05
      format:
        default: BAM
      buffer_size:
        default: 10000
    out:
      - peak_xls_file
      - narrow_peak_file
      - peak_summits_file
      - broad_peak_file
      - moder_r_file
      - gapped_peak_file
      - treat_pileup_bdg_file
      - control_lambda_bdg_file
      - macs_log

  macs_island_count:
    run: ../../tools/macs2-island-count.cwl
    in:
      input_file: macs2_callpeak/peak_xls_file
    out: [fragments, islands]

  macs2_callpeak_forced:
    run: ../../tools/macs2-callpeak.cwl
    in:
      trigger:
        source: [force_fragment_size, macs_island_count/fragments]
        valueFrom: |
          ${
            return !self[0] && parseInt(self[1]) < 80;
          }
      peak_xls_file_staged: macs2_callpeak/peak_xls_file
      narrow_peak_file_staged: macs2_callpeak/narrow_peak_file
      broad_peak_file_staged: macs2_callpeak/broad_peak_file
      gapped_peak_file_staged: macs2_callpeak/gapped_peak_file
      peak_summits_file_staged: macs2_callpeak/peak_summits_file
      moder_r_file_staged: macs2_callpeak/moder_r_file
      treat_pileup_bdg_file_staged: macs2_callpeak/treat_pileup_bdg_file
      control_lambda_bdg_file_staged: macs2_callpeak/control_lambda_bdg_file
      treatment: samtools_sort_index_after_rmdup/bam_bai_pair
      control: control_file
      nolambda:
        source: control_file
        valueFrom: |
          ${
            return !Boolean(self);
          }
      genome_size: genome_size
      mfold:
        default: "4 40"
      verbose:
        default: 3
      nomodel:
        default: true
      extsize: exp_fragment_size
      broad: broad_peak
      call_summits:
        source: broad_peak
        valueFrom: $(!self)
      keep_dup:
        default: auto
      q_value:
        default: 0.05
      format:
        default: BAM
      buffer_size:
        default: 10000
    out:
      - peak_xls_file
      - narrow_peak_file
      - peak_summits_file
      - broad_peak_file
      - moder_r_file
      - gapped_peak_file
      - treat_pileup_bdg_file
      - control_lambda_bdg_file
      - macs_log

  bamtools_stats:
    run: ../../tools/bamtools-stats.cwl
    in:
      input_files: samtools_sort_index_after_rmdup/bam_bai_pair
    out: [mappedreads]

  bam_to_bigwig:
    run: bam-genomecov-bigwig.cwl
    in:
      input: samtools_sort_index_after_rmdup/bam_bai_pair
      genomeFile: chrom_length
      mappedreads: bamtools_stats/mappedreads
    out: [outfile]

  get_stat:
      run: ../../tools/python-get-stat.cwl
      in:
        bowtie_log: bowtie_aligner/output_bowtie_log
        rmdup_log: samtools_rmdup/rmdup_log
      out: [output]

$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "run-dna-single-end"
s:downloadUrl: https://raw.githubusercontent.com/SciDAP/workflows/master/workflows/scidap/run-dna-single-end.cwl
s:codeRepository: https://github.com/SciDAP/workflows
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

doc: |
  Current workflow is used to run CHIP-Seq basic analysis with single-end input FASTQ file.
  In outputs it returns coordinate sorted BAM file alongside with index BAI file, quality
  statistics of the input FASTQ file, reads coverage in a form of BigWig file, peaks calling
  data in a form of narrowPeak or broadPeak files.

s:about: >
  Current workflow is used to run CHIP-Seq basic analysis with single-end input FASTQ file.
  In outputs it returns coordinate sorted BAM file alongside with index BAI file, quality
  statistics of the input FASTQ file, reads coverage in a form of BigWig file, peaks calling
  data in a form of narrowPeak or broadPeak files.
  Workflow starts with running fastx_quality_stats (Step fastx_quality_stats) from FASTX-Toolkit
  to calculate quality statistics for input FASTQ file. At the same time (Step bowtie_aligner)
  Bowtie is used to align reads from input FASTQ file to reference genome. The output of this step
  is unsorted SAM file which is being sorted and indexed by samtools sort and samtools index
  (Step samtools_sort_index). Depending on workflow’s input parameters indexed and sorted BAM file
  could be processed by samtools rmdup (Step samtools_rmdup) to remove all possible read duplicates.
  In a case when removing duplicates is not necessary the step returns original input BAM and BAI
  files without any processing. If the duplicates were removed the following step
  (Step samtools_sort_index_after_rmdup) reruns samtools sort and samtools index with BAM and BAI files,
  if not - the step returns original unchanged input files. Right after that (Step macs2_callpeak)
  macs2 callpeak performs peak calling. On the base of returned outputs the next step
  (Step macs_island_count) calculates the number and islands and estimated fragment size. If the last
  one is less that 80 (hardcoded in a workflow) macs2 callpeak is rerun again with forced fixed
  fragment size value (Step macs2_callpeak_forced). If at the very beginning it was set in workflow
  input parameters to run peak calling with fixed fragment size, this step is skipped and the
  original peak calling results are saved. The following two steps (Step bamtools_stats and bam_to_bigwig)
  are used to calculate coverage on the base of input BAM file and save it in BigWig format.
  For that purpose bamtools stats returns the number of mapped reads number which is then used as
  scaling factor by bedtools genomecov when it performs coverage calculation and saves
  it in BED format. The last one is then being sorted and converted to BigWig format by
  bedGraphToBigWig tool from UCSC utilities.
