cwlVersion: v1.0
class: CommandLineTool


requirements:
- class: ShellCommandRequirement
- class: InlineJavascriptRequirement
  expressionLib:
  - var get_output_prefix = function() {
        if (inputs.output_prefix) {
          return inputs.output_prefix;
        }
        var root = inputs.bowtie_alignment_report.basename.split('.').slice(0,-1).join('.');
        var suffix = "_collected_statistics_report";
        return (root == "")?inputs.bowtie_alignment_report.basename+suffix:root+suffix;
    };


hints:
- class: DockerRequirement
  dockerPull: biowardrobe2/scidap:v0.0.3


inputs:

  script:
    type: string?
    default: |
        #!/usr/bin/env python
        import os
        import sys
        import argparse
        import json


        def cut_int(s):
            return int(s.strip().split()[0])


        def cut_float(s):
            return float(s.strip().split()[0])


        TRIMGALORE = {
            "Trimming mode": {
                "alias": "trimming mode",
                "function": str,
                "pair_end_specific": False
            },
            "Quality Phred score cutoff": {
                "alias": "quality phred score cutoff",
                "function": int,
                "pair_end_specific": False
            },
            "Adapter sequence": {
                "alias": "adapter sequence",
                "function": str,
                "pair_end_specific": False
            },
            "Minimum required adapter overlap": {
                "alias": "minimum required adapter overlap",
                "function": cut_int,
                "pair_end_specific": False
            },
            "Maximum trimming error rate": {
                "alias": "maximum trimming error rate",
                "function": cut_float,
                "pair_end_specific": False
            },
            "Minimum required sequence length": {
                "alias": "minimum required read length",
                "function": cut_int,
                "pair_end_specific": False
            },
            "Total number of sequences analysed": {
                "alias": "number of reads/pairs analysed for length validation",
                "function": int,
                "pair_end_specific": False
            },
            "Number of sequence pairs removed": {
                "alias": "reads/pairs removed because of length cutoff",
                "function": cut_int,
                "pair_end_specific": False
            }
        }


        BOWTIE = {
            "reads processed": {
                "alias": "total reads/pairs processed",
                "function": int,
                "pair_end_specific": False
            },
            "reads with at least one reported alignment": {
                "alias": "reads/pairs with at least one reported alignment",
                "function": cut_int,
                "pair_end_specific": False
            },
            "reads that failed to align": {
                "alias": "reads/pairs unmapped",
                "function": cut_int,
                "pair_end_specific": False
            },
            "reads with alignments suppressed due to -m": {
                "alias": "reads/pairs suppressed due to multimapping",
                "function": cut_int,
                "pair_end_specific": False
            }
        }


        BAMSTATS = {
            "raw total sequences": {
                "alias": "total reads/pairs",
                "function": int,
                "pair_end_specific": True
            },
            "reads mapped": {
                "alias": "reads/pairs mapped",
                "function": int,
                "pair_end_specific": True
            },
            "reads unmapped": {
                "alias": "reads/pairs unmapped",
                "function": int,
                "pair_end_specific": True
            },
            "average length": {
                "alias": "reads average length",
                "function": float,
                "pair_end_specific": False
            },
            "maximum length": {
                "alias": "reads maximum length",
                "function": int,
                "pair_end_specific": False
            },
            "average quality": {
                "alias": "reads average quality",
                "function": float,
                "pair_end_specific": False
            },
            "insert size average": {
                "alias": "insert size average",
                "function": float,
                "pair_end_specific": False
            },
            "insert size standard deviation": {
                "alias": "insert size standard deviation",
                "function": float,
                "pair_end_specific": False
            }
        }


        MACS2 = {
            "total fragments in treatment": {
                "alias": "total fragments in treatment",
                "function": int,
                "pair_end_specific": False
            },
            "total tags in treatment": {
                "alias": "total tags in treatment",
                "function": int,
                "pair_end_specific": False
            },
            "fragments after filtering in treatment": {
                "alias": "fragments after filtering in treatment",
                "function": int,
                "pair_end_specific": False
            },
            "tags after filtering in treatment": {
                "alias": "tags after filtering in treatment",
                "function": int,
                "pair_end_specific": False
            },
            "Redundant rate in treatment": {
                "alias": "redundant rate in treatment",
                "function": float,
                "pair_end_specific": False
            }
        }


        def arg_parser():
            general_parser = argparse.ArgumentParser()
            general_parser.add_argument("--trim1",           help="Path to Trimgalore report file for FASTQ 1",           required=True)
            general_parser.add_argument("--trim2",           help="Path to Trimgalore report file for FASTQ 2. Optional", required=False)
            general_parser.add_argument("--bowtie",          help="Path to Bowtie report file",                           required=True)
            general_parser.add_argument("--bamstats",        help="Path to bam statistics report file",                   required=True)
            general_parser.add_argument("--bamstatsfilter",  help="Path to bam statistics report file after filtering",   required=True)
            general_parser.add_argument("--macs2",           help="Path to MACS2 called peaks xls file",                  required=True)
            general_parser.add_argument("--output",          help="Output filename prefix",                               required=True)
            return general_parser


        def normalize_args(args, skip_list=[]):
            normalized_args = {}
            for key,value in args.__dict__.items():
                if key not in skip_list:
                    normalized_args[key] = value if not value or os.path.isabs(value) else os.path.normpath(os.path.join(os.getcwd(), value))
                else:
                    normalized_args[key]=value
            return argparse.Namespace (**normalized_args)


        def open_file(filename):
            lines = []
            with open(filename, 'r') as infile:
                for line in infile:
                    if line.strip():
                        lines.append(line.strip())
            return lines


        def split_line(l, delimiter=":"):
            return [i.strip() for i in l.split(delimiter)]


        def get_correspondent_key(data, long_k):
            for short_k, v in data.items():
                if short_k in long_k:
                    return v["alias"], v["function"], v["pair_end_specific"]
            raise Exception


        def process_trimgalore_report(filepath, collected_results, header="adapter trimming statistics"):
            if not collected_results.get(header, None):
                collected_results[header] = {"fastq": []}
            fastq = {}
            for line in open_file(filepath):
                try:
                    key, value = split_line(line)
                    if "Total reads processed" in line:
                        fastq["total reads processed"] = int(value.strip().replace(",",""))
                    elif "Reads with adapters" in line:
                        fastq["reads with adapters"] = int(value.strip().replace(",","").split()[0])
                    else:
                        res_key, res_function, pair_end_specific = get_correspondent_key(TRIMGALORE, key)
                        if not collected_results[header].get(res_key, None):
                            collected_results[header][res_key] = res_function(value)
                except Exception:
                    pass
            collected_results[header]["fastq"].append(fastq)


        def process_custom_report(filepath, collected_results, header, key_dict, pair_end=False):
            if not collected_results.get(header, None):
                collected_results[header] = {}
            for line in open_file(filepath):
                try:
                    key, value = split_line(line)
                    res_key, res_function, pair_end_specific = get_correspondent_key(key_dict, key)
                    if not collected_results[header].get(res_key, None):
                        if pair_end_specific and pair_end:
                            collected_results[header][res_key] = res_function(res_function(value)/2)
                        else:
                            collected_results[header][res_key] = res_function(value)
                except Exception:
                    pass


        def process_macs2_xls(filepath, collected_results, header):
            if not collected_results.get(header, None):
                collected_results[header] = {}
            count, length, prev_start, prev_end, prev_chr = 0, 0, 0, 0, ""
            for line in open_file(filepath):
                if "#" in line or "start" in line:
                    continue
                line_list = [l.strip() for l in line.strip().split()]
                chr = line_list[0]
                start = int(line_list[1])
                end = int(line_list[2])
                if start == prev_start and end == prev_end and chr == prev_chr:
                    continue
                count = count + 1
                length = length + end - start
                prev_start, prev_end, prev_chr = start, end, chr
            collected_results[header]["number of peaks called"] = count
            collected_results[header]["mean peak size"] = round(float(length)/float(count), 2)


        def collect_stats(args):
            collected_results = {}
            process_trimgalore_report(args.trim1, collected_results)
            if args.trim2:
                process_trimgalore_report(args.trim2, collected_results)
            process_custom_report(args.bowtie, collected_results, "alignment statistics", BOWTIE)
            process_custom_report(args.bamstats, collected_results, "BAM statistics", BAMSTATS, bool(args.trim2))
            process_custom_report(args.bamstatsfilter, collected_results, "BAM statistics after filtering", BAMSTATS, bool(args.trim2))
            process_custom_report(args.macs2, collected_results, "Peak calling statistics", MACS2)
            process_macs2_xls(args.macs2, collected_results, "Peak calling statistics")    
            return (collected_results)


        def export_results_json(collected_data, filepath):
            with open(filepath+".json", 'w') as output_stream:
                output_stream.write(json.dumps(collected_data, indent=4))


        def export_results_table(collected_data, filepath):
            with open(filepath+".tsv", 'w') as output_stream:
                total = collected_data["alignment statistics"]["total reads/pairs processed"]
                mapped = collected_data["BAM statistics after filtering"]["reads/pairs mapped"]
                multimapped = collected_data["alignment statistics"].get("reads/pairs suppressed due to multimapping", 0)
                unmapped = collected_data["alignment statistics"]["reads/pairs unmapped"]
                filtered = collected_data["BAM statistics"]["reads/pairs mapped"] - collected_data["BAM statistics after filtering"]["reads/pairs mapped"]
                output_stream.write("Tags total\tMapped\tMulti-mapped\tUnmapped\tFiltered\n")
                output_stream.write("\t".join(str(l) for l in [total, mapped, multimapped, unmapped, filtered]))


        def main(argsl=None):
            if argsl is None:
                argsl = sys.argv[1:]
            args,_ = arg_parser().parse_known_args(argsl)
            args = normalize_args(args)
            export_results_json(collect_stats(args), args.output)
            export_results_table(collect_stats(args), args.output)


        if __name__ == "__main__":
            sys.exit(main(sys.argv[1:]))
    inputBinding:
      position: 5

  trimgalore_report_fastq_1:
    type: File
    inputBinding:
      position: 6
      prefix: "--trim1"

  trimgalore_report_fastq_2:
    type: File?
    inputBinding:
      position: 7
      prefix: "--trim2"

  bowtie_alignment_report:
    type: File
    inputBinding:
      position: 9
      prefix: "--bowtie"

  bam_statistics_report:
    type: File
    inputBinding:
      position: 8
      prefix: "--bamstats"

  bam_statistics_after_filtering_report:
    type: File
    inputBinding:
      position: 11
      prefix: "--bamstatsfilter"

  macs2_called_peaks:
    type: File
    inputBinding:
      position: 13
      prefix: "--macs2"

  output_prefix:
    type: string?
    inputBinding:
      position: 16
      prefix: "--output"
      valueFrom: $(get_output_prefix())
    default: ""


outputs:

  collected_statistics_json:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".json")

  collected_statistics_tsv:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".tsv")

  mapped_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_prefix()+".tsv")
      outputEval: $(parseInt(self[0].contents.split('\n')[1].split('\t')[1]))


baseCommand: [python, '-c']


$namespaces:
  s: http://schema.org/

$schemas:
- http://schema.org/docs/schema_org_rdfa.html

s:name: "collect-statistics-chip-seq-trim"
s:downloadUrl: https://raw.githubusercontent.com/Barski-lab/workflows/master/tools/collect-statistics-chip-seq-trim.cwl
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
  Tool processes and combines log files generated by Trimgalore, Bowtie, Samtools and MACS2.

  `get_output_prefix` function returns output file prefix equal to `output_prefix`+`_collected_statistics_report` (if this input is provided) or
  generated on the base of bowtie log basename with `_collected_statistics_report` extension.


s:about: |
  Runs python code from the script input
