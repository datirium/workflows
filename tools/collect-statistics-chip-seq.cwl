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
  dockerPull: biowardrobe2/sc-tools:v0.0.29


inputs:

  script:
    type: string?
    default: |
        #!/usr/bin/env python3
        import os
        import sys
        import argparse
        import pandas
        import yaml
        import math

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
            "Quality encoding type selected": {
                "alias": "quality encoding type",
                "function": str,
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
            },
            "Sequences removed because they became shorter": {
                "alias": "reads/pairs removed because of length cutoff",
                "function": cut_int,
                "pair_end_specific": False
            },
            "order": ["fastq",
                    "trimming mode",
                    "adapter sequence",
                    "number of reads/pairs analysed for length validation",
                    "reads/pairs removed because of length cutoff",
                    "minimum required read length",
                    "quality phred score cutoff",
                    "quality encoding type",
                    "minimum required adapter overlap",
                    "maximum trimming error rate"]
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
            },
            "order": ["total reads/pairs processed",
                    "reads/pairs with at least one reported alignment",
                    "reads/pairs suppressed due to multimapping",
                    "reads/pairs unmapped"]
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
            "reads duplicated": {
                "alias": "reads/pairs duplicated",
                "function": cut_int,
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
            },
            "order": ["total reads/pairs",
                    "reads/pairs mapped",
                    "reads/pairs unmapped",
                    "reads/pairs duplicated",
                    "reads average length",
                    "reads maximum length",
                    "reads average quality",
                    "insert size average",
                    "insert size standard deviation"]
        }


        MACS2 = {
            "total fragments in treatment": {
                "alias": "total reads/pairs in treatment",
                "function": int,
                "pair_end_specific": False
            },
            "total tags in treatment": {
                "alias": "total reads/pairs in treatment",
                "function": int,
                "pair_end_specific": False
            },
            "fragments after filtering in treatment": {
                "alias": "reads/pairs after filtering in treatment",
                "function": int,
                "pair_end_specific": False
            },
            "tags after filtering in treatment": {
                "alias": "reads/pairs after filtering in treatment",
                "function": int,
                "pair_end_specific": False
            },
            "Redundant rate in treatment": {
                "alias": "redundant rate in treatment",
                "function": float,
                "pair_end_specific": False
            },
            "order": ["number of peaks called",
                    "mean peak size",
                    "total reads/pairs in treatment",
                    "reads/pairs after filtering in treatment",
                    "redundant rate in treatment",
                    "fraction of reads in peaks"]
        }


        def arg_parser():
            general_parser = argparse.ArgumentParser()
            general_parser.add_argument("--trim1",           help="Path to Trimgalore report file for FASTQ 1. Optional", required=False)
            general_parser.add_argument("--trim2",           help="Path to Trimgalore report file for FASTQ 2. Optional", required=False)
            general_parser.add_argument("--bowtie",          help="Path to Bowtie report file",                           required=True)
            general_parser.add_argument("--bamstats",        help="Path to bam statistics report file",                   required=True)
            general_parser.add_argument("--bamstatsfilter",  help="Path to bam statistics report file after filtering",   required=True)
            general_parser.add_argument("--macs2",           help="Path to MACS2 called peaks xls file",                  required=True)
            general_parser.add_argument("--atdp",            help="Path to ATDP output TSV file",                         required=True)
            general_parser.add_argument("--preseq",          help="Path to Preseq output file",                           required=False)
            general_parser.add_argument("--paired",          help="Process as paired-end. Default: False",                action="store_true")
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
            if collected_results[header]["trimming mode"] == "single-end":
                collected_results[header]["number of reads/pairs analysed for length validation"] = fastq["total reads processed"]
            collected_results[header]["fastq"].append(fastq)
            if TRIMGALORE.get("order", None):
                collected_results[header] = {k: collected_results[header][k] for k in TRIMGALORE["order"] if k in collected_results[header]}


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
            if key_dict.get("order", None):
                collected_results[header] = {k: collected_results[header][k] for k in key_dict["order"] if k in collected_results[header]}


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
            in_treatment = collected_results[header]["reads/pairs after filtering in treatment"]
            mapped = collected_results["BAM statistics after filtering"]["reads/pairs mapped"]
            collected_results[header]["fraction of reads in peaks"] = round(float(in_treatment)/float(mapped),2)
            if MACS2.get("order", None):
                collected_results[header] = {k: collected_results[header][k] for k in MACS2["order"] if k in collected_results[header]}


        def process_atdp_results(filepath, collected_results, header):
            if not collected_results.get(header, None):
                collected_results[header] = {}
            collected_results[header]["maximum"] = str(pandas.read_csv(filepath, sep="\t")["Y"].max())


        def process_preseq_results(filepath, collected_results, header, threashold=0.001):
            px, py = 0, 0
            for line in open_file(filepath):
                if "TOTAL_READS" in line:
                    continue
                values = [float(l.strip()) for l in line.strip().split()]
                dx, dy = values[0]-px, values[1]-py
                px, py = values[0], values[1]
                if dx != 0:
                    angle = math.degrees(math.atan2(dy, dx))
                    if angle <= threashold:
                        collected_results[header] = {"maximum library diversity": values[0]}
                        break


        def collect_stats(args):
            collected_results = {}
            if args.trim1:
                process_trimgalore_report(args.trim1, collected_results)
            if args.trim2:
                process_trimgalore_report(args.trim2, collected_results)
            process_custom_report(args.bowtie, collected_results, "alignment statistics", BOWTIE)
            process_custom_report(args.bamstats, collected_results, "BAM statistics", BAMSTATS, bool(args.paired))
            process_custom_report(args.bamstatsfilter, collected_results, "BAM statistics after filtering", BAMSTATS, bool(args.paired))
            process_custom_report(args.macs2, collected_results, "peak calling statistics", MACS2)
            process_macs2_xls(args.macs2, collected_results, "peak calling statistics")
            process_atdp_results(args.atdp, collected_results, "average tag density")
            if args.preseq:
                process_preseq_results(args.preseq, collected_results, "library preparation")
            return (collected_results)


        def export_results_yaml(collected_data, filepath):
            with open(filepath+".yaml", 'w') as output_stream:
                output_stream.write(yaml.dump(collected_data, width=1000, sort_keys=False))


        def export_results_markdown(collected_data, filepath):
            with open(filepath+".md", 'w') as output_stream:
                for line in yaml.dump(collected_data, width=1000, sort_keys=False).split("\n"):
                    if not line.strip():
                        continue
                    if line.startswith("  - "):
                        output_stream.write(line+"\n")
                    elif line.startswith("    "):
                        output_stream.write("<br>"+line+"\n")
                    elif line.startswith("  "):
                        output_stream.write("- "+line+"\n")
                    else:
                        output_stream.write("### "+line+"\n")


        def export_results_table(collected_data, filepath):
            with open(filepath+".tsv", 'w') as output_stream:
                total = collected_data["alignment statistics"]["total reads/pairs processed"]
                mapped = collected_data["BAM statistics after filtering"]["reads/pairs mapped"]
                multimapped = collected_data["alignment statistics"].get("reads/pairs suppressed due to multimapping", 0)
                unmapped = collected_data["alignment statistics"]["reads/pairs unmapped"]
                filtered = collected_data["BAM statistics"]["reads/pairs mapped"] - collected_data["BAM statistics after filtering"]["reads/pairs mapped"]
                header = [
                            "Tags total",
                            "Mapped",
                            "Multi-mapped",
                            "Unmapped",
                            "Filtered",

                            "alignment statistics",
                            "total reads/pairs processed",
                            "reads/pairs with at least one reported alignment",
                            "reads/pairs suppressed due to multimapping",
                            "reads/pairs unmapped",
                            
                            "BAM statistics",
                            "total reads/pairs",
                            "reads/pairs mapped",
                            "reads/pairs unmapped",
                            "reads/pairs duplicated",
                            "insert size average",
                            "insert size standard deviation",
                            "reads average length",
                            "reads average quality",
                            "reads maximum length",

                            "BAM statistics after filtering",
                            "total reads/pairs",
                            "reads/pairs mapped",
                            "reads/pairs unmapped",
                            "reads/pairs duplicated",
                            "insert size average",
                            "insert size standard deviation",
                            "reads average length",
                            "reads average quality",
                            "reads maximum length",
                
                            "peak calling statistics",
                            "number of peaks called",
                            "mean peak size",
                            "total reads/pairs in treatment",
                            "reads/pairs after filtering in treatment",
                            "redundant rate in treatment",
                            "fraction of reads in peaks",

                            "average tag density",
                            "maximum"]

                if collected_data.get("adapter trimming statistics", None):
                    header.extend(["adapter trimming statistics",
                                "trimming mode",
                                "adapter sequence",
                                "quality phred score cutoff",
                                "minimum required adapter overlap",
                                "maximum trimming error rate",
                                "minimum required read length",
                                "number of reads/pairs analysed for length validation",
                                "reads/pairs removed because of length cutoff",
                                "fastq1: total reads processed",
                                "fastq1: reads with adapters"])
                    if len(collected_data["adapter trimming statistics"]["fastq"]) > 1:
                        header.extend(["fastq2: total reads processed", "fastq2: reads with adapters"])

                output_stream.write("\t".join(header)+"\n")

                data = [
                        total,
                        mapped,
                        multimapped,
                        unmapped,
                        filtered,

                        "",
                        collected_data["alignment statistics"]["total reads/pairs processed"],
                        collected_data["alignment statistics"]["reads/pairs with at least one reported alignment"],
                        collected_data["alignment statistics"].get("reads/pairs suppressed due to multimapping", 0),
                        collected_data["alignment statistics"]["reads/pairs unmapped"],
                        
                        "",
                        collected_data["BAM statistics"]["total reads/pairs"],
                        collected_data["BAM statistics"]["reads/pairs mapped"],
                        collected_data["BAM statistics"]["reads/pairs unmapped"],
                        collected_data["BAM statistics"]["reads/pairs duplicated"],
                        collected_data["BAM statistics"]["insert size average"],
                        collected_data["BAM statistics"]["insert size standard deviation"],
                        collected_data["BAM statistics"]["reads average length"],
                        collected_data["BAM statistics"]["reads average quality"],
                        collected_data["BAM statistics"]["reads maximum length"],
                        
                        "",
                        collected_data["BAM statistics after filtering"]["total reads/pairs"],
                        collected_data["BAM statistics after filtering"]["reads/pairs mapped"],
                        collected_data["BAM statistics after filtering"]["reads/pairs unmapped"],
                        collected_data["BAM statistics after filtering"]["reads/pairs duplicated"],
                        collected_data["BAM statistics after filtering"]["insert size average"],
                        collected_data["BAM statistics after filtering"]["insert size standard deviation"],
                        collected_data["BAM statistics after filtering"]["reads average length"],
                        collected_data["BAM statistics after filtering"]["reads average quality"],
                        collected_data["BAM statistics after filtering"]["reads maximum length"],
                        
                        "",
                        collected_data["peak calling statistics"]["number of peaks called"],
                        collected_data["peak calling statistics"]["mean peak size"],
                        collected_data["peak calling statistics"]["total reads/pairs in treatment"],
                        collected_data["peak calling statistics"]["reads/pairs after filtering in treatment"],
                        collected_data["peak calling statistics"]["redundant rate in treatment"],
                        collected_data["peak calling statistics"]["fraction of reads in peaks"],

                        "",
                        collected_data["average tag density"]["maximum"]]

                if collected_data.get("adapter trimming statistics", None):
                    data.extend(["",
                                collected_data["adapter trimming statistics"]["trimming mode"],
                                collected_data["adapter trimming statistics"]["adapter sequence"],
                                collected_data["adapter trimming statistics"]["quality phred score cutoff"],
                                collected_data["adapter trimming statistics"]["minimum required adapter overlap"],
                                collected_data["adapter trimming statistics"]["maximum trimming error rate"],
                                collected_data["adapter trimming statistics"]["minimum required read length"],
                                collected_data["adapter trimming statistics"]["number of reads/pairs analysed for length validation"],
                                collected_data["adapter trimming statistics"]["reads/pairs removed because of length cutoff"],
                                collected_data["adapter trimming statistics"]["fastq"][0]["total reads processed"],
                                collected_data["adapter trimming statistics"]["fastq"][0]["reads with adapters"]])
                    if len(collected_data["adapter trimming statistics"]["fastq"]) > 1:
                        data.extend([collected_data["adapter trimming statistics"]["fastq"][1]["total reads processed"],
                                    collected_data["adapter trimming statistics"]["fastq"][1]["reads with adapters"]])

                output_stream.write("\t".join(str(l) for l in data))


        def main(argsl=None):
            if argsl is None:
                argsl = sys.argv[1:]
            args,_ = arg_parser().parse_known_args(argsl)
            args = normalize_args(args, skip_list=["paired"])
            collected_data = collect_stats(args)
            export_results_yaml(collected_data, args.output)
            export_results_table(collected_data, args.output)
            export_results_markdown(collected_data, args.output)


        if __name__ == "__main__":
            sys.exit(main(sys.argv[1:]))
    inputBinding:
      position: 5

  trimgalore_report_fastq_1:
    type: File?
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
      position: 8
      prefix: "--bowtie"

  bam_statistics_report:
    type: File
    inputBinding:
      position: 9
      prefix: "--bamstats"

  bam_statistics_after_filtering_report:
    type: File
    inputBinding:
      position: 10
      prefix: "--bamstatsfilter"

  macs2_called_peaks:
    type: File
    inputBinding:
      position: 11
      prefix: "--macs2"

  atdp_results:
    type: File
    inputBinding:
      position: 12
      prefix: "--atdp"

  preseq_results:
    type: File?
    inputBinding:
      position: 13
      prefix: "--preseq"

  paired_end:
    type: boolean?
    inputBinding:
      position: 14
      prefix: "--paired"

  output_prefix:
    type: string?
    inputBinding:
      position: 15
      prefix: "--output"
      valueFrom: $(get_output_prefix())
    default: ""


outputs:

  collected_statistics_yaml:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".yaml")

  collected_statistics_tsv:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".tsv")

  collected_statistics_md:
    type: File
    outputBinding:
      glob: $(get_output_prefix()+".md")

  mapped_reads:
    type: int
    outputBinding:
      loadContents: true
      glob: $(get_output_prefix()+".tsv")
      outputEval: $(parseInt(self[0].contents.split('\n')[1].split('\t')[1]))


baseCommand: [python3, '-c']


$namespaces:
  s: http://schema.org/

$schemas:
- https://github.com/schemaorg/schemaorg/raw/main/data/releases/11.01/schemaorg-current-http.rdf

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
