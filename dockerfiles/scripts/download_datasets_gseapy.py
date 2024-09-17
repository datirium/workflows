#!/usr/bin/env python
import os
import sys
import argparse
from gseapy.utils import retry, mkdirs


# import json
# import requests
# def get_libraries():
#     lib_url='http://amp.pharm.mssm.edu/Enrichr/datasetStatistics'
#     response = requests.get(lib_url)
#     if not response.ok:
#         raise Exception("Error getting the Enrichr libraries")
#     libs_json = json.loads(response.text)
#     libs = [lib['libraryName'] for lib in libs_json['statistics']]
#     return sorted(libs)


def parse_lib(libname, suffix, location):
    tmpname = "enrichr." + suffix + ".gmt"
    tempath = os.path.join(location, tmpname)
    if os.path.isfile(tempath):
        print("Already downloaded", libname)
    else:
        return download_library(libname, location)


def download_library(libname, location):
    print("Downloading: ", libname)
    s = retry(5)
    ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary'
    query_string = '?mode=text&libraryName=%s'
    response = s.get( ENRICHR_URL + query_string % libname, timeout=None)
    if not response.ok:
        raise Exception('Error fetching enrichment results, check internet connection first.')
    mkdirs(location)
    genesets_dict = {}
    outname = "enrichr.%s.gmt"%libname
    gmtout = open(os.path.join(location, outname), "w")
    for line in response.iter_lines(chunk_size=1024, decode_unicode='utf-8'):
        line=line.strip()
        k = line.split("\t")[0]
        v = list(map(lambda x: x.split(",")[0], line.split("\t")[2:]))
        genesets_dict.update({ k: v})
        outline = "%s\t\t%s\n"%(k, "\t".join(v))
        gmtout.write(outline)
    gmtout.close()
    return genesets_dict


class ArgsParser():

    def __init__(self, args):
        self.args, _ = self.get_parser().parse_known_args(args)
        self.set_args_as_attributes()

    def set_args_as_attributes(self):
        for arg, value in vars(self.args).items():
            setattr(self, arg, value)

    def get_parser(self):
        general_parser = argparse.ArgumentParser()
        general_parser.add_argument(
            "--libraries",
            help="Name of GSEAPY genesets to select for downloading",
            type=str, required=True, nargs="+"
        )
        general_parser.add_argument(
            "--names",
            help="Names to save downloaded genesets",
            type=str, required=True, nargs="+"
        ),
        general_parser.add_argument(
            "--output",
            help="Path to the folder to save downloaded files",
            type=str, default="./"
        )
        return general_parser


def main(args=None):
    args = ArgsParser(sys.argv[1:] if args is None else args)
    for library, suffix in zip(args.libraries, args.names):
        try:
            print(f"Attempting to download {library} with the name {suffix}")
            parse_lib(library, suffix, args.output)
        except Exception as err:
            print(f"Failed to download {library} with the name {suffix}")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))