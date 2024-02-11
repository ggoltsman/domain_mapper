#!/usr/bin/env python


"""
domain_mapper.py
Author:  Gene Goltsman, 2024

This program maps protein coordinates to genomic coordinates for protein domains associated with a gene trascript.  

The inputs are:
1. A four-column (tab-separated) file containing information about the protein domains and the
associated genes. The first column is the gene symbol. The second column is the transcript
id. The third column is the genome build. The fourth column contains the protein domain(s)
found in that gene. The format of the fourth column is as follows

“[domain_name]:[aa_start]-[aa_end]"

Multiple domains should be separated by a semicolon, e.g.:

“[domain_name]:[aa_start]-[aa_end];[domain_name]:[aa_start]-[aa_end];"

Example:
"Zincfinger,RING-type:24-658;BRCA1,serine-richdomain:345-507;”


2. A two-column *list* of feature files corresponding to the genome builds referenced in the domain file. 
Each feature file must be in gff3 format and contain the genomic features for a specific Human genome build. 
The first column is the path to the gff3 file, the second column is the build name as it appears in the domain table. 
Example:

/path/to/Homo_sapiens.GRCh37.87.chr.gff3.gz GRCh37
/path/to/Homo_sapiens.GRCh38.111.chr.gff3.gz GRCh38


The output is a tab-delimited file with the following columns:

Column 1: The gene name.
Column 2: The genome build.
Column 3: The chromosome.
Column 4: The protein domain name
Column 5: The protein coordinates of the domain. (Format: start-end)
Column 6: The amino acid length of the domain.
Column 7: The genomic coordinates of the domain. (Format: start-end)
Column 8: The length of the domain in genomic space.
"""


import sys
import gffpandas.gffpandas as gffpd
from dataclasses import dataclass
import argparse


@dataclass
class Domain:
    name: str
    aa_start: int
    aa_end: int
    aa_len: int
    genomic_start = 0
    genomic_end = 0
    genomic_len = 0


class Transcript:
    def __init__(self, gene, trID, build, domains_str):
        self.gene = gene
        self.tr_ID = trID
        self.build = build
        self.domains = self.parse_domains_str(domains_str)

        self.chrom = ""
        self.strand = ""
        self.cds_coords = []

    def parse_domains_str(self, domains_str):
        """
        Parses the domains string and extracts aa coordinates on the protein
        Returns a dictionary of Domain objects containg the aa coords
        """

        domains = {}
        if domains_str[0] != '"' or domains_str[-1] != '"':
            raise ValueError(
                "Illegal domain info string; must be in double quotes", domains_str
            )
        domainInfos = domains_str.strip('"').rstrip(";").split(";")
        for dInfoStr in domainInfos:
            dName, dCoordsStr = dInfoStr.split(":")
            aa_coords = dCoordsStr.split("-")
            aa_coords = [int(i) for i in aa_coords]
            aa_len = aa_coords[1] - aa_coords[0] + 1
            domains[dName] = Domain(dName, aa_coords[0], aa_coords[1], aa_len)
        return domains

    def parse_cds(self, CDSBuilds):
        """
        Queries the cds dataframe for entries matching the transcript id.
        """

        if self.build not in CDSBuilds.keys():
            raise ValueError(
                "Build name in the input file not matiching any of the gffs loaded",
                self.build,
            )

        cds_df = CDSBuilds[self.build]
        if not cds_df.Parent.str.startswith('transcript:').all():
            raise ValueError(
                "At least one CDS feature has an invalid `Parent` attribute.  Please make sure all have the format 'transcript:[transcript id]' "
            )

        query_by_transcript = f'Parent == "transcript:{self.tr_ID}"'
        my_transcript_df = cds_df.query(query_by_transcript).reset_index()
        if not my_transcript_df.shape[0]:
            raise ValueError(
                "No match found in the gff file for transcript ", self.tr_ID, 
                "; Note: CDS 'Parent' attributes must have the format 'transcript:[transcript id]' "
            )

        self.chrom = my_transcript_df.seq_id[0]
        self.strand = my_transcript_df.strand[0]
        self.cds_coords = list(zip(my_transcript_df.start, my_transcript_df.end))

    def map_domains_to_genome(self):
        """
        Indexes all cds positions, then, after converting the aa coords to nt offsets relative to the transcript,
        accesses the domain positions by index to get the genomic coords.
        Updates the Domain objects with genomic coords
        """

        tr_idx = []
        for s, e in self.cds_coords:
            tr_idx.extend(list(range(s, e + 1)))
        assert (len(tr_idx) / 3).is_integer(), "Sum of CDS lengths must be divisible by three to be a valid transcript"

        if self.strand == "-":
            # if the transcript is on the reverse strand, we simply reverse the order of the cds positions so that tr_idx[0] is the last genomic position of the trancript
            tr_idx.reverse()

        for domain in self.domains.values():
            # domain offsets are converted from aa to nt
            domain_start_offset = (
                domain.aa_start * 3 - 2
            )  # we want the offest to the 1st nucleotide of the 1st codon of the domain
            domain_end_offset = domain.aa_end * 3

            # get genomic start,end. Regardless of the transcript strand, we report the positions in the nominal left-to-right orider
            gen_coords = [tr_idx[domain_start_offset - 1], 
                          tr_idx[domain_end_offset - 1]]
            domain.genomic_start, domain.genomic_end = sorted(gen_coords)
            domain.genomic_len = domain.genomic_end - domain.genomic_start + 1

    def print_mapped_domain_info(self, out_fh):
        """
        Writes transcript and domain info to the provided file handle
        header = ['#Gene','Build','Chrom','Dom_name', 'Dom_aa_coords','Dom_aa_len','Dom_genomic_coords',  'Dom_genomic_len']
        """
        for dm in self.domains.values():
            fields = [
                self.gene,
                self.build,
                self.chrom,
                dm.name,
                str(dm.aa_start) + "-" + str(dm.aa_end),
                str(dm.aa_len),
                str(dm.genomic_start) + "-" + str(dm.genomic_end),
                str(dm.genomic_len),
            ]
            
            # print('\t'.join(fields))
            out_fh.write("\t".join(fields) + "\n")


def load_transcripts(domain_file):
    transcripts = []
    with open(domain_file, "r", encoding="utf-8") as domf:
        for line in domf:
            if line.startswith("#"):
                continue
            ll = line.rstrip().split("\t")
            if len(ll) != 4 or any(v is None for v in ll):
                raise ValueError("Malformatted input entry:", line)

            gene, trID, build, domains_str = ll
            transcripts.append(Transcript(gene, trID, build, domains_str))

    return transcripts


def parse_gff(path):
    """
    Loads a gff3 file into a gffpd object,
    filters it to only store CDS entries and parses the attribute column into separate columns.
    Returns a pandas dataframe
    """

    gff = gffpd.read_gff3(path)
    cds_df = (
        gff.filter_feature_of_type(["CDS"])
        .attributes_to_columns()
        .astype({"seq_id": "str"})
    )
    assert cds_df.shape[0], "No CDS data could be extracted from the gff file"
    cds_df.drop(["attributes"], axis=1, inplace=True)
    return cds_df


def parse_gffs(gff_list):
    """
    Reads a *list* of gff files and build names
    Returns parsed and filtered dictionary of pandas dataframes, with the build names as key.
    """
    gff_cds_builds = {}
    with open(gff_list, "r", encoding="utf-8") as gl:
        for line in gl:
            path, build = line.rstrip().split("\t")
            df = parse_gff(path)
            gff_cds_builds[build] = df
    return gff_cds_builds


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d",
        "--domain_file",
        type=str,
        required=True,
        help="A file containing information about the protein domains and the associated genes and transcripts",
    )
    parser.add_argument(
        "-g",
        "--gff_list",
        type=str,
        required=True,
        help="A list of paths to feature files (in gff3 format), with corresponding Human genome build names specified",
    )
    parser.add_argument(
        "-o", "--output_file", type=str, required=True, help="Output file name"
    )

    if not len(sys.argv) > 1:
        parser.print_help()
        parser.exit()
    args = parser.parse_args()

    out = open(args.output_file, "w")
    header = [
        "#Gene",
        "Build",
        "Chrom",
        "Dom_name",
        "Dom_aa_coords",
        "Dom_aa_len",
        "Dom_genomic_coords",
        "Dom_genomic_len",
    ]
    out.write("\t".join(header) + "\n")

    transcripts = load_transcripts(args.domain_file)
    gffCDSBuilds = parse_gffs(args.gff_list)

    for t in transcripts:
        t.parse_cds(gffCDSBuilds)
        t.map_domains_to_genome()
        t.print_mapped_domain_info(out)

    print("Done")

    out.close()


if __name__ == "__main__":
    main()
