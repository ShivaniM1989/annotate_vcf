import subprocess
import vcf
import requests
from collections import namedtuple
import argparse
import numpy as np

VCLASS_COL = 13
VCONSEQUENCE_COL = 6
VINFO_COL = 7


class VcfAnnotate:
    """
    Class to annotate Vcf files

    Attributes
    ----------
    vcf : vcf file to annotate
    sitedict (dict): dictionary that stores variant class
                     and consequence information from VEP output
    output: output file with annotations for each variant

    Methods
    -------
    run_vep_and_parse_output: runs Ensembl VEP and returns dictionary
                              with relevant information about the site
    parse_vcf_and_extract_metrics: writes information about each variant
                                   to the output file specified
    """

    def __init__(self, vcf, output):
        """
        Constructs attributes for VcfAnnotate object

        Args:
             vcf: vcf file to annotate
             output: output file with annotations for each variant
        """
        self.vcf = vcf
        self.sitedict = self.run_vep_and_parse_output()
        self.output = output

    def run_vep_and_parse_output(self):
        """
        Runs Ensembl Variant Effect Predictor (VEP) to extract variant class
        and consequence information for each variant in the vcf

        Returns:
               site_dict (dict): dictionary storing variant class and
                                 consequence for each variant
        """
        info = namedtuple("variant_class", "consequence")
        site_dict = {}
        subprocess.check_call(
            [
                "vep",
                "--cache",
                "-i",
                self.vcf,
                "-o",
                "vep_output.txt",
                "--assembly",
                "GRCh37",
                "--variant_class",
                "--most_severe",
                "--force_overwrite",
            ]
        )
        with open("vep_output.txt") as f:
            for line in f:
                if not line.startswith("#"):
                    info.variant_class = line.split()[VCLASS_COL].split("=")[1]
                    info.consequence = line.split()[VCONSEQUENCE_COL]
                    site_dict[line.split()[0]] = info
        return site_dict

    def parse_vcf_and_extract_metrics(self):
        """
        Parses vcf file and writes information for each variant
        to the output file
        """
        input_vcf = vcf.Reader(open(self.vcf, "r"))
        with open(self.output, "w") as outfile:
            outfile.write(
                "\t".join(
                    [
                        "chrom",
                        "position",
                        "ref",
                        "alt",
                        "seq_depth_at_site",
                        "num_var_alleles",
                        "perc_var_versus_ref_alleles",
                        "variant_class",
                        "consequence",
                        "exac_af",
                    ]
                )
                + "\n"
            )
            for rec in input_vcf:
                metrics = get_info_per_site(rec, self.sitedict)
                outfile.write("\t".join(metrics) + "\n")


def extract_exac_af(site: str):
    """
    Extracts allele frequency of variant from ExAC API.
    For more information: http://exac.hms.harvard.edu/

    Args:
         site (str): information about variant formatted as
                     chrom-position-refallele-alternateallele

    Returns:
         af (float): allele frequency of the variant
    """
    url = "http://exac.hms.harvard.edu/rest/variant/variant/"
    af = "unavailable"
    variant_info = requests.get(f"{url}{site}")
    assert (
        variant_info.status_code == 200
    ), "Please check the format of the argument string"
    json_ = variant_info.json()
    if "allele_freq" in json_.keys():
        af = json_["allele_freq"]
    return af


def get_info_per_site(rec, site_dict: dict):
    """
    Extracts information for each variant in the vcf & VEP output file
    Args:
         rec (vcf.model._Record): record for a variant from parsing vcf
                                  file using vcf.Reader
         site_dict (dict): dictionary storing VEP variant class and
                           consequence information for each variant

    Returns:
         list with information about variant

    """
    chrom = rec.CHROM
    pos = rec.POS
    ref = rec.REF
    alt = rec.ALT
    depth = rec.INFO["DP"]
    num_var = np.sum(rec.INFO["AO"])
    perc_ref = 100 * rec.INFO["RO"] / (num_var + rec.INFO["RO"])
    perc_alt = 100 * num_var / (num_var + rec.INFO["RO"])
    af_list = []
    for aa in alt:
        af_list.append(extract_exac_af(f"{chrom}-{pos}-{ref}-{aa}"))
    sitekey = f"{chrom}_{str(pos)}_{ref}/{'/'.join(map(str,alt))}"
    if sitekey in site_dict.keys():
        consequence = site_dict[sitekey].consequence
        variant_class = site_dict[sitekey].variant_class
    else:
        consequence = "unknown"
        variant_class = "no_VEP_output:" + ",".join(map(str, rec.INFO["TYPE"]))
    return [
        str(chrom),
        str(pos),
        ref,
        ",".join(map(str, alt)),
        str(depth),
        str(num_var),
        f"{str(perc_alt)}:{str(perc_ref)}",
        variant_class,
        consequence,
        ",".join(map(str, af_list)),
    ]


def main(in_vcf, out_metrics):
    """
    Calls parse_vcf_and_extract_metrics for the VcfAnnotate object for
    a specific vcf and specified output file
    """
    VcfAnnotate(in_vcf, out_metrics).parse_vcf_and_extract_metrics()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_vcf", type=str, required=True)
    parser.add_argument("-o", "--output_metrics", type=str, required=True)
    args = parser.parse_args()
    main(args.input_vcf, args.output_metrics)
