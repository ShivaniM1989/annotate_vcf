import annotate_vcf
import pytest


def test_extract_exac_af_pass():
    """
    Checks if correct allele frequency
    reported from ExAC API
    """
    site = "14-21853913-T-C"
    af = annotate_vcf.extract_exac_af(site)
    assert af == 0.000046048996131884326


def test_extract_exac_af_fail():
    """
    Checks for assertion error when improperly
    formatted site provided
    """
    site = "14_21853913_T_C"
    with pytest.raises(AssertionError):
        annotate_vcf.extract_exac_af(site)
