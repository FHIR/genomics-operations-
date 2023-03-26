import pytest
from app import common
import datetime

# checking validate_subject method
def test_validate_subject():
    patient_id = "HG00XXX"

    with pytest.raises(Exception) as excinfo:
        common.validate_subject(patient_id)

    expected_msg = f"Patient ({patient_id}) not found."
    assert str(excinfo.value) == f"400 Bad Request: {expected_msg}"

# checking in_int method
def test_is_int():
    assert common.is_int('2') == True
    assert common.is_int('string') == False
    assert common.is_int('3.2') == False

# checking the date
def test_get_date():
    assert common.get_date("ge2022-01-01") == {'OPERATOR': '$gte','DATE': datetime.datetime(2022,1,1)}
    with pytest.raises(KeyError):
        common.get_date("2022-01-01")

 # checking get_build_and_chrom_by_ref_seq
def test_get_build_and_chrom_by_ref_seq():
    assert common.get_build_and_chrom_by_ref_seq('NC_012920.1') == {"chrom": 'chrM', "build": 'GRCh37'}
    assert common.get_build_and_chrom_by_ref_seq('NC_000014.9') == {"chrom": 'chr14', "build": 'GRCh38'}

#
def test_get_other_build():
    assert common.get_other_build('GRCh38') == "GRCh37"
    assert common.get_other_build('GRCh34') == "GRCh38"

