import hashlib
import pytest
from wigtools import wiggle

@pytest.mark.parametrize("line,expected", [
    ("fixedStep", True),
    ("variableStep", True),
    ("nonMeta", False)
])
def test_is_meta_line(line, expected):
    assert wiggle._is_meta_line(line) == expected

@pytest.mark.parametrize("line,meta", [
    ("fixedStep chrom=chr start=1 step=1",
     dict(is_fixed=True, chrom="chr", start=1, step=1)),
    ("variableStep chrom=chr start=1 step=1",
     dict(is_fixed=False, chrom="chr", start=1, step=1)),
])
def test_parse_meta_line(line, meta):
    assert wiggle._parse_meta_line(line) == meta

def test_parse_meta_line_error():

    with pytest.raises(wiggle.WiggleInvalidMetaLine):
        wiggle._parse_meta_line("fixedStep a=b=c")

@pytest.mark.parametrize("reg1,reg2,expected", [
    (("chr1", 1, 2),
     ("1", 1, 2),
      0),
    (("chr1", 2, 2),
     ("1", 0, 2),
     0),
    (("chr1", 1, 2),
     ("1", 1, 3),
     -1),
    (("chr1", 1, 2),
     ("chr11", 1, 3),
     -1),
    (("chr10", 1, 2),
     ("chr9", 1, 3),
     1),
])
def test_compare_regions(reg1, reg2, expected):
    assert wiggle._compare_regions(reg1, reg2) == expected

@pytest.mark.parametrize("reg1, reg2, base1, base2, expected", [
    (("chr1", 1, 10),
     ("1", 1, 3),
     1, 1, True),
    (("chr1", 1, 10),
     ("2", 1, 3),
     1, 1, False),
    (("chr1", 4, 10),
     ("chr1", 1, 3),
     1, 1, False),
    (("chr1", 4, 10),
     ("chr1", 11, 13),
     1, 1, False),
])
def test_intersect_regions(reg1, reg2, base1, base2, expected):
    assert wiggle._intersect_regions(reg1, reg2, base1, base2) == expected

@pytest.mark.parametrize("chrom,expected", [
    ("1", 1),
    ("X", 23),
    ("Y", 24),
    ("M", 25),
    ("MT", 26),
    ("27", 27),
    ("chr99", 99),
    ("abc", int(hashlib.sha1(b"abc").hexdigest(), 16) % (10 ** 8)),
])
def test_chrom_to_sortable(chrom, expected):
    assert wiggle._chrom_to_sortable(chrom) == expected