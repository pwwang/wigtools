from pathlib import Path
import pytest
from remotedata import remotedata
from wigtools.wiggle import Wiggle, WiggleUnsortedFile, WiggleReshapeError

@pytest.fixture
def here():
    return Path(__file__).parent.resolve()

@pytest.fixture
def rdata(here):
    return remotedata(dict(
        source = 'github',
        cachedir = str(here / 'remotedata'),
        ## if branch is not master: pwwang/remotedata/branch
        repos  = 'CRG-Barcelona/bwtool',
        ## optional, default is first part of repos
        # user = 'pwwang',
        ## github token, in case you have > 60 requests per hours to github API
        # token = 'xxx',
    ))

def test_init(tmp_path):
    wigfile = tmp_path / 'test_wiggle_wiggle_init.wig'
    wigfile.write_text("")
    wiggle = Wiggle(wigfile)
    assert wiggle.wigfile == wigfile
    assert wiggle.base == 1
    assert len(wiggle.blocks) == 0
    assert len(wiggle) == 0

def test_read(rdata):
    wigfile = rdata.get("tests/wigs/main.wig")
    wiggle = Wiggle(wigfile)
    assert len(wiggle) == 2

    assert wiggle.blocks['chr:1'].regions == list(i+1 for i in range(23))
    assert wiggle.blocks['chr:28'].regions == list(i+28 for i in range(9))

def test_stringify(rdata, tmp_path):
    wigfile = rdata.get("tests/wigs/main.wig")
    wiggle = Wiggle(wigfile)
    expected = """\
variableStep chrom=chr span=1
0	1.0
1	2.0
2	5.0
3	6.0
4	5.0
5	3.0
6	3.0
7	5.0
8	5.0
9	5.0
10	6.0
11	6.0
12	0.0
13	2.0
14	3.0
15	3.0
16	10.0
17	4.0
18	4.0
19	2.0
20	2.0
21	2.0
22	1.0
variableStep chrom=chr span=1
27	2.0
28	3.0
29	4.0
30	6.0
31	6.0
32	4.0
33	4.0
34	4.0
35	2.0
"""

    assert wiggle.stringify(base=0) == expected

    outfile = tmp_path / 'test_wiggle_wiggle_stringify.wig'
    assert wiggle.stringify(base=0, outfile=outfile) == ""
    assert outfile.read_text() == expected

def test_sort(rdata):
    wiggle = Wiggle(rdata.get("tests/wigs/main.wig"))
    assert list(wiggle.blocks) == ["chr:1", "chr:28"]
    # push the first block to the end
    wiggle.blocks["chr:1"] = wiggle.blocks.pop("chr:1")
    assert list(wiggle.blocks) == ["chr:28", "chr:1"]
    wiggle.sort()
    assert list(wiggle.blocks) == ["chr:1", "chr:28"]

def test_query(rdata):
    wiggle = Wiggle(rdata.get("tests/wigs/main.wig"))
    wiggle2 = wiggle.query([
        ("chr", 1, 4),
    ])
    assert len(wiggle2) == 1
    assert list(wiggle2.blocks) == ['chr:1']

    assert wiggle2.blocks['chr:1'].start == 1
    assert wiggle2.blocks['chr:1'].end == 23

def test_intersect_unsorted(rdata):

    wiggle = Wiggle(rdata.get("tests/wigs/main.wig"))
    # push the first block to the end
    wiggle.blocks["chr:1"] = wiggle.blocks.pop("chr:1")
    assert list(wiggle.blocks) == ["chr:28", "chr:1"]

    with pytest.raises(WiggleUnsortedFile):
        wiggle._intersect([
            ("chr", 130, 134),
        ])

    wiggle = Wiggle(rdata.get("tests/wigs/main.wig"))

    with pytest.raises(WiggleUnsortedFile):
        wiggle._intersect([
            ("chr", 3, 4),
            ("chr", 1, 2),
        ])

def test_reshape(rdata):
    wiggle = Wiggle(rdata.get("tests/wigs/main.wig"))
    # 1-23
    # 28-36
    reshaped = wiggle.reshape([
        ("chr", 2, 5),
        ("chr", 10, 12),
        ("chr", 35, 36)
    ])

    assert len(reshaped) == 3
    assert list(reshaped.blocks) == ["chr:2", "chr:10", "chr:35"]

    assert reshaped.blocks["chr:2"].start == 2
    assert reshaped.blocks["chr:2"].end == 5
    assert reshaped.blocks["chr:2"].regions == [2,3,4,5]
    assert reshaped.blocks["chr:2"].data == [2.0,5.0,6.0,5.0]

    assert reshaped.blocks["chr:10"].start == 10
    assert reshaped.blocks["chr:10"].end == 12
    assert reshaped.blocks["chr:10"].regions == [10,11,12]
    assert reshaped.blocks["chr:10"].data == [5.0,6.0,6.0]

    assert reshaped.blocks["chr:35"].start == 35
    assert reshaped.blocks["chr:35"].end == 36
    assert reshaped.blocks["chr:35"].regions == [35,36]
    assert reshaped.blocks["chr:35"].data == [4.0,2.0]

def test_reshape_partially_overlapping_block(rdata):
    wiggle = Wiggle(rdata.get("tests/wigs/main.wig"))
    # 1-23
    # 28-36
    # partially overlapping
    reshaped = wiggle.reshape([
        ("chr", 20, 25),
        ("chr", 26, 30)
    ])
    assert len(reshaped) == 2
    assert list(reshaped.blocks) == ["chr:20", "chr:26"]

    assert reshaped.blocks["chr:20"].start == 20
    assert reshaped.blocks["chr:20"].end == 23
    assert reshaped.blocks["chr:20"].regions == [20,21,22,23]
    assert reshaped.blocks["chr:20"].data == [2.0,2.0,2.0,1.0]

    assert reshaped.blocks["chr:26"].start == 28
    assert reshaped.blocks["chr:26"].end == 30
    assert reshaped.blocks["chr:26"].regions == [28,29,30]
    assert reshaped.blocks["chr:26"].data == [2.0,3.0,4.0]

def test_reshape_partially_overlapping_region(tmp_path):

    wiggle_str = """\
variableStep chrom=chr span=100
100	100.0
300 200.0
600	300.0
variableStep chrom=chr span=100
800	400.0
1200	500.0
1600	600.0
2000	700.0
fixedStep chrom=chr start=3000 step=100 span=200
1.0
2.0
"""
    wiggle_file = tmp_path / "test_wiggle_wiggle_reshape_partially_overlapping_region.wig"
    wiggle_file.write_text(wiggle_str)
    wiggle = Wiggle(wiggle_file)

    reshaped = wiggle.reshape([
        ("chr", 200, 1000),
        ("chr", 1650, 2049)
    ])

    assert len(reshaped) == 2
    assert list(reshaped.blocks) == ["chr:200", "chr:1650"]

    assert reshaped.blocks["chr:200"].start == 300
    assert reshaped.blocks["chr:200"].end == 899
    assert reshaped.blocks["chr:200"].regions == [300, 600, 800]
    assert reshaped.blocks["chr:200"].data == [200.0, 300.0, 400.0]

    assert reshaped.blocks["chr:1650"].start == 1650
    assert reshaped.blocks["chr:1650"].end == 2099
    assert reshaped.blocks["chr:1650"].regions == [1650, 2000]
    assert reshaped.blocks["chr:1650"].data == [300.0, 350.0]

    with pytest.raises(WiggleReshapeError):
        wiggle.reshape([
            ("chr", 1000, 4000)
        ])
