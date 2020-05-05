import pytest
from wigtools.wiggle import WiggleBlock, WiggleInvalidDataLine, \
    WiggleUnsupportedStringifyFormat

def test_init():
    block = WiggleBlock(is_fixed=True, chrom="chr")
    assert block.base == 1
    assert block.is_fixed == True
    assert block.chrom == "chr"
    assert block.start == None
    assert block.step == None
    assert block.span == 1
    assert block._end == None
    assert block.data == []
    assert block._regions == []

def test_regions():
    block = WiggleBlock(is_fixed=True, chrom="chr")
    block._regions = [1]
    assert block.regions == [1]

    block2 = WiggleBlock(is_fixed=True, chrom="chr", start=1, step=1)
    block2.data = [1.1, 2.1]
    assert block2.regions == [1, 2]

def test_end():
    block = WiggleBlock(is_fixed=True, chrom="chr")
    block._end = 10
    assert block.end == 10

    block2 = WiggleBlock(is_fixed=True, chrom="chr", start=1, step=1)
    block2.data = [1.1, 2.1]
    assert block2.end == 2

def test_block_id():
    block = WiggleBlock(is_fixed=True, chrom="chr", start=1, step=1)
    assert block.block_id == "chr:1"


def test_take():
    block = WiggleBlock(is_fixed=True, chrom="chr", start=1, step=1)
    block.take("1.1\n")
    block.take("2.1\n")

    assert block.data == [1.1, 2.1]
    assert block.regions == [1,2]
    assert block.end == 2

    with pytest.raises(WiggleInvalidDataLine):
        block.take("1\t2\n")

    block2 = WiggleBlock(is_fixed=False, chrom="chr")
    block2.take("10\t1.1\n")
    block2.take("20\t2.1\n")

    assert block2.data == [1.1, 2.1]
    assert block2.regions == [10, 20]
    assert block2.end == 20
    assert block2.start == 10

    with pytest.raises(WiggleInvalidDataLine):
        block2.take("3.1\n")

def test_stringify_to_wiggle(tmp_path):
    block = WiggleBlock(is_fixed=True, chrom="chr", start=1, step=1)
    block.take("1.1\n")
    block.take("2.1\n")
    block.take("3.1\n")
    block.take("4.1\n")
    block.take("5.1\n")

    expected = """\
fixedStep chrom=chr span=1 start=1 step=1
1.1
2.1
3.1
4.1
5.1
"""
    assert block._stringify_to_wiggle() == expected

    outfile = tmp_path / 'test_wiggle_block_stringify_to_wiggle.wig'
    with open(outfile, 'w') as fout:
        assert block._stringify_to_wiggle(writer=fout) == ""

    assert outfile.read_text() == expected

    block2 = WiggleBlock(is_fixed=False, chrom="chr")
    block2.take("1\t1.1\n")
    block2.take("2\t2.1\n")
    block2.take("3\t3.1\n")
    block2.take("4\t4.1\n")
    block2.take("5\t5.1\n")

    expected2 = """\
variableStep chrom=chr span=1
1	1.1
2	2.1
3	3.1
4	4.1
5	5.1
"""
    assert block2._stringify_to_wiggle() == expected2

    outfile2 = tmp_path / 'test_wiggle_block_stringify_to_wiggle2.wig'
    with open(outfile2, 'w') as fout2:
        assert block2.stringify(fmt='wiggle', writer=fout2) == ""

    assert outfile2.read_text() == expected2


def test_stringify_to_bedgraph(tmp_path):
    block = WiggleBlock(is_fixed=True, chrom="chr", start=1, step=1)
    block.take("1.1\n")
    block.take("2.1\n")
    block.take("3.1\n")
    block.take("4.1\n")
    block.take("5.1\n")

    expected = """\
chr	1	1	1.1
chr	2	2	2.1
chr	3	3	3.1
chr	4	4	4.1
chr	5	5	5.1
"""
    assert block._stringify_to_bedgraph() == expected

    outfile = tmp_path / 'test_wiggle_block_stringify_to_bedgraph.wig'
    with open(outfile, 'w') as fout:
        assert block.stringify(fmt='bedgraph', writer=fout) == ""

    assert outfile.read_text() == expected

def test_stringify_error():
    block = WiggleBlock(is_fixed=True, chrom="chr", start=1, step=1)
    block.take("1.1\n")

    with pytest.raises(WiggleUnsupportedStringifyFormat):
        block.stringify(fmt="xyz")


def test_subset():
    block = WiggleBlock(is_fixed=True, chrom="chr",
                        start=1, step=2, span=2)
    block.take("1.1\n") # 1, 2
    block.take("2.1\n") # 3, 4
    block.take("3.1\n") # 5, 6
    block.take("4.1\n") # 7, 8
    block.take("5.1\n") # 9, 10

    ssblock = block.subset(("chr", 2, 4))
    assert ssblock.regions == [2, 3]
    assert ssblock.data == [.55, 2.1]

    ssblock2 = block.subset(("chr", 2, 4), partial="whole")
    assert ssblock2.regions == [2, 3]
    assert ssblock2.data == [1.1, 2.1]

    ssblock3 = block.subset(("chr", 1, 4), qbase=0)
    assert ssblock3.regions == [2, 3]
    assert ssblock3.data == [.55, 2.1]

    ssblock4 = block.subset(("chr", 1, 4), qbase=0, partial="whole")
    assert ssblock4.regions == [2, 3]
    assert ssblock4.data == [1.1, 2.1]

    ssblock5 = block.subset(("chr", 1, 10))
    assert ssblock5.regions == [1, 3, 5, 7, 9]
    assert ssblock5.data == [1.1, 2.1, 3.1, 4.1, 5.1]

def test_stats():
    block = WiggleBlock(is_fixed=True, chrom="chr",
                        start=1, step=2, span=2)
    block.take("1.1\n") # 1, 2
    block.take("2.1\n") # 3, 4
    block.take("3.1\n") # 5, 6
    block.take("4.1\n") # 7, 8
    block.take("5.1\n") # 9, 10

    assert block.stats() == dict(
        min=1.1,
        max=5.1,
        mean=3.1,
        median=3.1,
        sum=15.5,
        count=5,
        bp=10
    )

    block.take("6.1\n")
    assert pytest.approx(block.stats("median")["median"]) == 3.6
