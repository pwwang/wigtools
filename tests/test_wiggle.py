import sys
import cmdy
import pytest

@pytest.fixture
def python():
    return cmdy.python.bake(_exe=sys.executable)

def test_switch_base(python):
    cmd = cmdy.echo("""\
variableStep chrom=chr span=1
1\t1
2\t2
variableStep chrom=chr span=1
5\t5
6\t6
""", _pipe=True) | python({"m": "wigtools"}, "switch-base", to=0)
    assert cmd.stdout == """\
variableStep chrom=chr span=1
0\t1.0
1\t2.0
variableStep chrom=chr span=1
4\t5.0
5\t6.0
"""


def test_sort(python):
    cmd = cmdy.echo("""\
variableStep chrom=chr span=1
5\t5
6\t6
variableStep chrom=chr span=1
1\t1
2\t2
""", _pipe=True) | python({"m": "wigtools"}, "sort")
    assert cmd.stdout == """\
variableStep chrom=chr span=1
1\t1.0
2\t2.0
variableStep chrom=chr span=1
5\t5.0
6\t6.0
"""

def test_stats(python):
    cmd = cmdy.echo("""\
variableStep chrom=chr span=1
5\t5
6\t6
variableStep chrom=chr span=1
1\t1
2\t2
""", _pipe=True) | python({"m": "wigtools"}, "stats")
    assert cmd.stdout == """\
Chrom\tStart\tEnd\tmin\tmax\tmean\tmedian\tsum\tcount\tbp
chr\t5\t6\t5.0\t6.0\t5.5\t5.5\t11.0\t2\t2
chr\t1\t2\t1.0\t2.0\t1.5\t1.5\t3.0\t2\t2
"""

def test_query(python, tmp_path):
    qfile = tmp_path / 'test_query.bed'
    qfile.write_text("""\
chr\t1\t2
""")
    cmd = cmdy.echo("""\
variableStep chrom=chr span=1
1\t1
2\t2
variableStep chrom=chr span=1
5\t5
6\t6
""", _pipe=True) | python({"m": "wigtools"}, "query", qfile=qfile)
    assert cmd.stdout == """\
variableStep chrom=chr span=1
1\t1.0
2\t2.0
"""

def test_reshape(python, tmp_path):
    qfile = tmp_path / 'test_query.bed'
    qfile.write_text("""\
# This is a bed file
chr\t1\t10
""")
    cmd = cmdy.echo("""\
variableStep chrom=chr span=1
1\t1
2\t2
variableStep chrom=chr span=1
5\t5
6\t6
""", _pipe=True) | python({"m": "wigtools"}, "reshape", qfile=qfile)
    assert cmd.stdout == """\
variableStep chrom=chr span=1
1\t1.0
2\t2.0
5\t5.0
6\t6.0
"""
