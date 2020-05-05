"""Implementation of functions for the tools"""
from typing import List, Iterable
from wigtools.wiggle import Wiggle

def _bed_to_regions(bedfile: str) -> Iterable:
    with open(bedfile) as fbed:
        for line in fbed:
            line = line.rstrip("\r\n")
            if line[:1] == '#':
                continue
            parts = line.split("\t")[:3]
            yield parts[0], int(parts[1]), int(parts[2])


def switch_base(infile: str, outfile: str, from_base: int, to_base: int):
    """Switch the coordinate base of a wiggle file"""
    wiggle = Wiggle(infile, base=from_base)
    wiggle.stringify(base=to_base, outfile=outfile)

def sort(infile: str, outfile: str):
    """Sort the blocks in a wiggle file by chrom and start. """
    wiggle = Wiggle(infile)
    wiggle.sort()
    wiggle.stringify(outfile=outfile)

def stats(infile: str, outfile: str, base: int,
          statistics: List[str], header: bool):
    "Statistics for data in a wiggle file for each block"
    wiggle = Wiggle(infile, base)

    with open(outfile, 'w') as fout:
        if header:
            fout.write("Chrom\tStart\tEnd\t{}\n".format('\t'.join(statistics)))
        for block in wiggle.blocks.values():
            bstats = block.stats(statistics)
            stats_str = "\t".join(str(bstats[stat]) for stat in statistics)
            fout.write(f"{block.chrom}\t{block.start}\t"
                       f"{block.end}\t{stats_str}\n")

def reshape(infile: str, # pylint: disable=too-many-arguments
            outfile: str,
            base: int,
            qfile: str,
            qbase: int,
            partial: str):
    """Summarize data in a wiggle file for the regions in given region file"""
    wiggle = Wiggle(infile, base)
    regions = _bed_to_regions(qfile)

    wiggle = wiggle.reshape(regions, qbase, partial)

    wiggle.stringify(outfile=outfile)


def query(infile: str,
          outfile: str,
          base: int,
          qfile: str,
          qbase: int):
    """Summarize data in a wiggle file for the regions in given region file"""
    wiggle = Wiggle(infile, base)
    regions = _bed_to_regions(qfile)

    wiggle = wiggle.query(regions, qbase)

    wiggle.stringify(outfile=outfile)
