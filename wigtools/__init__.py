"""A set of tools for wiggle file"""
from diot import Diot
from pyparam import commands
from wigtools import functional

__version__ = "0.0.1"

commands._desc = __doc__

# switch base
SWITCH_BASE_COMMAND = commands['switch-base']
SWITCH_BASE_COMMAND._desc = "Switch the coordinate base of a wiggle file."
SWITCH_BASE_COMMAND.to.required = True
SWITCH_BASE_COMMAND.to.type = int
SWITCH_BASE_COMMAND.to.desc = ("Either 0 or 1. "
                               "Switch the coordinate base to `<to>`, "
                               "implying the original base `1-<to>`")
SWITCH_BASE_COMMAND.i = "/dev/stdin"
SWITCH_BASE_COMMAND.i.desc = "The input wiggle file"
SWITCH_BASE_COMMAND.o = "/dev/stdout"
SWITCH_BASE_COMMAND.o.desc = "The output wiggle file"

# sort the wiggle file by chrom and start
commands.sort = ("Sort the blocks in a wiggle file by chrom and start. "
                 "Chromosomes will be sorted the way `sort -V` does.")
commands.sort.i = SWITCH_BASE_COMMAND.i
commands.sort.o = SWITCH_BASE_COMMAND.o
commands.sort._hbald = False

commands.stats = ("Statistics for data in a wiggle file for each block")
commands.stats.i = SWITCH_BASE_COMMAND.i
commands.stats.o = SWITCH_BASE_COMMAND.o
commands.stats.base = 1
commands.stats.base.desc = "The coordinate base of the input and output file"
commands.stats.stats = []
commands.stats.stats.desc = ("The data stats for each region. "
                             "Default: ['min', 'max', 'mean', "
                             "'median', 'sum', 'count', 'bp']")
commands.stats.stats.callback = lambda opt: (
    opt.set_value(['min', 'max', 'mean', 'median', 'sum', 'count', 'bp'])
    if not opt.value else None
)
commands.stats.nohead = False
commands.stats.nohead.desc = "Don't put a header for output file."
commands.stats._hbald = False

# reshape: generate a new wiggle file in the query regions,
commands.reshape = ("Generate a new wiggle file and reshape the blocks "
                    "to the query regions")
commands.reshape.i = SWITCH_BASE_COMMAND.i
commands.reshape.o = SWITCH_BASE_COMMAND.o
commands.reshape.base = commands.stats.base
commands.reshape.qfile.required = True
commands.reshape.qfile.desc = "The query file in BED format."
commands.reshape.qbase.desc = "The coordinate base of `qfile`"
commands.reshape.qbase.callback = lambda opt, ps: (
    opt.set_value(ps.base.value) if opt.value is None else None
)
commands.reshape.partial = "fraction"
commands.reshape.partial.desc = [
    "How to assign the data for partially overlapping regions",
    "- `fraction`: proportional to the overlapping length",
    "- `whole`: using the whole data"
]

# query
commands.query = "Find the blocks that intersect with the query regions"
commands.query.i = SWITCH_BASE_COMMAND.i
commands.query.o = SWITCH_BASE_COMMAND.o
commands.query.base = commands.stats.base
commands.query.qfile = commands.query.qfile
commands.query.qbase = commands.query.qbase

# window: make blocks with given window

def switch_base(opts):
    """Switch the coordinate base of a wiggle file"""
    functional.switch_base(opts.i, opts.o, from_base=1-opts.to, to_base=opts.to)

def sort(opts):
    """Sort the blocks in a wiggle file by chrom and start."""
    functional.sort(opts.i, opts.o)

def reshape(opts):
    """Summarize data in a wiggle file for the regions in given region file"""
    functional.reshape(opts.i, opts.o, base=opts.base,
                       qfile=opts.qfile, qbase=opts.qbase,
                       partial=opts.partial)

def stats(opts):
    """Statistics for data in a wiggle file for each block"""
    functional.stats(opts.i, opts.o, opts.base, opts.stats, not opts.nohead)

def query(opts):
    """Find the blocks that intersect with the query regions"""
    functional.query(opts.i, opts.o, opts.base,
                     qfile=opts.qfile, qbase=opts.qbase)

def main():
    """Main entry"""
    command, opts, _ = commands._parse(dict_wrapper=Diot)
    globals()[command.replace('-', '_')](opts)
