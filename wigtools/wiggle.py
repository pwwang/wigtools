"""Classes for wigtools"""
import hashlib
from typing import Tuple
import attr
from diot import OrderedDiot

def _is_meta_line(line):
    """Check if a line is a meta line or a data line"""
    return line[:9] == "fixedStep" or line[:12] == "variableStep"

def _parse_meta_line(line) -> dict:
    """Parse the meta line"""
    ret = {}
    parts = line.rstrip("\r\n").split()
    blocktype = parts.pop(0)
    ret["is_fixed"] = blocktype == "fixedStep"

    for part in parts:
        items = part.split('=')
        if len(items) != 2:
            raise WiggleInvalidMetaLine("Invalid value assignment in meta line")
        # Check the meta name
        ret[items[0]] = (int(items[1])
                         if items[0] in ("start", "step", "span")
                         else items[1])
    return ret

def _compare_regions(region1: Tuple, region2: Tuple) -> int:
    """Compare two regions using the chrom and start
    A region is considered smaller if:
    1. chromsosome is smaller
    2. end position is smaller than the end position of the other region
    No need to consider base, since 0-based coordinate doesn't include end
    and 1-based does"""
    chr1 = _chrom_to_sortable(region1[0])
    chr2 = _chrom_to_sortable(region2[0])
    end2 = region2[2]

    return (0 if (chr1, region1[2]) == (chr2, end2)
            else 1 if (chr1, region1[2]) > (chr2, end2)
            else -1)

def _intersect_regions(region1: Tuple, region2: Tuple,
                       base1: int = None, base2: int = None) -> bool:
    """Tell if two regions are overlapping"""
    chr1 = _chrom_to_sortable(region1[0])
    chr2 = _chrom_to_sortable(region2[0])
    if chr1 != chr2:
        return False
    # convert all to 1-based
    start1 = region1[1] if base1 is None else region1[1] + 1 - base1
    start2 = region2[1] if base2 is None else region2[1] + 1 - base2
    if region1[2] < start2 or region2[2] < start1:
        return False
    return True

def _chrom_to_sortable(chrom):
    """Convert chromosomes to numbers that are sorted like version sort"""
    non_number_chrom_mappings = {
        "X": 23,
        "Y": 24,
        "M": 25,
        "MT": 26
    }
    chrom = chrom[3:] if chrom[:3] == "chr" else chrom
    chrom = non_number_chrom_mappings.get(chrom, chrom)
    if not str(chrom).isdigit():
        chrom = int(hashlib.sha1(chrom.encode()).hexdigest(),
                    16) % (10 ** 8)
    return int(chrom)

class WiggleInvalidDataLine(Exception):
    """When the format of a data line is invalid"""

class WiggleInvalidMetaLine(Exception):
    """When the format of a meta line is invalid"""

class WiggleUnsupportedStringifyFormat(Exception):
    """When the format of stringifying a wiggle block is not supported"""

class WiggleUnsortedFile(Exception):
    """When the input file is not sorted"""

class WiggleReshapeError(Exception):
    """When trying to merge blocks with different spans"""

@attr.s(kw_only=True, slots=True)
class WiggleBlock: # pylint: disable=too-many-instance-attributes
    """A wiggle block that marked by variableStep or fixedStep
    in a wiggle file"""

    # the coordinate base
    base = attr.ib(default=1)
    is_fixed = attr.ib()
    chrom = attr.ib()
    start = attr.ib(default=None)
    step = attr.ib(default=None)
    span = attr.ib(default=1)
    # infer the end of the block to check overlaps
    _end = attr.ib(init=False, repr=False, default=None)

    data = attr.ib(init=False, repr=False, default=attr.Factory(list))
    # the start positions of each regions
    # this will be inferred for fixedStep blocks after intaking is done
    _regions = attr.ib(init=False, repr=False, default=attr.Factory(list))

    @property
    def end(self):
        """Get the end position of block"""
        if self._end:
            return self._end
        return self.regions[-1] + self.span - self.base

    @property
    def regions(self):
        """Get the starts of regions"""
        if self._regions:
            return self._regions
        return [self.start + i * self.step for i in range(len(self.data))]

    @property
    def block_id(self):
        """Get the id of the block"""
        return f"{self.chrom}:{self.start}"

    def take(self, line: str):
        """Take in a data line"""
        parts = line.strip().split()

        if self.is_fixed and len(parts) != 1:
            raise WiggleInvalidDataLine("Wrong columns in data line "
                                        "for a fixedStep block")
        if not self.is_fixed and len(parts) != 2:
            raise WiggleInvalidDataLine("Wrong columns in data line "
                                        "for a variableStep block")

        if self.is_fixed:
            self.data.append(float(parts[0]))
        else:
            if not self.start:
                self.start = int(parts[0])
            self.data.append(float(parts[1]))
            self._regions.append(int(parts[0]))

    def _stringify_to_wiggle(self, base=None, writer=None):
        """Stringify the block to wiggle format"""
        base = self.base if base is None else base
        meta = [
            "fixedStep" if self.is_fixed else "variableStep",
            f"chrom={self.chrom}",
            f"span={self.span}"
        ]
        if self.is_fixed:
            meta.extend([
                f"start={self.start + base - self.base}",
                f"step={self.step}"
            ])
        ret = " ".join(meta) + "\n"

        if writer:
            writer.write(ret)
            # clear buffer
            ret = ""

        for i, dat in enumerate(self.data):
            if self.is_fixed:
                ret += f"{dat}\n"
            else:
                ret += f"{self._regions[i] + base - self.base}\t{dat}\n"
            if writer:
                writer.write(ret)
                ret = ""

        return ret

    def _stringify_to_bedgraph(self, base=None, writer=None):
        """Stringify the block to bedGraph"""
        base = self.base if base is None else base
        ret = ""
        for start, dat in zip(self.regions, self.data):
            start = start + base - self.base
            ret += (f"{self.chrom}\t{start}\t"
                    f"{start+self.span-self.base}\t{dat}\n")
            if writer:
                writer.write(ret)
                ret = ""

        return ret

    def stringify(self, base=None, fmt='wiggle', writer=None):
        """Stringify the block"""
        base = self.base if base is None else base
        if fmt == 'wiggle':
            return self._stringify_to_wiggle(base, writer)
        if fmt.lower() == 'bedgraph':
            return self._stringify_to_bedgraph(base, writer)
        raise WiggleUnsupportedStringifyFormat(
            "Unsupported stringifying format"
        )

    def subset(self, query: Tuple,
               qbase: int = None,
               partial: str = "fraction"):
        """Subset this block by query region"""
        qbase = self.base if qbase is None else qbase
        qstart = query[1] - qbase + self.base
        qend = query[2]
        ret = WiggleBlock(is_fixed=False, chrom=self.chrom,
                          span=self.span, base=self.base)
        for i, region in enumerate(self.regions):
            regend = region + self.span - self.base
            if not _intersect_regions((self.chrom, region, regend),
                                      query, self.base, qbase):
                continue
            intersect_start = max(qstart, region)
            intersect_end = min(qend, regend)
            ret._regions.append(intersect_start)
            if partial == "fraction":
                ret.data.append(self.data[i] *
                                (intersect_end - intersect_start + self.base) /
                                self.span)
            else:
                ret.data.append(self.data[i])
        return ret

    def stats(self, what=None):
        """Calculate stats of this block"""
        what = what or ['min', 'max', 'mean', 'median', 'sum', 'count', 'bp']
        if not isinstance(what, list):
            what = [what]
        ret = {}
        lendata = len(self.data)
        if 'min' in what:
            ret['min'] = min(self.data)
        if 'max' in what:
            ret['max'] = max(self.data)
        if 'mean' in what:
            ret['mean'] = sum(self.data) / lendata
        if 'median' in what:
            if lendata % 2 == 1:
                ret['median'] = self.data[lendata // 2]
            else:
                ret['median'] = (self.data[lendata // 2 - 1] +
                                 self.data[lendata // 2]) / 2.
        if 'sum' in what:
            ret['sum'] = sum(self.data)
        if 'count' in what:
            ret['count'] = lendata
        if 'bp' in what:
            ret['bp'] = lendata * self.span
        return ret

@attr.s(slots=True)
class Wiggle:
    """A wiggle file"""

    wigfile = attr.ib()
    base = attr.ib(default=1)

    blocks = attr.ib(init=False, default=attr.Factory(OrderedDiot),
                     repr=False)

    def __attrs_post_init__(self):
        if self.wigfile:
            # only read the data while a wiggle file is provided
            self._read()

    def __len__(self):
        """Get the number of blocks"""
        return len(self.blocks)

    def _read(self):
        """Read the wiggle file"""
        current_block = None
        with open(self.wigfile, 'r') as fwig:
            for line in fwig:
                if _is_meta_line(line):
                    meta = _parse_meta_line(line)
                    current_block = WiggleBlock(**meta, base=self.base)
                    # start cannot be calculated for variableStep
                elif line.rstrip("\r\n") and current_block:
                    current_block.take(line)
                    if current_block.block_id not in self.blocks:
                        self.blocks[current_block.block_id] = current_block

    def stringify(self, fmt='wiggle', base=None, outfile=None):
        """Stringify the object.
        Only to it for small file, otherwise there may be memory issues"""
        base = self.base if base is None else base

        fout = open(outfile, 'w') if outfile else None
        ret = ""
        for block in self.blocks.values():
            ret += block.stringify(base, fmt, fout)
        if fout:
            fout.close()
        return ret

    def sort(self):
        """Sort the blocks in a wiggle file by chrom and start. """
        block_ids = sorted(
            self.blocks.keys(),
            key=lambda block: (_chrom_to_sortable(self.blocks[block].chrom),
                               self.blocks[block].start)
        )
        orig_blocks = self.blocks
        self.blocks = OrderedDiot()
        for block_id in block_ids:
            self.blocks[block_id] = orig_blocks[block_id]
        del orig_blocks

    def _region(self, block_id):
        """Make a block id a region to compare"""
        return (self.blocks[block_id].chrom,
                self.blocks[block_id].start,
                self.blocks[block_id].end)

    def _intersect(self, qreg, qbase=None,
                   reshape=False, partial="fraction"):
        """Get the blocks that have intersect with qreg. Blocks will be subset
        with the query regions"""
        # pylint: disable=too-many-branches
        qbase = self.base if qbase is None else qbase
        iter_self = iter(self.blocks)
        iter_query = iter(qreg)
        curr_self = curr_query = None
        prev_self = prev_query = None

        ret = Wiggle(None, self.base)
        block = None
        while True:
            if curr_query is None:
                try:
                    curr_query = next(iter_query)
                except StopIteration:
                    break
                if (prev_query and _compare_regions(prev_query,
                                                    curr_query) != -1):
                    raise WiggleUnsortedFile(
                        "Query file is not sorted. "
                        "Region {} appears after {}".format(curr_query,
                                                            prev_query)
                    )
                if reshape:
                    # We are doing reshape, we should generate new blocks
                    # We can't do fixedStep, since we don't know if the
                    # coming blocks are fixedStep or not
                    if block and block.regions and not block.start:
                        block.start = block.regions[0]
                    # have to update the span
                    block = WiggleBlock(base=qbase,
                                        is_fixed=False,
                                        chrom=curr_query[0],
                                        span=None)
                    ret.blocks[f"{curr_query[0]}:{curr_query[1]}"] = block
            if curr_self is None:
                try:
                    curr_self = next(iter_self)
                except StopIteration:
                    break
                if (prev_self and
                        _compare_regions(self._region(prev_self),
                                         self._region(curr_self)) != -1):
                    raise WiggleUnsortedFile(
                        "Current wiggle file is not sorted. "
                        "Region {} appears after {}".format(curr_self,
                                                            prev_self)
                    )

            if _intersect_regions(self._region(curr_self),
                                  curr_query, self.base, qbase):

                if not reshape:
                    ret.blocks[curr_self] = self.blocks[curr_self]
                else:
                    ssblock = self.blocks[curr_self].subset(curr_query,
                                                            qbase,
                                                            partial)
                    if block.span and block.span != self.blocks[curr_self].span:
                        raise WiggleReshapeError(
                            "Cannot merge blocks with different spans "
                            f"({block.span}, {self.blocks[curr_self].span}) "
                            f"that intersect with region: {curr_query}"
                        )
                    block.span = self.blocks[curr_self].span
                    block._regions.extend(ssblock._regions)
                    block.data.extend(ssblock.data)
                    block._end = None
            comp = _compare_regions(self._region(curr_self), curr_query)
            if comp <= 0:
                prev_self = curr_self
                curr_self = None
            else:
                prev_query = curr_query
                curr_query = None

        if reshape and block and block.regions and not block.start:
            block.start = block.regions[0]

        return ret

    def query(self, query, qbase=None):
        """Query the blocks that have intersection with query"""
        return self._intersect(query, qbase)

    def reshape(self, query, qbase=None, partial="fraction"):
        """Reshape the blocks in the query regions"""
        return self._intersect(query, qbase, reshape=True, partial=partial)
