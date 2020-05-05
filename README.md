# wigtools
A set of tools for wiggle file

## Installation
```
pip install wigtools
```

## Usage
```bash console
> wigtools

Description:
  A set of tools for wiggle file

Usage:
  wigtools <command> [OPTIONS]

Global optional options:
  -h, -H, --help      - Show help message and exit.

Available commands:
  switch-base         - Switch the coordinate base of a wiggle file.
  sort                - Sort the blocks in a wiggle file by chrom and start. Chromosomes will be \
                        sorted the way  sort -V  does.
  stats               - Statistics for data in a wiggle file for each block
  reshape             - Generate a new wiggle file and reshape the blocks to the query regions
  query               - Find the blocks that intersect with the query regions
  help [COMMAND]      - Print help message for the command and exit.
```

### Switch coordinate base for a wiggle file

```bash console
> cat test.wig
variableStep chrom=chr
1	1.0
2	2.0

> cat test.wig | wigtools switch-base --to 0
variableStep chrom=chr span=1
0	1.0
1	2.0
```

### Sort a wiggle file

```bash console
> cat test-unsorted.wig
variableStep chrom=chr
5	1.0
6	2.0
variableStep chrom=chr
1	1.0
2	2.0

> cat test.wig-unsorted.wig | wigtools sort
variableStep chrom=chr span=1
1	1.0
2	2.0
variableStep chrom=chr span=1
5	1.0
6	2.0
```

### Calculate the statistics of each block

```bash console
> cat test-unsorted.wig | wigtools sort | wigtools stats
Chrom   Start   End     min     max     mean    median  sum     count   bp
chr     1	2	1.0     2.0     1.5     1.5     3.0     2	2
chr     5	6	1.0     2.0     1.5     1.5     3.0     2	2

> cat test-unsorted.wig | wigtools sort | wigtools stats --stats mean count --nohead
chr     1	2	1.5     2
chr     5	6	1.5     2
```

### Query a wiggle file to find blocks

```bash console
> cat query.bed
chr	2	3

> wigtools query -i test-unsorted.wig --qfile query.bed
variableStep chrom=chr span=1
1	1.0
2	2.0

> wigtools query -i test-unsorted.wig --qfile query.bed --qbase 0
# No overlapping blocks
```

### Reshape the blocks in query regions

```bash console
> cat reshape.bed
chr	1	8

> cat test-unsorted.wig | wigtools sort | wigtools reshape --qfile reshape.bed
variableStep chrom=chr span=1
1	1.0
2	2.0
5	1.0
6	2.0
```