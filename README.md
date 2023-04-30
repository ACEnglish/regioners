# regione_rust
A rust implementation of [regioneR](https://academic.oup.com/bioinformatics/article/32/2/289/1744157) 
for interval overlap permutation testing.


## Install

```bash
git clone https://github.com/ACEnglish/regione_rust
cd regione_rust
cargo build --release
# executable in ./target/release/regione_rust
```

## Quick Start

### Full Params
```bash
Usage: regione_rust [OPTIONS] --genome <GENOME> -A <BED_A> -B <BED_B> --output <OUTPUT>

Options:
  -g, --genome <GENOME>        chromosome lengths (chrom<tab>length)
  -A <BED_A>                   bed file of regions (chrom<tab>start<tab>end)
  -B <BED_B>                   bed file of regions (chrom<tab>start<tab>end)
  -n, --num-times <NUM_TIMES>  number of permutations to perform [default: 100]
  -o, --output <OUTPUT>        output json file
  -t, --threads <THREADS>      number of threads to use [default: 1]
      --mask <MASK>            bed file of genome regions to mask (chrom<tab>start<tab>end)
      --no-swap                do not swap A and B
      --any                    count any overlap instead of number of overlaps
      --per-chrom              randomize regions within each chromosome
      --circle                 use circularization for randomization
  -h, --help                   Print help
  -V, --version                Print version
```

## Introduction

`regione_rust` performs a perumutation on the intersection of two bed files. It first counts the number of intersections
between the two bed files and then will randomly shuffle one of the bed files and count the number of intersections.
This permutation is repeated `--num-times`. The mean and standard deviation of the permutations is compared to the
original intersection and a p-value is computed.

### Parameters
There are a number of options for controlling how `regione_rust` runs. 

#### `--genome`
A two column file with `chrom\tsize`. This becomes the space over which we can shuffle regions. If there are any regions
in the bed files on chromosomes not inside the `--genome` file, those regions will not be loaded.

#### `-A` and `-B`
Bed files with genomic regions to test.

#### `--threads`
Permutations can be calculated by multiple `--threads` in order to speed up computation.

#### `--mask`
The regions may have spans of the genome on which they should not be placed (e.g. reference gaps). Use `--mask`
to hide unused genome spans.

#### `--no-swap`
By default, `regione_rust` will swap `-A` and `-B` if `-A` contains fewer intervals.
Because `-A` is randomly shuffled and then intersected with `-B`, if `-A` is smaller than `-B` we can 
speed up runtime by performing the swap. However, one may have a reason to shuffle `-A` regardless of size. 
To accomplish this, simply specify `--no-swap`

#### `--any`
By default, `regione_rust` calculates intersections as the number of overlaps. For example, if one `-A` region hits two
`-B` regions, that counts as two intersections. By using `--any`, we count that there is any intersection between the
two bed files. So for our example, `--any` would count it as one intersection.

#### `--per-chrom`
By default, `regione_rust` will allow regions to be placed anywhere on the genome. With `--per-chrom` regions are
shuffled into positions on the same chromosome.

#### `--circle`
By default, `regione_rust` will randomly shuffle each region. With `--circle`, all regions are shifted by a set amount
such that their spatial distances are preserved. For example, two regions `(x, y)`, `(x+s, y+s)` by default will be
shuffled to `(x±r1, y±r1)` and `(x±r2, y±r2)`. However, `--circle` shuffles to `(x±r1, y±r1), (x±r1, y±r1)`.

## Things to Note

Is there any reason to believe your regions aren't randomly distributed across chromosomes? use `--per-chrom`. For
example, occur at different densities per-chromosome (e.g. chr1 has more genes than expected whereas chrY has fewer).

Do your regions have a spacial relationship? Use `--circle`. For example, genes are typically grouped.

Are your regions being merged by `regione_rust`? Yes. as of now, we're going for speed and `rust_lapper` works more
quickly with merged overlaps by seek-ing. We're working on making this optional as well as implementing regioneR's 
`allow.overlaps=FALSE` during shuffling..

Cites:
- [non-random genes](https://pubmed.ncbi.nlm.nih.gov/20642358/#:~:text=Genes%20are%20nonrandomly%20distributed%20in,genes%20with%20similar%20expression%20profiles.)

## Performance Test

Days are pseudo-versioning based on end of development day testing.

Day2 test of 1,000 permutations on 29,598 (epd promoters) regions intersection with 1,784,804 (TRs) on 4 cores.
- regione_rust : Im guessing 133.417s
- regione_rust --per-chrom : 121.620s 2m11s
- regione_rust --per-chrom --circle : 1m2.301s
- regione_rust --circle : 34.439s

Day2 test of 1,000 permutations on 29,598 (epd promoters) regions intersection with 1,784,804 (TRs) on 4 cores.
- regione_rust : 169.400s
- regione_rust --per-chrom : 158.159s
- regione_rust --per-chrom --circle : 55.629s
- regione_rust --circle : 44.529s

Day1 test of 100 permutations on 29,598 regions intersection with 1,784,804 on 4 cores.
- regioneR : 1292.313s
- regione_rust (no params) : 641.454s

## Output

The output is a json with structure:
- alt : (l)ess or (g)reater alternate hypothesis used
- n : number of permutations performed
- obs : observed number of intersections
- perm_mu : permutations' mean
- perm_sd : permutations' standard deviation
- perms : permutations' number of intersections
- pval : permutation test's p-value
- zscore : permutation test's zscore

## Plotting

Using python with seaborn:
```python
data = json.load(open("regione_rust_output.json"))
# 
p = sb.histplot(data=data, x="perms",
		color='gray', edgecolor='gray', kde=False, stat='density')
p = sb.kdeplot(data=data, x="perms",
		color='black', ax=p)

x = data['obs']
props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
plt.axvline(x, color='blue')
plt.text(x, y, 'observed intersections',rotation=90, bbox=props, ma='center')
p.set(xlabel="Intersection Count", ylabel="Permutation Density")

p = sb.histplot(data=data['perms'])
x = data['obs']
y = p.get_ylim()[1] // 2 - 20
plt.axvline(x)
plt.text(x - 150, y, 'observed intersections',rotation=90)
```

<img src="img_girl.jpg" alt="Girl in a jacket" style="width:250px;">


## ToDos:

- implement `--no-overlaps` (will require `--max-retry`)
- can save memory by making sending Lappers as read-only to threads? [src](https://stackoverflow.com/questions/68908091/how-do-i-send-read-only-data-to-other-threads-without-copying)
- gzip file reading
- local z-score
