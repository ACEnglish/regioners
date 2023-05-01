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
      --random <RANDOM>        randomization strategy [default: shuffle] [possible values: shuffle, circle]
      --count <COUNT>          overlap counting strategy [default: all] [possible values: all, any]
      --mask <MASK>            bed file of genome regions to mask (chrom<tab>start<tab>end)
      --per-chrom              randomize regions within each chromosome
      --merge-overlaps         merge inputs' overlaps before processing
      --no-swap                do not swap A and B
  -h, --help                   Print help
  -V, --version                Print version
```

## Introduction

`regione_rust` performs a perumutation on the intersection of two bed files. It first counts the number of intersections
between the two bed files and then will randomly shuffle one of the bed files and count the number of intersections.
This permutation is repeated `--num-times`. The mean and standard deviation of the permutations is compared to the
original intersection and a p-value is computed.

### Parameter details
There are a number of options for controlling how `regione_rust` runs. 

#### `--genome`
A two column file with `chrom\tsize`. This becomes the space over which we can shuffle regions. If there are any regions
in the bed files on chromosomes not inside the `--genome` file, those regions will not be loaded.

#### `-A` and `-B`
Bed files with genomic regions to test. They must be sorted and every `start < stop`.

#### `--num-times`
Number of permutations to perform. See [this](https://stats.stackexchange.com/questions/80025/required-number-of-permutations-for-a-permutation-based-p-value) 
for help on selecting a value.

#### `--random`
Randomization strategy. 

##### shuffle
By default, `regione_rust` will randomly shuffle each region. For example, two regions 
`(x, y)`, `(x+s, y+s)` will be shuffled to `(x±r1, y±r1)` and `(x±r2, y±r2)`

##### circle
With `--circle`, all regions are shifted by a set amount such that their spatial distances 
are preserved. i.e. `(x±r1, y±r1), (x±r1, y±r1)`

#### `--count`
Counting strategy.

##### all 
Calculates intersections as the number of overlaps. For example, if one `-A` region hits two `-B` regions, that counts as two intersections. 

##### any
Count that there is any intersection of an interval. So our example above would count a single intersection.

#### `--mask`
The regions may have spans of the genome on which they should not be placed (e.g. reference gaps). Use `--mask`
to hide unused genome spans.

#### `--per-chrom`
By default, `regione_rust` randomization strategies will allow regions to be placed anywhere on the genome. 
With `--per-chrom` regions are placed into positions on the same chromosome.

#### `--merge-overlaps`
Merge overlapping intervals in `-A` and `-B` before processing. This is mainly a convenience function.

#### `--no-swap`
By default, `regione_rust` will swap `-A` and `-B` if `-A` contains fewer intervals.
Because `-A` is randomly shuffled and then intersected with `-B`, if `-A` is smaller than `-B` we can 
speed up runtime by performing the swap. However, one may have a reason to shuffle `-A` regardless of any size
difference. To accomplish this, simply specify `--no-swap`


## Things to Note

Is there any reason to believe your regions aren't randomly distributed across chromosomes? use `--per-chrom`. For
example, occur at different densities per-chromosome (e.g. chr1 has more genes than expected whereas chrY has fewer).

Do your regions have a spacial relationship? Use `--circle`. For example, genes are typically grouped.

Are your regions being merged by `regione_rust`? Yes. as of now. We're working on making this optional as well 
as implementing regioneR's `allow.overlaps=FALSE` during shuffling..

Cites:
- [non-random genes](https://pubmed.ncbi.nlm.nih.gov/20642358/#:~:text=Genes%20are%20nonrandomly%20distributed%20in,genes%20with%20similar%20expression%20profiles.)

## Performance Test

Test of 1,000 permutations on 29,598 epd promoters regions intersection with 1,784,804 TRs using 4 cores.

Docker file running ubuntu:latest
- regione_rust (defaults) : 3.887s
- regione_rust --per-chrom : 3.702s
- regione_rust --per-chrom --circle : 3.107s
- regione_rust --circle : 3.327s


Mac build target x86_64-apple-darwin
- regione_rust (defaults) : 3.210s
- regione_rust --per-chrom : 3.547s
- regione_rust --per-chrom --circle : 3.210s
- regione_rust --circle : 3.413s


regioneR test of 100 permutations on above data in Rstudio docker.
- regioneR : 1292.313s

## Output

The output is a json with structure:
- A_cnt : number of entries in `-A` (note may be swapped from original paramter)
- B_cnt : number of entries in `-A` (note may be swapped from original paramter)
- alt : alternate hypothesis used for p-value - 'l'ess or 'g'reater
- any : value of `--any` parameter
- circle : value -f `--circle` parameter
- merged : value of `merge-overlaps` parameter
- n : number of permutations performed
- obs : observed number of intersections
- per_chrom : value of `--per-chrom` parameter
- perm_mu : permutations' mean
- perm_sd : permutations' standard deviation
- perms : list of permutations' number of intersections
- pval : permutation test's p-value
- swapped : were `-A` and `-B` swapped
- zscore : permutation test's zscore

## Plotting

Using python with seaborn:
```python
data = json.load(open("regione_rust_output.json"))
p = sb.histplot(data=data, x="perms",
		color='gray', edgecolor='gray', kde=False, stat='density')
p = sb.kdeplot(data=data, x="perms",
		color='black', ax=p)

x = data['obs']
props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
plt.axvline(x, color='blue')
plt.text(x, y, 'observed intersections',rotation=90, bbox=props, ma='center')
p.set(xlabel="Intersection Count", ylabel="Permutation Density")
```

<img src="https://raw.githubusercontent.com/ACEnglish/regione_rust/main/figs/example_plot.png" alt="Girl in a jacket" style="width:250px;">


## ToDos:

- gzip file reading
- implement `--no-overlaps` (will require `--max-retry`)
- local z-score
- can save memory by making sending Lappers as read-only to threads?
  [src](https://stackoverflow.com/questions/68908091/how-do-i-send-read-only-data-to-other-threads-without-copying)
  though memory isn't much of a problem.
