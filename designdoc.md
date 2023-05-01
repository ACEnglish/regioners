Day 4 Log:
==========

CLI cleaning:

The parameters are getting a little wonky. So I'm going to clean them up.

--random shuffle, circle, novl
--count all, any

no-overlap randomizer

build_gaps(lapper, genome)

no_overlap_shuffle is its own thing

So I gotta make it known that there's 
	<default> shuffle with overlaps
	--circle preserves spacing and moves
	--no-overlap shuffle but without allowing overlaps

Day 3 Log:
==========

Didn't keep tight notes.
Main thing was that I realized how I was performing the was incorrect (needed to use std::mem::swap instead of `let
(a,b) = (b.clone(), a.clone())`)

Speaking of weird syntax- turns out I had a bunch of it. But I found `rustup component add clippy` and `rustfmt`. The
latter seemed to help with indenting. The former was able to tell I was using poor design patterns. There were two main
errors:

I was manually implementing map with:
```rust
let mask = match args.mask { Some(p) => Some(io::read_mask(&p)), None => None }
```
But I should have been doing
```rust
let mask = args.mask.map(|p| io::read_mask(&p));
```
I don't know how long it would have taken me to realize what `.map` was doing without clippy calling it out.

Similarly, I found that my lines:
```rust
if let Ok(lines) = read_lines(file) {
    for line in lines {
        if let Ok(cur_line) = line {
```
Should have been written
```rust
if let Ok(lines) = read_lines(file) {
    for line in lines.flatten() {
```
That third `if let` was the only ... thing ... inside the loop. (things like `Some` and `None`).
So because I was only handling `Ok`, I could just use `flatten`. I'll have to look up if flatten is only for `Ok` and if
its only on iterables.

This helped me figure out how to use the match/map/collect a little better and lead to the refactoring of the
overlap counters in main (which I think might be a smidge faster)

I started making `--no-overlaps`. There's commented out code that will run it. However, it's brute force and will fail
often. I have an idea for considering the gaps between regions to help shuffle the intervals. But I just realized it'll
only work if we enforce `--merge-overlaps`.

Day 2 Log:
=========
Turns out that build without `--release` is much slower. But I can do about 1 iteration every 0.7s compared to regioneR 1.5s

[done ] : Data - I need to move everything regioneRust..test_files to here.
[done ] : And document the getters
[done ] : chrom check: only load -A/-B that have a key in GenomeShift
[done ]--any : store_true - count as having any overlap instead of number of overlaps
[done ]--no-swap : do not swap A and B
[done ]--overlaps : store_true - allow overlapping entries during randomization (maybe make switch to no-overlaps)
[done ]--mask : easy enough after changing the GenomeShift
[done ]--per-chrom : easy enough after changing the GenomeShift
[done ]--circle and I'm pretty much done... 
[done ] zeros checkd for zsocre
[most ] write documentation
[half ]--max-retry : 1000 - when using shuffle, maximum number of attempts before quitting
	Currently won't get an infinite loop, but will silently give a corrupt result.

[later ] :
gzip filereading: Would be nice to not require uncompressed bed files as they can get pretty large
no_clone : 
tests : Turn CpGi and prom into toy data

After all of that, I *might* work on profiling. Maybe I can make some more improvements to speed.
later : local z-score?
later : I would love to figure out how to dedup the `io::read_*` .. maybe same as --circle

https://stackoverflow.com/questions/36390665/how-do-you-pass-a-rust-function-as-a-parameter

Day 1 Log:
=====
```
-A               bed file of regions
-B               bed file of regions
[-g|--genome]    genome length definition.. chrom\tlength
[-m|--mask]      bed file of genome regions to mask
[-n|--num-times] number of permutations to perform
[-o|--output]    output - json file to save each permutation's intersect counts as well as stats: {p-val, etc}
```

Build into the ArgParser a 'validate' that will ensure all the files exist as a pre-flight check
including a check on num-times being high enough to calculate a meaningful p-value


Internals
=====

File parsing :
1 - open A/B/genome/mask (read/print a few lines to debug, type casing will be hard)
Wasn't too bad .. but now I want to move it to a submodule.. and also done, cool


rust-lapper
========
for the first version, I don't want to deal with a HashMap of chromosomes.

So, step 1 is to parse the --genome and build a `HashMap<str, u32>`
let mut cur_start:u32 = 0
for each line in genome:
	hashmap[line[chrom]] = cur_start
	cur_start += line[chrom_length]

Then, when reading the bed file, we push the boundaries up by the hashmap[chrom] shift
That way we only have one giant Vec<Iv> to Lapper against.

build a per-chromosome (from genome) rust-lapper
populate lines from files into the rust-lapper

Now, do I have to randomly place across (or within) chromosomes? Would it suffice to turn
the genome into a single super long sequence, translate 
intervals into a single 'super long' sequence

randomizer
----------
I will probably have to build the region randomizer

I will need to figure out how to use the mask. Like, does it exclude regions which -A/-B overlap?
I'm pretty sure it will do the exclude during the randomization i.e. if random place is in mask, retry.

And then I'll just have to copy the 


Advanced features
--
--force-A-random : first version will just randomize A, second version will randomize whichever has fewer regions A or B but
with a flag to force A to be randomized

Commands for the different randomizations / overlaps


```

      #Compute the p-value
      if (alt == "less") {
        pval <- (sum(orig.ev >= rand.ev, na.rm=TRUE) + 1) / (num.valid.values + 1)
      } else { #alt == "greater"
        pval <- (sum(orig.ev <= rand.ev, na.rm=TRUE) + 1) / (num.valid.values + 1)
      }
      #if the original alternative was not the best one, suggest the user to change it
      if(alternative=="greater" & orig.ev<mean(rand.ev,na.rm=TRUE)) message("Alternative is greater and the observed statistic is less than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
      if(alternative=="less" & orig.ev>mean(rand.ev,na.rm=TRUE)) message("Alternative is less and the observed statistic is greater than the permuted statistic mean. Maybe you want to use recomputePermTest to change the alternative hypothesis.")
      #Compute the z-score
      if(orig.ev == 0 & all(rand.ev == 0)){ #If everything is 0, warning and "empty" results
        warning(paste0("All permuted values and the original evaluation value are equal to 0. Z-score cannot be computed."))
        pval <- 1
        zscore <- NA
      } else{
        zscore <- round((orig.ev - mean(rand.ev, na.rm=TRUE)) / stats::sd(rand.ev, na.rm=TRUE), 4)
      }
    } else {
      pval <- NA
      zscore <- NA
      alt <- alternative
    }
```

chunk_size = num_elements / threads
for i in range(0..threads)
    for i in range(0..chunk_size)
        span a thread and do the work
