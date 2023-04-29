
Day 2 Log:
Turns out that build without `--release` is much slower. But I can do about 1 iteration every 0.7s compared to regioneR 1.5s

[done ] : Data - I need to move everything regioneRust..test_files to here.
[done ] : And document the getters
[done ] : chrom check: only load -A/-B that have a key in GenomeShift
[done ]--any : store_true - count as having any overlap instead of number of overlaps
[done ]--no-swap : do not swap A and B
[done ]--overlaps : store_true - allow overlapping entries during randomization (maybe make switch to no-overlaps)
[stop ]--max-retry : 1000 - when using shuffle, maximum number of attempts before quitting
	Won't get an infinite loop, but will silently give a corrupt result.

[idea ]--mask : I have the option, but I'm not using it.
	Make a new mask object that holds m_lapper = HashMap<String, Lapper<u64, u64>>
	As well as HashMap<String, maskedBaseCount: u64>
	When loading the genome, shorten each chromosome by maskedBaseCount

	When loading bed files entries, 
		exclude entries overlapping the mask
		calculate left_shift = m_lapper[chrom][0 : start].sum()
		set bed entry's m_start = start + m_shift - left_shift
	
[idea ]--per-chrom : .. okay ..
	So, I think there's logic for subsetting the new random position if I hold on to what the master interval tree
	hits.
	But this breaks the eventual --circle

	And I might get performance increases if I just do the work to split by chromosome
	However, I don't know how to translate split-by-chromosome design back to ignoring chromosome.

	So randomizers take a randomize_range object. It will hold the range each hit is allowed to be mapped to.
	For not --per-chrom, it will hold genome_size and just pick a number
	For --per-chrom, it will hold chromosome interals,
		for each input interval (which are unchanged)
			map to randomize_range to figure out my bounds.
	For --circle and not --per-chrom
		it will hold genome_size, pick a number to shift everything
		for each interval - move that shift while 'wrapping' anything that moves off genome
	For --circle and --per-chrom
		for each randomize_range interval, pick a number within its range.
		for each input interval (which are unchanged)
			map to the random number's made above, 'wrapping' anything that moves off chromosome to the
			interval's start.

	I'll also probably need this.
	https://stackoverflow.com/questions/36390665/how-do-you-pass-a-rust-function-as-a-parameter

[later ] :
zeros check: if there's no overlaps observed/permuted, that's an edge case for the math
gzip filereading: Would be nice to not require uncompressed bed files as they can get pretty large
no_clone : https://stackoverflow.com/questions/68908091/how-do-i-send-read-only-data-to-other-threads-without-copying
tests : Turn CpGi and prom into toy data

After all of that, I *might* work on profiling. Maybe I can make some more improvements to speed.
later : local z-score?
later : I would love to figure out how to dedup the `io::read_*` .. maybe same as --circle

Day 1 Log:
CLI only that takes four parameters
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
