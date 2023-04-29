
Day 1 Summary:
I got it going and I think it works.

Turns out that build without `--release` is much slower. But I can do about 1 iteration every 0.7s compared to regioneR 1.5s

There's a few more options I want to implement
--overlaps : store_true - allow overlapping entries during randomization (maybe make switch to no-overlaps)
--circularize : store_true - use circularize instead of shuffle (just add a random number to all, rotate end to start)
--mask : I have the option, but I'm not using it.
	Now, the easy implementation is to just retry if a random hits
	But I bet I can shorten each chrom during GenomeShift building, and then just mask beds on input by
	if intersects mask, exclude... 
--max-retry : 1000 - when using shuffle, maximum number of attempts before quitting

zeros check: if there's no overlaps observed/permuted, that's an edge case
chrom check: only load -A/-B that have a key in GenomeShift
gzip filereading: Would be nice to not require uncompressed bed files as they can get pretty large

After all of that, I *might* work on profiling. Maybe I can make some more improvements to speed.

local z-score?

And I'm probably not doing per-chrom

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
