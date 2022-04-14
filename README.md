# VariaNNNts

`VariaNNNts` is a command-line interface tool that generates synthetic reads with specified variants (insertions, deletions, and substitutions). 

To run `VariaNNNts`, provide a table of variants and the reference genome to use. 


# Usage

```
usage: VariaNNNts [-h] [-v] [-g] [-o] [--overwrite] [--read_count_bbt]
                  [--seed_bbt] [--q_bbt] [--qv_bbt] [--adderrors_bbt]

  -h, --help            show this help message and exit
  -v , --variant_table
                        CSV table with desired variants. [./variant_table.csv]
  -g , --genome_directory
                        directory containing genomes to generate variants
                        from. [./]
  -o , --output_directory
                        directory to output files. [./VariaNNNts_output]
  --overwrite           flag to overwrite all files in the output_directory.
                        [False]
  --read_count_bbt      number of synthetic read pairs to generate. [500000]
  --seed_bbt            random number generator seed to make synthetic reads.
                        [5]
  --q_bbt               average base quality value. use -1 for a random seed.
                        [40]
  --qv_bbt              standard deviation of variation around q_bbt [0]
  --adderrors_bbt       add substitution errors based on quality values [f]
  ```

 # Citation 

[VariaNNNts - Alex McFarland - github.com/agmcfarland/VariaNNNts](https://github.com/agmcfarland/VariaNNNts)

[BBMap - Bushnell B. - sourceforge.net/projects/bbmap/](https://sourceforge.net/projects/bbmap/)

[Follow me on twitter! @alexmcfarland_](https://twitter.com/alexmcfarland_)