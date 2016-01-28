#!/bin/bash
# all taxa models
# switch to just the fauna models
FILES=/home/psmits/preserve/data/data_dump/fauna*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    /home/psmits/perserve/stan/survival_full sample num_samples=10000 num_warmup=10000 thin=10 \
      init=0 \
      id=$i \
      data file=$f \
      output file=/home/psmits/preserve/data/mcmc_out/faun_weib_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    /home/psmits/preserve/stan/survival_nosample sample num_samples=10000 num_warmup=10000 thin=10 \
      init=0 \
      id=$i \
      data file=$f \
      output file=/home/psmits/perserve/data/mcmc_out/faun_nosamp_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    peter/preserve/stan/survival_constantalpha sample num_samples=10000 num_warmup=10000 thin=10 \
      init=0 \
      id=$i \
      data file=$f \
      output file=peter/preserve/data/mcmc_out/faun_conalph_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    peter/preserve/stan/survival_conalp_nosamp sample num_samples=10000 num_warmup=10000 thin=10 \
      init=0 \
      id=$i \
      data file=$f \
      output file=peter/preserve/data/mcmc_out/faun_conalp_nosamp_${i}.csv &
  done
  wait
done
