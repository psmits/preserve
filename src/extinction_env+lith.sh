#!/bin/bash
# all taxa models
# switch to just the fauna models
FILES=../data/data_dump/fauna*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ./survival_wei_lith sample num_samples=50000 num_warmup=50000 thin=25 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/faun_expo_lith_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    ./survival_wei_lith sample num_samples=50000 num_warmup=50000 thin=25 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/faun_weib_lith_${i}.csv &
  done
  wait
done
