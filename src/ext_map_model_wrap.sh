#!/bin/bash
# all taxa models
# switch to just the fauna models
FILES=../data/data_dump/fauna*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ./survival_exp_map sample num_samples=7500 num_warmup=7500 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/faun_expo_map_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    ./survival_wei_map sample num_samples=7500 num_warmup=7500 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/faun_weib_map_${i}.csv &
  done
  wait
done
