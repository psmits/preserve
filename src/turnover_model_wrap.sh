#!/bin/bash
# all taxa models
# switch to just the fauna models
FILES=../data/data_dump/sight*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ./turnover_model sample num_samples=1000 num_warmup=1000 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/state_space_${i}.csv &
  done
  wait
done

