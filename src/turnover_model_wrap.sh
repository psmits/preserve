#!/bin/bash
# all taxa models
# switch to just the fauna models
FILES=../data/data_dump/sight_info*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ./turnover_model sample num_samples=5000 num_warmup=5000 \
      random seed=420 \
      id=$i \
      data file=$f \
      init=1 \
      output file=../data/mcmc_out/state_space_${i}.csv &
  done
  wait
done

