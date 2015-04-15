#!/bin/bash
# all taxa models
# switch to just the fauna models
FILES=../data/data_dump/fauna*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ./survival_exponential sample num_samples=5000 num_warmup=5000 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/faun_expo_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    ./survival_weibull sample num_samples=5000 num_warmup=5000 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/faun_weib_${i}.csv &
  done
  wait
done
