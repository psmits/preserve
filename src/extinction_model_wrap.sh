#!/bin/bash
# all taxa models
# switch to just the fauna models
FILES=../data/data_dump/fauna*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ./weibull_review sample num_samples=10000 num_warmup=10000 thin=10 \
      init=0 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/faun_weib_${i}.csv &
  done
  wait
done
