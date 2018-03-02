#!/bin/bash
# all taxa models
# switch to just the fauna models
FILES=../data/data_dump/impute*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ../stan/survival_geometric \
      sample num_samples=1000 num_warmup=1000 thin=1 \
      id=$i \
      init=0 \
      data file=$f \
      output file=../data/mcmc_out/faun_geometric_${i}.csv &
  done
  wait
done
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ../stan/survival_disc_weibull \
      sample num_samples=1000 num_warmup=1000 thin=5 \
      id=$i \
      init=0 \
      data file=$f \
      output file=../data/mcmc_out/faun_dweibull_${i}.csv &
  done
  wait
done
