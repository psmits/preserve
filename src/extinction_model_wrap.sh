#!/bin/bash
# all taxa models
FILES=../data/data_dump/survival*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ./survival_model sample num_samples=10000 num_warmup=10000 thin=10 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/surv_${i}.csv &
  done
  wait
done
# switch to just the fauna models
FILES=../data/data_dump/fauna*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ./fauna_surv_model sample num_samples=10000 num_warmup=10000 thin=10 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/faun_surv_${i}.csv &
  done
  wait
done
