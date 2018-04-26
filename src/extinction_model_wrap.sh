#!/bin/bash
FILES=../data/data_dump/impute*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ../stan/surv_cweib_base sample \
      adapt delta=0.9 \
      num_samples=2000 num_warmup=2000 thin=2 \
      algorithm=hmc engine=nuts max_depth=10 stepsize=1 \
      id=$i \
      init=0 \
      data file=$f \
      output file=../data/mcmc_out/surv_cweib_base_${i}.csv &
  done
  wait
done
FILES=../data/data_dump/impute*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ../stan/surv_cweib_cens sample \
      adapt delta=0.9 \
      num_samples=2000 num_warmup=2000 thin=2 \
      algorithm=hmc engine=nuts max_depth=10 stepsize=1 \
      id=$i \
      init=0 \
      data file=$f \
      output file=../data/mcmc_out/surv_cweib_cens_${i}.csv &
  done
  wait
done
