#!/bin/bash
# all taxa models
# switch to just the fauna models
FILES=../data/data_dump/impute*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ../stan/surv_cweib_base sample \
      adapt delta=0.999 \
      num_samples=1000 num_warmup=1000 thin=1 \
      algorithm=hmc engine=nuts max_depth=10 stepsize=1 \
      id=$i \
      init=0 \
      data file=$f \
      output file=../data/mcmc_out/surv_cweib_base_${i}.csv &
  done
  wait
done
