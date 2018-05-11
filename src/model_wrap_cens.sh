#!/bin/bash
FILES=../data/data_dump/impute*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    ../stan/surv_cweib_cens sample \
      adapt delta=0.9999 \
      num_samples=30000 num_warmup=30000 thin=30 \
      algorithm=hmc engine=nuts max_depth=15 stepsize=1 \
      id=$i \
      init=0 \
      data file=$f \
      output file=../data/mcmc_out/surv_cweib_cens_${i}.csv &
  done
  wait
done

