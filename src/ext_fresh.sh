#!/bin/bash
for i in `seq 1 4`;
do
  ../stan/surv_fresh \
    sample num_samples=1000 num_warmup=1000 \
    id=$i \
    init=0 \
    data file=../data/data_dump/clean_info.data.R \
    output file=../data/mcmc_out/clean_${i}.csv &
done
