#!/bin/bash
FILES=../data/data_dump/survival*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    # this is the poisson model
    ./survival_model sample num_samples=10000 num_warmup=10000 thin=10 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/surv_${i}.csv &
  done
  wait
  grep lp__ ../data/mcmc_out/surv_1.csv > \
    ../data/mcmc_out/combined_surv.csv
  sed '/^[#1]/d' ../data/mcmc_out/surv_*.csv >> \
    ../data/mcmc_out/combined_surv.csv
done
FILES=../data/data_dump/fauna*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    # this is the poisson model
    ./survival_model sample num_samples=10000 num_warmup=10000 thin=10 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/faun_surv_${i}.csv &
  done
  wait
  grep lp__ ../data/mcmc_out/faun_surv_1.csv > \
    ../data/mcmc_out/combined_faun_surv.csv
  sed '/^[#1]/d' ../data/mcmc_out/faun_surv_*.csv >> \
    ../data/mcmc_out/combined_faun_surv.csv
done
