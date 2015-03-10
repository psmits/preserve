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
  grep lp__ ../data/mcmc_out/surv_1.csv > \
    ../data/mcmc_out/combined_surv.csv
  sed '/^[#1]/d' ../data/mcmc_out/surv_*.csv >> \
    ../data/mcmc_out/combined_surv.csv
  for i in `seq 1 4`;
  do
    ./complex_surv_mod sample num_samples=10000 num_warmup=10000 thin=10 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/comp_surv_${i}.csv &
  done
  wait
  grep lp__ ../data/mcmc_out/comp_surv_1.csv > \
    ../data/mcmc_out/combined_comp_surv.csv
  sed '/^[#1]/d' ../data/mcmc_out/comp_surv_*.csv >> \
    ../data/mcmc_out/combined_comp_surv.csv
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
  grep lp__ ../data/mcmc_out/faun_surv_1.csv > \
    ../data/mcmc_out/combined_faun_surv.csv
  sed '/^[#1]/d' ../data/mcmc_out/faun_surv_*.csv >> \
    ../data/mcmc_out/combined_faun_surv.csv
  for i in `seq 1 4`;
  do
    ./complex_fauna_surv sample num_samples=10000 num_warmup=10000 thin=10 \
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/comp_faun_surv_${i}.csv &
  done
  wait
  grep lp__ ../data/mcmc_out/comp_faun_surv_1.csv > \
    ../data/mcmc_out/combined_comp_faun_surv.csv
  sed '/^[#1]/d' ../data/mcmc_out/comp_faun_surv_*.csv >> \
    ../data/mcmc_out/combined_comp_faun_surv.csv
done
