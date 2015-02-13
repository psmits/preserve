#!/bin/bash
FILES=../data/data_dump/*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    # this is the poisson model
    ./neg_bin_mod sample num_samples=10000 num_warmup=10000 thin = 10 \\
      random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/nbn_samp_${i}.csv &
  done
  wait
  grep lp__ ../data/mcmc_out/nbn_samp_1.csv > \
    ../data/mcmc_out/combined_nbn.csv
  sed '/^[#1]/d' ../data/mcmc_out/nbn_samp_*.csv >> \
    ../data/mcmc_out/combined_nbn.csv
done
