#!/bin/bash
FILES=../data/data_dump/*
for f in $FILES;
do
  n=${f//[^0-9]/};
  for i in `seq 1 4`;
  do
    # this is the poisson model
    ./neg_bin_mod sample random seed=420 \
      id=$i \
      data file=$f \
      output file=../data/mcmc_out/nbn_samp_${n}_${i}.csv &
  done
  wait
  grep lp__ ../data/mcmc_out/nbn_samp_${n}_1.csv > \
    ../data/mcmc_out/combined_nbn_${n}.csv
  sed '/^[#1]/d' ../data/mcmc_out/nbn_samp_{$n}_*.csv >> \
    ../data/mcmc_out/combined_nbn_${n}.csv
done
