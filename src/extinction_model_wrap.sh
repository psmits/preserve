#!/bin/bash
# all taxa models
# switch to just the fauna models
FILES=/home/psmits/preserve/data/data_dump/impute*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    /home/psmits/preserve/stan/survival_edit \
      sample num_samples=50000 num_warmup=50000 thin=50 \
      init=0 \
      id=$i \
      data file=$f \
      output file=/home/psmits/preserve/data/mcmc_out/faun_weib_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    /home/psmits/preserve/stan/survival_relab \
      sample num_samples=50000 num_warmup=50000 thin=50 \
      init=0 \
      id=$i \
      data file=$f \
      output file=/home/psmits/preserve/data/mcmc_out/faun_relab_${i}.csv &
  done
  wait
  for i in `seq 1 4`;
  do
    /home/psmits/preserve/stan/survival_impute \
      sample num_samples=50000 num_warmup=50000 thin=50 \
      init=0 \
      id=$i \
      data file=$f \
      output file=/home/psmits/preserve/data/mcmc_out/faun_impute_${i}.csv &
  done
  wait
done
FILES=/home/psmits/preserve/data/data_dump/high*
for f in $FILES;
do
  for i in `seq 1 4`;
  do
    /home/psmits/preserve/stan/survival_highgrade \
      sample num_samples=50000 num_warmup=50000 thin=50 \
      init=0 \
      id=$i \
      data file=$f \
      output file=/home/psmits/preserve/data/mcmc_out/faun_highgrade_${i}.csv &
  done
  wait
done
