library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(rstan)
library(survival)
library(stringr)

counts <- read_rdump('../data/data_dump/count_info.data.R')

orders <- counts$order[counts$genus]
count.df <- with(counts, data.frame(count = count, genus = factor(genus), 
                                    orders = factor(orders)))

sum.stat <- ddply(count.df, .(orders), summarize, 
                  x.bar = mean(count), 
                  n.count = length(count),
                  n.zed = sum(count == 0),
                  n.gen = length(unique(genus)))

