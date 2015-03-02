Chance of fossil preservation
-----------------------------

Taxon ranges through T stages. Each stage has X occurrences of taxon.
X ~ Pois(\lambda).

Set this up hierarchically. Each taxon is in a group (genus in order). Each
group is part of a larger group (order in class). Final level is the super group
(class in hard-part inverts).

Overdispersed Poisson is Negative Binomial (Gamma mixed with Poisson).


Genus duration of marine invertebrates based on class, origination cohort, and background regime
-----------------------------------------------------------------------------

Taxon duration in stages. Survival model with intercept + two hierarchical
terms. Class and cohort within regime.

Idea is estimate how different periods of background are in terms of basic
extinction risk and how/if that might vary within a background regime.

How do the different classes differ in "fundamental extinction risk?" 
