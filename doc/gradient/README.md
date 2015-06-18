Gradients
=========

Data
----

Paleozoic brachiopod occurrence information. Sourced from Foote/Miller work.

Split based on geographic provinces. Currently N Tropical, N Temperate,
S Tropical and S Temperate (defined at 22degree). Would prefer better structure
to the provinces so that tropical/temperate can be higher-level binary
predictor. 


Approach
--------

Hierarchical Jolly-Seber capture-mark-recapture model. Both within and between
geographic provinces.

Need to include taxon effect as taxa can be in more than one province which
violates exchangability.

Need to include binary predictor of "seen else where?"
