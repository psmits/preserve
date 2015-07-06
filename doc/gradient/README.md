Gradients
=========

Questions/theoretical framework
-------------------------------

Latitudinal diversity gradients, "normal" or "reversed."

Higher extinction and origination at tropics.

Hypotheses based on diversification

-  Museum vs Cradle
  -  Preserve versus create diversity.
-  Out of the Tropics / In to the Tropics
  - Diversify in one and then move into other.

These hypotheses have distinct predictions about the relative birth/death
probabilities, and migration/extirpation, between the environmental
types.

Greater recruitment and/or loss in tropical or temperate zones?

Correlation between recruitment and loss?


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

Model is focused on first entry and last exit, not re-colonization dynamics.
That would require a more traditional HMM, not one with an absorbing state
(death/extinction).

Do I want to distinguish between northern and southern hemispheres?
