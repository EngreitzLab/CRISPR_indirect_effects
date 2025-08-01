## Estimating indirect effects in CRISPR datasets

This analysis workflow empirically computes indirect effects rates from CRISPR datasets analyzed
using sceptre and uses these to compute direct effect rates and the probability of CRISPR E-G pairs 
being the result of a direct effect (vs. indirect, trans-acting effects).

### Required inputs
The workflow computes indirect and direct effect rates from the output of the following workflows:
- https://github.com/EngreitzLab/DC_TAP_Paper
- https://github.com/EngreitzLab/ENCODE_Test_Dataset_Analysis

For integrating with the K562 ENCODE-rE2G training data, it uses indirect effect rates
computed by the workflow generating this combined dataset:
- https://github.com/argschwind/ENCODE_CRISPR_data
