# dghrank

This package is intended for usage by California Prenatal Quality Care Collective (CPQCC) and Vermont Oxford Network (VON) for institutional quality analysis based on performance indicators. The primary function, fitBabyMonitor(), takes in a matrix of individual data and returns data.frames of institution and (if desired) subset summary statistics and quality scores.


# Usage

The source code for *dghrank* may be found [here](https://github.com/dhelkey/babymonitor) and can be loaded into R using the *devtools* pacakge by:

```
devtools::install_github('dhelkey/babymonitor', force = TRUE)
library('babymonitor')
```

The function *babyMonitorSummary* may be used to construct standardized scores for a given quality indicator.
Composite scores are created by generaating standardized components using *fitBabyMonitor* and combing the results with *scoreComposite*.

# Model





# Standardizing



# Composite Score







