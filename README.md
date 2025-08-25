# microseisms_LOPS
## purpose
This repository is designed to transfer old matlab codes that were used for papers such as Ardhuin (GRL 2018) "Large-Scale Forces Under Surface Gravity Waves at a Wavy
Bottom: A Mechanism for the Generation of Primary Microseisms", first focusing on primary microseisms and hum. Note that the scripts used to generate the plots in the 2018 paper had a error: there was an extra alpha factor, and the seismic response was overestimated by a factor about 10 : this means that it actually takes a 10 times larger rocky area to generate the measured microseisms. 

This repository is not designed to be as user friendly as WMSAN (https://github.com/lystom/WMSAN) but that other package does not yet include primary microseisms. 

## general principle
Using information about: 
- ocean waves (wind-generated waves) : ideally we need the full frequency-direction spectrum E(f,theta). Assuming some directional distribution one could also use E(f), making the manipulation of maps of sources a bit easier. Note that maps of E(f,theta) can be obtained from ECMWF (either for the ERA5 reanalysis or the operational forecasts, more on this later... ) 
- bottom topography                  : reprensented by a bottom elevation PSD (see Ardhuin 2018). 
- solid Earth structure              : represented by a C2 (as in Longuet-Higgins 1950) coupling coefficient, here based on a solid half-space (the water layer can be neglected due to the interaction in relatively shallow water), and attenuation coeffients (Q) and a propagation term meant to reprensent source and receiver site effects. 

a synthetic spectrogram is generated and can be compared to measurements.

## example script
comp_primary_from_spec.ipynb : compares synthetic microseisms at Geoscope station EDA for year 2023 with measurements. This uses an equivalent area dA for the bottom-wave interactions over which both the bottom and wave spectra are assumed homogeneions 

## installation
git clone --quiet https://github.com/ardhuin/microseisms_LOPS

## more later
