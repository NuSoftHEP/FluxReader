## Overview:

In a nutshell, FluxReader reads `Dk2Nu` flux files and outputs a set of spectra. 
For a given set of spectra, each spectrum in the set will plot the same variable(s). 
However, each spectrum will make the distribution for a specific neutrino flavor, parent species, 
applied cross section, and detector location. The user can configure almost everything, 
including (but not necessarily limited to): what flavors, parents, cross sections, detectors, binning, 
labels, variables, and weights to use to create output spectra.

To create its output, FluxReader reprojects each entry in the flux files, or neutrino ray, 
to the particular detector using a function within Dk2Nu itself that also calculates the 
new energy and importance weight. 
The projection is generally done to the same point on the detector, 
typically the center of the front face. For straight flux distributions, 
or when no cross section is applied, the weights of each neutrino ray are summed to output a flux, 
and the units are scaled to be neutrinos per square cm.

For event rate distibutions, or when cross sections are applied, all necessary cross section splines 
are generated from a root file created by the `genie_xsec` product. 
Splines are created for each neutrino flavor, target nucleus, and interaction current. 
The new energy from each reprojected neutrino ray is given to the appropriate spline, 
and the cross section for the given energy is returned. 
The event rate is calculated as `flux * cross-section * number-of-targets`, 
where number of targets is calculated based on standard molar mass and 
an estimate that each detector is a single material. 
The result for each neutrino ray is summed to output an event rate, 
and the units are scaled to be events per kton.

The neutrino fluxes and event rates are not scaled to any POT value; 
however, FluxReader stores the number of POT used to generate its distributions. 
The POT value is stored in a histogram called TotalPOT with a single bin set to the number of POT. 
This means that the user can easily scale to the desired POT by the simple ratio:
```
(user-defined-scale)/TotalPOT->GetBinContent(1)
```

