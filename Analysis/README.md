# CONTROL PLOTS

## Using controlPlots.py to create control plots during different stages of cut-flow

The procedure of `controlPlots.py` script is to read:
* NanoAOD files for data and MC 
* or read pkl files with the desired branches that were already produced by this same script

to generated control plots from selected observables.

Firts of all, execute:
```bash
. /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc11-opt/setup.sh
```

To generate only pkl files with the desired branches, turn on the flag [OnlyGeneratePKLFiles](https://github.com/cmunozdi/XConeReclustering/blob/a02bb5ca9d64fdbfd85cb4728cabfe21cdeaeabf/Analysis/controlPlots.py#L19). 

Additionaly, use the [flags](https://github.com/cmunozdi/XConeReclustering/blob/a02bb5ca9d64fdbfd85cb4728cabfe21cdeaeabf/Analysis/controlPlots.py#L9-L17) for the different MC contributions to read from NanoAOD (on) or from preexisting pkl (off).

[Indicate](https://github.com/cmunozdi/XConeReclustering/blob/a02bb5ca9d64fdbfd85cb4728cabfe21cdeaeabf/Analysis/controlPlots.py#L21) the output directory where the plots will come out.

Include the branches you are interested in for [MC](https://github.com/cmunozdi/XConeReclustering/blob/a02bb5ca9d64fdbfd85cb4728cabfe21cdeaeabf/Analysis/controlPlots.py#L635) and [data](https://github.com/cmunozdi/XConeReclustering/blob/a02bb5ca9d64fdbfd85cb4728cabfe21cdeaeabf/Analysis/controlPlots.py#L639)

At the [end](https://github.com/cmunozdi/XConeReclustering/blob/a02bb5ca9d64fdbfd85cb4728cabfe21cdeaeabf/Analysis/controlPlots.py#L1640-L1670) of the script, call the function `plot_variable` as in the examples to plot the different branches. You can also move them out of the cutflow for-loop, and select the datasets at the specific stage you are interested in.