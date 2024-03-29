# OpponencyModels
Code for the models for direction opponency from Badwan, B.A., Creamer, M.S., Zavatone-Veth, J.A., and Clark, D.A. Dynamic nonlinearities enable direction opponency in *Drosophila* elementary motion detectors. *Nature Neuroscience* 22, 1318–1326 (2019) doi: [10.1038/s41593-019-0443-y](https://doi.org/10.1038/s41593-019-0443-y)


## Feedforward models for direction-opponency

The `RunFigure6Models` function runs all feedforward models, while the `RunFigure6ModelParameterSweep` function runs parameter sweeps for all models with tunable parameters.

## Natural scene velocity discriminability

The natural scene database from Meyer _et al_. 2014 may be downloaded [here (click for link)](https://pub.uni-bielefeld.de/data/2689637). The database is provided as a set of .rar archives, each containing a set of .mat files. Extract the .mat files from the archives using a tool such as [7-Zip](https://www.7-zip.org/), and copy the resulting .mat files into a folder named `imageData` within a root directory, which we refer to as `localPath`. Then, run the `NaturalScenesToContrast` function with an appropriate `localPath` input (the default is the working directory). With that step complete, the `PlotOpponencyInformation` and `SpatialCorrFigure` functions may be run.
