# Quantile Mapping
This repository contains tools to derive and apply corrections based on quantile mapping.

The quantile mapping is defined via two cumulative distribution functions (CDF), i.e. the *source* CDF derived from the 'wrong' distribution and the *target* CDF derived from the 'correct' distribution. Usually the source is based on MC and the target on data.

## Derivation
The CDFs are modelled as splines (TSpline3) here in order to fit deliberate shapes and to produce continuous distributions at the same time. Histograms of the source and the target distributions (not CDFs) need to be provided in a ROOT file. Dealing with interpolating splines here, those histograms should not suffer from large statistical uncertainties, which should be considered in the choice of the binning.

The script `create_CDF_splines.py` calculates the CDF splines from a deliberate number of input distributions and stores them in a separate ROOT file which will later be used for the application. At this stage it is not specified which CDFs are sources and which ones are targets since they are technically equal structures and the transformations could be done in both directions.

The derivatives of the splines at the endpoints are set to zero as boundary conditions. If the considered distributions do not go to zero at the edges (e.g. due to applied cuts), this is suboptimal and should be adapted. Alternative boundary conditions can either concern the first or the second derivative (see options of TSpline3).

## Application
This repository provides a standalone Python script `apply_quantile_mapping.py` which can be used to apply corrections to plain ROOT trees. It appends a branch with the shifted values to the tree. The actual shifting procedure of a single value is performed by instances of the class `QuantileShifter` implemented in `QuantileShifter.py`. It can therefore be imported and used in more complex applications. Instanciating this class requires the path to the ROOT file containg CDF splines and the names of the intended source and target splines, which defines the correction performed by this instance.

If this package is not embedded into a CMSSW setup, `<parent_path>/quantile_mapping` should be added to your PYTHONPATH.

The `QuantileShifter` class is also available as c++ module in `quantile_shifter.h` and `quantile_shifter.cc`.
In addition to the constructor asking for the CDFs, as used in the Python class, an additional function `init` with corresponding arguments and a default constructor are available, which enables to define the mapping independently from the declaration of the QuantileShifter instance.

Events that are out of the range of the source splines are not shifted.

### Inversion method
The necessary inversion of the target CDF spline is done numerically. By default, the spline is narrowed by linear interpolations between the knots in a first step. The result for a given quantile is then corrected by the difference to the spline scaled by the inverse derivative of the spline at this point.

In regions with small derivatives, which are usually the margins, this method can lead to large deviations and bad results. For this reason, a different inversion method is applied if the default method includes a correction step larger than half of the distance between the surrounding spline knots or a correction step passing one of the knots: In this case, the interval between two knots that contains the solution is bisected in five iterations (assuming that the CDF spline is monotonic) and in the final interval the spline is linearly interpolated to determine the solution.

The latter inversion method runs a bit slower but is more robust. Instances of the `QuantileShifter` can be created with an additional flag set to `true` in order to use this method for all events.
