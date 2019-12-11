Code tested on lxplus7.cern.ch. In general it shall work in any reasonably modern Linux with ROOT installed.


# I. Input data collection

The sources are in `data/TOTEM_<energy>` directories. They contain the original data sources as well as pre-processed input `data.root` including
  * `g_dsdt`: differential cross-section data points and statistical uncertainties
  * `m_dsdt_cov_syst_t_dep`: covariance matrix with the t-dependent contributions to the systematic uncertainty
  * `v_dsdt_syst_t_indep`: vector of t-independent (normalisation) contributions to the systematic uncertainty



# II. Definition of t-ranges

Edit `run_multiple` such that the model "bootstrap/bootstrap" is enabled.

Then go to the `t_ranges` sub-directory and run
```
make && ./t_investigation_new
```
The output can then be used to update `datasets.h`.


# III. Fits and extrapolation

Run
```
./run_multiple
```
to execute a series of fits. Edit the file to select the t-ranges, fit models and uncertainty models of your interest. The output is saved in directories
```
fits/<t_range>/<fit model>/<uncertainty>
```
There are two output files:
  * `do_fits.root`: output of dsigma/dt fit at each energy
  * `s_extrapolation.root`: output of extrapolation to D0 energy
