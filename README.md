Code tested on lxplus7.cern.ch. In general it shall work in any reasonably moder Linux with ROOT installed.


# I. Input data collection

The sources are in `data/TOTEM_<energy>` directories. They contain the orginal data sources as well as pre-process input `data.root` inlcuding
  * `g_dsdt`: differential cross-section data points and statistical uncertainties
  * `m_dsdt_cov_syst_t_dep`: covariance matrix with the t-dependent contributions
  * `v_dsdt_syst_t_indep`: vector of t-independent (normalisation) contributions



# II. Definition of t-ranges

The relevant code is in `t_ranges` sub-directory. Run
```
make && ./t_investigation
```
to apply the rules defined in `t_investigation.cc`.



# III. Fits and extrapolation

The script `run_multiple` can be used to execute a series of fits. Edit the file to select the t-ranges, fit models and uncertainty models. The output is saved in directories
`fits/<t_range>/<fit model>/<uncertainty>`. There are two output files:
  * `do_fits.root`: output of dsigma/dt fits at each energy
  * `s_extrapolation.root`: output of extrapolation to D0 energy
