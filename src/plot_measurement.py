from info import basic, paths, config, zbins, plotting, source_nofz_pars, sysmaps
from ggl import GGL, Measurement, ResponsesScale, ResponsesProjection, TestStars, TestPSF, TestSizeSNR, TestSysMaps


print 'Starting measurement class, plotting...'
print 'zbinning used:', zbins
measurement = Measurement(basic, config, paths, zbins, plotting)

measurement.plot_boostfactors()
measurement.plot_boost_factors_cov()
measurement.plot_randoms()
measurement.plot_gammax()	
measurement.compute_sn_ratio('gt')
measurement.compute_sn_ratio('gt_boosted')
measurement.plot_from_twopointfile('gt')
measurement.plot_from_twopointfile('gt_boosted')

