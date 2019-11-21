from info import basic, paths, config, zbins, alt_zbins, plotting, source_nofz_pars, sysmaps
from ggl import GGL, Measurement, ResponsesScale, ResponsesProjection, TestStars, TestPSF, TestSizeSNR, TestSysMaps

T = True
F = False


run_measurement = T
run_responses_nk = F
run_responses_ng = F
run_stars = F
run_psf = F
run_size_snr = F
run_sysmaps = F


if 'combined_sample_fid' in config['lens_v']:
    zbinning = zbins
if 'maglim' in config['lens_v']:
    zbinning = alt_zbins

if run_measurement:
    print 'Starting measurement class...'
    print 'zbinning used:', zbinning
    measurement = Measurement(basic, config, paths, zbinning, plotting)

    if not basic['plot_blinded']:
        measurement.run()
        measurement.save_2pointfile('gt')
        measurement.save_2pointfile('gt_boosted')
        measurement.save_2pointfile('boost_factor')
        if not basic['blind']:
            measurement.plot()
        measurement.plot_boostfactors()
        measurement.plot_randoms()
        measurement.plot_gammax()
	
    if basic['blind'] and basic['plot_blinded']:
        measurement.compute_sn_ratio('gt')
        measurement.compute_sn_ratio('gt_boosted')
        measurement.plot_from_twopointfile('gt')
        measurement.plot_from_twopointfile('gt_boosted')

if run_responses_nk:
    responses = ResponsesScale(basic, config, paths, zbinning, plotting)
    responses.run()
    responses.plot('lens', mask_scales =False)
    responses.plot('lens', mask_scales =True)
    #responses.plot_sigmas('lens', mask_scales =False)
    #responses.plot('random')

if run_responses_ng:
    responses = ResponsesProjection(basic, config, paths, zbinning, plotting)
    responses.run()

if run_stars:
    stars = TestStars(basic, config, paths, zbinning, plotting)
    stars.run('bright')
    stars.run('faint')
    stars.plot()

if run_psf:
    ggl_int = init_ggl_class(basic,config, paths, zbins)
    testpsf = TestPSF(ggl_int, config,paths, zbins, plotting)
    bands = ['r', 'i', 'z']
    testpsf.run_y3(bands)
    testpsf.plot(bands)


if run_size_snr:
    size_snr = TestSizeSNR(basic, config, paths, zbinning, plotting, source_nofz_pars)
    size_snr.run('size')
    #size_snr.run('snr')
    size_snr.plot()

if run_sysmaps:
    ggl_int = init_ggl_class(basic, config, paths, zbins)
    testsys = TestSysMaps(ggl_int, config, paths, zbins, plotting, source_nofz_pars, sysmaps)
    bands = ['r', 'i', 'z']
    testsys.run_y3(['airmass', 'fwhm', 'maglimit', 'skybrite'], bands)
    testsys.plot_y3(bands)

