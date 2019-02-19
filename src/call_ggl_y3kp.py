from info import basic, paths, config, zbins, plotting, source_nofz_pars, sysmaps
from ggl import GGL, Measurement, Responses, TestStars, TestPSF, TestSizeSNR, TestSysMaps

T = True
F = False

run_measurement = T
run_responses_nk = F
run_stars = F
run_psf = F
run_size_snr = F
run_sysmaps = F


if run_measurement:
    print 'Starting measurement class...'
    #gglensing = GGL(basic, config, paths)
    measurement = Measurement(basic, config, paths, zbins, plotting)
    if not basic['plot_blinded']:
        measurement.run()
        #measurement.save_boostfactors_2pointfile() #deprecated, without errors
        measurement.save_2pointfile('gt')
        measurement.save_2pointfile('boost_factor')
        if not basic['blind']:
            measurement.plot()
        measurement.plot_boostfactors()
        measurement.plot_randoms()
        measurement.plot_gammax()

    if basic['blind'] and basic['plot_blinded']:
        measurement.compute_sn_ratio()
        measurement.plot_from_twopointfile()

if run_responses_nk:
    responses = Responses(basic, config, paths, zbins, plotting)
    responses.run()
    #responses.plot('lens')
    #responses.plot_sigmas('lens')
    #responses.plot('random')

if run_stars:
    stars = TestStars(basic, config, paths, zbins, plotting)
    stars.run('bright')
    stars.run('faint')
    stars.plot()

if run_psf:
    psf = TestPSF(basic, config, paths, zbins, plotting)
    ra_lims = (-1, 361)
    dec_lims = (-90, -35)
    psf.save_psf_residuals_y3(ra_lims, dec_lims)
    # psf.run()
    # psf.plot()

if run_size_snr:
    size_snr = TestSizeSNR(basic, config, paths, zbins, plotting, source_nofz_pars)
    size_snr.run('size')
    size_snr.run('snr')
    size_snr.plot()

if run_sysmaps:
    sysmaps = TestSysMaps(basic, config, paths, zbins, plotting, source_nofz_pars, sysmaps)
    sysmaps.run(['airmass', 'fwhm', 'maglimit', 'skybrite'], ['r'])
    sysmaps.plot()
