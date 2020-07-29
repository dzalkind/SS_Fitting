from SS_FitTools import SSFit_Excitation, WAMIT_Out, idDOF, FDI_Fitting, TimeDomain_Fit


if __name__=='__main__':

    ss_fitExt = SSFit_Excitation()

    ss_fitExt.InputFile = '/Users/dzalkind/Tools/SS_Fitting/test_cases/OC3_Spar/Excitation_Options.inp'

    # read SS_Fitting inputs
    ss_fitExt.readInputs()

    # get WAMIT info
    wam = WAMIT_Out(HydroFile=ss_fitExt.InputVars['HydroFile'], Type=ss_fitExt.WamType)
    wam.readFRFs(ss_fitExt.InputVars['WaveDir'],ss_fitExt.InputVars['Grav'],ss_fitExt.InputVars['WtrDen'])   # read output file

    # find causal impulse response function (IRF)
    ss_fitExt.causalIRF(wam,Plot=False)

    # fit state space model based using causal IRF
    tdf = TimeDomain_Fit(ss_fitExt.Kt,ss_fitExt.tt,ss_fitExt.InputVars['FitReq'])
    tdf.fit()

    # write output matrices
    ss_fitExt.writeMats(tdf.K_fit)


    # print(idDOF(5))

    # # FDI Fitting
    # fdi = FDI_Fitting()
    # fdi.fit(ss_fitExt,wam)

    # # Visualize Fit
    # fdi.visualizeFits(wam)

    # ss_fitExt.outputMats(fdi)
    # print('here')
    
