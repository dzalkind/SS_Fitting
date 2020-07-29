from SS_FitTools import SSFit_Radiation, WAMIT_Out, idDOF, FDI_Fitting


if __name__=='__main__':

    ss_fitRad = SSFit_Radiation()

    ss_fitRad.InputFile = '/Users/dzalkind/Tools/SS_Fitting/test_cases/UMaine_Semi/SS_Fitting_Options.inp'

    # read SS_Fitting inputs
    ss_fitRad.readInputs()

    # get WAMIT info
    wam = WAMIT_Out(HydroFile=ss_fitRad.InputVars['HydroFile'], Type=1)
    wam.readFRFs()   # read output file

    # set weights
    ss_fitRad.setWeights(wam.omega)

    print(idDOF(5))

    # FDI Fitting
    fdi = FDI_Fitting()
    fdi.fit(ss_fitRad,wam)

    # Visualize Fit
    fdi.visualizeFits(wam)

    ss_fitRad.outputMats(fdi)
    print('here')
    
