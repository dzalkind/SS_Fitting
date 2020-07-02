import numpy as np
import numpy.matlib as ml
import matplotlib.pyplot as plt 
from scipy import interpolate

class SSFit_Radiation(object):

    def __init__(self, **kwargs):

        self.InputFile = ''
        self.InputVars = dict()

        # Optional population of class attributes from key word arguments
        for (k, w) in kwargs.items():
            try:
                setattr(self, k, w)
            except:
                pass

        super(SSFit_Radiation, self).__init__()

    
    def readInputs(self):

        with open(self.InputFile,'r') as f:

            f.readline()
            self.InputVars['HydroFile']      = f.readline().split('%')[0].strip()
            self.InputVars['DOFs']          = np.fromstring(f.readline().split('%')[0].strip()[1:-1],dtype=int,sep = ' ')
            self.InputVars['FreqRange']     = np.fromstring(f.readline().split('%')[0].strip()[1:-1],dtype=float,sep = ' ')
            self.InputVars['WeightFact']    = float(f.readline().split('%')[0])
            self.InputVars['IDMethod']      = int(f.readline().split('%')[0])
            self.InputVars['FitReq']        = float(f.readline().split('%')[0])
            self.InputVars['PlotFit']       = int(f.readline().split('%')[0])
            self.InputVars['ManReduction']  = int(f.readline().split('%')[0])

    def setWeights(self,omega):
        
        omega_Range             = self.InputVars['FreqRange']

        weight_Ind              = np.bitwise_and(omega >= omega_Range[0],omega <= omega_Range[1])

        weight                  = np.zeros([len(omega),1])
        weight[weight_Ind]      = self.InputVars['WeightFact']
        weight[~weight_Ind]     = 1 - self.InputVars['WeightFact']

        #normalize & save
        self.weights            = weight/sum(weight)


class WAMIT_Out(object):
    def __init__(self,**kwargs):
        self.HydroFile   = ''
        self.Type       = None

        # Optional population of class attributes from key word arguments
        for (k, w) in kwargs.items():
            try:
                setattr(self, k, w)
            except:
                pass

        super(WAMIT_Out, self).__init__()


    def readFRFs(self):

        if self.Type == 1:

            with open(self.HydroFile+'.1','r') as f:
                all_data = [x.split() for x in f.readlines()]
                
                per     = [float(x[0]) for x in all_data]
                I       = [int(x[1]) for x in all_data]
                J       = [int(x[2]) for x in all_data]
                Abar_ij = [float(x[3]) for x in all_data]
                Bbar_ij = np.zeros([1,len(Abar_ij)]).tolist()[0]

                for iPer, p in enumerate(per):
                    if p > 0:
                        Bbar_ij[iPer] = all_data[iPer][4]

                # check that all the lengths are the same
                if len(per) != len(I) != len(J) != len(Abar_ij) != len(Bbar_ij):
                    print('The lengths of per, I, J, Abar, and Bbar are not the same!')


                # put Abar_ij, Bbar_ij into 6x6x(nFreq) matrices
                Abar = np.zeros([6,6,len(np.unique(per))])
                Bbar = np.zeros([6,6,len(np.unique(per))])

                Period = np.unique(per)

                # loop through because I don't know a cleverer way
                # note that the index will be -1 compared to matlab and the normal way of writing them
                # A(i,j,per)
                for k, p in enumerate(per):
                    Abar[I[k]-1,J[k]-1,np.where(Period==p)[0]] = Abar_ij[k]
                    Bbar[I[k]-1,J[k]-1,np.where(Period==p)[0]] = Bbar_ij[k]
                
                
                # Break Abar into A_inf and Abar(\omega)
                Abar_inf        = np.squeeze(Abar[:,:,Period==0])
                Abar            = Abar[:,:,Period!=0]
                Bbar            = Bbar[:,:,Period!=0]
                omega           = 2 * np.pi / Period[Period!=0]
                omega[omega<0]  = 0
                
                # Sort based on frequency because
                freqIndSort     = np.argsort(omega)
                omega           = omega[freqIndSort]
                Abar            = Abar[:,:,freqIndSort]
                Bbar            = Bbar[:,:,freqIndSort]



                if False:
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    ax.plot(omega,Bbar[5,5,:])
                    plt.show()
                

                ### Dimensionalize
                # Coefficients
                rho     = 1000. # water density as set in matlab script
                L_ref   = 1.    # reference length as set in WAMIT
                
                # scale factor depends on dimension
                # Abar_ij = A_ij / (\rho L^k), Bbar_ij = B_ij / (\rho L^k \omega)
                # k = 3 for i,j <= 3, k = 5 for i,j > 3, and k = 4 for cross terms, kind of moot since L = 1
                kk      = np.array([(3,3,3,4,4,4),\
                                    (3,3,3,4,4,4),\
                                    (3,3,3,4,4,4),\
                                    (4,4,4,5,5,5),\
                                    (4,4,4,5,5,5),\
                                    (4,4,4,5,5,5)], dtype=int)

                # scaling matrices
                scaleA  = np.reshape(np.tile(rho * L_ref**kk,[1,1,len(omega)]),np.shape(Abar))
                scaleB  = [rho * L_ref**kk * om for om in omega ]  #close! but we need to rearange indices
                scaleB  = np.transpose(scaleB,[1,2,0])

                # scale, finally
                A       = Abar * scaleA
                A_inf    = Abar_inf * scaleA[:,:,0]
                B       = Bbar * scaleB

                # Form K
                K = B + 1j * np.tile(omega,[6,6,1]) * (A - np.transpose(np.tile(A_inf,[len(omega),1,1]),[1,2,0]))



                # Interpolate to get evenly spaced matrices
                exp     = 8
                om_samp = np.linspace(0,max(omega),2**exp)

                fA      = interpolate.interp1d(omega,A)
                fB      = interpolate.interp1d(omega,B)
                fK      = interpolate.interp1d(omega,K)
                
                A_samp   = fA(om_samp)
                B_samp   = fB(om_samp)
                K_samp   = fK(om_samp)




                if False:
                    plt.plot(omega,K.imag[3,3,:],'o',om_samp,K_samp.imag[3,3,:],'-')
                    plt.show()


                # Set insignificant entries of K to zero
                # Doing this the same way as matlab for now

                # Find maximum diagonal element of damping (B) matrix
                B_max = np.diag(B.max(2)).max()

                # if any diagonal elements are less than 0.0001 * B_max  -> 0
                # if any non-diagonal elements are less than 0.0005 * B_max -> 0
                for i in range(0,6):
                    for j in range(0,6):
                        if i == j:
                            if max(abs(B[i,j,:])) < 0.0001 * B_max:
                                K[i,j,:] = 0
                        else:
                            if max(abs(B[i,j,:])) < 0.0005 * B_max:
                                K[i,j,:] = 0

                
                # store values
                self.omega      = om_samp
                self.K          = K_samp
                self.A          = A_samp
                self.A_inf      = A_inf
                self.B          = B_samp
                self.B_max      = B_max
                
                print('here')



        return





if __name__=='__main__':

    ss_fitRad = SSFit_Radiation()

    ss_fitRad.InputFile = '/Users/dzalkind/Tools/SS_Fitting/test_cases/Haliade/SS_Fitting_Options.inp'

    # read SS_Fitting inputs
    ss_fitRad.readInputs()

    # get WAMIT info
    wam = WAMIT_Out(HydroFile=ss_fitRad.InputVars['HydroFile'], Type=1)
    wam.readFRFs()   # read output file

    # set weights
    ss_fitRad.setWeights(wam.omega)

    print('here')
    
