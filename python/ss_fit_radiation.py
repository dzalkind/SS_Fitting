import numpy as np
import numpy.matlib as ml
import matplotlib.pyplot as plt 
from invfreqs import invfreqs
from scipy import interpolate
from scipy.signal import freqs, bode, tf2ss, ss2tf
from scipy.linalg import block_diag, eig
from datetime import datetime

# from control import tf2ss, ss, ss2tf
# from control.matlab import bode

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

    def outputMats(self,sysID):
        
        # Find number of states in each sys
        states = [np.shape(S[0])[0] for S in sysID.sys]

        # Find number of states per dof (first index)
        statesPerDoF = np.zeros(6, dtype=int)

        for i in range(0,6):
            for k, DoF in enumerate(sysID.sysDOF):
                if DoF[0] == i:
                    statesPerDoF[i] += states[k]


        # states = np.zeros(len(sysID.sys))
        # for i,S in enumerate(sysID.sys):
        #     states[i] = np.shape(A.A)[1]


        # Make block diagonal matrices
        AA = np.empty(0)
        for S in sysID.sys:
            AA = block_diag(AA,S[0])

        # Remove first row because of initialization 
        AA = AA[1:]
        print('here')

        # Input BB and output CC matrix, first initialize: 6 inputs and outputs
        BB = np.zeros([1,6])
        CC = np.zeros([6,1])

        for k, (S, dof) in enumerate(zip(sysID.sys,sysID.sysDOF)):
            B_k = np.squeeze(np.zeros((states[k],6)))
            B_k[:,dof[0]]= np.squeeze(np.asarray(S[1]))
            BB = np.concatenate((BB,B_k),axis=0)

            C_k = np.squeeze(np.zeros((6,states[k])))
            C_k[dof[1],:] = -np.squeeze(np.asarray(S[2]))
            CC = np.concatenate((CC,C_k),axis=1)
            # print('here')
            
        BB = BB[1:,:] 
        CC = CC[:,1:]

        # check stability
        tfAll   = ss2tf(AA,BB,CC,np.zeros([6,6]))
        ssAll   = tf2ss(tfAll[0],tfAll[1])

        # via the roots of denominator of tfAll
        if any(np.roots(tfAll[1]) > 0):
            print('Warning: SS system unstable.  Try with a lower R^2 value or check inputs and be careful using some DoFs')

        # check if proper
        #  order of numerator       order of denominator
        if np.shape(tfAll[0])[1] > np.shape(tfAll[1])[0]:
            print('Warning: SS system is not proper')

        # check if passive: TODO: DO
        # started, but went back and forth on whether to use control or scipy.signal
        # print('here')
        # for i in range(0,6):
        #     mag = bode(sysAll[i,i],np.linspace(0,5))

        # write output files
        with open(self.InputVars['HydroFile']+'_py.ss','w') as f:


            # header info
            now = datetime.now()
            print('here')
            f.write('{}\n'.format('SS_Fitting v1.00.01: State-Spaces Matrices obtained using FDI method, on ' + now.strftime('%b-%d-%Y %H:%M:%S')))
            f.write('{}\t\t\t{}\n'.format(np.array_str(self.InputVars['DOFs'])[1:-1],'%Enabled DoFs'))
            f.write('{}\t\t\t\t\t{}\n'.format(sum(states), '%Radiation States'))
            f.write('{}\t\t\t{}\n'.format(np.array_str(statesPerDoF)[1:-1], '%Radiation States per DoF'))

            np.savetxt(f,AA,fmt='%.6e')
            np.savetxt(f,BB,fmt='%.6e')
            np.savetxt(f,CC,fmt='%.6e')


        print('here')


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
                Bbar_max = np.diag(Bbar.max(2)).max()

                # if any diagonal elements are less than 0.0001 * Bbar_max  -> 0
                # if any non-diagonal elements are less than 0.0005 * Bbar_max -> 0
                for i in range(0,6):
                    for j in range(0,6):
                        if i == j:
                            if max(abs(Bbar[i,j,:])) < 0.0001 * Bbar_max:
                                K_samp[i,j,:] = 0
                                print(idDOF(i) + ' - Unimportant DoF, not fit')
                        else:
                            if max(abs(Bbar[i,j,:])) < 0.0005 * Bbar_max:
                                K_samp[i,j,:] = 0

                
                # store values
                self.omega      = om_samp
                self.K          = K_samp
                self.A          = A_samp
                self.A_inf      = A_inf
                self.B          = B_samp
                self.Bbar_max   = Bbar_max
                
                print('here')



        return

class FDI_Fitting(object):

    def __init__(self,**kwargs):

            #Input option structure
        self.OrdMax = 7           # - Maximum order of transfer function to be considered. Typical value 20.
        self.AinfFlag = 1          # - if set to 1, the algorithm uses Ainf (3D hydrodynamic data),  #if set to 0, the algorithm estimates Ainf (2D hydrodynamic data)
        self.Method = 2            # - There are 3 parameter estimation methods (1,2,3). See help of fit_siso_fresp. Recomended method 2 (best trade off between accuracy and speed)
        self.Iterations = 20       # - Related to parameter estimation methods 2 and 3. See help of fit_siso_fresp. Typical value 20.
        self.PlotFlag = 0          # - If set to 1 the function plots the results of each iteration of the automatic order detection in a separate figure.
        self.LogLin = 1            #  - logarithmic or linear frequency scale for plotting.
        self.wsFactor = 0.1        # - Sample faster than the Hydrodynamic code for plotting. Typical value 0.1.
        self.wminFactor = 0.1      # - The minimum frequency to be used in the plot is self.wminFactor*Wmin, where  #Wmin is the minimum frequency of the dataset used for identification.%Typical value 0.1.
        self.wmaxFactor = 5        #- the maximum frequency to be used in the plot is self.wmaxFactor*Wmax, where Wmax is the maximum frequency of the dataset used for identification. Typical value 5.
        
        # Optional population of class attributes from key word arguments
        for (k, w) in kwargs.items():
            try:
                setattr(self, k, w)
            except:
                pass

        super(FDI_Fitting, self).__init__()

    def fit(self,ss_fit,wam):
        # inputs: ss_fit  - state space fitting options
        #         wam     - wamit(-like) hydrodynamic coefficients
        #

        # unpack mass and damping matrices
        A       = wam.A
        B       = wam.B
        # K       = wam.K
        A_inf   = wam.A_inf

        sys     = []   #initialize state space control
        sysDOF  = []
        Khat    = []

        for i in range(0,6):
            for j in range(0,6):
                if max(abs(wam.K[i,j,:])) > 0:
                    
                    # one index pair at a time
                    A_ij = A[i,j,:]
                    B_ij = B[i,j,:]

                    # use only weighted frequency indices
                    useFreq = np.bitwise_and(wam.omega > ss_fit.InputVars['FreqRange'][0], wam.omega < ss_fit.InputVars['FreqRange'][1])
                    A_ij    = A_ij[useFreq]
                    B_ij    = B_ij[useFreq]
                    om      = wam.omega[useFreq]

                    # Compute Memory Function Frequency Response K(j\omega)
                    self.K_ij       = B_ij + 1j * om * (A_ij - A_inf[i,j])

                    # Initial order
                    self.order = 2

                    K_Num, K_Den, Khat_ij = self.ident_memory(om,Plot=False)
                    

                    # compute error
                    r2b, r2a = computeError(self.K_ij,Khat_ij,om,A_inf[i,j])

                    # Increase order of transfer function until R^2 satisfied, unless maximum order reached
                    # Currently, only A or B needs to match, as it is written in matlab, could update
                    while (self.order < self.OrdMax) and  \
                         ((r2b < ss_fit.InputVars['FitReq']) and (r2a < ss_fit.InputVars['FitReq'])):

                         self.order += 1
                         K_Num, K_Den, Khat_ij = self.ident_memory(om,Plot=False)

                         r2b, r2a = computeError(self.K_ij,Khat_ij,om,A_inf[i,j])

                    # Convert to state-space system
                    sys.append(tf2ss(K_Num,K_Den))
                    sysDOF.append([i,j])

                    _, Khat_samp  = freqs(K_Num,K_Den,worN=wam.omega)

                    Khat.append(Khat_samp)

                    # Check Stability
                    if any(np.real(np.roots(K_Den)) > 0):
                        print('WARNING: The system representing ' + idDOF(i)+'->'+idDOF(j) + ' is UNSTABLE')


        # Save models
        self.sys    = sys
        self.sysDOF = sysDOF
        self.Khat   = Khat


    def ident_memory(self,om,Plot=False):

        ### Start of ident_retardation_FD.m
        # Scale Data for ID
        K = self.K_ij

        scale = max(abs(K))

        K_bar = K/scale
        F_resp  = K_bar / (1j * om)      # remove zero at s=0 from the data

        # Frequency response fitting
        ord_den = self.order
        ord_num = self.order - 2


        ### Start of fit_siso_resp
        # Using iterative method with frequency dependent weighting

        Weight = np.ones(np.shape(om))

        for iter in range(0,20):
            b,a = invfreqs(F_resp, om, ord_num, ord_den, wf=Weight)
            # a   = makeStable(a)    # DZ: I think this is potentially problematic, can lead to divergent solutions
            
            Weight = 1/abs(np.polyval(a,1j*om))**2

        ### End of fit_siso_resp

        _ , F_hat = freqs(b,a,worN=om)


        # rescale and incorporate zero
        b   = scale * np.concatenate([b,np.zeros(1)])

        _,K_hat     = freqs(b,a,worN=om)

        if Plot:
            

            fig = plt.figure()
            ax1 = fig.add_subplot(211)
            ax1.plot(om,abs(K),'o',om,abs(K_hat))
            # ax1.title(idDOF)
            

            ax2 = fig.add_subplot(212)
            ax2.plot(om,np.angle(K),'o',om,np.angle(K_hat))

            plt.show()

        return b, a, K_hat


    def visualizeFits(self,wamData):

        fig = []

        for q, P in enumerate(self.sys):

            normalIdx = np.array(self.sysDOF[q]) + 1
            sub = str(normalIdx[0]) + str(normalIdx[1])

            plt.figure()
            plt.plot(wamData.omega,np.real(wamData.K[self.sysDOF[q][0],self.sysDOF[q][1],:]),'o',label='K_'+sub)
            plt.plot(wamData.omega,np.real(self.Khat[q]),label= 'Khat_'+sub)

            plt.title(idDOF(self.sysDOF[q][0])+'->'+idDOF(self.sysDOF[q][1])+ ' Transfer Function')
            

            # plt.rc('text', usetex=True)
            plt.legend()

        plt.show()

        print('here')


def computeError(K,K_hat,om,A_inf):

    B_ij        = np.real(K)
    A_ij        = np.imag(K) / om + A_inf

    Bhat_ij     = np.real(K_hat)
    Ahat_ij     = np.imag(K_hat) / om + A_inf

    # using FDIRadMod.m notation to compute R^2 for both A and B
    sseb        = np.matmul(np.transpose((B_ij - Bhat_ij)),(B_ij - Bhat_ij))
    sstb        = np.matmul(np.transpose((B_ij - np.mean(B_ij))),(B_ij - np.mean(B_ij)))
    r2b         = 1- sseb/sstb

    ssea        = np.matmul(np.transpose((A_ij - Ahat_ij)),(A_ij - Ahat_ij))
    ssta        = np.matmul(np.transpose((A_ij - np.mean(A_ij))),(A_ij - np.mean(A_ij)))
    r2a         = 1- ssea/ssta

    return r2b, r2a

def makeStable(p):
    # check that polynomial has roots in the RHP, if not, flip
    r = np.roots(p)
    r = [-root if np.real(root) > 0 else root for root in r]

    p = np.poly(r)

    return p

                
def idDOF(dof):
    DOFs    = []
    DOFs.append('Surge')
    DOFs.append('Sway')
    DOFs.append('Heave')
    DOFs.append('Roll')
    DOFs.append('Pitch')
    DOFs.append('Yaw')

    return DOFs[dof]



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
    
