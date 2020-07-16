'''
Testing Script for comparing matlab and python SS_Fit tools
'''

import numpy as np

class SS_FitModel(object):
    '''
    The full block diagonal SS model input to fast
    '''

    def __init__(self):
        self.A = []
        self.B = []
        self.C = []

        self.DoFs = []
        self.States = 0
        self.StatesPerDoF = np.zeros(6)


    def readInput(self):
        with open(self.ss_file,'r') as f:
            # Header and Options
            f.readline()
            self.DoFs     = np.fromstring(f.readline().split('%')[0].strip()[1:-1],dtype=int,sep = ' ')
            self.States          = int(f.readline().split('%')[0])
            self.StatesPerDoF     = np.fromstring(f.readline().split('%')[0].strip()[1:-1],dtype=int,sep = ' ')


            all_data = [np.array(x.split()) for x in f.readlines()]
            self.A = np.squeeze([a_dat.astype(np.float)  for a_dat in all_data[:self.States]])
            self.B = np.squeeze([a_dat.astype(np.float)  for a_dat in all_data[self.States:2*self.States]])
            self.C = np.squeeze([a_dat.astype(np.float)  for a_dat in all_data[2*self.States:]])


                
if __name__ == '__main__':

    # python .ss
    ss_py = SS_FitModel()


    ss_py.ss_file = '/Users/dzalkind/Tools/SS_Fitting/test_cases/Haliade/pelastar_py.ss'
    # ml_file = '/Users/dzalkind/Tools/SS_Fitting/test_cases/Haliade/pelastar.ss'

    ss_py.readInput()
