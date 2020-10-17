import numpy as np
from Pygrad.Pygrad import Pygrad

#Data type to store results with uncertainties
# Initialized with val and % error
class PGRes:
    '''Class used to store the output of Pygrad with a val and error.'''
    err=0
    '''The error variable.'''
    val=0
    '''The value variable.'''
    def __init__(self, a=0,b=0):
        self.val=a
        self.err=b*a/100.
    def __str__(self):
        return "value: " + str(self.val)+ "; error: "+ str(self.err)

# PygradRun helper class
class PygradRun:
    '''Class used to be the wrapper object of Pygrad.'''
    #Default settings for running Pygrad
    PGSettings   ={'Gases'                  :['CF4','Ar'],
                   'Fractions'              :[90,10],
                   'EField_Vcm'             :100, 
                   'InitialElectronEnergy':5900.0,
                   'ThermalCut'             :2.0,
                   'EnergyCut'              :0.0,
                   'Temperature_C'          :23,
                   'Pressure_Torr'          :750.062,
                   'BField_Tesla'           :0,
                   'BField_angle'           :0,
                   'Enable_penning'         :False,
                   'Imip'                   :3,
                   'BeamDirection'          :-1,
                   'Seed'                   :0,
                   'OutputVerbosity'        :2,
                   'DetectorEfficiency'     :100.0,
                   'ExcitationWeight'       :0.5,
                   'kgas'                   :1,
                   'lgas'                   :1,
                   'lcmp'                   :True,
                   'lray'                   :True,
                   'lpap'                   :True,
                   'lbrm'                   :True,                   'iecasc'                 :True
                   }
    '''Dictionary used to store the inputs/settings for the Pygrad simulation.'''
    # Available Gases
    Gases = [np.nan, 'CF4', 'ARGON', 'HELIUM4', 'HELIUM3', 'NEON', 'KRYPTON', 'XENON', 'CH4', 'ETHANE', 'PROPANE'
         , 'ISOBUTANE', 'CO2', np.nan, 'H2O', 'OXYGEN', 'NITROGEN', np.nan, np.nan, np.nan, np.nan
         , 'HYDROGEN', 'DEUTERIUM', np.nan, np.nan, 'DME']
    '''Array of gases in Pygrad.'''

    # Print list of available gases
    def ListGases(self):
        '''Function used to print all the gases names in Pygrad.'''
        for g in self.Gases:
            if(type(g)==str):
                print(g,self.GasCode(g))

    # Convert GasName into Degrad GasCode            
    def GasCode(self,GasName):
        '''Function used to get the ID of the gas. The ID is simply the index of that gas in that array.'''
        return self.Gases.index(GasName)

    # Convert Degrad GasCode in GasName
    def GasName(self,Code):
        '''Function used to return the name of the Gas ID given.'''
        return Gases[Code]
    
    # Load Input Dictionary into Degrad object
    def ProcessInputs(self,DGObject, Inputs):
        '''Function used to setup the Pygrad Object with the given inputs in the PGSettings dictionary.'''
        for key in self.PGSettings.keys():
            if not key in Inputs:
                print("Input "+str(key)+ " not set, using default "+str(self.PGSettings[key]))
                Inputs[key]=self.PGSettings[key]
        if(len(Inputs['Gases'])!=len(Inputs['Fractions'])):
            print("Error! Gas and fraction lists not the same length")
            return False
        if(len(Inputs['Gases'])>6):
            print("Error! Too many gases. Max is 6.")
            return False
        if(abs(sum(Inputs['Fractions'])-100)>1e-6):
            print("Error! Gas fractions don't add to 100%")
            return False
        DGObject.EField=Inputs['EField_Vcm']
        DGObject.NumberOfGases=len(Inputs['Gases'])
        GasIDs=np.zeros(6,dtype='int')
        GasFractions=np.zeros(6,dtype='float')
        for i in range(DGObject.NumberOfGases):
            GasIDs[i] = self.GasCode(Inputs['Gases'][i])
            GasFractions[i]  = Inputs['Fractions'][i]
        if Inputs['Seed']:
            DGObject.Random_Seed = Inputs['Seed']
        DGObject.GasIDs  = GasIDs
        DGObject.GasFractions   = GasFractions
        DGObject.Enable_Penning   = Inputs['Enable_penning']
        DGObject.TemperatureCentigrade  = Inputs['Temperature_C']
        DGObject.Pressure_Torr   = Inputs['Pressure_Torr']
        DGObject.BField_Mag   = Inputs['BField_Tesla']
        DGObject.BField_Angle = Inputs['BField_angle']
        DGObject.ThermalCut = Inputs['ThermalCut']
        DGObject.EnergyCut = Inputs['EnergyCut']
        DGObject.nDelta = Inputs['nDelta']
        DGObject.imip = Inputs['Imip']
        DGObject.BeamDirection = Inputs['BeamDirection']
        DGObject.OutputVerbosity = Inputs['OutputVerbosity']
        DGObject.DetectorEfficiency = Inputs['DetectorEfficiency']
        DGObject.ExcitationWeight = Inputs['ExcitationWeight']
        DGObject.kgas = Inputs['kgas']
        DGObject.lgas = Inputs['lgas']
        DGObject.lcmp = Inputs['lcmp']
        DGObject.lray = Inputs['lray']
        DGObject.lpap = Inputs['lpap']
        DGObject.lbrm = Inputs['lbrm']
        DGObject.IECascade = Inputs['iecasc']
        return True

    # Extract Outputs into Output Dictionary
    def ProcessOutputs(self, DGObject):
        '''Function used to fill the Outputs dictionary with the results from the Pygrad object.'''
        Outputs={}
        return Outputs

    # Run Pygrad with chosen settings
    def Run(self,MySettings):
        '''Function used to run the Pygrad simulation. Note that the PGSettings dictionary needs to be set up.'''
        DGObject = Pygrad()
        Status=self.ProcessInputs(DGObject,MySettings)
        if(Status):
            DGObject.Start()
        Outputs = self.ProcessOutputs(DGObject)
        return Outputs
