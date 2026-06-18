import numpy as np
from dataclasses import dataclass

@dataclass
class DataExporter():
    """class to export data from experimental files as txt"""

    # define type of data to be exported 
    abs_spec : bool = None
    em_spec :  bool = None

    def export(self):
        """exports data"""

        # loop through files
        for i in range(len(self.files)):

            ### Absorption ###
            if self.abs_spec:
                if self.IR:
                    if self.norm[i] or self.devide[i]:
                        head = r'wavenlength / µm,    wavenumber / cm-1,  norm. absorbance'
                    elif self.c:
                        head = r'wavenlength / µm,    wavenumber / cm-1,  extinction coeffcient / M-1 cm-1'
                    else:
                        head = r'wavenlength / µm,    wavenumber / cm-1,  absorbance'
                    wn = self.x[i]
                    wl = 1e4/wn
                    A = self.y[i]
                else:
                    if self.norm[i] or self.devide[i]:
                        head = r'wavenlength / nm,    wavenumber / cm-1,  norm. absorbance'
                    elif self.c:
                        head = r'wavenlength / nm,    wavenumber / cm-1,  extinction coeffcient / M-1 cm-1'
                    else:
                        head = r'wavenlength / nm,    wavenumber / cm-1,  absorbance' 
                    if self.wn:
                        wn = self.x[i]
                        wl = 1e4/wn
                        A = self.y[i] 
                    else:
                        wl = self.x[i]
                        wn = 1e4/wl
                        A = self.y[i]  

                print("exported absorption data")
                np.savetxt('%s.txt'%(self.files[i][:self.files[i].find('.')]), 
                            np.column_stack([wl, wn, A]),
                            header=head, delimiter=',')

            ### emission ###
            if self.em_spec:
                if self.norm[i] or self.devide[i]:
                    head = r'wavenlength / nm,  wavenumber / cm-1,  norm. intensity for wavelength, norm. intensity for wavenumber'
                else:
                    head = r'wavenlength / nm,  wavenumber / cm-1,  intensity for wavelength / a.u., intensity for wavenumber / a.u.'
                if self.wn:
                    wn = self.x[i]
                    wl = 1e4/wn
                    I_wn = self.y[i]
                    I_wl = I_wn / wl**2
                else:
                    wl = self.x[i]
                    wn = 1e4/wn
                    I_wl = self.y[i]
                    I_wn = I_wl * wl**2  

                print("exported emission data")
                np.savetxt('%s.txt'%(self.files[i][:self.files[i].find('.')]), 
                            np.column_stack([wl, wn, I_wl, I_wn]),
                            header=head, delimiter=',') 

            ### transient absorption ###    
            if self.two_dim:
                if self.wavelength:
                    if self.norm[i] or self.devide[i]:
                        head_dA = r'norm. transient absorption'
                    else:
                        head_dA = r'norm. '                                        

        