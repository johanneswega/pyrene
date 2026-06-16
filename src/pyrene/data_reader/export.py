import numpy as np
from dataclasses import dataclass

@dataclass
class DataExporter():
    """class to export data from experimental files as txt"""

    def export(self):
        """exports data"""

        # loop through files
        for i in range(len(self.files)):

            ### Absorption ###
            if self.IR:
                if self.norm[i]:
                    head = r'wavenlength / µm,    wavenumber / cm-1,  norm. absorbance'
                elif self.c:
                    head = r'wavenlength / µm,    wavenumber / cm-1,  extinction coeffcient / M-1 cm-1'
                else:
                    head = r'wavenlength / µm,    wavenumber / cm-1,  absorbance'
                wn = self.x[i]
                wl = 1e4/wn
                A = self.y[i]
            else:
                if self.norm[i]:
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

            print("exportet data")
            np.savetxt('%s.txt'%(self.files[i][:self.files[i].find('.')]), 
                        np.column_stack([wl, wn, A]),
                        header=head, delimiter=',')                                 

        