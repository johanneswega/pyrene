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
                        head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  norm. absorbance'
                    elif self.c:
                        head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  extinction coeffcient / M-1 cm-1'
                    else:
                        head = r'wavenlength / nm,    wavenumber / 10^3 cm-1,  absorbance' 
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
                    head = r'wavenlength / nm,  wavenumber / 10^3 cm-1,  norm. intensity for wavelength, norm. intensity for wavenumber'
                else:
                    head = r'wavenlength / nm,  wavenumber / 10^3 cm-1,  intensity for wavelength / a.u., intensity for wavenumber / a.u.'
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
                # make folder structure
                import os
                if not os.path.exists('exported'):
                    os.makedirs('exported')
                if not os.path.exists('exported/' + self.files[i][:self.files[i].find('.')]):
                    os.makedirs('exported/' + self.files[i][:self.files[i].find('.')])
                if self.wavelength:
                    if not os.path.exists('exported/' + self.files[i][:self.files[i].find('.')] + '/kinetics'):
                        os.makedirs('exported/' + self.files[i][:self.files[i].find('.')] + '/kinetics')    
                if self.delay:
                    if not os.path.exists('exported/' + self.files[i][:self.files[i].find('.')] + '/spectra'):
                        os.makedirs('exported/' + self.files[i][:self.files[i].find('.')] + '/spectra') 
                if self.contour:
                    if not os.path.exists('exported/' + self.files[i][:self.files[i].find('.')] + '/contour'):
                        os.makedirs('exported/' + self.files[i][:self.files[i].find('.')] + '/contour') 

                if self.experiment=='femto':
                    head_dt = 'delay / ps'
                else:
                    head_dt = 'delay / ns'
                if self.norm[i] or self.devide[i]:
                    head_dA = r'norm. transient absorption'
                else:
                    head_dA = r'transient absorption / mOD'

                # export kinetics
                if self.wavelength:
                    if self.IR:
                        name = 'exported/' + self.files[i][:self.files[i].find('.')] + '/kinetics/' + '%s_cm-1.txt'%(self.wavelength[i])
                    else:
                        name = 'exported/' + self.files[i][:self.files[i].find('.')] + '/kinetics/' + '%s_nm.txt'%(self.wavelength[i])
                    np.savetxt(name, np.column_stack([self.x[i], self.y[i]]), header=head_dt + ' , ' + head_dA, delimiter=',')
                    print("exported kinetics data")

                # export spectra 
                if self.delay:
                    if self.IR or self.wn:
                        wn = self.x[i]
                        wl = 1e4/self.x[i]
                    else:
                        wl = self.x[i]
                        wn = 1e4/self.x[i]
                    if self.IR: 
                        head = r'wavenlength / µm,    wavenumber / cm-1, ' + head_dA
                    else:
                        head = r'wavenlength / nm,    wavenumber / 10^3 cm-1, ' + head_dA
                    dt = self.get_delay_labels(self.delay)[i]
                    dt = dt.replace(' ', '_')
                    name = 'exported/' + self.files[i][:self.files[i].find('.')] + '/spectra/' + dt + '.txt'
                    np.savetxt(name, np.column_stack([wl, wn, self.y[i]]), header=head, delimiter=',')
                    print("exported TA spectra")
                        
                # export 2D map / contour 
                if self.contour:
                    name = 'exported/' + self.files[i][:self.files[i].find('.')] + '/contour/'
                    if self.wn:
                        np.savetxt(name + 'wavelength.txt', 1e4/self.x[i], header='wavelength / nm')
                    else:
                        np.savetxt(name + 'wavelength.txt', 1e4/self.x[i], header='wavelength / nm')
                    np.savetxt(name + 'time.txt', self.y[i], header=head_dt)
                    np.savetxt(name + 'TA.txt', self.z[i], header=head_dA)
                    print("exported TA contour")
                    
                                                            

        