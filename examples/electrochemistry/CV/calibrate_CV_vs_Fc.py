from pyrene.electrochemistry import CV

c = CV(files=['data/raw_CV_data.csv', 'data/Fc_calibration.csv'], 
       labels=['raw data', 'Fc measurement'], colors=['r', 'orange'])
# determine E1/2 of Fc+/Fc
c.get_half_wave_potential(file_index=1, E_range=[0.0, 0.5], plot=True)
c.show()

# use it calibrate the CV measurement
c = CV(files=['data/raw_CV_data.csv'], Fc=[0.2149],
       labels=['calibrated data'], colors=['r'])
c.show()

# you can also plot it against SCE then
c = CV(files=['data/raw_CV_data.csv'], Fc=[0.2149], SCE=True, figsize=[5, 3.5], savefig='halfwave.png',
       labels=['calibrated data'], colors=['r'])
c.get_half_wave_potential(file_index=0, E_range=[0.0, 0.5], plot=True)
c.show()