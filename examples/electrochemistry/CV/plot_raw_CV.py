from pyrene.electrochemistry import CV

c = CV(files=['data/raw_CV_data.csv', 'data/background_current.csv'], 
       labels=['raw data', 'background'], colors=['r', 'gray'])
c.show()

# for US convention #
c = CV(files=['data/raw_CV_data.csv', 'data/background_current.csv'], US=True,
       labels=['raw data', 'background'], colors=['r', 'gray'])
c.show()