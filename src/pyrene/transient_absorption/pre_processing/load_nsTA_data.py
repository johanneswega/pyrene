import numpy as np

def load_nsTA_file(file, error=False):
    fh = open(file)
    # initialize array for data
    data = []
    a = []
    b = []
    # list for nsamples in each delay bin
    nsamples = []
    for line in fh:
        # split line
        line = line.split()
        # if line not empty, i.e. not blank
        if len(line)!=0:
            # reject header
            if line[0]!='%':
                    # append time delay
                    a.append(float(line[0]))
                    # append TA signal
                    if not error:
                        a.append(float(line[2]))
                    else:
                        a.append(float(line[3]))
                    # append samples in bin
                    b.append(float(line[-1]))
        # if line blank --> next delay
        if len(line)==0:
            # add a-list to data 
            data.append(a)
            # add n samples to b list
            nsamples.append(b[0])
            # reset place holder lists
            a = []
            b = []
    fh.close()

    # the data array is now organized as:
    # [delay0, TA@px0, delay0, TA@px1, delay0, ...,delay0, TA@px524]
    # [delay1, TA@px0, delay1, TA@px1, delay1, ...,delay1, TA@px524]
    # and so on until the last delay

    # make empty list for delays
    t = []
    for i in range(len(data)):
        if len(data[i])!=0:
            t.append(data[i][0]*10**9)

    # number of pixel of CCD
    npx = 524
    # initialize empty dA-array
    dA = np.zeros((len(t), npx))
    # add values into dA-array
    for i in range(len(t)):
        spec = np.array(data[i])
        # only take the TA@px values 
        # i.e. reject every second entry or where data==delay_i
        spec = spec[spec!=data[i][0]]
        dA[i,:] = spec

    return np.array(t), dA, np.array(nsamples)