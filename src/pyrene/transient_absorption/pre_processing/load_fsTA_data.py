import numpy as np

def load_fsTA_file(file, error=False):
    # get the number of scans 
    with open(file, 'r') as f:
        for line in f: 
            if r'scans' in line.split():
                n_scans = int(line.split()[-1])

    # get shots
    nshots = np.loadtxt(file, usecols=[1], skiprows=20)[0]

    # get time axis 
    data = np.loadtxt(file, skiprows=20, usecols=[0])
    dt = np.split(data, n_scans)[0]

    # get data frame
    # there are 520 pixels/wavelengths
    # the data is organized like this:
    # delay, n_shots, TA@px1, err_TA@px1, TA@px2, err_TA@px2, ...
    
    # get column indices for TA data 
    if not error:
        indices = [2 + 2*n for n in range(520)]
    else:
        indices = [3 + 2*n for n in range(520)]
    
    # load correct coloumns
    df = np.loadtxt(file, skiprows=20, usecols=indices)
    # convert from µOD to mOD
    df = df/1000
    # convert RMS to error 
    if error:
        df = df/np.sqrt(nshots)
    # expected number of rows per scan
    n_delays = len(dt)

    # number of complete scans actually present
    n_complete_scans = len(df) // n_delays

    # discard incomplete tail
    df = df[:n_complete_scans * n_delays]

    # split array into n scans 
    df = np.vsplit(df, n_scans)
    return dt*10**(12), df