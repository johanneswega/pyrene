def main():
    import numpy as np
    import os

    # Read eps and n data
    file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "epsilon.txt")
    data = np.loadtxt(file, dtype=str, delimiter='\t')

    # 1 - Solvent name
    # 2 - Dipole moment in D
    # 4 - Diectric constant at 25°C
    # 6 - -1000*dln(eps)/dT in K-1
    # 8 - Refractive index at 25°C
    # 9 - -1000 dn/dT

    # Read viscosity data
    file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "viscosity.txt")
    visd = np.loadtxt(file, dtype=str, delimiter='\t')

    # 1 - Solvent name
    # 4 - Viscosity in cP
    # 6 - -100*dln(eta)/dT in K-1

    # ideal gas constant
    R = 8.31446261815324

    # input solvent 
    print("")
    solv = input('Solvent : ')

    if len(solv)==1:
        print("")
        for i in range(len(data[:,0])):
            solv_name = data[i,1]
            if solv_name[0]==solv:
                print(solv_name)
        for i in range(len(data[:,0])):
            solv_name = data[i,1]
            if '-' in solv_name:
                index = solv_name.find('-')
                if solv_name[index+1]==solv:
                    print(solv_name)
            if '"'==solv_name[0]:
                if solv_name[1]==solv:
                    print(solv_name)                
        print("")
    elif 'search' in solv:
        print("")
        if 'eps' in solv:
            substring = 'eps '
            index = solv.index(substring) + len(substring)
            eps = float(solv[index:])
            # go through all solvents and print solvents within +- 10 % of eps
            for i in range(len(data[:,0])):
                if data[i,4]=='':
                    continue
                if 0.9*eps <= float(data[i,4]) <= 1.1*eps:
                    print('%s,  eps = %.4g'%(data[i,1], float(data[i,4])))
        if 'eta' in solv:
            substring = 'eta '
            index = solv.index(substring) + len(substring)
            eta = float(solv[index:])
            # go through all solvents and print solvents within +- 10 % of eps
            for i in range(len(visd[:,0])):
                if visd[i,4]=='':
                    continue
                if 0.9*eta <= float(visd[i,4]) <= 1.1*eta:
                    print('%s,  eta = %.4g'%(visd[i,1], float(visd[i,4])))
        if 'df' in solv:
            substring = 'df '
            index = solv.index(substring) + len(substring)
            df = float(solv[index:])        
            # go through all solvents and print solvents within +- 0.05 of df
            for i in range(len(data[:,0])):
                if data[i,4]=='':
                    continue
                if data[i,8]=='':
                    continue
                else:
                    eps = float(data[i,4])
                    n =float(data[i,8])
                    ons = (2*(eps - 1)/(2*eps+1)) - (2*(n**2 - 1)/(2*n**2+1))
                    if df-0.05 <= ons <= df+0.05:
                        print('%s,  df = %.3g'%(data[i,1], ons))
    else:
        print("")
        k = 0 
        if not '-T' in solv:
            # do all calculations at 25°C
            T = 298.15
            print('T = %.5g K'%T)
            for i in range(len(data[:,0])):
                if solv==data[i,1]:
                    # get dipole moment
                    mu = float(data[i,2])
                    print("mu = %.4g D"%mu)
                    # get dielectric constant
                    eps = float(data[i,4])
                    print("eps = %.4g"%eps)
                    # get refractive index
                    n = float(data[i,8])
                    print("n = %.6g"%n)
                    # onsager function
                    df = (2*(eps - 1)/(2*eps+1)) - (2*(n**2 - 1)/(2*n**2+1))
                    print("df = %.4g"%df)
                    break
                k+=1
            for i in range(len(visd[:,0])):
                if solv==visd[i,1]:
                    eta = float(visd[i,4])
                    print("eta = %.4g cP"%eta)
                    kd = (8*R*1000*T)/(3*eta/1000)
                    print("k_diff = %.3g M-1 s-1"%kd)
        else: 
            # do custom calculations at a different temperature
            T = float(solv[solv.find('T')+1:]) + 273.15
            solv = solv[:solv.find('T')-2]
            print('T = %.5g K'%T)
            for i in range(len(data[:,0])):
                if solv==data[i,1]:
                    # get dielectric constant at 25°C
                    eps25 = float(data[i,4])
                    # get dipole moment
                    mu = float(data[i,2])
                    print("mu = %.4g D"%mu)
                    # get conversion factor
                    if data[i,6]=='':
                        print('no conv. factor found, using data at 25°C')
                        d = 0
                    else:
                        d = -1*float(data[i,6])/1000
                    # convert to different temperature
                    eps = eps25*np.exp(d*(T - 298.15))
                    print("eps = %.4g"%eps)

                    # get refractive index at 25°C
                    n25 = float(data[i,8])
                    # get conversion factor
                    if data[i,9]=='':
                        print('no conv. factor found, using data at 25°C')
                        g = 0
                    else:              
                        g = -float(data[i,9])/1000
                    # convert to different temperature
                    n = n25 + g*(T - 298.15)   
                    print("n = %.6g"%n) 
                    df = (2*(eps - 1)/(2*eps+1)) - (2*(n**2 - 1)/(2*n**2+1))
                    print("df = %.4g"%df)              
                    break
                k+=1
            for i in range(len(visd[:,0])):
                if solv==visd[i,1]:
                    # get viscosity at 25°C
                    eta25 = float(visd[i,4])
                    # get conversion factor
                    if visd[i,4]=='':
                        print('no conv. factor found, using data at 25°C')
                        h = 0
                    else:               
                        h = -1*float(visd[i,4])/100
                    # convert to different temperature
                    eta = eta25*np.exp(h*(T - 298.15))
                    print("eta = %.4g cP"%eta)
                    kd = (8*R*1000*T)/(3*eta/1000)
                    print("k_diff = %.3g M-1 s-1"%kd)        

        if k==len(data[:,0]):
            print(""" '%s' not found in database! Perhaps check the spelling.
                if you'd like to list all solvents beginning with a certain letter or character, rerun the program and type e.g. 'a' and all solvents
                starting with a in the database are listed."""%solv)
            
    print("")