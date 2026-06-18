def main():
    import numpy as np 
    import scipy.constants as sc

    print('')
    E = input(' Input energy gap [nm, kK, cm-1, eV, Hz, kJ/mol, kcal/mol, kT, Eh]: ')

    E = E.split()
    value = float(E[0])
    unit = E[1]

    if unit=='nm':
        print('')
        value = value*10**-9
        print(' %.4g kK'%((1/value)*10**-2/1000))
        print(' %.4g cm-1'%((1/value)*10**-2))
        print(' %.4g eV'%(((sc.h*sc.c)/(value*sc.e))))
        print(' %.4g Hz'%(((sc.c)/(value))))
        print(' %.4g kJ/mol'%(((sc.h*sc.c*sc.Avogadro*10**-3)/(value))))
        print(' %.4g kcal/mol'%(((sc.h*sc.c*sc.Avogadro*10**-3)/(value*4.1868))))
        print(' %.4g kT'%(((sc.h*sc.c)/(value))/ (sc.Boltzmann * 298.15)))
        print(' %.4g Eh'%(((sc.h*sc.c)/(value*sc.e))*0.0367493))
        print('')

    if unit=='cm-1' or unit=='kk' or unit=='kK':
        print('')
        if unit=='cm-1':
            value = (1/value)*10**-2
        else: 
            value = (1/value)*10**-2 / 1000
        print(' %.4g nm'%((value*10**9)))
        print(' %.4g eV'%(((sc.h*sc.c)/(value*sc.e))))
        print(' %.4g Hz'%(((sc.c)/(value))))
        print(' %.4g kJ/mol'%(((sc.h*sc.c*sc.Avogadro*10**-3)/(value))))
        print(' %.4g kcal/mol'%(((sc.h*sc.c*sc.Avogadro*10**-3)/(value*4.1868))))
        print(' %.4g kT'%(((sc.h*sc.c)/(value))/ (sc.Boltzmann * 298.15)))
        print(' %.4g Eh'%(((sc.h*sc.c)/(value*sc.e))*0.0367493))
        print('')

    if unit=='eV':
        print('')
        value = (sc.h*sc.c)/(value*sc.e)
        print(' %.4g nm'%((value*10**9)))
        print(' %.4g kK'%((1/value)*10**-2/1000))
        print(' %.4g cm-1'%((1/value)*10**-2))
        print(' %.4g Hz'%(((sc.c)/(value))))
        print(' %.4g kJ/mol'%(((sc.h*sc.c*sc.Avogadro*10**-3)/(value))))
        print(' %.4g kcal/mol'%(((sc.h*sc.c*sc.Avogadro*10**-3)/(value*4.1868))))
        print(' %.4g kT'%(((sc.h*sc.c)/(value))/ (sc.Boltzmann * 298.15)))
        print(' %.4g Eh'%(((sc.h*sc.c)/(value*sc.e))*0.0367493))
        print('')

    if unit=='Hz':
        print('')
        value = sc.c/value
        print(' %.4g nm'%((value*10**9)))
        print(' %.4g kK'%((1/value)*10**-2/1000))
        print(' %.4g cm-1'%((1/value)*10**-2))
        print(' %.4g eV'%(((sc.h*sc.c)/(value*sc.e))))
        print(' %.4g kJ/mol'%(((sc.h*sc.c*sc.Avogadro*10**-3)/(value))))
        print(' %.4g kcal/mol'%(((sc.h*sc.c*sc.Avogadro*10**-3)/(value*4.1868))))
        print(' %.4g kT'%(((sc.h*sc.c)/(value))/ (sc.Boltzmann * 298.15)))
        print(' %.4g Eh'%(((sc.h*sc.c)/(value*sc.e))*0.0367493))
        print('')

    if unit=='kJ/mol':
        print('')
        value = (sc.c*sc.h*sc.Avogadro)/(value*10**3)
        print(' %.4g nm'%((value*10**9)))
        print(' %.4g kK'%((1/value)*10**-2/1000))
        print(' %.4g cm-1'%((1/value)*10**-2))
        print(' %.4g eV'%(((sc.h*sc.c)/(value*sc.e))))
        print(' %.4g Hz'%(((sc.c)/(value))))
        print(' %.4g kcal/mol'%(((sc.h*sc.c*sc.Avogadro*10**-3)/(value*4.1868))))
        print(' %.4g kT'%(((sc.h*sc.c)/(value))/ (sc.Boltzmann * 298.15)))
        print(' %.4g Eh'%(((sc.h*sc.c)/(value*sc.e))*0.0367493))
        print('')

    if unit=='kcal/mol':
        print('')
        value = (sc.c*sc.h*sc.Avogadro)/(value*10**3*4.1868)
        print(' %.4g nm'%((value*10**9)))
        print(' %.4g kK'%((1/value)*10**-2/1000))
        print(' %.4g cm-1'%((1/value)*10**-2))
        print(' %.4g eV'%(((sc.h*sc.c)/(value*sc.e))))
        print(' %.4g Hz'%(((sc.c)/(value))))
        print(' %.4g kJ/mol'%(((sc.h*sc.c*sc.Avogadro*10**-3)/(value))))
        print(' %.4g kT'%(((sc.h*sc.c)/(value))/ (sc.Boltzmann * 298.15)))
        print(' %.4g Eh'%(((sc.h*sc.c)/(value*sc.e))*0.0367493))
        print('')

    if unit=='kT':
        print('')
        value = (sc.h*sc.c)/(value*sc.Boltzmann*298.15)
        print(' %.3g nm'%((value*10**9)))
        print(' %.3g kK'%((1/value)*10**-2/1000))
        print(' %.3g cm-1'%((1/value)*10**-2))
        print(' %.3g eV'%(((sc.h*sc.c)/(value*sc.e))))
        print(' %.3g Hz'%(((sc.c)/(value))))
        print(' %.3g kJ/mol'%(((sc.h*sc.c*sc.Avogadro*10**-3)/(value))))
        print(' %.3g Eh'%(((sc.h*sc.c)/(value*sc.e))*0.0367493))
        print('')

    if unit=='Eh':
        print('')
        value = (sc.h*sc.c)/((value/0.0367493)*sc.e)
        print(' %.4g nm'%((value*10**9)))
        print(' %.4g kK'%((1/value)*10**-2/1000))
        print(' %.4g cm-1'%((1/value)*10**-2))
        print(' %.4g eV'%(((sc.h*sc.c)/(value*sc.e))))
        print(' %.4g Hz'%(((sc.c)/(value))))
        print(' %.4g kJ/mol'%(((sc.h*sc.c*sc.Avogadro*10**-3)/(value))))
        print(' %.4g kT'%(((sc.h*sc.c)/(value))/ (sc.Boltzmann * 298.15)))
        print('')

