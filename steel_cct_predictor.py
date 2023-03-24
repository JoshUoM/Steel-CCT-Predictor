# import packages
import numpy as np
import math
from scipy import integrate
from scipy.integrate import quad
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import warnings
%matplotlib inline
warnings.filterwarnings("ignore", category=RuntimeWarning) 
warnings.filterwarnings("ignore", category=UserWarning) 

# define sub-functions to be used in primary functions
def Ae3_Calc(comp):
    required_elm = ['C','Si','Mn','Ni','Cr']
    for e in required_elm:
        if e not in comp:
            comp[e] = 0
    return int(((1570-(323*comp['C'])-(25*comp['Mn'])+(80*comp['Si'])-(3*comp['Cr'])-(32*comp['Ni']))-32)*(5/9))

def Ae1_Calc(comp):
    required_elm = ['Si','Mn','Ni','Cr']
    for e in required_elm:
        if e not in comp:
            comp[e] = 0
    return int(((1333-(25*comp['Mn'])+(40*comp['Si'])+(42*comp['Cr'])-(26*comp['Ni']))-32)*(5/9))

def Acm_Calc(comp):
    return int(224.4 + (992.4*comp['C']) - (465.1*(comp['C']**2)) + (46.7*comp['Cr']) + (19*comp['C']*comp['Cr']) - (6.1*(comp['Cr']**2)) + (7.6*comp['Mn']) + (10*comp['Mo']) - (6.8*comp['Cr']*comp['Mo']) - (6.9*comp['Ni']) + (3.7*comp['C']*comp['Ni']) - (2.7*comp['Cr']*comp['Ni']) + (0.8*(comp['Ni']**2)) + (16.7*comp['Si']))

def Ae3C_Calc(comp,T):
    required_elm = ['Si','Mn','Ni','Cr']
    for e in required_elm:
        if e not in comp:
            comp[e] = 0
    return (((T*(9/5))+32)-1570+(25*comp['Mn'])-(80*comp['Si'])+(3*comp['Cr'])+(32*comp['Ni']))/(-323)

def densityf(T,comp,sf):
    return 7875.96-(0.297*T)-((5.62*(10**(-5)))*(T**2))+(sf)*(-206.35+(0.00778*T)+((1.472*(10**(-6)))*(T**2))) - (36.86*comp['Si']) - (7.24*comp['Mn'])

def densitya(T,comp):
    return 8099.79 - (0.5060*T)+((-118.26+(0.00739*T))*(comp['C'])) - (68.24*comp['Si']) - (6.01*comp['Mn'])

def maxf(comp,T,sf):
    Wf = (Ae3C_Calc(comp,T) - comp['C'])/(Ae3C_Calc(comp,T) - sf)
    return (Wf/densityf(T,comp,sf))/( (Wf/densityf(T,comp,sf)) + ((1-Wf)/densitya(T,comp)) )

def Gl(l):
    return -3.2877-(6.6439*math.log10(l*10**(-3)))

def Gd(d):
    l = ((math.pi/4)*(d**2))**(1/2)
    return -3.2877-(6.6439*math.log10(l*10**(-3)))

def mass2mole(mass):
    molweight = {'Fe':55.85,'C':12.01,'Si':28.09,'Mn':54.94,'Ni':58.69,'Cr':52,'Mo':95.94,'W':183.85,'Co':58.93,'V':50.94,'Nb':92.91,'Cu':63.55,'Al':26.98,'Ti':47.88,'O':16,'N':14.01,'B':10.81,'P':30.97,'S':32.06,'As':74.92}
    if 'Fe' in mass.keys():
        pass
    else:
        mass['Fe'] = 100 - sum(mass.values())
    a = []
    for elm in mass.keys():
        mole = mass[elm]/molweight[elm]
        a.append(mole)
    tot_moles = sum(a)
    moles = {}
    n = 0
    for mol in a:
        moles[list(mass.keys())[n]] = mol/tot_moles
        n += 1
    return moles

def Bs_Calc(comp):
    return int(637-(58*comp['C'])-(35*comp['Mn'])-(15*comp['Ni'])-(34*comp['Cr'])-(41*comp['Mo']))

def Ms_Calc(comp):
    return int(539-(423*comp['C'])-(30.4*comp['Mn'])-(17.7*comp['Ni'])-(12.1*comp['Cr'])-(7.5*comp['Mo'])-(7.5*comp['Si'])+(10*comp['Co']))

def FC_Calc(comp):
    return np.exp(1+(6.31*comp['C'])+(1.78*comp['Mn'])+(0.31*comp['Si'])+(1.12*comp['Ni'])+(2.7*comp['Cr'])+(4.06*comp['Mo']))

def PC_Calc(comp):
    return np.exp(-4.25+(4.12*comp['C'])+(4.36*comp['Mn'])+(0.44*comp['Si'])+(1.71*comp['Ni'])+(3.33*comp['Cr'])+(5.19*np.sqrt(comp['Mo'])))

def BC_Calc(comp):
    return np.exp(-10.23+(10.18*comp['C'])+(0.85*comp['Mn'])+(0.55*comp['Ni'])+(0.9*comp['Cr'])+(0.36*comp['Mo']))

def T0_Calc(comp):
    dTM = {'Si':-3,'Mn':-37.5,'Ni':-6,'Mo':-26,'Cr':-19,'V':-44,'Co':19.5,'Al':8,'Cu':4.5}
    dTNM = {'Si':0,'Mn':-39.5,'Ni':-18,'Mo':-17,'Cr':-18,'V':-32,'Co':16,'Al':15,'Cu':-11.5}
    To = 970 - (80*mass2mole(comp)['C']*100) - 273.15
    for e in mass2mole(comp).keys():
        if e in dTM.keys():
            dTsub = mass2mole(comp)[e]*100*((7*dTNM[e])+(-1*dTM[e]))/(7-1)
            To = To - dTsub
            return To 
        
def t_theta(comp,T):
    c_moles = mass2mole(comp)['C']
    E = 0.01
    return (((math.log(1-E))/(-4.07*(10**4)*(c_moles**0.635)*math.exp(-33598/(8.314*(T+273.15)))))**(1/0.62))*3600

def DC_Lee(C,comp,T):
    R = 8.314*10**(-3)
    sum1, sum2 = 0, 0
    k1 = {'Mn':-0.0315,'Si':0.0509,'Ni':-0.0085,'Cr':0,'Mo':0.3031,'Al':-0.0520}
    k2 = {'Mn':-4.3663,'Si':4.0507,'Ni':-1.2407,'Cr':7.7260,'Mo':12.1266,'Al':-6.7886}
    for e in k1.keys():
        sum1 = sum1 + k1[e]*comp[e]
        sum2 = sum2 + k2[e]*comp[e]
    D = (0.146-(0.036*C*(1-(1.075*comp['Cr'])))+(sum1))*math.exp(-(144.3-(15*C)+(0.37*(C**2))+sum2)/(R*(T+273.15)))
    return D*10**(-4)

def t_diff(comp,T):
    w = 0.2*10**(-6)
    Dc = []
    XS = np.linspace(comp['C'],Ae3C_Calc(comp,T),20)
    for x in XS:
        Dc.append(DC_Lee(x,comp,T)/(Ae3C_Calc(comp,T)-comp['C']))
    D = integrate.trapz(Dc,XS)
    return ((w**2)*math.pi*(comp['C']-0.03)**2)/(16*D*((Ae3C_Calc(comp,T)-comp['C'])**2))

def roundup(x):
    return x if x % 1000 == 0 else x + 1000 - x % 1000

def SX_eq(x):
    return (1/(x**(0.4*(1-x))*(1-x)**(0.4*x)))

# primary functions to be run by users

def CCT_Plotter(Ts,comp,rates):
    colors = {'f':'steelblue','p':'mediumorchid','b':'chocolate','bu':'sandybrown','bl':'chocolate','m':'mediumseagreen'}
    
    plt.figure(figsize=(12,8))

    Ae3i, Ae1i = Ae3_Calc(comp), Ae1_Calc(comp)
    
    maxX = roundup(Ae3_Calc(comp)/min(rates))

    plt.plot([0.1,maxX],[Ae3i,Ae3i], color = 'k', linestyle = '--', linewidth = 1.5)
    plt.text(maxX*0.1, Ae3i-10, '$Ae_3$='+str(Ae3i)+'\N{DEGREE SIGN}C', bbox={'facecolor': 'white'},fontsize=14)  
    plt.plot([0.1,maxX],[Ae1i,Ae1i], color = 'dimgray', linestyle = '--', linewidth = 1.5)
    plt.text(maxX*0.2, Ae1i-10, '$Ae_1$='+str(Ae1i)+'\N{DEGREE SIGN}C', bbox={'facecolor': 'white'},fontsize=14)
    
    dT = 1
    temp = np.linspace(Ae3i,0,(Ae3i+1)*int(dT**(-1)))
    
    for r in rates:
        times = []
        for T in temp[1:]:
            times.append((Ae3i-T)/r)
        if r in [0.01,0.1,1,10,100]:
            plt.plot(times,temp[1:],color='red',linewidth=0.75)
            plt.text((Ae3i-250)/r, 25, str(r)+'\N{DEGREE SIGN}C/s', bbox={'facecolor': 'white'},fontsize=12)  
        else:
            plt.plot(times,temp[1:],color='red',linewidth=0.25)

        for phase in ['f','p','bu','bl','m']:
            for T in Ts[r][phase]:
                plt.scatter((Ae3i-T)/r,T,color=colors[phase],marker='o')

    for phase in ['f','p','b','m']:
        ti,Ti = [],[]
        for r in rates:
            try:
                ti.append((Ae3i-Ts[r][phase][0])/r)
                Ti.append(Ts[r][phase][0])
            except IndexError:
                ti.append(np.nan)
                Ti.append(np.nan)
        plt.plot(ti,Ti,color=colors[phase])#,marker='o')

    handle = []
    for phase in ['f','p','bu','bl','m']:
        handle.append(mpatches.Patch(color=colors[phase],label=phase))
    plt.legend(handles=handle,loc='lower left',
              fancybox=True, shadow=True, ncol=1, prop={'size': 16})

    plt.xlabel('Time (s)', fontsize = 18)
    plt.ylabel('Temperature (\N{DEGREE SIGN}C)', fontsize=18)
    plt.xscale('log')
    plt.xlim(1,maxX)
    plt.ylim(0,Ae3i+25)
    
    plt.tick_params(axis='x', labelsize=16)
    plt.tick_params(axis='y', labelsize=16)
    return

def CCT_Fractions(Ts,rates):
    df = pd.DataFrame({'Rate':rates})
    for phase in ['f','p','bu','bl','m']:
        X = []
        for r in rates:
            X.append(len(Ts[r][phase])/100)
        df['X'+phase] = X
    X = []
    for n in range(len(df)):
        X.append(abs(round(0.99-sum(df.iloc[n][1:]),2)))
    df['Xa'] = X
    return df

# CCT calulator function

def Steel_CCT_Calculator(comp,G,rates):
    
    elements = ['C', 'Si', 'Mn', 'Ni', 'Cr', 'Mo', 'V', 'Cu', 'Al', 'Ti', 'Co', 'W', 'P', 'As']
    for e in elements:
        if e in comp.keys():
            pass
        else:
            comp[e] = 0
    
    SX = {}
    for X in np.linspace(0.01,0.99,99):
        SX[round(X,2)] = quad(SX_eq, 0, round(X,2))[0]

    Ae3i, Ae1i, PC = Ae3_Calc(comp), Ae1_Calc(comp), PC_Calc(comp)
    Ts, dT, Tfinish = {}, 1, 0

    for r in rates:
        
        latest_trans = Ae3i

        compi = comp.copy()

        Ts[r] = {}
        for phase in ['f','p','b','bu','bl','m','a']:
            Ts[r][phase]=[]

        Xa,Xf,Xp,Xb,Xbu,Xbl,Xm = 0.99,0.01,0.01,0.01,0.01,0.01,0.01
        XaF,XaU,XaL = 0.99,0.99,0.99
        dt,Cbu,Cbl,sf,sbu,sbl = dT/r,compi['C'],compi['C'],0.03,0.03,0.27

        while Xa > 0:
            Xa_current = Xa

            temp = np.linspace(Ae3i,0,(Ae3i+1)*int(dT**(-1)))

            if Xf == 0.01:
                Ae3, FC = Ae3_Calc(compi), FC_Calc(compi)
            if Xb == 0.01:
                Bs, BC = Bs_Calc(compi), BC_Calc(compi) 
            Ms = Ms_Calc(compi)
            To = T0_Calc(compi)

            rf, rp, rb = [], [], []

            for T in temp:

                if Xa != Xa_current:
                    break

                # FERRITE:
                if T <= Ae3 and Xp == 0.01 and Xb == 0.01 and Xm == 0.01 and Xf <= maxf(comp,T,sf) and Xa > 0:

                    tF = (FC/((2**(0.41*G)*((Ae3-T)**3)*np.exp(-27500/(1.987*(T+273.15))))))*SX[round(Xf,2)]
                    rf.append(dt/tF)

                    if sum(rf) >= 1.00:

                        Ts[r]['f'].append(round(T,1))
                        latest_trans = T

                        compi['C'] = (comp['C']+((Xf)*(comp['C']-sf))/((1-Xf)*XaF))
                        Cbu, Cbl = compi['C'], compi['C']
                        
                        Xf += 0.01
                        Xa -= 0.01
                        XaU, XaL = Xa, Xa

                # PEARLITE:
                if T <= Ae1i and Xb == 0.01 and Xm == 0.01 and Xa > 0:

                    tP = (PC/((2**(0.32*G)*((Ae1i-T)**3)*np.exp(-27500/(1.987*(T+273.15))))))*SX[round(Xp,2)]
                    rp.append(dt/tP)

                    if sum(rp) >= 1.00 and T <= latest_trans and T <= min(Acm_Calc(compi),Ae3_Calc(compi)):

                        Ts[r]['p'].append(round(T,1))
                        latest_trans = T

                        Xp += 0.01
                        Xa -= 0.01

                # BAINITE:
                if T <= Bs and T < To and Xm == 0.01 and Xa > 0: 

                    tB = (BC/((2**(0.29*G)*((Bs-T)**2)*np.exp(-27500/(1.987*(T+273.15))))))*SX[round(Xb,2)]
                    rb.append(dt/tB)

                    if sum(rb) >= 1.00 and T <= latest_trans:

                        Ts[r]['b'].append(round(T,1))
                        latest_trans = T
                         
                        if t_theta(compi,T) >= t_diff(compi,T):
                            Xa0 = 1
                            C_aust_u = (Cbu+(Xbu*(Cbu-sbu)/((1-Xbu)*XaU)))
                            Cbl = C_aust_u
                            XaL = Xa - 0.01
                            
                            compi['C'] = Cbl
                            Ts[r]['bu'].append(round(T,1))
                            Xbu += 0.01
                        
                        elif t_theta(compi,T) < t_diff(compi,T):
                            C_aust_l = (Cbl+(Xbl*(Cbl-sbl)/((1-Xbl)*XaL)))
                            compi['C'] = C_aust_l
                            Ts[r]['bl'].append(round(T,1))
                            Xbl += 0.01

                        Xb += 0.01
                        Xa -= 0.01
                        
                # MARTENSITE:   
                if T <= Ms:

                    if Xm == 0.01:

                        Xa0 = Xa
                        dTm = 215*Xa0

                        Ts[r]['m'].append(round(T,1))

                        Xm += 0.01
                        Xa -= 0.01

                    elif Xm != 0.01:

                        k = -np.log(0.01)/(dTm)
                        TM = Ms+(1/k)*np.log(1-(Xm/(Xa0+0.01)))

                        if round(T,0) == round(TM,0):

                            Ts[r]['m'].append(round(T,1))

                            Xm += 0.01
                            Xa -= 0.01

                if T <= Tfinish and Xa > 0:
                    Xa = 0
                    Ts[r]['a'].append(round(T,1))
                    break
    return Ts
