#Importing the necessities
import numpy as np
import matplotlib.pyplot as plt

"""
Masses of particles in units of kg. The "Mev" is used to convert energies
to/from MeV and Joules (J). We include them outside our function due to
calculations done on the outside.
"""
Mev = 1.602176634*10**(-19)*10**6 #J
mu = 1.6605*10**(-27) #kg. This is used to convert the atomic weight into units
#of kg
#The following masses are for the particles used in the calculations.
mH = 1.00782503207*mu #kg
mDe = 2.01410177785*mu #kg
m3He = 3.01602931914*mu #kg
m4He = 4.00260325415*mu #kg
m7Be = 7.016929828*mu #kg
m7Li = 7.016004548*mu #kg
m8Be = 8.005305103*mu #kg
m8Bor = 8.024607233*mu #kg
m12C = 12*mu #kg
m13N = 13.005738609*mu #kg
m13C = 13.00335483778*mu #kg
m14N = 14.00307400478*mu #kg
m15O = 15.99491461956*mu #kg
m15N = 15.00010889823*mu #kg

"""
Constants used. The rho_core and temp_core are values for our Sun.
T9 is the core temperature in units of 10^9K
"""
rho_core = 1.62*10**5 #kg m**-3
temp_core = 1.57*10**7 #K
T9 = temp_core /10**9

"""
This is our function for calculating the energy production. Its variables are T9
which is the temperature in the units of K/10^9. rho_core is the core density
and sanity_check is where we decide if we want to run it or not.
To run the sanity check, put "run" where it says sanity_check, or else just put
"no" or "stop" etc.
"""
def energy_production(T9, rho_core, sanity_check):
    """
    Constant used in the calculations. c is the speed of light, u is the atomic
    unit, Mev is used to convert energies from MeV to Joule. NA ia Avogrados
    number.
    """
    c = 299792458 #m/s
    u = 1.66053904*10**(-27) #kg
    Mev = 1.602176634*10**(-19)*10**6 #J
    NA = 6.0221*10**(23) #1/mol

    """
    T9_str1 and T9_str2 are the two T9* in Table 3.1 in the Lecture Note.
    """
    T9_str1 = T9/(1+4.95*10**(-2)*T9) #This
    T9_str2 = T9/(1+0.759*T9)

    #These are our mass fractions
    X = 0.7
    Y_3He = 10**(-10)
    Y = 0.29
    Z_7Li = 10**(-7)
    Z_7Be = 10**(-7)
    Z_14N = 10**(-11)

    #These are the masses used in the energy production. They are the same as
    #the ones on the outside
    mu =  1.6605*10**(-27) #kg
    mH = 1.00782503207*mu#kg
    mDe = 2.01410177785*mu #kg
    m3He = 3.01602931914*mu #kg
    m4He = 4.00260325415*mu #kg
    m7Be = 7.016929828*mu #kg
    m7Li = 7.016004548*mu #kg
    m8Be = 8.005305103*mu #kg
    m8Bor = 8.024607233*mu #kg
    m12C = 12*mu#kg
    m13N = 13.005738609*mu#kg
    m13C = 13.00335483778*mu#kg
    m14N = 14.00307400478*mu#kg
    m15O = 15.003065617*mu#kg
    m15N = 15.00010889823*mu#kg

    """
    Here are the energy output of all the reactions in PPI, PPII, PPIII and CNO
    """
    Q1H_1H = (2*mH-mDe)*c**2 - 0.265*Mev #J
    QDe_1H = (mDe + mH - m3He)*c**2
    #PP1
    Q3He_3He = (2*m3He - m4He - 2*mH)*c**2
    #PP2
    Q3He_4He = (m3He + m4He - m7Be)*c**2
    Q7Be_e = (m7Be - m7Li)*c**2 - 0.81*Mev
    Q7Li_1H = (m7Li + mH - 2*m4He)*c**2
    #PP3
    Q3He_4He = (m3He + m4He - m7Be)*c**2
    Q7Be_1H = (m7Be + mH - m8Bor)*c**2
    Q8Bor = (m8Bor - m8Be)*c**2 - 6.711*Mev
    Q8Be = (m8Be - 2*m4He)*c**2
    #This is our combination of decay energy production that we use instead of
    #each of them seperately
    Q_decay = Q8Bor + Q8Be
    #CNO
    Q12C_1H = (m12C+mH-m13N)*c**2
    Q13N = (m13N - m13C)*c**2 - 0.707*Mev
    Q13C_1H = (m13C + mH - m14N)*c**2
    Q14N_1H = (m14N + mH - m15O)*c**2
    Q15O = (m15O - m15N)*c**2 -0.997*Mev
    Q15N_1H = (m15N + mH - m12C -m4He)*c**2

    """
    This is the combinated energy output of the CNO cycle, which is the sum of
    the individual energy outputs.
    """
    QCNO = Q12C_1H + Q13N + Q13C_1H + Q14N_1H + Q15O + Q15N_1H

    #Making a list with all the energy outputs, so I can return this and print
    #them nicely
    Q_list = np.array([Q1H_1H,QDe_1H,Q3He_3He,Q3He_4He,Q7Be_e,Q7Li_1H,Q7Be_1H,Q8Bor,Q8Be,Q_decay,Q12C_1H,Q13N,Q13C_1H,Q14N_1H,Q15O,Q15N_1H])

    """
    Here we calculate the total energy output in each branch, including the
    neutrino energy so we can calculate the energy loss
    """
    Q_PP1_tot = (Q3He_3He + 2*(Q1H_1H + QDe_1H)) + 2*0.265*Mev
    Q_PP2_tot = (Q3He_4He + Q1H_1H + QDe_1H + 0.265*Mev) + Q7Be_e + Q7Li_1H + 0.81*Mev
    Q_PP3_tot = (Q3He_4He + Q1H_1H + QDe_1H + 0.265*Mev) + (Q7Be_1H + (Q8Bor + Q8Be)) + 6.711*Mev
    Q_CNO_tot = QCNO + 0.707*Mev + 0.997*Mev
    #Making an array of the total energy output so we can print this nicely
    Q_tot = np.array([Q_PP1_tot,Q_PP2_tot,Q_PP3_tot,Q_CNO_tot])

    """
    Here we calculate the neutrino energy loss percentages
    """
    neutrino_loss_pp1 = 100*2*0.265*Mev/Q_PP1_tot
    neutrino_loss_pp2 = 100*(0.81+0.265)*Mev/Q_PP2_tot
    neutrino_loss_pp3 = 100*(6.711+0.265)*Mev/Q_PP3_tot
    neutrino_loss_cno = 100*(0.707*Mev + 0.997*Mev)/Q_CNO_tot
    neutrino_loss = np.array([neutrino_loss_pp1,neutrino_loss_pp2,neutrino_loss_pp3,neutrino_loss_cno])

    """
    Number densities of the particles
    """
    n_1H = rho_core*X/mH #1/m^3
    n_3He = rho_core*Y_3He/(m3He)#1/m^3
    n_4He = rho_core*Y/(m4He)#1/m^3
    n_7Be = rho_core*Z_7Be/(m7Be)#1/m^3
    n_7Li = rho_core*Z_7Li/(m7Li)#1/m^3
    n_14N = rho_core*Z_14N/(m14N)#1/m^3
    #Number density of electron based on the other number density of the other
    #particles
    ne = 1*n_1H+2*n_3He+2*n_4He+3*n_7Li+4*n_7Be+7*n_14N#/(19*mu)

    """
    Nuclear reaction rates based on table 3.1 in Lecture Note
    """
    NA_lamPP = (4.01*10**(-15)*T9**(-2/3)*np.exp(-3.380*T9**(-1/3))*(1+0.123*T9**(1/3) + 1.09*T9**(2/3) + 0.938*T9))/100**3
    NA_lam33 = (6.04*10**(10)*T9**(-2/3)*np.exp(-12.276*T9**(-1/3))*(1+0.034*T9**(1/3) - 0.522*T9**(2/3)-0.124*T9 + 0.353*T9**(4/3) + 0.213*T9**(5/3)))/100**3
    NA_lam34 = (5.61*10**6*T9_str1**(5/6)*T9**(-3/2)*np.exp(-12.826*T9_str1**(-1/3)))/100**3
    NA_lam_m_17 = (1.096*10**9*T9**(-2/3)*np.exp(-8.472*T9**(-1/3))-4.830*10**8*T9_str2**(5/6)*T9**(-3/2)*np.exp(-8.472*T9_str2**(-1/3))+1.06*10**(10)*T9**(-3/2)*np.exp(-30.442*T9**(-1)))/100**3
    NA_lam17 = (3.11*10**5*T9**(-2/3)*np.exp(-10.262*T9**(-1/3))+2.53*10**3*T9**(-3/2)*np.exp(-7.306*T9**(-1)))/100**3
    NA_lamP14 = (4.90*10**7*T9**(-2/3)*np.exp(-15.228*T9**(-1/3) -0.092*T9**2)*(1+0.027*T9**(1/3)-0.778*T9**(2/3) - 0.149*T9 + 0.261*T9**(4/3) + 0.127*T9**(5/3))+2.37*10**3*T9**(-3/2)*np.exp(-3.011*T9**(-1))+2.19*10**4*np.exp(-12.53*T9**(-1)))/100**3

    """
    Here we have added an if-test for the electron capture, making sure it has
    the values for the upper limits when temperature is below T = 10^6K
    """
    if (10**(6-9)> T9):#.all():
        NA_lame7 = (1.57*10**(-7))/ne/100**3
    else:
        NA_lame7 = (1.34*10**(-10)*T9**(-1/2)*(1 - 0.537*T9**(1/3) + 3.86*T9**(2/3) + 0.0027*T9**(-1)*np.exp(2.515*10**(-3)*T9**(-1))))/100**3

    """
    Here we write the reaction rates for all our reactions, starting with the
    the common steps and the PPI reaction
    """
    r1H_1H = n_1H**2/(rho_core*2)*NA_lamPP/NA
    r3He_3He = n_3He**2/(rho_core*2)*NA_lam33/NA
    r3He_4He = n_3He*n_4He/(rho_core)*NA_lam34/NA

    """
    In order to prevent over consuming the helium produced in the common steps
    we make a scale that we use to make sure the PP branches use only what is
    being produced.
    """
    if (2*r3He_3He + r3He_4He) > r1H_1H:
        R = (r1H_1H/(2*r3He_3He + r3He_4He))
        r3He_3He = r3He_3He*R
        r3He_4He = r3He_4He*R

    r7Be_e = n_7Be*ne /(rho_core)*NA_lame7/NA
    r7Be_1H = n_7Be*n_1H/(rho_core)*NA_lam17/NA

    """
    Here we make sure that the rates of the first reactions in PPII and PPIII
    do not consume more than what PPI produces.
    """
    if (r7Be_e + r7Be_1H) > r3He_4He:
        R = (r3He_4He/(r7Be_e + r7Be_1H))
        r7Be_e = r7Be_e*R
        r7Be_1H = r7Be_1H*R

    r7Li_1H = n_7Li*n_1H/(rho_core)*NA_lam_m_17/NA

    """
    Our final scaling is to make sure the final step in PPII does noe consume
    more than what the previous reaction produces.
    """
    if r7Li_1H > r7Be_e:
        r7Li_1H = r7Li_1H*(r7Be_e/r7Li_1H)

    r14N_1H = n_14N*n_1H/(rho_core)*NA_lamP14/NA

    """
    This is our sanity check. To run it, use "run" when calling the function
    It will print out the check for the values of core temperature and density
    that is placed in the function.
    """
    if sanity_check == "run":
        print("---------------------SANITY CHECK-------------------------")
        print(f"T_core = {temp_core:.3f} 10^9 K and rho_core = {rho_core} kg m^-3")
        print(f"r1H_1H*(Q1H_1H + QDe_1H)*rho_core   = {r1H_1H*(Q1H_1H + QDe_1H)*rho_core:.3g} J m^-3 s^-1")
        print(f"r3He_3He*Q3He_3He*rho_core          = {r3He_3He*Q3He_3He*rho_core:.3g} J m^-3 s^-1")
        print(f"r3He_4He*Q3He_4He*rho_core          = {r3He_4He*Q3He_4He*rho_core:.3g} J m^-3 s^-1")
        print(f"r7Be_e*Q7Be_e*rho_core              = {r7Be_e*Q7Be_e*rho_core:.3g} J m^-3 s^-1")
        print(f"r7Li_1H*Q7Li_1H*rho_core            = {r7Li_1H*Q7Li_1H*rho_core:.3g} J m^-3 s^-1")
        print(f"r7Be_1H*(Q7Be_1H + QDecay)*rho_core = {r7Be_1H*(Q7Be_1H + Q8Bor + Q8Be)*rho_core:.3g} J m^-3 s^-1")
        print(f"r14N_1H*QCNO*rho_core               = {r14N_1H*QCNO*rho_core:.3g} J m^-3 s^-1")

    """
    Here we calculate the energy production for PPI, PPII, PPIII, CNO and the
    total combined energy output.
    """
    eps_PP1 = r3He_3He*(Q3He_3He + 2*(Q1H_1H + QDe_1H))
    eps_PP2 = r3He_4He*(Q3He_4He + Q1H_1H + QDe_1H) + r7Be_e*Q7Be_e + r7Li_1H*Q7Li_1H
    eps_PP3 = r3He_4He*(Q3He_4He + Q1H_1H + QDe_1H) + r7Be_1H*(Q7Be_1H + (Q8Bor + Q8Be))
    eps_CNO = r14N_1H*QCNO
    eps_tot = eps_PP1 + eps_PP2 + eps_PP3 + eps_CNO
    """
    In order to print everything we need, we return the energy production in an
    array, the _list with all the energy outputs, the neutrino loss percentages
    and the total energy output Q_tot
    """
    return np.array([eps_PP1, eps_PP2, eps_PP3, eps_CNO, eps_tot]), Q_list, neutrino_loss, Q_tot

"""
This is where we specify what core temperature and density we want ot use.
For both sanity check and all of our calculations we use the rho_core from the
Sun. For the temperature we have firs the solar core temperature, which is the
one we used for our calculations and which is the one that will be the default
when handin in this code. The other temperature, which is currently #-out is
the other temperature for the sanity check that can be tested. This means the
Sanity check can test one value of the temperature at a time.
"""
rho_core = 1.62*10**5 #kg m**-3
temp_core = 1.57*10**7/10**9 #K/10^9
#temp_core = 10**8/10**9

"""
Here we call on the function to make it easier to use further
"""
energy, Q_list, neutrino_loss, Q_tot = energy_production(temp_core, rho_core, "run")
print()

"""
Here we make arrays with the name and the acutal values of the energy outputs.
The actual values are in the units of J, since this is what we used in the
calculations
"""
Q_list_name = np.array(["Q1H_1H","QDe_1H","Q3He_3He","Q3He_4He","Q7Be_e","Q7Li_1H","Q7Be_1H","Q8Bor","Q8Be","Q_decay","Q12C_1H","Q13N","Q13C_1H","Q14N_1H","Q15O","Q15N_1H"])

Q_real = np.array([1.177,5.494,12.860,1.586,0.049, 17.346,0.137,8.367,2.995,8.367+2.995,1.944,1.513,7.551,7.297,1.757,4.966])*Mev

"""
Now we start the printing of the results in the project.
The first result is the energy outputs. Here we print the name, calculated value
in units of J, then in MeV and then the percent error between the acutal values
and what we calculated.
"""
print("----------ENERGY OUTPUT FROM PP CHAIN AND CNO CYCLE-------")
for i in range(len(Q_list_name)):
    print(f"{Q_list_name[i]:10s}","=", f"{Q_list[i]:.4g}", "J", "=", f"{Q_list[i]/Mev:.4g}", "MeV" ,f", Percent error % = {100* np.abs(Q_list[i] - Q_real[i])/Q_real[i]:.2g}" )
print()

"""
Now we print the percentage of energies lost to neutrinos, as well as the Total
energy output in each of the branches.
"""
print("------ENERGY LOST TO NEUTRINOS------")
neutrino_name = ["Percentage lost PPI","Percentage lost PPII","Percentage lost PPIII","Percentage lost CNO"]
for i in range(len(neutrino_loss)):
    print(f"{neutrino_name[i]:22s}","=", f"{neutrino_loss[i]:.2f}%", f"Total energy output Q in branch: {Q_tot[i]/Mev:.4f} MeV")
print()

"""
Here we print the energy production for each reaction branch
"""
print("----------------ENERGY PRODUCTION----------------")
print(f"Energy production PPI   = {energy[0]:.3g} [J kg^-1 s^-1]")
print(f"Energy production PPII  = {energy[1]:.3g} [J kg^-1 s^-1]")
print(f"Energy production PPIII = {energy[2]:.3g} [J kg^-1 s^-1]")
print(f"Energy production CNO   = {energy[3]:.3g} [J kg^-1 s^-1]")
print(f"Total energy production = {energy[4]:.3g} [J kg^-1 s^-1]")

"""
For plotting the energy production we used logaritms, which is why we redefine
the temperature as such. Here we use the temperature intervall specified in the
task.
"""
temp = np.linspace(4-9, 9-9, 1000)
temp = 10**temp

"""
To easier plot the energy production we made arrays with zeros with the same
lenght as the temperature temp.
"""
eps_PP1_L = np.zeros(len(temp))
eps_PP2_L = np.zeros(len(temp))
eps_PP3_L = np.zeros(len(temp))
eps_CNO_L = np.zeros(len(temp))
eps_tot_L = np.zeros(len(temp))

"""
We make arrays of the energy productions that we can plot
"""
for i in range(len(temp)):
    eps_PP1_L[i], eps_PP2_L[i], eps_PP3_L[i], eps_CNO_L[i], eps_tot_L[i] = energy_production(temp[i], rho_core, "no")[0]

"""
Here we plot the energy production
"""
plt.plot(np.log10(temp), eps_PP1_L/eps_tot_L, label="εPP1/εtot")
plt.plot(np.log10(temp), eps_PP2_L/eps_tot_L, label="εPP2/εtot")
plt.plot(np.log10(temp), eps_PP3_L/eps_tot_L, label="εPP3/εtot")
plt.plot(np.log10(temp), eps_CNO_L/eps_tot_L, label="εCNO/εtot")
plt.legend()
plt.xlabel("log10(T[K]/10^9)")
plt.ylabel("ε/ε_tot")
plt.title("Relative energy production")
plt.show()

"""
This final part is the Gamow peak. We start by defining the Gamow peak as a
function. Its variables are the temperature of the core, an energy interval, the
masses of the two reaction particles and their atomic number.
"""

def Gamow(T, E, m1, m2, Z1, Z2):
    kB = 1.3806*10**(-23)
    e = 1.6022*10**(-19)
    h = 6.6261*10**(-34)
    eps = 8.8542*10**(-12)
    m = m1*m2/(m1+m2)
    return np.exp(-E/(kB*T))*np.exp(-np.sqrt(m/(2*E))*Z1*Z2*e**2*np.pi/(eps*h))

"""
Here we define the variables in Gamow. The masses are listed further up, as well
as mass fractions. The temperature is the core temperature of the Sun. The
energy is converted into logaritmic scale.
"""
T = 1.57*10**7
E = np.linspace((-17), (-13), 10000) #J
E = 10**E

"""
Here we plot the Gamow peaks for all the reactions we wanted. To get the
relative probability we divided the Gamow peak by the sum for each of the
reactions
"""
plt.plot(np.log10(E), Gamow(T, E, mH, mH, 1, 1)/np.sum(Gamow(T, E, mH, mH, 1, 1)), label="λpp")
plt.plot(np.log10(E), Gamow(T, E, mH, mDe, 1, 1)/np.sum(Gamow(T, E, mH, mDe, 1, 1)), label="λpd")
plt.plot(np.log10(E), Gamow(T, E, m3He, m3He, 2, 2)/np.sum(Gamow(T, E, m3He, m3He, 2, 2)), label="λ33")
plt.plot(np.log10(E), Gamow(T, E, m3He, m4He, 2, 2)/np.sum(Gamow(T, E, m3He, m4He, 2, 2)), label="λ34")
#plt.plot(np.log10(E), Gamow(T, E, m7Be, me)/np.mean(Gamow(T, E, m7Be, me)), label="λe7")
plt.plot(np.log10(E), Gamow(T, E, m7Li, mH, 3, 1)/np.sum(Gamow(T, E, m7Li, mH, 3, 1)), label="λ'17")
plt.plot(np.log10(E), Gamow(T, E, m7Be, mH, 4, 1)/np.sum(Gamow(T, E, m7Be, mH, 4, 1)), label="λ17")
plt.plot(np.log10(E), Gamow(T, E, m12C, mH, 6, 1)/np.sum(Gamow(T, E, m12C, mH, 6, 1)), label="λp12")
plt.plot(np.log10(E), Gamow(T, E, m13C, mH, 6, 1)/np.sum(Gamow(T, E, m13C, mH, 6, 1)), label="λp13")
plt.plot(np.log10(E), Gamow(T, E, m14N, mH, 7, 1)/np.sum(Gamow(T, E, m14N, mH, 7, 1)), label="λp14")
plt.plot(np.log10(E), Gamow(T, E, m15N, mH, 7, 1)/np.sum(Gamow(T, E, m15N, mH, 7, 1)), label="λp15")
plt.xlim(-16, -14)
plt.title("Relative probability of Gamow peak")
plt.ylabel("Relative probability")
plt.xlabel("log10(Energy [J])]")
plt.legend()
plt.show()

"""
Here we print the energies for the Gamow peaks
"""
print()
print("----------------GAMOW PEAKS----------------")
print(f"Gamow peak", "λpp :","E =", f"{E[np.where(Gamow(T, E, mH, mH, 1, 1) == np.max(Gamow(T, E, mH, mH, 1, 1)))[0][0]]:.4g}", "J")
print(f"Gamow peak", "λpd :","E =", f"{E[np.where(Gamow(T, E, mH, mDe, 1, 1) == np.max(Gamow(T, E, mH, mDe, 1, 1)))[0][0]]:.4g}", "J")
print(f"Gamow peak", "λ33 :","E =", f"{E[np.where(Gamow(T, E, m3He, m3He, 2, 2) == np.max(Gamow(T, E, m3He, m3He, 2, 2)))[0][0]]:.4g}", "J")
print(f"Gamow peak", "λ34 :","E =", f"{E[np.where(Gamow(T, E, m3He, m4He, 2, 2) == np.max(Gamow(T, E, m3He, m4He, 2, 2)))[0][0]]:.4g}", "J")
print(f"Gamow peak", "λ'17:","E =", f"{E[np.where(Gamow(T, E, m7Li, mH, 3, 1) == np.max(Gamow(T, E, m7Li, mH, 3, 1)))[0][0]]:.4g}", "J")
print(f"Gamow peak", "λ17 :","E =", f"{E[np.where(Gamow(T, E, m7Be, mH, 4, 1) == np.max(Gamow(T, E, m7Be, mH, 4, 1)))[0][0]]:.4g}", "J")
print(f"Gamow peak", "λp12:","E =", f"{E[np.where(Gamow(T, E, m12C, mH, 6, 1) == np.max(Gamow(T, E, m12C, mH, 6, 1)))[0][0]]:.4g}", "J")
print(f"Gamow peak", "λp13:","E =", f"{E[np.where(Gamow(T, E, m13C, mH, 6, 1) == np.max(Gamow(T, E, m13C, mH, 6, 1)))[0][0]]:.4g}", "J")
print(f"Gamow peak", "λp14:","E =", f"{E[np.where(Gamow(T, E, m14N, mH, 7, 1) == np.max(Gamow(T, E, m14N, mH, 7, 1)))[0][0]]:.4g}", "J")
print(f"Gamow peak", "λp15:","E =", f"{E[np.where(Gamow(T, E, m15N, mH, 7, 1) == np.max(Gamow(T, E, m15N, mH, 7, 1)))[0][0]]:.4g}", "J")
