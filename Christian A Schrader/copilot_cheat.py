"""
This code is for modeling the energy production
in the solar core, 1, 10 min to make
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as c
import tabulate as tab

class EnergyProduction:
    """
    This class is for modeling the energy production
    by the fusion chains in the solar core
    """
    def __init__(self, temperature, density):
        self.temperature = temperature
        self.density = density

        # Mass fractions
        self.X = 0.7
        self.Y32He = 1e-10
        self.Y = 0.29
        self.Z73Li = 1e-7
        self.Z74Be = 1e-7
        self.Z147N = 1e-11

        self.convert_MeV2J = 1.60218e-13    # MeV to J
        self.cubic_cm2cubic_m = 1e-6        # cm^3 to m^3
        self.u = c.atomic_mass              # Atomic mass unit
        self.N_A = c.Avogadro               # Avogadro's number
        self.k_B = c.Boltzmann              # Boltzmann's constant

    def number_densities(self, of_type, number):
        """
        Calculates the number density of a species
        """
        ans = (of_type*self.density)/(number*self.u)
        return ans
    
    def reaction_rate(self, n_i, n_k, lambda_ik):
        """
        Calculates the reaction rate
        """
        if n_i == n_k:
            """
            Kroenicker delta function
            """
            den = 2*self.density
        else:
            den = self.density
        
        ans = n_i*n_k*lambda_ik/den
        return ans
    
    def all_reaction_rates(self, temperature=None):
        """
        Calculates all the reaction rates
        for the PP and CNO from the method reaction_rate
        and uses the values from the tables in the project

        """
        if temperature is None:
            T = self.temperature
        else:
            T = temperature
            
        # Get the mass fractions from global variables
        X = self.X
        Y = self.Y
        Y32He = self.Y32He
        Z73Li = self.Z73Li
        Z74Be = self.Z74Be
        Z147N = self.Z147N

        # Sets the special temperatures for lambda calculations
        T9 = T*1e-9
        T9_34 = T9/(1 + 4.95e-2*T9)
        T9_17 = T9/(1 + 0.759*T9)
        
        # Calculates the number densities
        #nx = self.number_densities(X, 1)
        #ny3 = self.number_densities(Y32He, 3)
        #ny = self.number_densities(Y, 1)
        #nz = self.number_densities(Z147N, 7)
        #nz3 = self.number_densities(Z73Li, 7)
        #nz4 = self.number_densities(Z74Be, 7)
        """
        The above was autocalculated by github
        """
        # Calculate the number densites of the different species
        # Could you do that for the mass fractions provided?
        nx = self.number_densities(X, 1)
        ny3 = self.number_densities(Y32He, 3)
        ny = self.number_densities(Y, 4)
        nz14 = self.number_densities(Z147N, 14)
        nz3 = self.number_densities(Z73Li, 3)
        nz7 = self.number_densities(Z74Be, 7)

        nz = nz14 + nz7 + nz3     # Total number density of Z
        ne = nx + 2*ny + nz         # Electron number density
        
        # Calculate the lambda values needed for the reaction rates
        # lambda for 1H + 1H -> 2H + e+ + nu_e
        lambda_pp = (4.01e-15*T9**(-2/3)*np.exp(-3.380*T9**(-1/3))*(1 + 0.123*T9**(1/3) 
                    + 1.09*T9**(2/3) + 0.938*T9 + 0.797*T9**(4/3)))
        lambda_pp = lambda_pp*self.cubic_cm2cubic_m/self.N_A
        
        # lambda for 3He + 3He -> 4He + 2 1H
        lambda_33 = (6.04e10*T9**(-2/3)*np.exp(-12.276*T9**(-1/3))*(1 + 0.034*T9**(1/3) 
                    - 0.522*T9**(2/3) - 0.124*T9 + 0.353*T9**(4/3) + 0.213*T9**(5/3)))
        lambda_33 = lambda_33*self.cubic_cm2cubic_m/self.N_A

        # lambda for 3He + 4He -> 7Be + gamma
        #lambda_34 = 5.61e6*T9_34**(-2/3)*np.exp(-12.826*T9_34**(-1/3))*(1 + 0.034*T9_34**(1/3) - 0.522*T9_34**(2/3) - 0.124*T9_34 + 0.353*T9_34**(4/3) + 0.213*T9_34**(5/3))
        #That was wrong, but keeping it for reference
        lambda_34 = 5.61e6*T9_34**(5/6)*np.exp(-12.826*T9_34**(-1/3))*T9**(-3/2)
        lambda_34 = lambda_34*self.cubic_cm2cubic_m/self.N_A

        # lambda for 7Be + e- -> 7Li + nu_e
        lambda_e7 = (1.34e-10*T9**(-1/2)*(1 - 0.537*T9**(1/3) 
                    + 3.86*T9**(2/3) + 0.0027*T9**(-1)*np.exp(2.515e-3*T9_17**(-1))))
        lambda_e7 = lambda_e7*self.cubic_cm2cubic_m/self.N_A
        # Here it used T9_17, but I it should be T9, see confusion

        # Need to check for upper limit for lambda_e7
        if T < 1e6:
            # Does not happend unless temperature is below 1e6
            """
            if lambda_e7 > 1.34e-10:
                lambda_e7 = 1.34e-10
            """
            if lambda_e7 > 1.57e-7/(ne*self.N_A):
                lambda_e7 = 1.57e-7/(ne*self.N_A)
        
        if np.abs(lambda_e7 - 1.34e-10) < 1e-10:
            print("Lambda_e7 is equal to 1.34e-10")

        # Need to double check this one with raw values
        # Are the two above the same?
        # lambda for 7Be + 1H -> 8B + gamma
        lambda_p7 = 3.11e5*T9**(-2/3)*np.exp(-10.262*T9**(-1/3)) + 2.53e3*T9**(-3/2)*np.exp(-7.306*T9**(-1))
        lambda_p7 = lambda_p7*self.cubic_cm2cubic_m/self.N_A

        # lambda 7Li + 1H -> 2 4He
        #lambda_p7_prime = 1.096e9*T9**(-2/3)*np.exp(-8.472*T9**(-1/3)) - 4.83e8**T9_17**(5/6)*np.exp(-8.472*T9_17**(-1/3)) + 1.06e10*T9**(-3/2)*np.exp(-30.442*T9**(-1))
        #lambda_p7_prime = lambda_p7_prime*self.cubic_cm2cubic_m/self.N_A
        lambda_p7_prime = self.cubic_cm2cubic_m*(1.096e9*T9**(-2/3)*np.exp(-8.472*T9**(-1/3)) -
                                4.830e8**T9_17**(5/6)*T9**(-3/2) *
                                np.exp(-8.472*T9_17**(-1/3)) +
                                1.06e10*T9**(-3/2)*np.exp(-30.442/T9))/self.N_A

        # lambda for 14N + 1H -> 15O + gamma
        lambda_p14 = 4.9e7*T9**(-2/3)*np.exp(-15.228*T9**(-1/3) - 0.092*T9**2) * (1 + 0.027*T9**(1/3) - 0.778*T9**(2/3) - 0.149*T9 + 0.261*T9**(4/3) + 0.127*T9**(5/3)) + 2.37e3*T9**(-3/2)*np.exp(-3.011*T9**(-1)) + 2.19e4*np.exp(-12.53*T9**(-1))
        lambda_p14 = lambda_p14*self.cubic_cm2cubic_m/self.N_A


        # reaction rates
        # reaction rate H + H -> 2H + e+ + nu_e
        r_pp = self.reaction_rate(nx, nx, lambda_pp1)

        # reaction rate 3He + 3He -> 4He + 2 1H
        r_33 = self.reaction_rate(ny3, ny3, lambda_33)

        # reaction rate 3He + 4He -> 7Be + gamma
        r_34 = self.reaction_rate(ny3, ny, lambda_34)

        # reaction rate 7Be + e- -> 7Li + nu_e
        r_e7 = self.reaction_rate(ne, nz7, lambda_e7)

        # reaction rate 7Be + 1H -> 8B + gamma
        r_p7 = self.reaction_rate(nx, nz7, lambda_p7)

        # reaction rate 7Li + 1H -> 2 4He
        r_p7_prime = self.reaction_rate(nx, nz3, lambda_p7_prime)

        # reaction rate 14N + 1H -> 15O + gamma
        r_p14 = self.reaction_rate(nx, nz14, lambda_p14)

        """
        Check element abundance
        """
        if (r_33 + r_34) >r_pp:
            r_33 = (r_33*r_pp)/(r_33 + r_34)
            r_34 = (r_34*r_pp)/(r_33 + r_34)

        if (r_e7 + r_p7) > r_34:
            r_e7 = (r_e7*r_34)/(r_e7 + r_p7)
            r_p7 = (r_p7*r_34)/(r_e7 + r_p7)

        if r_p7_prime > r_e7:
            """
            Upper limit of r_p7_prime is r_e7
            """
            r_p7_prime = r_e7

        r = [r_pp, r_33, r_34, r_e7, r_p7, r_p7_prime, r_p14]

        return r
    
    def fusion_reaction(self, temperature=None):
        """
        Method to calculate the fusion reactions
        using the Q values from the paper
        """
        if temperature is None:
            temperature = self.temperature
        else:
            temperature = temperature
        MeV2J = self.convert_MeV2J

        r = self.all_reaction_rates(temperature)
        r_pp = r[0]
        r_33 = r[1]
        r_34 = r[2]
        r_e7 = r[3]
        r_p7 = r[4]
        r_p7_prime = r[5]
        r_p14 = r[6]

        # Q values
        Q_pp = 1.177*MeV2J
        Q_33 = 12.86*MeV2J
        Q_34 = 1.586*MeV2J
        Q_e7 = 0.049*MeV2J
        Q_p7 = 0.137*MeV2J
        Q_p7_prime = 17.346*MeV2J
        Q_p14 = 7.297*MeV2J

        Q_pd = 5.494*MeV2J
        Q_8 = 8.367*MeV2J
        Q_8_prime = 2.995*MeV2J
        Q_decay = Q_8 + Q_8_prime

        Q_CNO = (1.944 + 4.966 + 1.757 + 1.513 + 7.551 + 7.297)*MeV2J
        Q = [Q_pp, Q_33, Q_34, Q_e7, Q_p7, Q_p7_prime, Q_p14, Q_pd, Q_decay, Q_CNO]

        # The energy generation rates
        epsilon_ppI = r_33*(Q_33 + 2*(Q_pp + Q_pd))
        epsilon_ppII = r_34*(Q_34 + Q_pp + Q_pd) + r_e7*Q_e7 + r_p7_prime*Q_p7_prime
        epsilon_ppIII = (r_34*(Q_34 + Q_pp + Q_pd) + r_p7*(Q_p7 + Q_decay))
        epsilon_CNO = r_p14*Q_CNO

        epsilon = [epsilon_ppI, epsilon_ppII, epsilon_ppIII, epsilon_CNO]
        
        if temperature is None:
            self.r, self.Q, self.epsilon = r, Q, epsilon
        return r, Q, epsilon
    

    def reduced_mass(self, lambda_val):
        """
        Method to calculate the reduced mass
        """
        
        zi, zk = lambda_val
        mi = zi*self.u
        mk = zk*self.u

        mu = (mi*mk)/(mi + mk)
        return mu
    
    def gamow_peak(self, lambda_val, values=[-17, -13]):
        """
        Method to calculate the Gamow peak
        """
        N = 1000
        E = np.logspace(values[0], values[1], N)

        # Constants
        beta = 1/(self.k_B*self.temperature)
        eps0 = c.epsilon_0
        h = c.h
        e = c.e
        m = self.reduced_mass(lambda_val)
        
        # The gamow peak
        zi, zk = lambda_val
        prod1 = np.exp(-E*beta)
        arg = np.sqrt(m/(2*E))*zi*zk*e**2*np.pi/(eps0*h)
        prod2 = np.exp(-arg)
        P_g = prod1*prod2

        P = P_g/np.max(P_g)
        return E, P
    
    def plot_gamow_peak(self, values=[-17, -13]):
        """
        Plots the Gamow peak for all the lambdas used in the 
        fusion_reaction from the gamow_peak method
        """
        lambda_pp = [1, 1]
        lambda_33 = [2, 2]
        lambda_34 = [2, 2]
        lambda_p7_prime = [3, 1]
        lambda_p7 = [4, 1]
        lambda_e7 = [7, 1]
        lambda_p14 = [7, 1]

        # list of all the lambdas
        lambda_list = [lambda_pp, lambda_33, lambda_34, lambda_e7, lambda_p7, lambda_p7_prime, lambda_p14]

        # list of all the names
        names = [r'$\lambda_{pp}$', r'$\lambda_{33}$', r'$\lambda_{34}$', r'$\lambda_{e7}$', r'$\lambda_{p7}$', r'$\lambda_{p7\prime}$', r'$\lambda_{p14}$']

        N = len(lambda_list)
        E = np.zeros((N, 1000))
        P = np.zeros((N, 1000))
        for i in range(N):
            E[i], P[i] = self.gamow_peak(lambda_list[i], values=values)
            plt.semilogx(E[i], P[i], label=names[i])

        plt.xlabel(r'$E$ [MeV]')
        plt.ylabel(r'$P$')
        plt.legend()
        plt.show()
        
    def sanity_check(self):
        """
        Method to check that the energy generation rates
        are equal to the energy loss rates
        """
        temp = self.temperature
        rho = self.density

        expected = [4.04e2, 8.68e-9, 4.68e-5, 1.49e-6,
                    5.29e-4, 1.63e-6, 9.18e-8]
        
        r, Q, epsilon = self.fusion_reaction(temp)
        r_pp, r_33, r_34, r_e7, r_p7, r_p7_prime, r_p14 = r
        Q_pp, Q_33, Q_34, Q_e7, Q_p7, Q_p7_prime, Q_p14, Q_pd, Q_decay, Q_CNO = Q
        epsilon_ppI, epsilon_ppII, epsilon_ppIII, epsilon_CNO = epsilon

        # list to store all the values for sanity check
        tmp = []
        tmp.append(r_pp*(Q_pp + 2*Q_pd)*rho)
        tmp.append(r_33*Q_33*rho)
        tmp.append(r_34*Q_34*rho)
        tmp.append(r_e7*Q_e7*rho)
        tmp.append(r_p7_prime*Q_p7_prime*rho)
        tmp.append(r_p7*(Q_p7 + Q_decay)*rho)
        tmp.append(r_p14*Q_CNO*rho)

        # print the difference between the expected and calculated

        for i in range(len(tmp)):
            diff = np.abs(tmp[i] - expected[i])
            # print the value in table format with the difference
            print(tab.tabulate([[tmp[i], expected[i], diff]], headers=['Calculated', 'Expected', 'Difference'], tablefmt='latex'))

    def print_energy(self):
        """
        Prints the table of energy
        in a nicely printed table
        """
        # Verify that the energy is calculated
        try:
            epsilon = self.epsilon
        except AttributeError:
            epsilon = self.fusion_reaction(self.temperature)[-1]
        text = ['PP I', 'PP II', 'PP III', 'CNO']
        zipped = zip(text, epsilon)
        print(tab.tabulate(zipped, headers=['Reaction', 'Energy [J/s/kg]'], tablefmt='latex'))

    def temperature_plot(self, N=1000):
        """
        Method to plot the temperature dependence
        """
        T = np.logspace(4, 9, N)
        PPI = np.zeros(N)
        PPII = np.zeros(N)
        PPIII = np.zeros(N)
        CNO = np.zeros(N)
        total = np.zeros(N)

        for i in range(len(T)):
            eps = self.fusion_reaction(T[i])[-1]
            total[i] = np.sum(eps)
            PPI[i] = eps[0]/total[i]
            PPII[i] = eps[1]/total[i]
            PPIII[i] = eps[2]/total[i]
            CNO[i] = eps[3]/total[i]

        # plt.title('Energy production plot for T$\in[10^{4}, 10^{9}]$')
        plt.xscale('log')
        plt.plot(T, PPI, label='$\epsilon_{PPI}/\epsilon_{tot}$')
        plt.plot(T, PPII, label='$\epsilon_{PPII}/\epsilon_{tot}$')
        plt.plot(T, PPIII, label='$\epsilon_{PPIII}/\epsilon_{tot}$')
        plt.plot(T, CNO, label='$\epsilon_{CNO}/\epsilon_{tot}$')
        plt.xlabel('T[K]')
        plt.ylabel('$\epsilon/\epsilon_{tot}$')
        plt.legend()
        plt.show()
        

if __name__ == "__main__":
    T = 1.57e7      # K
    rho = 1.62e5    # kg/m^3
    test = EnergyProduction(T, rho)
    test.temperature_plot()
    test.plot_gamow_peak()
    test.print_energy()
    test.sanity_check()