import matplotlib.pyplot as plt
import numpy as np
class StarCore:
    def __init__(self, T, rho, X, Y3, Y, Z7Li, Z7Be, Z14N):
        # Store parameters
        self.T = T
        self.rho = rho
        self.X = X
        self.Y3 = Y3
        self.Y = Y
        self.Z7Li = Z7Li
        self.Z7Be = Z7Be
        self.Z14N = Z14N

        # Conversion factor MeV to J
        self.MeV_to_J = 1.60218e-13
# Masses of particles in kg
        self.mass_H = 1.67e-27
        self.mass_He3 = 5.01e-27
        self.mass_He4 = 6.64e-27
        self.mass_Li7 = 1.16e-26
        self.mass_Be7 = 1.17e-26
        self.mass_N14 = 2.33e-26

        # Number densities
        self.n_H = self.X * self.rho / self.mass_H
        self.n_He3 = self.Y3 * self.rho / self.mass_He3
        self.n_He4 = self.Y * self.rho / self.mass_He4
        self.n_Li7 = self.Z7Li * self.rho / self.mass_Li7
        self.n_Be7 = self.Z7Be * self.rho / self.mass_Be7
        self.n_N14 = self.Z14N * self.rho / self.mass_N14

    def reaction_rate(self, n_i, n_j, lambda_ij, i, j):
        # Calculate the reaction rate using the provided formula
        delta_ij = 1 if i == j else 0
        return n_i * n_j * lambda_ij / (self.rho * (1 + delta_ij))

    def lambda_pp(self):
        T9 = self.T / 1e9
        return 4.01e-15 * T9**(-2/3) * np.exp(-3.380 * T9**(-1/3))

    def lambda_33(self):
        T9 = self.T / 1e9
        return 6.04e10 * T9**(-2/3) * np.exp(-12.276 * T9**(-1/3))

    def lambda_34(self):
        T9 = self.T / 1e9
        T9_star = T9 / (1 + 4.95e-2 * T9)
        return 5.61e6 * T9_star**(5/6) * T9**(-3/2) * np.exp(-12.826 * T9_star**(-1/3))

    def lambda_e7(self):
        T9 = self.T / 1e9
        return 1.34e-10 * T9**(-1/2) * (1 - 0.537 * T9**(1/3) + 3.86 * T9**(2/3) + 0.0027 / T9 * np.exp(2.515e-3 / T9))

    def lambda_17_prime(self):
        T9 = self.T / 1e9
        T9_star = T9 / (1 + 0.759 * T9)
        return (1.096e9 * T9**(-2/3) * np.exp(-8.472 * T9**(-1/3)) -
                4.83e8 * T9_star**(5/6) * T9**(-3/2) * np.exp(-8.472 * T9_star**(-1/3)) +
                1.06e10 * T9**(-3/2) * np.exp(-30.442 / T9))

    def lambda_17(self):
        T9 = self.T / 1e9
        return 3.11e5 * T9**(-2/3) * np.exp(-10.262 * T9**(-1/3))

    def lambda_p14(self):
        T9 = self.T / 1e9
        return (4.9e7 * T9**(-2/3) * np.exp(-15.228 * T9**(-1/3) - 0.092 * T9**2) *
               (1 + 0.027 * T9**(1/3) - 0.778 * T9**(2/3) - 0.149 * T9 + 0.261 * T9**(4/3) + 0.127 * T9**(5/3)) +
               2.37e3 * T9**(-3/2) * np.exp(-3.011 / T9) + 2.19e4 * np.exp(-12.53 / T9))

    def r_pp(self):
        return self.reaction_rate(self.n_H, self.n_H, self.lambda_pp(), 1, 1)

    def r_33(self):
        return self.reaction_rate(self.n_He3, self.n_He3, self.lambda_33(), 2, 2)

    def r_34(self):
        return self.reaction_rate(self.n_He3, self.n_He4, self.lambda_34(), 3, 4)

    def r_e7(self):
        return self.reaction_rate(self.n_Be7, self.n_Be7, self.lambda_e7(), 7, 7)

    def r_p7_prime(self):
        return self.reaction_rate(self.n_Li7, self.n_H, self.lambda_17_prime(), 3, 1)

    def r_p7(self):
        return self.reaction_rate(self.n_Be7, self.n_H, self.lambda_17(), 4, 1)

    def r_p14(self):
        return self.reaction_rate(self.n_N14, self.n_H, self.lambda_p14(), 7, 1)

    def energy_production(self):
        # Reaction rates
        r_33 = self.r_33()
        r_34 = self.r_34()
        r_e7 = self.r_e7()
        r_p7 = self.r_p7()
        r_p7_prime = self.r_p7_prime()
        r_p14 = self.r_p14()

        # Energy released per reaction (Q values) in MeV
        Q_pp = 1.442
        Q_pd = 5.494
        Q_33 = 12.860
        Q_34 = 1.586
        Q_e7 = 0.861
        Q_p7 = 17.346
        Q_decay = 0.478
        Q_p7_prime = 17.346
        Q_CNO = 25.033

        # Convert Q values to Joules
        Q_pp *= self.MeV_to_J
        Q_pd *= self.MeV_to_J
        Q_33 *= self.MeV_to_J
        Q_34 *= self.MeV_to_J
        Q_e7 *= self.MeV_to_J
        Q_p7 *= self.MeV_to_J
        Q_decay *= self.MeV_to_J
        Q_p7_prime *= self.MeV_to_J
        Q_CNO *= self.MeV_to_J

        # Calculate the energy production in each chain
        PPI = r_33*(Q_33 + 2*(Q_pp + Q_pd))
        PPII = (r_34*(Q_34 + Q_pp + Q_pd) + r_e7*Q_e7 + r_p7_prime*Q_p7_prime)
        PPIII = (r_34*(Q_34 + Q_pp + Q_pd) + r_p7*(Q_p7 + Q_decay))
        CNO = r_p14*Q_CNO

        # Calculate the total energy production
        total_energy = PPI + PPII + PPIII + CNO

        return PPI, PPII, PPIII, CNO, total_energy

# Define the temperature range.
T = np.logspace(4, 9, 1000)
RHO = 1.62e5
total_energy_production = []
relative_energy_production_pp1 = []
relative_energy_production_pp2 = []
relative_energy_production_pp3 = []
relative_energy_production_cno = []
for temp in T:
    calc = StarCore(temp, RHO, 0.7, 1e-10, 0.29, 1e-7, 1e-7, 1e-11)
    energy_production_pp1, energy_production_pp2, energy_production_pp3, energy_production_cno, total_energy = calc.energy_production()

    # Normalize the energy productions so they sum to 1 at each temperature.
    total_energy_production.append(total_energy)
    relative_energy_production_pp1.append(energy_production_pp1 / total_energy)
    relative_energy_production_pp2.append(energy_production_pp2 / total_energy)
    relative_energy_production_pp3.append(energy_production_pp3 / total_energy)
    relative_energy_production_cno.append(energy_production_cno / total_energy)

# Plot the relative energy production of each process.
plt.figure(figsize=(10, 6))
plt.loglog(T, relative_energy_production_pp1, label="PPI")
plt.loglog(T, relative_energy_production_pp2, label="PPII")
plt.loglog(T, relative_energy_production_pp3, label="PPIII")
plt.loglog(T, relative_energy_production_cno, label="CNO")
plt.xlabel('Temperature (K)')
plt.ylabel('Relative energy production')
plt.title('Relative energy production as a function of temperature')
plt.legend()
plt.grid(True)
plt.show()

def plot_gamow_peak(E_range):
    # replace with your actual function for calculating relative probability
    relative_probability = np.random.random(len(E_range))

    plt.figure(figsize=(10, 6))
    plt.semilogx(E_range, relative_probability, label="Relative probability")
    plt.xlabel('Energy (J)')
    plt.ylabel('Relative Probability')
    plt.title('Gamow Peak')
    plt.legend()
    plt.grid(True)
    plt.show()

# Specify the energy range E.
E_range = np.logspace(-17, -13, 1000)
plot_gamow_peak(E_range)
