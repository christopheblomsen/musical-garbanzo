import numpy as np
import matplotlib.pyplot as plt

class EnergyProductionCalculator:
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

    def calculate_energy_production(self):
        # Calculate energy production for each fusion chain
        energy_ppi = self.calculate_energy_ppi()
        energy_ppii = self.calculate_energy_ppii()
        energy_ppiii = self.calculate_energy_ppiii()
        energy_cno = self.calculate_energy_cno()

        # Total energy production
        total_energy_production = energy_ppi + energy_ppii + energy_ppiii + energy_cno

        return total_energy_production

    def calculate_energy_ppi(self):
        # Constants
        Q01 = 1.177 * 1e-5  # Energy release per reaction (J)
        rho = self.density  # Density (kg/m^3)

        # Reaction rates
        lambda01 = 4.01 * 1e-15 * self.X**2 * self.Y**2 * (self.temperature / 1e9)**(-2/3) * \
                    np.exp(-3.380 * (self.temperature / 1e9)**(-1/3))  # Reaction rate (cm^3/s/mol)

        # Energy production rate
        r11 = Q01 * lambda01 * rho  # Energy production (J/m^3/s)

        return r11

    def calculate_energy_ppii(self):
        # Constants
        Q23 = 5.494 * 1e-3  # Energy release per reaction (J)
        rho = self.density  # Density (kg/m^3)

        # Reaction rates
        lambda23 = 5.61 * 1e-7 * self.X * self.Y**3 * (self.temperature / 1e9)**(-2/3) * \
                    np.exp(-12.276 * (self.temperature / 1e9)**(-1/3))  # Reaction rate (cm^3/s/mol)

        # Energy production rate
        r32 = Q23 * lambda23 * rho  # Energy production (J/m^3/s)

        # Check if step 1 can produce enough helium
        if r32 > self.calculate_energy_ppi():
            r32 = self.calculate_energy_ppi()  # Limit the energy production to what step 1 can produce

        return r32

    def calculate_energy_ppiii(self):
        # Constants
        Q34 = 1.098 * 1e-2  # Energy release per reaction (J)
        rho = self.density  # Density (kg/m^3)

        # Reaction rates
        lambda34 = 1.65 * 1e-7 * self.X * self.Y**2 * (self.temperature / 1e9)**(-2/3) * \
                    np.exp(-12.826 * (self.temperature / 1e9)**(-1/3))  # Reaction rate (cm^3/s/mol)

        # Energy production rate
        r32 = Q34 * lambda34 * rho  # Energy production (J/m^3/s)

        # Check if step 2 can produce enough helium
        if r32 > self.calculate_energy_ppii():
            r32 = self.calculate_energy_ppii()  # Limit the energy production to what step 2 can produce

        return r32

    def calculate_energy_cno(self):
        # Constants
        Qcnocycle = 8.705 * 1e-4  # Energy release per reaction (J)
        rho = self.density  # Density (kg/m^3)

        # Reaction rate
        lambda_cno = 8.24 * 1e-26 * self.X * self.Z147N * (self.temperature / 1e9)**(-2/3) * \
                     np.exp(-15.228 * (self.temperature / 1e9)**(-1/3))  # Reaction rate (cm^3/s/mol)

        # Energy production rate
        r14 = Qcnocycle * lambda_cno * rho  # Energy production (J/m^3/s)

        # Check if step 3 can produce enough beryllium
        if r14 > self.calculate_energy_ppiii():
            r14 = self.calculate_energy_ppiii()  # Limit the energy production to what step 3 can produce

        return r14


# Sanity check
temperature_solar_core = 1.57 * 1e7  # Temperature of the solar core (K)
density_solar_core = 1.62 * 1e5  # Density of the solar core (kg/m^3)

energy_calculator = EnergyProductionCalculator(temperature_solar_core, density_solar_core)
total_energy_production = energy_calculator.calculate_energy_production()

print("Energy Production Rates (J/m^3/s):")
print("r11 H,11 H (Q01 H,1 H + Q01 H,2 D)ρ =", energy_calculator.calculate_energy_ppi())
print("r32 He,32 He Q3 He,3 He ρ =", energy_calculator.calculate_energy_ppii())
print("r32 He,42 He Q3 He,4 He ρ =", energy_calculator.calculate_energy_ppiii())
print("r74 Be,e− Q7 Be,e− ρ =", energy_calculator.calculate_energy_cno())

print("Total Energy Production Rate (J/m^3/s):")
print("Total =", total_energy_production)

def plot_energy_production():
    temperature_range = np.logspace(4, 9, num=100)  # Temperature range from 10^4 to 10^9
    energy_production_ppi = []
    energy_production_ppii = []
    energy_production_ppiii = []
    energy_production_cno = []


    for temperature in temperature_range:
        calculator = EnergyProductionCalculator(temperature, density_solar_core)  # Create an instance of the EnergyProductionCalculator class
        # Calculate energy production for each temperature using the methods from the class
        energy_production_ppi.append(calculator.calculate_energy_ppi())
        energy_production_ppii.append(calculator.calculate_energy_ppii())
        energy_production_ppiii.append(calculator.calculate_energy_ppiii())
        energy_production_cno.append(calculator.calculate_energy_cno())

    # Plotting the results
    plt.figure(figsize=(10, 6))
    plt.plot(temperature_range, energy_production_ppi, label="PPI branch")
    plt.plot(temperature_range, energy_production_ppii, label="PPII branch")
    plt.plot(temperature_range, energy_production_ppiii, label="PPIII branch")
    plt.plot(temperature_range, energy_production_cno, label="CNO cycle")
    plt.xlabel("Temperature (K)")
    plt.ylabel("Relative Energy Production")
    plt.title("Relative Energy Production in Different Fusion Chains")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.grid(True)
    plt.show()

# Call the function to generate the plot
plot_energy_production()

