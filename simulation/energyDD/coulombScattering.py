import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class Material:
    name: str
    density: float      # g/cm^3
    X0: float           # cm
    I: float            # eV
    alpha: float        # range-energy parameter
    p: float            # range-energy exponent
    birks: float = 0.0  # optional


@dataclass
class Particle:
    name: str
    mass: float         # MeV/c^2
    charge: int         # elementary charge units
    energy: float       # MeV

def RangeFromEnergy(energy, material):
    return material.alpha*energy**material.p

def EnergyFromRange(range, material):
    return (range/material.alpha)**(1/material.p)

def ScatteringAngle(material, particle, LayerThickness):
    return 13.6 * particle.charge * np.sqrt(LayerThickness / material.X0) *(particle.mass+particle.energy) / (particle.energy**2+2*particle.mass * particle.energy)*(1+0.038*np.log(LayerThickness/material.X0))

def CumulativeScatteringAngle(material, particle, InitialEnergy, LayerThickness):
    """
    Calculate cumulative scattering angle vs depth.
    Uses constant depth increments through the material.
    """
    # Total range at initial energy
    total_range = RangeFromEnergy(InitialEnergy, material)
    
    # Create constant depth steps
    num_steps = int(total_range / LayerThickness)
    depths = np.linspace(0, total_range, num_steps)
    
    depths_list = []
    energies_list = []
    scatteringAngles = []
    individual_angles = []  # Track individual layer angles too
    cumulativeAngle = []
    total_angle_sq = 0
    
    for i in range(len(depths)-1):  # Step through each layer
        depth_start = depths[i]
        depth_end = depths[i+1]
        
        # Energy at START of this layer (before it scatters)
        remaining_range_start = total_range - depth_start
        
        if remaining_range_start > 0:
            energy_at_start = EnergyFromRange(remaining_range_start, material)
        else:
            break
        
        particle.energy = energy_at_start
        
        # Scattering angle for this layer (constant thickness)
        angle = ScatteringAngle(material, particle, LayerThickness)
        scatteringAngles.append(angle)
        individual_angles.append(angle)
        
        # Accumulate in quadrature
        total_angle_sq += angle**2
        cumulativeAngle.append(np.sqrt(total_angle_sq))
        
        depths_list.append(depth_start)
        energies_list.append(energy_at_start)
    
    return np.array(depths_list), np.array(cumulativeAngle), np.array(energies_list), np.array(individual_angles)


pbwo4 = Material(
    "PbWO4",
    density=8.28,
    X0=0.89,
    I=823,
    alpha=6.741e-3,
    p=1.707
)

h2o = Material(
    "H2O",
    density=1.0,
    X0=36.08,
    I=75,
    alpha=2.369e-2,
    p=1.757
)

proton = Particle(
    "Proton",
    mass=938.2720813,
    charge=1,
    energy=0.0
)

LayerThickness = 0.01  # 0.01 cm depth increments
InitialEnergy = 220.1   # MeV

depths, scatteringAngles, energies_at_depth, individual_angles = CumulativeScatteringAngle(h2o, proton, InitialEnergy, LayerThickness)

sigma_r = depths[-1]/10*(1-np.cos(scatteringAngles[-1])) 
print(f"Depth of material: {depths[-1]:.4f} cm")
print(f"Final cumulative scattering angle: {scatteringAngles[-1]:.4f} mrad")

print(f"Final range: {sigma_r:.4f} cm")

# Create figure with multiple plots
fig, axes = plt.subplots(2, 3, figsize=(15, 8))

# Plot 1: Cumulative angle vs depth
axes[0, 0].plot(depths, scatteringAngles, 'b-', linewidth=2)
axes[0, 0].set_xlabel('Depth in Material [cm]')
axes[0, 0].set_ylabel('Cumulative Scattering Angle [mrad]')
axes[0, 0].set_title('Cumulative Scattering Angle vs Depth')
axes[0, 0].grid(True, alpha=0.3)

# Plot 2: Cumulative angle vs energy (for reference)
axes[0, 1].plot(energies_at_depth, scatteringAngles, 'r-', linewidth=2)
axes[0, 1].set_xlabel('Kinetic Energy [MeV]')
axes[0, 1].set_ylabel('Cumulative Scattering Angle [mrad]')
axes[0, 1].set_title('Cumulative Scattering Angle vs Energy')
axes[0, 1].grid(True, alpha=0.3)

# Plot 3: Energy vs Depth
axes[0, 2].plot(depths, energies_at_depth, 'g-', linewidth=2)
axes[0, 2].set_xlabel('Depth in Material [cm]')
axes[0, 2].set_ylabel('Kinetic Energy [MeV]')
axes[0, 2].set_title('Proton Energy Loss vs Depth')
axes[0, 2].grid(True, alpha=0.3)

# Plot 4: Individual layer scattering angles
axes[1, 0].semilogy(depths, individual_angles, 'purple', linewidth=1.5)
axes[1, 0].set_xlabel('Depth in Material [cm]')
axes[1, 0].set_ylabel('Angle per Layer [mrad] (log scale)')
axes[1, 0].set_title('Individual Layer Scattering Angles')
axes[1, 0].grid(True, alpha=0.3, which='both')

# Plot 5: Individual angles vs energy
axes[1, 1].semilogy(energies_at_depth, individual_angles, 'orange', linewidth=1.5)
axes[1, 1].set_xlabel('Kinetic Energy [MeV]')
axes[1, 1].set_ylabel('Angle per Layer [mrad] (log scale)')
axes[1, 1].set_title('Individual Layer Angles vs Energy')
axes[1, 1].grid(True, alpha=0.3, which='both')

# Plot 6: Derivative of cumulative angle (instantaneous contribution rate)
if len(scatteringAngles) > 1:
    derivative = np.gradient(scatteringAngles, depths)
    axes[1, 2].semilogy(depths, derivative, 'brown', linewidth=1.5)
    axes[1, 2].set_xlabel('Depth in Material [cm]')
    axes[1, 2].set_ylabel('d(Cumulative Angle)/d(Depth)')
    axes[1, 2].set_title('Rate of Angle Change')
    axes[1, 2].grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.show()