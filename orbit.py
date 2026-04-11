import numpy as np
import matplotlib.pyplot as plt

# constants
G = 6.674e-11
M = 5.972e24  # Earth mass (kg)

# initial conditions
r = np.array([7e6, 0.0])
v = np.array([0.0, 7500.0])

dt = 0.1
steps = 100000

positions = []
energies = []

dist = np.linalg.norm(r)
a = -G * M * r / dist**3

for _ in range(steps):
    r = r + v * dt + 0.5 * a * dt**2

    dist_new = np.linalg.norm(r)
    a_new = -G * M * r / dist_new**3

    v = v + 0.5 * (a + a_new) * dt
    a = a_new

    positions.append(r.copy())

    kinetic = 0.5 * np.linalg.norm(v)**2
    potential = -G * M / dist_new
    total_energy = kinetic + potential
    energies.append(total_energy)

positions = np.array(positions)

plt.figure(figsize=(6,6))
plt.plot(positions[:,0], positions[:,1])
plt.scatter(0, 0, label="Earth")
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.title("Orbital Motion Simulation")
plt.legend()
plt.axis("equal")
plt.grid()

plt.figure()
plt.plot(energies)
plt.xlabel("Time step")
plt.ylabel("Total Energy")
plt.title("Energy Conservation Check")
plt.grid()
plt.savefig("orbit_trajectory.png", dpi=300)
plt.savefig("energy_plot.png", dpi=300)

plt.show()
