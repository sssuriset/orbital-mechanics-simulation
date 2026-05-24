import os
import time
import csv
import numpy as np
import matplotlib.pyplot as plt

# orbital integrator stability study
# compares euler, semi-implicit euler, and velocity verlet

# constants
G = 6.674e-11
M = 5.972e24
mu = G * M

# starting values
rzero = np.array([7.0e6, 0.0])
vzero = np.array([0.0, 7500.0])

# simulation settings
dtvals = [1.0, 5.0, 10.0, 30.0, 60.0, 120.0, 300.0]
totaltime = 7 * 24 * 60 * 60

methods = ["Euler", "Semi-Implicit Euler", "Velocity Verlet"]

os.makedirs("results", exist_ok=True)


def accel(r):
    dist = np.linalg.norm(r)
    return -mu * r / dist**3


def energy(r, v):
    kin = 0.5 * np.linalg.norm(v)**2
    pot = -mu / np.linalg.norm(r)
    return kin + pot


def angmom(r, v):
    return abs(r[0] * v[1] - r[1] * v[0])


def trueperiod():
    rmag = np.linalg.norm(rzero)
    vmag = np.linalg.norm(vzero)

    en = 0.5 * vmag**2 - mu / rmag
    axis = -mu / (2 * en)

    return 2 * np.pi * np.sqrt(axis**3 / mu)


def findperiod(pos, dt):
    y = pos[:, 1]
    x = pos[:, 0]

    crosses = []

    for i in range(1, len(y)):
        if y[i - 1] < 0 and y[i] >= 0:
            frac = -y[i - 1] / (y[i] - y[i - 1])
            tcross = (i - 1 + frac) * dt
            xcross = x[i - 1] + frac * (x[i] - x[i - 1])

            if xcross > 0:
                crosses.append(tcross)

    if len(crosses) == 0:
        return None

    if len(crosses) == 1:
        return crosses[0]

    periods = np.diff(crosses)
    return np.mean(periods)


def runsim(method, dt):
    steps = int(totaltime / dt)

    r = rzero.copy()
    v = vzero.copy()
    a = accel(r)

    pos = [r.copy()]
    en = [energy(r, v)]
    am = [angmom(r, v)]

    start = time.time()

    for _ in range(steps):
        if method == "Euler":
            aold = accel(r)
            r = r + v * dt
            v = v + aold * dt

        elif method == "Semi-Implicit Euler":
            v = v + accel(r) * dt
            r = r + v * dt

        elif method == "Velocity Verlet":
            rnew = r + v * dt + 0.5 * a * dt**2
            anew = accel(rnew)
            vnew = v + 0.5 * (a + anew) * dt

            r = rnew
            v = vnew
            a = anew

        pos.append(r.copy())
        en.append(energy(r, v))
        am.append(angmom(r, v))

    runtime = time.time() - start

    pos = np.array(pos)
    en = np.array(en)
    am = np.array(am)

    edrift = ((en[-1] - en[0]) / abs(en[0])) * 100
    ldrift = ((am[-1] - am[0]) / abs(am[0])) * 100

    mperiod = findperiod(pos, dt)
    tperiod = trueperiod()

    if mperiod is None:
        perr = None
    else:
        perr = ((mperiod - tperiod) / tperiod) * 100

    return {
        "pos": pos,
        "en": en,
        "am": am,
        "edrift": edrift,
        "ldrift": ldrift,
        "runtime": runtime,
        "mperiod": mperiod,
        "perr": perr
    }


# reference simulation
ref_dt = 0.5
ref = runsim("Velocity Verlet", ref_dt)
ref_final_pos = ref["pos"][-1]


# run simulations
data = {}

for dt in dtvals:
    data[dt] = {}

    for method in methods:
        result = runsim(method, dt)
        result["position_error"] = np.linalg.norm(result["pos"][-1] - ref_final_pos)
        data[dt][method] = result


# print table
print("\nOrbital Integrator Results")
print("Method | dt (s) | Energy Drift (%) | Angular Momentum Drift (%) | Final Position Error (m) | Period Error (%) | Runtime (s)")
print("-" * 145)

for dt in dtvals:
    for method in methods:
        row = data[dt][method]

        ptext = "N/A" if row["perr"] is None else f"{row['perr']:.6f}"

        print(
            f"{method} | {dt} | "
            f"{row['edrift']:.6f} | "
            f"{row['ldrift']:.6f} | "
            f"{row['position_error']:.2f} | "
            f"{ptext} | "
            f"{row['runtime']:.4f}"
        )


# save csv
csvpath = "results/results.csv"

with open(csvpath, "w", newline="") as file:
    writer = csv.writer(file)

    writer.writerow([
        "method",
        "dt",
        "edrift",
        "ldrift",
        "position_error",
        "mperiod",
        "perr",
        "runtime"
    ])

    for dt in dtvals:
        for method in methods:
            row = data[dt][method]

            writer.writerow([
                method,
                dt,
                row["edrift"],
                row["ldrift"],
                row["position_error"],
                row["mperiod"],
                row["perr"],
                row["runtime"]
            ])


# figure 1: orbit comparison
dtplot = 300.0

plt.figure(figsize=(7, 7))

for method in methods:
    pos = data[dtplot][method]["pos"]
    plt.plot(pos[:, 0], pos[:, 1], label=method)

earth = plt.Circle((0, 0), 6.371e6, color="blue", alpha=0.3)
plt.gca().add_patch(earth)
plt.xlabel("x position (m)")
plt.ylabel("y position (m)")
plt.title(f"Orbital Trajectory Comparison, dt = {dtplot} s")
plt.legend()
plt.axis("equal")
plt.grid()
plt.savefig("results/orbitcomparison.png", dpi=300)


# figure 2: energy conservation
plt.figure(figsize=(8, 5))

for method in methods:
    en = data[dtplot][method]["en"]
    echange = ((en - en[0]) / abs(en[0])) * 100
    t = np.arange(len(en)) * dtplot / 3600

    plt.plot(t, echange, label=method)

plt.xlabel("Time (hours)")
plt.ylabel("Energy Change (%)")
plt.title(f"Energy Conservation Comparison, dt = {dtplot} s")
plt.legend()
plt.grid()
plt.savefig("results/energycomparison.png", dpi=300)


# figure 3: angular momentum conservation
plt.figure(figsize=(8, 5))

for method in methods:
    am = data[dtplot][method]["am"]
    lchange = ((am - am[0]) / abs(am[0])) * 100
    t = np.arange(len(am)) * dtplot / 3600

    plt.plot(t, lchange, label=method)

plt.xlabel("Time (hours)")
plt.ylabel("Angular Momentum Change (%)")
plt.title(f"Angular Momentum Conservation Comparison, dt = {dtplot} s")
plt.legend()
plt.grid()
plt.savefig("results/angmomcomparison.png", dpi=300)


# figure 4: timestep sensitivity
plt.figure(figsize=(8, 5))

for method in methods:
    driftvals = [abs(data[dt][method]["edrift"]) for dt in dtvals]
    plt.plot(dtvals, driftvals, marker="o", label=method)

plt.xlabel("Timestep, dt (s)")
plt.ylabel("Absolute Energy Drift (%)")
plt.title("Timestep Sensitivity of Numerical Integrators")
plt.yscale("log")
plt.legend()
plt.grid()
plt.savefig("results/timestepsensitivity.png", dpi=300)


# figure 5: period error
plt.figure(figsize=(8, 5))

for method in methods:
    perrs = []

    for dt in dtvals:
        err = data[dt][method]["perr"]
        perrs.append(np.nan if err is None else abs(err))

    plt.plot(dtvals, perrs, marker="o", label=method)

plt.xlabel("Timestep, dt (s)")
plt.ylabel("Absolute Orbital Period Error (%)")
plt.title("Orbital Period Error by Timestep")
plt.legend()
plt.grid()
plt.savefig("results/perioderror.png", dpi=300)


# figure 6: final position error
plt.figure(figsize=(8, 5))

for method in methods:
    errors = [data[dt][method]["position_error"] for dt in dtvals]
    plt.plot(dtvals, errors, marker="o", label=method)

plt.xlabel("Timestep, dt (s)")
plt.ylabel("Final Position Error vs Reference (m)")
plt.title("Final Position Error Compared with High-Accuracy Reference")
plt.yscale("log")
plt.legend()
plt.grid()
plt.savefig("results/positionerror.png", dpi=300)

plt.show()