def main():
    import numpy as np

    eps0 = 8.8541878128e-12  # F/m
    c = 299792458.0          # m/s

    print("\nLaser quantity examples: 10 uJ, 2 mW, 5 uJ/cm2")
    inp = input(
        "Laser quantity [average power / energy per pulse / fluence]: "
    ).strip().replace("µ", "u").replace("μ", "u")

    # Beam geometry: diameter in um, area in cm^2
    d_um = float(input("Beam diameter [um]: "))
    if d_um <= 0:
        raise ValueError("Beam diameter must be positive.")

    d_cm = d_um * 1e-4
    area_cm2 = np.pi * (d_cm / 2)**2
    area_m2 = area_cm2 * 1e-4

    # Repetition rate
    rep_khz = float(input("Rep-rate [kHz]: "))
    if rep_khz <= 0:
        raise ValueError("Rep-rate must be positive.")
    rep_hz = rep_khz * 1e3

    # Parse numerical value and units
    import re
    match = re.fullmatch(
        r"\s*([0-9]*\.?[0-9]+(?:[eE][+-]?\d+)?)\s*([a-zA-Z0-9/^]+)\s*",
        inp
    )
    if not match:
        raise ValueError(
            "Enter a value such as 10uJ, 2mW, or 5uJ/cm2."
        )

    value = float(match.group(1))
    unit = match.group(2).lower()

    # Convert input to pulse energy E [J]
    if "/cm2" in unit:
        fluence_units = {
            "uj/cm2": 1e-6,
            "mj/cm2": 1e-3,
            "j/cm2": 1.0,
        }
        if unit not in fluence_units:
            raise ValueError("Fluence units: uJ/cm2, mJ/cm2, or J/cm2.")

        fluence_J_cm2 = value * fluence_units[unit]
        E = fluence_J_cm2 * area_cm2

    elif unit.endswith("w"):
        power_units = {"uw": 1e-6, "mw": 1e-3, "w": 1.0}
        if unit not in power_units:
            raise ValueError("Power units: uW, mW, or W.")

        P_avg = value * power_units[unit]
        E = P_avg / rep_hz

    elif unit.endswith("j"):
        energy_units = {"uj": 1e-6, "mj": 1e-3, "j": 1.0}
        if unit not in energy_units:
            raise ValueError("Energy units: uJ, mJ, or J.")

        E = value * energy_units[unit]

    else:
        raise ValueError("Unrecognized units.")

    if E <= 0:
        raise ValueError("Input energy must be positive.")

    # Average power and fluence
    P_avg = E * rep_hz
    fluence = E / area_cm2  # J/cm^2

     # Pulse duration [fs]
    tau_fs = float(input("Pulse duration [fs]: "))
    if tau_fs <= 0:
        raise ValueError("Pulse duration must be positive.")

    # Convert fs to seconds
    tau_s = tau_fs * 1e-15

    # Peak power [W], using E/t
    P_peak = E / tau_s

    # Peak irradiance [W/cm^2]
    I_peak = E / (tau_s * area_cm2)

    # Convert irradiance to W/m^2 for electric-field calculation
    I_SI = I_peak * 1e4
    E_field = np.sqrt(2 * I_SI / (eps0 * c))

    print("\n--- Results ---")
    print(f"Pulse energy                      = {E:.4g} J = {E*1e6:.4g} uJ")
    print(f"Average power                     = {P_avg:.4g} W = {P_avg*1e3:.4g} mW")
    print(f"Rep-rate                          = {rep_khz:.4g} kHz")
    print(f"Beam area (circle)                = {area_cm2:.4g} cm^2")
    print(f"Fluence                           = {fluence*1e3:.4g} mJ/cm^2 = {fluence*1e6:.4g} uJ/cm^2")
    print("Pulse duration                    = %.3g fs" % tau_fs)
    print("Peak power                        = %.3g W" % P_peak)
    print("Peak irradiance (flat-top)        = %.3g W/cm^2 (%.3g V/m)"%(I_peak, E_field))
    print("")