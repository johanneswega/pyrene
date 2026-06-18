
Dear all,

Some people ask me to share my solvent database program with them. I thought, it might be of interest for you.

In brief, I managed to digitize solvent data (dipole moment, dielectric constant, refractive index, viscosity at 25°C but the program allows also to extract these parameters at a different temperature) of around 260 solvents (all from the same source namely: The properties of solvents by Y. Marcus, 1998).

I can quickly (should take less than 5 min) show how to use the program during our group meeting tomorrow .

If you want to try it yourself here are some written instructions:

    Requirements: python 3 + numpy library

    Usage : 

    1.) Open the terminal and type "python3 solvents.py". (You have to be in the same folder with terminal as the files)

    2.) A prompt will appear asking for the solvent name. If you want to list all solvents starting with a certain letter, e.g., "a", just type the letter and press Enter. All solvents starting with "a" in the database will be printed. Make sure to use the exact spelling of the solvent name in the database for the properties to be listed.

    3.) Type the desired solvent, and the following parameters at 25°C will be printed:
       - Temperature (T) in K
       - Dipole moment (mu) in D
       - Dielectric constant (eps)
       - Refractive index (n)
       - Onsager function (df)
       - Viscosity (eta) in cP
       - Diffusion rate constant (k_diff) in M-1 s-1

    4.) If you want to see the parameters at a different temperature, you can write the following in the prompt: "solvent_name -T xx", where "xx" is the temperature in °C. For example, to see the properties of water at 80°C, you would type: "water -T 80".

    5.) The program also has a search function. Instead of the solvent name, you can write, for instance: "search eps 20". The program will then print all solvents that have a dielectric constant close to eps = 20(±10%). The same applies to the viscosity with the following format: "search eta 2.5" as well as for the Onsager function
    e.g. "search df 0.3" 

Best regards,
Johannes