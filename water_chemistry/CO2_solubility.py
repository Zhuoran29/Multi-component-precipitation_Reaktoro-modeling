from reaktoro import *


# Initialize a thermodynamic database
db = PhreeqcDatabase("pitzer.dat")
# db = PhreeqcDatabase("minteq.v4.dat")

# Create an aqueous phase automatically selecting all species with provided elements
# aqueousphase = AqueousPhase(speciate("H O C Na Cl"))
# aqueousphase.set(ActivityModelPitzer(db))

aqueousphase = AqueousPhase(speciate("H O C Na Cl Ca Mg K S Si"))
aqueousphase.set(ActivityModelPitzer())

# Create a gaseous phase with CO2(g)
gaseousphase = GaseousPhase("CO2(g)")
gaseousphase.set(ActivityModelPengRobinson())

# Create the chemical system
system = ChemicalSystem(db, aqueousphase, gaseousphase)

# Create the equilibrium solver
solver = EquilibriumSolver(system)

import numpy as np
import pandas as pd

temperatures = np.arange(25.0, 90.0, 5.0)
mgsNaCl  = np.array([ 100.0, 200.0, 300.0])
molsNaCl = np.array([ 0, 1, 2, 4])
P = 1.01325

df = pd.DataFrame(columns=["T", "amountNaCl", "amountCaq"])

for mgNaCl in molsNaCl:
    for T in temperatures:

        # Initial amount of the CO2 gas
        n0CO2g = 10.0

        # Define initial chemical state corresponding to the NaCl-brine of the given concentration
        state = ChemicalState(system)
        state.setTemperature(T, "celsius")
        state.setPressure(P, "bar")
        state.set("H2O"   , 1.0   , "kg")
        state.set("CO2(g)", n0CO2g, "mol")
        # state.set("Na+"   , mgNaCl*22.99/58.44 , "mg")
        # state.set("Cl-"   , mgNaCl*35.45/58.44 , "mg")
        state.set("Na+"   , mgNaCl , "mol")
        state.set("Cl-"   , mgNaCl, "mol")

        # Calculate equilibrium state
        res = solver.solve(state)
        # Stop if the equilibration did not converge or failed
        if not res.succeeded(): continue

        # Fetch resulting aqueous properties of the chemical state
        aqprops = AqueousProps(state)
        props = ChemicalProps(state)

        # Update value ["T", "amountNaCl", "amountCaq"] in the dataframe
        # df.loc[len(df)] = [T, mgNaCl, float(aqprops.elementMolality("C"))]
        df.loc[len(df)] = [T, mgNaCl, float(props.speciesAmount("CO3-2"))]


from reaktplot import * 

fig = Figure()
fig.title("SOLUBILITY OF CO2 IN NACL BRINES")
fig.xaxisTitle('TEMPERATURE [°C]')
fig.yaxisTitle('AMOUNT OF DISSOLVED CO2 [mol/kgw]')

for mgNaCl in molsNaCl:
    df_NaCl = df[df['amountNaCl'] == mgNaCl]
    fig.drawLineWithMarkers(df_NaCl["T"], df_NaCl["amountCaq"], name=f'{mgNaCl} mol of NaCl')

fig.show()