#%%
from reaktoro import *
from reaktplot import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
#%% INSPECTING DATABASE
db_thermofun = ThermoFunDatabase("cemdata18")
db_phreeqc_pitzer = PhreeqcDatabase("pitzer.dat")
db_phreeqc_frezchem = PhreeqcDatabase("frezchem.dat")
db_phreeqc_minteq = PhreeqcDatabase("minteq.v4.dat")

print(f"{'Species':<20}{'Formula':<60}{'Molar Mass (kg/mol)':<20}")
for species in db_phreeqc_minteq.species():
    print(f"{species.name():<20}{species.formula().str():<60}{species.molarMass():<20.6f}")


# %% GET SPECIES PROPERTIES
hydromagnesite = db_phreeqc_minteq.species().get('Hydromagnesite')
hydromagnesite.props(298.15, 101325)
hydromagnesite.elements()


#%%  SOLUBILITY of HYDROMAGNESITE
import pandas as pd
# Create the database

# db = SupcrtDatabase("supcrt07")
db = PhreeqcDatabase("thermoddem-v1.10.dat")
# db = PhreeqcDatabase("minteq.v4.dat")
# db = PhreeqcDatabase("pitzer.dat")

# Create an aqueous phase automatically selecting all species with provided elements
aqueousphase = AqueousPhase(speciate("H O C Ca Mg K Cl Na S N"))
aqueousphase.set(ActivityModelPitzer())

# Create a gaseous phase
gaseousphase = GaseousPhase("CO2(g)")
gaseousphase.set(ActivityModelPengRobinson())

# Create a mineral phase
mineral = MineralPhase("Hydromagnesite")
mineral3 = MineralPhase("Magnesite(Synth)")
mineral2 = MineralPhase("Calcite")

# Create the chemical system
system = ChemicalSystem(db, aqueousphase, gaseousphase, mineral, mineral2, mineral3)

# Define equilibrium specs
specs = EquilibriumSpecs (system)
specs.temperature()
specs.pressure()
specs.charge()
specs.openTo("Cl-")

# Define conditions to be satisfied at the chemical equilibrium state
conditions = EquilibriumConditions(specs)
conditions.charge(0.0, "mol") # to make sure the mixture is charge neutral

# Define the equilibrium solver
solver = EquilibriumSolver(specs)

# Define aqueous properties
aprops = AqueousProps(system)

def water():

    state = ChemicalState(system)
    state.add("H2O(aq)", 1.0, "kg")

    return state

def seawater(CO2=0, NaOH=0):

    state = ChemicalState(system)
    # Seawater composition
    state.setTemperature(25, "celsius")
    state.setPressure(1.0, "bar")
    state.add("H2O",         1.0, "kg")
    state.add("Ca+2"   ,   412.3, "mg")
    state.add("Mg+2"   ,  1290.0, "mg")
    state.add("Na+"    , 10768.0, "mg")
    state.add("K+"     ,   399.1, "mg")
    state.add("Cl-"    , 19353.0, "mg")
    state.add("HCO3-"  ,   141.7, "mg")
    state.add("SO4-2"  ,  2712.0, "mg")


    # state.add("CO2(g)", 0.36, "mol")  
    # state.add("Na+", 0.66, "mol")  
    # state.add("OH-", 0.66, "mol")  

    return state


df_seawater = pd.DataFrame(columns=["CO2", "NaOH", "pH",  "Hydromagnesite", "Magnesite"])

# Initial amount of calcite
n0Calcite = 0
n0Hydromagnesite = 0
n0Magnesite = 0

def solubility_of_calcite( T=25, P=1.01325, CO2 = 0, NaOH = 0,state=None, tag=None):


    state.add("CO2(g)", CO2, "mol")  
    state.add("Na+", NaOH, "mol")  
    state.add("OH-", NaOH, "mol")  

    conditions.temperature(T, "celsius")
    conditions.pressure(P, "bar")

    # Equilibrate the solution with given initial chemical state and desired conditions at equilibrium
    try:
        res = solver.solve(state, conditions)
    except:
        df_seawater.loc[len(df_seawater)] = [CO2, NaOH, 'fail', 'fail', 'fail']

    # Check calculation succeeded
    assert res.succeeded()

    # Update aqueous properties
    aprops.update(state)

    # print('pH = ', float(aprops.pH()), 'T = ', T)
    # Fetch the amount of final calcite in the equilibrium state
    nCalcite = float(state.speciesAmount("Calcite"))
    nHydromagnesite = float(state.speciesAmount("Hydromagnesite"))
    nMagnesite = float(state.speciesAmount("Magnesite(Synth)"))

    if tag == "seawater":
        df_seawater.loc[len(df_seawater)] = [CO2, NaOH, float(aprops.pH()), nHydromagnesite, nMagnesite]

    else:
        return

import numpy as np
temperatures = np.arange(20.0, 91.0, 5.0)   # in celsius
pressures = np.array([1, 10, 100])          # in bar


CO2s = np.arange(0,5,0.25)   # in celsius
NaOHs = np.arange(0,5,0.25)        

state = seawater()
state.set("Magnesite(Synth)", n0Magnesite, "mol")
state.set("Hydromagnesite", n0Hydromagnesite, "mol")
# [solubility_of_calcite(state, T, P, "seawater") for P in pressures for T in temperatures]
[solubility_of_calcite(state=state, CO2=CO2, NaOH=NaOH, tag="seawater") for CO2 in CO2s for NaOH in NaOHs]

from reaktplot import *

df_seawater_P1 = df_seawater[df_seawater['P'] == 1.0]

fig = Figure()

fig.xaxisTitle('Temperature [°C]')
fig.yaxisTitle('Hydromagnesite precipitation [mol/kgw]')

fig.drawLineWithMarkers(df_seawater_P1["T"], df_seawater_P1["Hydromagnesite"], name="seawater")
fig.drawLineWithMarkers(df_seawater_P1["T"], df_seawater_P1["Magnesite"], name="seawater")

fig.show()



# %% Solubility of all species

from pathlib import Path
from reaktoro import (PhreeqcDatabase,
                      AqueousPhase, 
                      ActivityModelPitzer, 
                      ActivityModelPengRobinson,
                      speciate, 
                      GaseousPhase,
                      ChemicalSystem,
                      EquilibriumSpecs,
                      EquilibriumConditions, 
                      EquilibriumSolver,
                      AqueousProps, 
                      ChemicalState,
                      )

db = PhreeqcDatabase("thermoddem-v1.10.dat")
# db = PhreeqcDatabase("minteq.v4.dat")

# Create an aqueous phase automatically selecting all species with provided elements
aqueousphase = AqueousPhase(speciate("H O C Ca Mg K Cl Na S N"))
aqueousphase.set(ActivityModelPitzer())

# Create a gaseous phase
gaseousphase = GaseousPhase("CO2(g)")
gaseousphase.set(ActivityModelPengRobinson())

# Create a mineral phase
mineral = MineralPhase("Hydromagnesite Calcite Polyhalite Glauberite Halite Gypsum Anhydrite Magnesite(Synth)")

Hydromagnesite = MineralPhase("Magnesite(Synth)")
Magnesite = MineralPhase("Magnesite(Synth)")


# Create the chemical system
system = ChemicalSystem(db, aqueousphase, gaseousphase, mineral)

# Define equilibrium specs
specs = EquilibriumSpecs (system)
specs.temperature()
specs.pressure()
specs.charge()
specs.openTo("Cl-")

# Define conditions to be satisfied at the chemical equilibrium state
conditions = EquilibriumConditions(specs)
conditions.charge(0.0, "mol") # to make sure the mixture is charge neutral

# Define the equilibrium solver
solver = EquilibriumSolver(specs)

# Define aqueous properties
aprops = AqueousProps(system)

def seawater(H2O=1):

    state = ChemicalState(system)
    # Seawater composition
    state.setTemperature(25, "celsius")
    state.setPressure(1.0, "bar")
    state.add("H2O",     H2O, "kg")
    state.add("Ca+2"   ,   412.3, "mg")
    state.add("Mg+2"   ,  1290.0, "mg")
    state.add("Na+"    , 10768.0, "mg")
    state.add("K+"     ,   399.1, "mg")
    state.add("Cl-"    , 19353.0, "mg")
    state.add("HCO3-"  ,   141.7, "mg")
    state.add("SO4-2"  ,  2712.0, "mg")


    state.add("CO2(g)", 0.36, "mol")  
    state.add("Na+", 0.66, "mol")  
    state.add("OH-", 0.66, "mol")  

    return state


df_seawater = pd.DataFrame(columns=["T", "H2O", "pH",  
                                    "Calcite",
                                    "Magnesite",
                                    "Hydromagnesite",
                                    "Polyhalite", 
                                    "Glauberite", 
                                    "Halite", 
                                    "Gypsum", 
                                    "Anhydrite",
                                    # "dissolved MgCO3"
                                    ])

def solubility(T, H2O, tag):
    # Initial species amount
    n0Hydromagnesite = 0
    n0Polyhalite = 0
    n0Glauberite = 0
    n0Calcite = 0
    n0Magnesite = 0
    n0Halite = 0
    n0Gypsum = 0
    n0Anhydrite = 0

    state = seawater(H2O = H2O)
    state.set("Hydromagnesite", n0Hydromagnesite, "mol")
    state.set("Polyhalite", n0Polyhalite, "mol")
    state.set("Glauberite", n0Glauberite, "mol")
    state.set("Calcite", n0Calcite, "mol")
    state.set("Magnesite(Synth)", n0Magnesite, "mol")
    state.set("Halite", n0Halite, "mol")
    state.set("Gypsum", n0Gypsum, "mol")
    state.set("Anhydrite", n0Anhydrite, "mol")

    conditions.temperature(T, "celsius")
    conditions.pressure(1, "bar")

    # Equilibrate the solution with given initial chemical state and desired conditions at equilibrium
    res = solver.solve(state, conditions)

    # Check calculation succeeded
    assert res.succeeded()

    # Update aqueous properties
    aprops.update(state)

    # Fetch the amount of final calcite in the equilibrium state
    nCalcite = float(state.speciesAmount("Calcite"))
    nMagnesite = float(state.speciesAmount("Magnesite(Synth)"))
    nHydromagnesite = float(state.speciesAmount("Hydromagnesite"))
    nPolyhalite = float(state.speciesAmount("Polyhalite"))
    nGlauberite = float(state.speciesAmount("Glauberite"))
    nHalite = float(state.speciesAmount("Halite"))
    nGypsum = float(state.speciesAmount("Gypsum"))
    nAnhydrite = float(state.speciesAmount("Anhydrite"))

    # nMgCO3 = float(state.speciesAmount("MgCO3"))

    if tag == "seawater":
        df_seawater.loc[len(df_seawater)] = [T, H2O, float(aprops.pH()), 
                                             nCalcite,
                                             nMagnesite,
                                             nHydromagnesite, 
                                             nPolyhalite, 
                                             nGlauberite,
                                             nHalite,
                                             nGypsum,
                                             nAnhydrite,
                                            #  nMgCO3
                                             ]
        
    else:
        return
    
# Sweeping T and concentration

temperatures = np.arange(20.0, 91.0, 5.0)
H2Os = np.arange(1, 0.01, -0.05)

for H2O in H2Os:
    try:
        solubility(25, H2O, "seawater")
    except:
        print(H2O)
    
# [solubility(state, T, H2O, "seawater") for H2O in H2Os for T in temperatures]
display(df_seawater)
# filepath = Path('/Users/zhuoranzhang/Documents/Crystallization_paper/reaktoro_seawater.csv')  
# df_seawater.to_csv(filepath)

# %% Figure
# fig = Figure()
# fig.xaxisTitle('H2O [kg]')
# fig.yaxisTitle('Precipitation [mol/kgw]')
# fig.drawLineWithMarkers(df_water_P1["T"], df_water_P1["Calcite"], name="water")
# fig.drawLineWithMarkers(df_rainwater_P1["T"], df_rainwater_P1["Calcite"], name="rainwater")
# fig.drawLine(df_seawater["P"], df_seawater["Magnesite"], name="Magnesite")
# fig.drawLine(df_seawater["P"], df_seawater["Halite"], name="Halite")


# fig.show()

# %% Solubility of all species
import pandas as pd
from pathlib import Path
from reaktoro import (PhreeqcDatabase,
                      AqueousPhase, 
                      ActivityModelPitzer, 
                      ActivityModelPengRobinson,
                      speciate, 
                      GaseousPhase,
                      ChemicalSystem,
                      EquilibriumSpecs,
                      EquilibriumConditions, 
                      EquilibriumSolver,
                      AqueousProps, 
                      ChemicalState,
                      MineralPhase,
                      ChemicalProps
                      )

db = PhreeqcDatabase("thermoddem-v1.10.dat")
# db = PhreeqcDatabase("minteq.v4.dat")
# db = PhreeqcDatabase("pitzer.dat")

# Create an aqueous phase automatically selecting all species with provided elements
aqueousphase = AqueousPhase(speciate("H O C Ca Mg K Cl Na S N"))
aqueousphase.set(ActivityModelPitzer())

# Create a gaseous phase
gaseousphase = GaseousPhase("CO2(g)")
gaseousphase.set(ActivityModelPengRobinson())

# Create a mineral phase
mineral = MineralPhase("Hydromagnesite Calcite Polyhalite Glauberite Halite Gypsum Anhydrite Magnesite(Synth)")

Hydromagnesite = MineralPhase("Hydromagnesite")
Magnesite = MineralPhase("Magnesite(Synth)")
Calcite = MineralPhase("Calcite")
Polyhalite = MineralPhase("Polyhalite")
Glauberite = MineralPhase("Glauberite")
Halite = MineralPhase("Halite")
Gypsum = MineralPhase("Gypsum")
Anhydrite = MineralPhase("Anhydrite")


# Create the chemical system
# system = ChemicalSystem(db, aqueousphase, gaseousphase, mineral)
system = ChemicalSystem(db, 
                        aqueousphase, 
                        # gaseousphase, 
                        Calcite,
                        # Hydromagnesite, 
                        Magnesite, 
                        Polyhalite,
                        Glauberite,
                        Halite,
                        Anhydrite,
                        )

# Define equilibrium specs
specs = EquilibriumSpecs (system)
specs.temperature()
specs.pressure()
# specs.pH()
specs.charge()
specs.openTo("Cl-")

# Define conditions to be satisfied at the chemical equilibrium state
conditions = EquilibriumConditions(specs)
conditions.charge(0.0, "mol") # to make sure the mixture is charge neutral

# Define the equilibrium solver
solver = EquilibriumSolver(specs)

# Define aqueous properties
aprops = AqueousProps(system)

def seawater(H2O=1):

    state = ChemicalState(system)
    # Seawater composition
    state.setTemperature(118, "celsius")
    state.setPressure(1.0, "bar")
    state.add("H2O",         H2O, "kg")
    # state.add("Ca+2"   ,   21100, "mg")
    # state.add("Mg+2"   ,   65940, "mg")
    # state.add("Na+"    ,  130680, "mg")
    # state.add("K+"     ,   20240, "mg")
    # state.add("Cl-"    ,  351138, "mg")
    # state.add("HCO3-"  ,   7490, "mg")
    # state.add("SO4-2"  ,  137250, "mg")

    state.add("Ca+2"   ,   412.3, "mg")
    state.add("Mg+2"   ,  1290.0, "mg")
    state.add("Na+"    , 10768.0, "mg")
    state.add("K+"     ,   399.1, "mg")
    state.add("Cl-"    , 19353.0, "mg")
    state.add("HCO3-"  ,   141.7, "mg")
    state.add("SO4-2"  ,  2712.0, "mg")
 
    # state.add("CO2(g)", 6.9, "mol")  
    # state.add("Na+", 10, "mol")  
    # state.add("OH-", 10, "mol")  

    return state


df_seawater = pd.DataFrame(columns=["T", "H2O", "pH",  
                                    "Calcite",
                                    "Magnesite",
                                    # "Hydromagnesite",
                                    "Polyhalite", 
                                    "Glauberite", 
                                    "Halite", 
                                    # "Gypsum", 
                                    "Anhydrite",
                                    # "dissolved MgCO3"
                                    ])

def solubility(T, H2O, tag):
    # Initial species amount
    n0Hydromagnesite = 0
    n0Polyhalite = 0
    n0Glauberite = 0
    n0Calcite = 0
    n0Magnesite = 0
    n0Halite = 0
    n0Gypsum = 0
    n0Anhydrite = 0

    state = seawater(H2O = H2O)

    state.set("Calcite", n0Calcite, "mol")
    state.set("Magnesite(Synth)", n0Magnesite, "mol")
    # state.set("Hydromagnesite", n0Hydromagnesite, "mol")
    state.set("Polyhalite", n0Polyhalite, "mol")
    state.set("Glauberite", n0Glauberite, "mol")
    state.set("Halite", n0Halite, "mol")
    # state.set("Gypsum", n0Gypsum, "mol")
    state.set("Anhydrite", n0Anhydrite, "mol")

    props = ChemicalProps(state)

    conditions.temperature(T, "celsius")
    conditions.pressure(1, "bar")
    # conditions.pH(8.22)

    # Equilibrate the solution with given initial chemical state and desired conditions at equilibrium
    res = solver.solve(state, conditions)

    # Check calculation succeeded
    assert res.succeeded()

    # Update aqueous properties
    aprops.update(state)

    print(props)
    # print(aprops)
    # print(state)
    print('pH = ', float(aprops.pH()))

    # Fetch the amount of final calcite in the equilibrium state
    nCalcite = float(state.speciesAmount("Calcite"))
    nMagnesite = float(state.speciesAmount("Magnesite(Synth)"))
    # nHydromagnesite = float(state.speciesAmount("Hydromagnesite"))
    nPolyhalite = float(state.speciesAmount("Polyhalite"))
    nGlauberite = float(state.speciesAmount("Glauberite"))
    nHalite = float(state.speciesAmount("Halite"))
    # nGypsum = float(state.speciesAmount("Gypsum"))
    nAnhydrite = float(state.speciesAmount("Anhydrite"))

    # nMgCO3 = float(state.speciesAmount("MgCO3"))

    if tag == "seawater":
        df_seawater.loc[len(df_seawater)] = [T, H2O, float(aprops.pH()), 
                                             nCalcite,
                                             nMagnesite,
                                            #  nHydromagnesite, 
                                             nPolyhalite, 
                                             nGlauberite,
                                             nHalite,
                                            #  nGypsum,
                                             nAnhydrite,
                                            #  nMgCO3
                                             ]
        return props

    else:
        return
    
props = solubility(25, 1, "seawater")





# %% Solubility at defined state
from reaktoro import * 

# db = PhreeqcDatabase("minteq.v4.dat")
# db = PhreeqcDatabase("thermoddem-v1.10.dat")
db = PhreeqcDatabase("pitzer.dat")
# db = SupcrtDatabase("supcrt07")

solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S"))
solution.set(ActivityModelPitzer())

# Create a gaseous phase
gaseousphase = GaseousPhase("H2O(g)")
gaseousphase.set(ActivityModelPengRobinson())

# Hydromagnesite = MineralPhase("Hydromagnesite")
Magnesite = MineralPhase("Magnesite")
Calcite = MineralPhase("Calcite")
Polyhalite = MineralPhase("Polyhalite")
Glauberite = MineralPhase("Glauberite")
Halite = MineralPhase("Halite")
Gypsum = MineralPhase("Gypsum")
Anhydrite = MineralPhase("Anhydrite")
Epsomite = MineralPhase("Epsomite")

system = ChemicalSystem(db, 
                        gaseousphase,
                        solution,
                        Magnesite,
                        Calcite,
                        Polyhalite,
                        Glauberite,
                        Halite,
                        Anhydrite,
                        Gypsum,
                        Epsomite
                        )

specs = EquilibriumSpecs(system)
# specs.temperature()
specs.pressure()
specs.phaseAmount("GaseousPhase")
# specs.pH()
specs.charge()
specs.openTo("Cl-")

# solver = EquilibriumSolver(specs)
solver = SmartEquilibriumSolver(specs)
 

state = ChemicalState(system)
state.temperature(90, "celsius")
# state.pressure(1.0, "atm")
state.set("H2O", 1.0, "kg")
# state.add("Ca+2"   ,   2670, "mg")
# state.add("Mg+2"   ,  8330, "mg")
# state.add("Na+"    ,  91740, "mg")
# state.add("K+"     ,   2560, "mg")
# state.add("Cl-"    , 14160, "mg")
# state.add("HCO3-"  ,   141.7, "mg")
# state.add("SO4-2"  ,  17350, "mg")

state.add("Cl-"    ,  351138, "mg")
state.add("Na+"    ,  130680, "mg")
state.add("Mg+2"   ,   65940, "mg")
state.add("K+"     ,   20240, "mg")
state.add("SO4-2"  ,  137250, "mg")
state.add("Ca+2"   ,   21100, "mg")
# state.add("HCO3-"  ,   7490, "mg")

# state.add("H2O(g)"  ,  1, "mol")

# equilibrate(state)

conditions = EquilibriumConditions(specs)
# conditions.temperature(118, "celsius")
conditions.pressure(0.99, "bar")

conditions.phaseAmount("GaseousPhase", 0.0001, "umol")
# conditions.setLowerBoundPressure(0.01, "bar")
# conditions.setUpperBoundPressure(1000.0, "bar")
conditions.setLowerBoundTemperature(25, "celsius")
conditions.setUpperBoundTemperature(200, "celsius")
# conditions.pH(2.0)
conditions.charge(0.0)

result = solver.solve(state, conditions)


# print(state)
 
props = ChemicalProps(state)
aprops = AqueousProps(system)
aprops.update(state)

# print(props)
print(state)
print('density = ', props.phaseProps("AqueousPhase").density())
print('pH = ', float(aprops.pH()))
# print(props.phaseProps("AqueousPhase").speciesAmounts())
props = ChemicalProps(state)
print(props)

#%% Useful functions


# Print saturation ratio
for i in aprops.saturationSpecies():
    print(i.name())
aprops.saturationRatios()

#%% Fix boiling temperature
from reaktoro import *

db = PhreeqcDatabase("pitzer.dat")

solution = AqueousPhase(speciate("H O"))
# solution = AqueousPhase(speciate("H O Na Cl C"))
solution.set(ActivityModelPitzer())

gases = GaseousPhase("H2O(g)")
# gases = GaseousPhase("CO2(g) H2O(g)")
gases.set(ActivityModelPengRobinsonPhreeqcOriginal())

system = ChemicalSystem(db, solution, gases)


# Initial state
state = ChemicalState(system)
state.temperature(100.0, "celsius")
# state.pressure(1, 'bar')
state.set("H2O", 1.0, "kg")
# state.set("Na+", 1.0, "mol")
# state.set("Cl-", 1.0, "mol")
# state.set("CO2", 1.0, "mol")

# print("INITIAL STATE")
# print(state)


# Specify sonsrtaints
specs = EquilibriumSpecs(system)
specs.temperature()
# specs.pressure()
specs.phaseAmount("GaseousPhase")

solver = EquilibriumSolver(specs)

# Make constraint
conditions = EquilibriumConditions(specs)
conditions.temperature(100, "celsius")
# conditions.pressure(1, 'bar')
conditions.phaseAmount("GaseousPhase", 1.0, "umol")  # umol = 1e-6 moles
conditions.setLowerBoundPressure(0.1, "bar")
conditions.setUpperBoundPressure(1000.0, "bar")
# conditions.setLowerBoundTemperature(50, "celsius")
# conditions.setUpperBoundTemperature(200, "celsius")

solver.solve(state, conditions)

print("FINAL STATE")
print(state)


#%% Fix boiling pressure

from reaktoro import *

db = PhreeqcDatabase("pitzer.dat")

solution = AqueousPhase(speciate("H O"))
# solution = AqueousPhase(speciate("H O Na Cl C"))
solution.set(ActivityModelPitzer())

gases = GaseousPhase("H2O(g)")
# gases = GaseousPhase("CO2(g) H2O(g)")
gases.set(ActivityModelPengRobinsonPhreeqcOriginal())

system = ChemicalSystem(db, solution, gases)


# Initial state
state = ChemicalState(system)
# state.temperature(100.0, "celsius")
state.pressure(1, 'bar')
state.set("H2O", 1.0, "kg")
# state.set("Na+", 1.0, "mol")
# state.set("Cl-", 1.0, "mol")
# state.set("CO2", 1.0, "mol")

# print("INITIAL STATE")
# print(state)


# Specify sonsrtaints
specs = EquilibriumSpecs(system)
# specs.temperature()
specs.pressure()
specs.phaseAmount("GaseousPhase")

solver = EquilibriumSolver(specs)

# Make constraint
conditions = EquilibriumConditions(specs)
# conditions.temperature(100, "celsius")
conditions.pressure(1, 'bar')
conditions.phaseAmount("GaseousPhase", 0.001, "umol")  # umol = 1e-6 moles
# conditions.setLowerBoundPressure(0.1, "bar")
# conditions.setUpperBoundPressure(1000.0, "bar")
conditions.setLowerBoundTemperature(25, "celsius")
conditions.setUpperBoundTemperature(150, "celsius")

solver.solve(state, conditions)

print("FINAL STATE")
print(state)
#%% Kinetics
from reaktoro import *

db = SupcrtDatabase("supcrtbl")  # let's use the SUPCRTBL database for this simulation

params = Params.embedded("PalandriKharaka.yaml")

system = ChemicalSystem(db,
    AqueousPhase().set(ActivityModelPitzer()),
    GaseousPhase("CO2(g)").set(ActivityModelPengRobinsonPhreeqcOriginal()),
    MineralPhases("Calcite Magnesite Dolomite"),
    MineralReaction("Calcite").setRateModel(ReactionRateModelPalandriKharaka(params)),
    MineralReaction("Magnesite").setRateModel(ReactionRateModelPalandriKharaka(params)),
    MineralReaction("Dolomite").setRateModel(ReactionRateModelPalandriKharaka(params)),
    MineralSurface("Calcite", 6.0, "cm2/cm3"),
    MineralSurface("Magnesite", 6.0, "cm2/cm3"),
    MineralSurface("Dolomite", 6.0, "cm2/cm3"),
) 

state = ChemicalState(system)
state.set("H2O(aq)"  , 1.0, "kg")
state.set("CO2(g)"   , 5.0, "mol")
state.set("Calcite"  , 1.0, "mol")
state.set("Magnesite", 1.0, "mol")
state.set("O2(aq)"   , 1.0, "umol")  # Important: when O2(aq) or H2(aq) are in the system, add an insignificant tiny amount of one of them to avoid numerical problems due to floating-point rounding errors! 
#%% Parametric at equilibrium state
def seawater(desal_rr = 0.8,
             pressure = 0.99, # bar
             Na2CO3 = 0, # mg
             CO2 = 0, # mg
             temp = 35 + 273.15, 
    ): 
    db = PhreeqcDatabase("pitzer.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S"))
    solution.set(ActivityModelPitzer())

    # Create a gaseous phase
    # gaseousphase = GaseousPhase("H2O(g)")
    # gaseousphase.set(ActivityModelPengRobinson())

    # Hydromagnesite = MineralPhase("Hydromagnesite")
    Magnesite = MineralPhase("Magnesite")
    Calcite = MineralPhase("Calcite")
    Polyhalite = MineralPhase("Polyhalite")
    Glauberite = MineralPhase("Glauberite")
    Halite = MineralPhase("Halite")
    Gypsum = MineralPhase("Gypsum")
    Anhydrite = MineralPhase("Anhydrite")
    Epsomite = MineralPhase("Epsomite")

    system = ChemicalSystem(db, 
                            # gaseousphase,
                            solution,
                            Magnesite,
                            Calcite,
                            Polyhalite,
                            Glauberite,
                            Halite,
                            Anhydrite,
                            Gypsum,
                            Epsomite
                            )

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    # specs.phaseAmount("GaseousPhase")
    # specs.pH()
    specs.charge()
    specs.openTo("Cl-")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)
    

    state = ChemicalState(system)
    state.temperature(90, "celsius")

    # Seawater composition (mg/kg)
    Na  = 10752
    Cl  = 19245
    SO4 = 2701
    Mg = 1295
    Ca = 416
    K  = 390
    HCO3= 145

    cf = 1 / (1-desal_rr)

    # state.pressure(1.0, "atm")
    state.set("H2O", 1.0, "kg")
    state.add("Ca+2"   ,   Ca * cf, "mg")
    state.add("Mg+2"   ,  Mg * cf, "mg")
    state.add("Na+"    ,  Na * cf, "mg")
    state.add("K+"     ,   K * cf, "mg")
    state.add("Cl-"    , Cl * cf, "mg")
    state.add("HCO3-"  ,   HCO3 * cf, "mg")
    state.add("SO4-2"  ,  SO4 * cf, "mg")

    # Additional chemicals
    added_Na = Na2CO3 / 106 * 46
    added_CO3 = Na2CO3 / 106 * 60

    state.add("CO2"  ,  CO2, "mg")
    state.add("Na+"  ,  added_Na, "mg")
    state.add("CO3-2"  ,  added_CO3, "mg")



    # equilibrate(state)

    conditions = EquilibriumConditions(specs)
    conditions.temperature(temp, "celsius")
    conditions.pressure(pressure, "bar")

    # conditions.phaseAmount("GaseousPhase", 0.0001, "umol")
    # conditions.setLowerBoundPressure(0.01, "bar")
    # conditions.setUpperBoundPressure(1000.0, "bar")
    conditions.setLowerBoundTemperature(25, "celsius")
    conditions.setUpperBoundTemperature(200, "celsius")
    # conditions.pH(2.0)
    conditions.charge(0.0)

    result = solver.solve(state, conditions)


    # print(state)
    
    props = ChemicalProps(state)
    aprops = AqueousProps(system)
    aprops.update(state)

    # print(props)
    # print(state)
    # print('density = ', props.phaseProps("AqueousPhase").density())
    # print('pH = ', float(aprops.pH()))
    # print(props.phaseProps("AqueousPhase").speciesAmounts())
    props = ChemicalProps(state)
    # print(props)
    output = {'pH': float(aprops.pH()),
              'density': props.phaseProps("AqueousPhase").density(),
              'temperature': float(aprops.temperature()) - 273.15,
              'pressure': float(aprops.pressure()),
              'Amount': state.speciesAmounts().asarray()
              }

    return output, props

rr = [0, 0.1,  0.3,  0.5,  0.7, 0.8, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95]
NaCls = []
Anhydrites = []
Magnesites = []
Calcites = []
Polyhalites = []
Glauberites = []
Gypsums = []
temperatures = []
pressures = []
pHs = []

rrs = []

#%% parametric on temp, pres
for i in rr:
    output,props = seawater(desal_rr=i,
                      pressure = 1.01325,
                      temp = 35)
    rrs.append(i)

    NaCls.append(output['Amount'][-4])
    Anhydrites.append(output['Amount'][-3])
    Magnesites.append(output['Amount'][-8])
    Calcites.append(output['Amount'][-7])
    Polyhalites.append(output['Amount'][-6])
    Glauberites.append(output['Amount'][-5])

    temperatures.append(output['temperature'])
    pressures.append(output['pressure'])
    pHs.append(output['pH'])

# Solid precipitation

fig = Figure()
fig.xaxisTitle('Desal recovery ratio')
fig.yaxisTitle('Precipitation [mol/kgw]')
fig.drawLineWithMarkers(rrs, NaCls, name="NaCl")
fig.drawLineWithMarkers(rrs , Anhydrites, name="Anhydrite")
fig.show()

# Temp and pressure

#%% parametric on NaOH, CO2

# rrs = []

CO2s = [i * 500 for i in range(31)]
NaOHs = [i * 500 for i in range(31)]

NaCls = pd.DataFrame(index = CO2s, columns=NaOHs)
Anhydrites = pd.DataFrame(index = CO2s, columns=NaOHs)
Magnesites = pd.DataFrame(index = CO2s, columns=NaOHs)
Calcites = pd.DataFrame(index = CO2s, columns=NaOHs)
Polyhalites = pd.DataFrame(index = CO2s, columns=NaOHs)
Glauberites = pd.DataFrame(index = CO2s, columns=NaOHs)
Gypsums = pd.DataFrame(index = CO2s, columns=NaOHs)
temperatures = pd.DataFrame(index = CO2s, columns=NaOHs)
pressures = pd.DataFrame(index = CO2s, columns=NaOHs)
pHs = pd.DataFrame(index = CO2s, columns=NaOHs)

Ca_rr = pd.DataFrame(index = CO2s, columns=NaOHs)
Mg_rr = pd.DataFrame(index = CO2s, columns=NaOHs)

for i in CO2s:
    for j in NaOHs:
        output,props = seawater(
                        desal_rr=0.85,
                        pressure = 1.01325,
                        temp = 35,
                        # CO2 = i,
                        # NaOH = j,
                        )
        # rrs.append(i)

        NaCls.loc[i, j] = output['Amount'][-4]
        Anhydrites.loc[i, j] = output['Amount'][-3]
        Magnesites.loc[i, j] = output['Amount'][-8]
        Calcites.loc[i, j] = output['Amount'][-7]
        Polyhalites.loc[i, j] = output['Amount'][-6]
        Glauberites.loc[i, j] = output['Amount'][-5]

        temperatures.loc[i, j] = output['temperature']
        pressures.loc[i, j] = output['pressure']
        pHs.loc[i, j] = output['pH']

        Ca_rr.loc[i,j] = output['Amount'][-7] * 40.087 / (0.416 / 0.25)
        Mg_rr.loc[i,j] = output['Amount'][-8] * 24.305 / (1.295 / 0.25)

# import pandas as pd

# results = np.array(
#     [
#     CO2s,
#     NaOHs,
#     pHs,
#     NaCls,
#     Anhydrites,
#     Magnesites,
#     Calcites,
#     Polyhalites,
#     Glauberites
#     ]
# )

# data_table = pd.DataFrame(
#     data=results,
#     index=[
#     "CO2",
#     "NaOH",
#     "pH",
#     "NaCl",
#     "Anhydrite",
#     "Magnesite",
#     "Calcite",
#     "Polyhalite",
#     "Glauberite",
#     ],
# )

# Solid precipitation

# fig = Figure()
# fig.xaxisTitle('Desal recovery ratio')
# fig.yaxisTitle('Precipitation [mol/kgw]')
# fig.drawLineWithMarkers(rrs, NaCls, name="NaCl")
# fig.drawLineWithMarkers(rrs , Anhydrites, name="Anhydrite")
# fig.show()

# %%
# Mg_rr.to_numpy()


Mg = Mg_rr.to_numpy().astype(float)
Ca = Ca_rr.to_numpy().astype(float)
pH = pHs.to_numpy().astype(float)

extent = np.min(NaOHs), np.max(NaOHs), np.min(CO2s), np.max(CO2s)

fig1, ax1 = plt.subplots(1,1)

im1 = ax1.imshow(np.flip(Mg,0), cmap='Blues', interpolation='nearest',extent = extent,
                  vmin = 0, vmax = 1)


fig1.colorbar(im1).set_label('Recovery rate of Mg')
# plt.show()
ax1.set_xlabel('NaOH addition (mg)')

ax1.set_ylabel('CO2 addition (mg)')

# fig2
fig2, ax2 = plt.subplots(1,1)

im2 = ax2.imshow(np.flip(Ca,0), cmap='Blues', interpolation='nearest',extent = extent,
                  vmin = 0, vmax = 1)


fig2.colorbar(im2).set_label('Recovery rate of Ca')

ax2.set_xlabel('NaOH addition (mg)')
ax2.set_ylabel('CO2 addition (mg)')

# fig 3
fig3, ax3 = plt.subplots(1,1)
im3 = ax3.imshow(np.flip(pH,0), interpolation='nearest',extent = extent,
                  vmin = 0, vmax = 14)

fig3.colorbar(im3).set_label('pH')
ax3.set_xlabel('NaOH addition (mg)')
ax3.set_ylabel('CO2 addition (mg)')

fig1.savefig('Mg_85.png', dpi=300)
fig2.savefig('Ca_85.png', dpi=300)
# fig3.savefig('pH_85.png', dpi=300)

# %% Replicate of paper, add NaOH
def brine_NaOH(desal_rr = 0.8,
          pressure = 1.01325, # bar
        NaOH = 0, # mg
        CO2 = 0, # mg
        temp = 25, 
    ): 
    db = PhreeqcDatabase("pitzer.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S"))
    solution.set(ActivityModelPitzer())

    # Create a gaseous phase
    # gaseousphase = GaseousPhase("H2O(g)")
    # gaseousphase.set(ActivityModelPengRobinson())

    # Hydromagnesite = MineralPhase("Hydromagnesite")
    Magnesite = MineralPhase("Magnesite")
    Calcite = MineralPhase("Calcite")
    Polyhalite = MineralPhase("Polyhalite")
    Glauberite = MineralPhase("Glauberite")
    Halite = MineralPhase("Halite")
    Gypsum = MineralPhase("Gypsum")
    Anhydrite = MineralPhase("Anhydrite")
    Epsomite = MineralPhase("Epsomite")
    Brucite = MineralPhase("Brucite")
    Portlandite = MineralPhase("Portlandite")
    

    system = ChemicalSystem(db, 
                            # gaseousphase,
                            solution,
                            Magnesite,
                            Calcite,
                            Polyhalite,
                            Glauberite,
                            Halite,
                            Anhydrite,
                            Gypsum,
                            Epsomite,
                            Brucite,
                            Portlandite,
                            )

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    # specs.phaseAmount("GaseousPhase")
    # specs.pH()
    specs.charge()
    specs.openTo("Cl-")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)
    

    state = ChemicalState(system)
    state.temperature(90, "celsius")

    # Seawater composition (mg/kg)
    Na  = 17987
    Cl  = 32220
    SO4 = 4766
    Mg = 2226
    Ca = 714
    K  = 858
    # HCO3= 300


    # state.pressure(1.0, "atm")
    state.set("H2O", 1.0, "kg")
    state.add("Ca+2"   ,   Ca , "mg")
    state.add("Mg+2"   ,  Mg , "mg")
    state.add("Na+"    ,  Na , "mg")
    state.add("K+"     ,   K , "mg")
    state.add("Cl-"    , Cl , "mg")
    state.add("SO4-2"  ,  SO4 , "mg")
    # state.add("HCO3-"  ,  HCO3 , "mg")


    # Additional chemicals (1 M NaOH)
    added_NaOH = NaOH * 0.039
    added_H2O = NaOH * 0.961

    state.add("CO2"  ,  CO2, "mg")
    state.add("Na+"  ,  added_NaOH / 40 * 23, "mg")
    state.add("OH-"  ,  added_NaOH / 40 * 17, "mg")
    state.add("H2O"  ,  added_H2O, "mg")



    # equilibrate(state)

    conditions = EquilibriumConditions(specs)
    conditions.temperature(temp, "celsius")
    conditions.pressure(pressure, "bar")

    # conditions.phaseAmount("GaseousPhase", 0.0001, "umol")
    # conditions.setLowerBoundPressure(0.01, "bar")
    # conditions.setUpperBoundPressure(1000.0, "bar")
    conditions.setLowerBoundTemperature(25, "celsius")
    conditions.setUpperBoundTemperature(200, "celsius")
    # conditions.pH(2.0)
    conditions.charge(0.0)

    result = solver.solve(state, conditions)


    # print(state)
    
    props = ChemicalProps(state)
    aprops = AqueousProps(system)
    aprops.update(state)

    # print(props)
    # print(state)
    # print('density = ', props.phaseProps("AqueousPhase").density())
    # print('pH = ', float(aprops.pH()))
    # print(props.phaseProps("AqueousPhase").speciesAmounts())
    props = ChemicalProps(state)
    # print(props)
    output = {'pH': float(aprops.pH()),
              'density': props.phaseProps("AqueousPhase").density(),
              'temperature': float(aprops.temperature()) - 273.15,
              'pressure': float(aprops.pressure()),
              'Amount': state.speciesAmounts().asarray()
              }

    return output, props

test_output,  test_props = brine_NaOH()
print('pH: ', test_output['pH'])

#%% Plotting
NaCls = []
Anhydrites = []
Magnesites = []
Calcites = []
Polyhalites = []
Glauberites = []
Gypsums = []
temperatures = []
pressures = []
pHs = []
Brucites = []
Portlandites = []

Nas = []
Mgs = []
Cas = []


# CO2s = [i * 500 for i in range(31)]
NaOHs = [i * 5000 for i in range(91)]

# for i in CO2s:
for j in NaOHs:
    output,props = brine_NaOH(
                    pressure = 1.01325,
                    temp = 25,
                    # CO2 = i,
                    NaOH = j,
                    )
    # rrs.append(i)

    NaCls.append(output['Amount'][-6])
    Anhydrites.append(output['Amount'][-5])
    Magnesites.append(output['Amount'][-10])
    Calcites.append(output['Amount'][-9])
    Polyhalites.append(output['Amount'][-8])
    Glauberites.append(output['Amount'][-7])
    Brucites.append(output['Amount'][-2])
    Portlandites.append(output['Amount'][-1])

    temperatures.append(output['temperature'])
    pressures.append(output['pressure'])
    pHs.append(output['pH'])

    Nas.append(output['Amount'][-12] * 23 / (output['Amount'][1] * 18 / 1000) * 1000)
    Cas.append(output['Amount'][4] * 40 / (output['Amount'][1] * 18 / 1000) * 1000)
    Mgs.append(output['Amount'][9] * 24 / (output['Amount'][1] * 18 / 1000) * 1000)


fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(pHs, Cas, 'b-', label='Ca')
ax1.plot(pHs, Mgs, 'r-', label='Mg')

ax2.plot(pHs, Nas, 'k-', label='Na')

ax1.set_xlabel('brine pH')
ax1.set_ylabel('Mg, Ca (mg/L)', color='k')
ax2.set_ylabel('Na (mg/L)', color='k')

ax1.set_ylim(-100, 2200)
ax2.set_ylim(17000, 20000)
ax1.set_xlim(7.8,13.2)
ax1.set_yticks([i * 200 for i in range(13)])

ax1.legend(loc='upper center')
ax2.legend(loc = 'best')
plt.show()


#%% replicate of paper, add CO2
def brine_CO2(desal_rr = 0.8,
        pressure = 1.01325, # bar
        NaOH = 5000 * 34, # mg (of 1M solution)
        CO2 = 0, # mg
        temp = 25, 
    ): 
    db = PhreeqcDatabase("pitzer.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S"))
    solution.set(ActivityModelPitzer())

    # Create a gaseous phase
    # gaseousphase = GaseousPhase("H2O(g)")
    # gaseousphase.set(ActivityModelPengRobinson())

    # Hydromagnesite = MineralPhase("Hydromagnesite")
    Magnesite = MineralPhase("Magnesite")
    Calcite = MineralPhase("Calcite")
    Polyhalite = MineralPhase("Polyhalite")
    Glauberite = MineralPhase("Glauberite")
    Halite = MineralPhase("Halite")
    Gypsum = MineralPhase("Gypsum")
    Anhydrite = MineralPhase("Anhydrite")
    Epsomite = MineralPhase("Epsomite")
    Brucite = MineralPhase("Brucite")
    Portlandite = MineralPhase("Portlandite")
    

    system = ChemicalSystem(db, 
                            # gaseousphase,
                            solution,
                            Magnesite,
                            Calcite,
                            Polyhalite,
                            Glauberite,
                            Halite,
                            Anhydrite,
                            Gypsum,
                            Epsomite,
                            Brucite,
                            Portlandite,
                            )

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    # specs.phaseAmount("GaseousPhase")
    # specs.pH()
    specs.charge()
    specs.openTo("Cl-")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)
    

    state = ChemicalState(system)
    state.temperature(90, "celsius")

    # Seawater composition (mg/kg)
    Na  = 17987
    Cl  = 32220
    SO4 = 4766
    Mg = 2226
    Ca = 714
    K  = 858
    # HCO3= 300


    # state.pressure(1.0, "atm")
    state.set("H2O", 1.0, "kg")
    state.add("Ca+2"   ,   Ca , "mg")
    state.add("Mg+2"   ,  Mg , "mg")
    state.add("Na+"    ,  Na , "mg")
    state.add("K+"     ,   K , "mg")
    state.add("Cl-"    , Cl , "mg")
    state.add("SO4-2"  ,  SO4 , "mg")
    # state.add("HCO3-"  ,  HCO3 , "mg")


    # Additional chemicals (1 M NaOH)

    added_NaOH = NaOH * 0.039
    added_H2O = NaOH * 0.961

    state.add("CO2"  ,  CO2, "mg")
    state.add("Na+"  ,  added_NaOH/40*23, "mg")
    state.add("OH-"  ,  added_NaOH/40*17, "mg")
    state.add("H2O"  ,  added_H2O, "mg")



    # equilibrate(state)

    conditions = EquilibriumConditions(specs)
    conditions.temperature(temp, "celsius")
    conditions.pressure(pressure, "bar")

    # conditions.phaseAmount("GaseousPhase", 0.0001, "umol")
    # conditions.setLowerBoundPressure(0.01, "bar")
    # conditions.setUpperBoundPressure(1000.0, "bar")
    conditions.setLowerBoundTemperature(25, "celsius")
    conditions.setUpperBoundTemperature(200, "celsius")
    # conditions.pH(2.0)
    conditions.charge(0.0)

    result = solver.solve(state, conditions)


    # print(state)
    
    props = ChemicalProps(state)
    aprops = AqueousProps(system)
    aprops.update(state)

    # print(props)
    # print(state)
    # print('density = ', props.phaseProps("AqueousPhase").density())
    # print('pH = ', float(aprops.pH()))
    # print(props.phaseProps("AqueousPhase").speciesAmounts())
    props = ChemicalProps(state)
    # print(props)
    output = {'pH': float(aprops.pH()),
              'density': props.phaseProps("AqueousPhase").density(),
              'temperature': float(aprops.temperature()) - 273.15,
              'pressure': float(aprops.pressure()),
              'Amount': state.speciesAmounts().asarray()
              }

    return output, props

test_output,  test_props = brine_CO2()
print('pH: ', test_output['pH'])


#%% Plotting
NaCls = []
Anhydrites = []
Magnesites = []
Calcites = []
Polyhalites = []
Glauberites = []
Gypsums = []
temperatures = []
pressures = []
pHs = []
Brucites = []
Portlandites = []

Nas = []
Mgs = []
Cas = []


CO2s = [i * 50 for i in range(75)]
# NaOHs = [i * 5000 for i in range(91)]

for i in CO2s:
    output,props = brine_CO2(
                    pressure = 1.01325,
                    temp = 25,
                    CO2 = i,
                    # NaOH = j,
                    )
    # rrs.append(i)
    if output['pH']>10:
        continue

    NaCls.append(output['Amount'][-6])
    Anhydrites.append(output['Amount'][-5])
    Magnesites.append(output['Amount'][-10])
    Calcites.append(output['Amount'][-9])
    Polyhalites.append(output['Amount'][-8])
    Glauberites.append(output['Amount'][-7])
    Brucites.append(output['Amount'][-2])
    Portlandites.append(output['Amount'][-1])

    temperatures.append(output['temperature'])
    pressures.append(output['pressure'])
    pHs.append(output['pH'])

    Nas.append(output['Amount'][-12] * 23 / (output['Amount'][1] * 18 / 1000) * 1000)
    Cas.append(output['Amount'][4] * 40 / (output['Amount'][1] * 18 / 1000) * 1000)
    Mgs.append(output['Amount'][9] * 24 / (output['Amount'][1] * 18 / 1000) * 1000)


fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(pHs, Cas, 'b-')
ax1.plot(pHs, Mgs, 'r-')

ax2.plot(pHs, Nas, 'k-')

ax1.set_xlabel('brine pH')
ax1.set_ylabel('Mg, Ca (mg/L)', color='k')
ax2.set_ylabel('Na (mg/L)', color='k')


# ax1.set_ylim(-100, 2200)
ax2.set_ylim(15000, 20000)
# ax1.set_xlim(7.8,13.2)
ax1.set_yticks([i * 200 for i in range(13)])
plt.gca().invert_xaxis()
plt.show()



# %%
from reaktoro import *

# Initialize the chemical system with Pitzer model
db = PhreeqcDatabase("pitzer.dat")

# Define the aqueous phase and solid phase
aqueous = AqueousPhase(speciate("H2O Na+ K+ Mg+2 Ca+2 Cl- SO4-2 OH- H+"))
aqueous.setActivityModel("Pitzer")
solid = SolidPhase("Mg(OH)2")  # Define Mg(OH)2 as a solid phase

# Create the chemical system
system = ChemicalSystem([aqueous, solid])

# Define the initial chemical state
state = ChemicalState(system)
state.setTemperature(25, "celsius")
state.setPressure(1, "atm")
state.setSpeciesMass("Na+", 17987, "mg/L")
state.setSpeciesMass("K+", 858, "mg/L")
state.setSpeciesMass("Mg+2", 2226, "mg/L")
state.setSpeciesMass("Ca+2", 714, "mg/L")
state.setSpeciesMass("Cl-", 32220, "mg/L")
state.setSpeciesMass("SO4-2", 4766, "mg/L")

# Define the kinetic reaction for Mg(OH)2 precipitation
kinetic_reaction = KineticReaction("Mg(OH)2(s) = Mg+2 + 2OH-")

# Define rate law parameters
kf = 1.0e-3  # Forward rate constant (mol/m^3/s)
kr = 1.0e-6  # Reverse rate constant (mol/m^3/s)
kinetic_reaction.set_rate_expression("kf * activity(Mg+2) * activity(OH-)^2 - kr")
kinetic_reaction.set_parameters({"kf": kf, "kr": kr})

# Add kinetic reaction to the system
kinetics = Kinetics(system)
kinetics.add(kinetic_reaction)

# Initialize time and solver settings
time = 0.0  # Start time (s)
time_step = 1.0  # Time step (s)
total_time = 600.0  # Total simulation time (s)
solver = KineticsSolver(system)
solver.set_time_step(time_step)
solver.set_tolerance(1e-8)

# Store results
results = []

# Simulate the kinetic process
while time < total_time:
    solver.solve(state, kinetics)
    pH = state.pH()
    mg_concentration = state.speciesAmount("Mg+2")
    results.append((time, pH, mg_concentration))
    time += time_step

# Output results
for t, pH, mg in results:
    print(f"Time: {t:.1f} s, pH: {pH:.2f}, Mg2+ concentration: {mg:.6f} mol/L")
# %%
