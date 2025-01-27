#%%
from reaktoro import *
from reaktplot import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

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
    # Magnesites.append(output['Amount'][-10])
    Calcites.append(output['Amount'][-9])
    Magnesites.append(float(props.speciesAmount("Magnesite")))
    Polyhalites.append(output['Amount'][-8])
    Glauberites.append(output['Amount'][-7])
    # Brucites.append(output['Amount'][-2])
    Brucites.append(float(props.speciesAmount("Brucite")))
    Portlandites.append(output['Amount'][-1])

    temperatures.append(output['temperature'])
    pressures.append(output['pressure'])
    pHs.append(output['pH'])

    Nas.append(output['Amount'][-12] * 23 / (output['Amount'][1] * 18 / 1000) * 1000)
    Cas.append(output['Amount'][4] * 40 / (output['Amount'][1] * 18 / 1000) * 1000)
    Mgs.append(output['Amount'][9] * 24 / (output['Amount'][1] * 18 / 1000) * 1000)


fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(pHs, Magnesites, 'b-', label='MgCO3')
ax1.plot(pHs, Brucites, 'r-', label='Mg(OH)2')

ax2.plot(pHs, Nas, 'k-', label='Na')

ax1.set_xlabel('brine pH')
ax1.set_ylabel('Mg, Ca (mg/L)', color='k')
ax2.set_ylabel('Na (mg/L)', color='k')

# ax1.set_ylim(-100, 2200)
ax2.set_ylim(17000, 20000)
ax1.set_xlim(7.8,13.2)
# ax1.set_yticks([i * 200 for i in range(13)])

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
    Mg = 2226 - 2226
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


CO2s = [i * 50 for i in range(150)]
# NaOHs = [i * 5000 for i in range(91)]

for i in CO2s:
    output,props = brine_CO2(
                    pressure = 1.01325,
                    temp = 25,
                    CO2 = i,
                    # NaOH = j,
                    )
    # rrs.append(i)
    # if output['pH']>10:
    #     continue

    # NaCls.append(output['Amount'][-6])
    # Anhydrites.append(output['Amount'][-5])
    # Magnesites.append(output['Amount'][-10])
    # Calcites.append(output['Amount'][-9])
    # Polyhalites.append(output['Amount'][-8])
    # Glauberites.append(output['Amount'][-7])
    # Brucites.append(output['Amount'][-2])
    # Portlandites.append(output['Amount'][-1])

    Calcites.append(float(props.speciesAmount("Calcite")))
    Cas.append(float(props.speciesAmount("Ca+2")))

    temperatures.append(output['temperature'])
    pressures.append(output['pressure'])
    pHs.append(output['pH'])

    # Nas.append(output['Amount'][-12] * 23 / (output['Amount'][1] * 18 / 1000) * 1000)
    # Cas.append(output['Amount'][4] * 40 / (output['Amount'][1] * 18 / 1000) * 1000)
    # Mgs.append(output['Amount'][9] * 24 / (output['Amount'][1] * 18 / 1000) * 1000)


fig, ax1 = plt.subplots()

ax2 = ax1.twinx()
ax1.plot(CO2s, Calcites, 'b-', label='CaCO3')
ax1.plot(CO2s, Cas, 'r-',label='Ca+2')

ax2.plot(CO2s, pHs, 'k-', label = 'pH')

ax1.set_xlabel('CO2 add (mg)')
ax1.set_ylabel('Mg, Ca (mg/L)', color='k')
ax2.set_ylabel('pH', color='k')


ax1.legend(loc='best')
ax2.legend(loc = 'upper center')
# ax1.set_ylim(-100, 2200)
# ax2.set_ylim(15000, 20000)
# ax1.set_xlim(7.8,13.2)
# ax1.set_yticks([i * 200 for i in range(13)])
# plt.gca().invert_xaxis()
plt.show()


# %%
