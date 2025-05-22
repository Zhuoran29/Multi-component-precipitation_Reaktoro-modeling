#%% imports
from reaktoro import *
from reaktplot import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import csv

__all__ = [
    "NaOH",
    "permian_brine",
    "permian_step_1",
    "plot_step_1_permian",
    "permian_step_2_alt",
    "plot_step_2_alt",
    "permian_step_3_alt",
    "plot_step_3_alt",
]

global molar_mass
molar_mass = {
    "H2O": 18.01528,
    "H2O(aq)": 18.01528,
    "Na+": 22.989769,
    "Mg+2": 24.305, 
    "Cl-": 35.453,
    "Ca+2": 40.078,
    "K+": 39.0983,
    "SO4-2": 96.06,
    "OH-": 17.008,
    "HCO3-": 61.0168,
    "CO3-2": 60.008,
    "Brucite":58.3197,
    "Calcite": 100.0869,
    "Magnesite": 84.3139,
    "CaSO4(aq)": 136.14,
    "CaOH+": 57,
    "MgOH+": 41.305,
    "MgSO4(aq)": 120.366,
    "NaSO4-": 22.99+96.06,
}

permian_feed = {
    "Na+"  : 40896,
    "Cl-"  : 78648,
    "SO4-2" : 496.3,
    "Mg+2" : 745,
    "Ca+2" : 3821,
}

#%% NaOH solution 
def NaOH(conc = 5, # mol/L
         ):
    # db = PhreeqcDatabase("pitzer.dat")
    db = SupcrtDatabase("supcrtbl-organics")

    solution = AqueousPhase(speciate("H O Na"))
    solution.set(ActivityModelPitzer())

    system = ChemicalSystem(db, 
                            solution)
    
    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.volume()
    specs.charge()
    specs.openTo("Na+")
    specs.openTo("H2O")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)
    

    state = ChemicalState(system)
    state.temperature(25, "celsius")

    state.set("H2O(aq)", 1.0, "kg")
    state.add("Na+"  ,  conc , "mol")
    state.add("OH-"  ,  conc , "mol")

    conditions = EquilibriumConditions(specs)
    conditions.temperature(25, "celsius")
    conditions.pressure(1.01325, "bar")
    conditions.charge(0.0)
    conditions.volume(0.001, "m3")

    result = solver.solve(state, conditions)


    # print(state)
    props = ChemicalProps(state)
    aprops = AqueousProps(system)

    density = float(props.phaseProps("AqueousPhase").density())
    H2O_mass = float(props.speciesAmount("H2O(aq)")) * molar_mass["H2O"] # g/L
    Na_mass = float(props.speciesAmount("Na+")) * molar_mass["Na+"] # g/L
    OH_mass = float(props.speciesAmount("OH-")) * molar_mass["OH-"]    # g/L 

    return props, [H2O_mass, Na_mass, OH_mass]

#%% Original solution (Permian basin)
def permian_brine(
        feed_props = {    
            "Na+"  : 40896,
            "Cl-"  : 78648,
            "SO4-2" : 496.3,
            "Mg+2" : 745,
            "Ca+2" : 3821,
            # "H4SiO4":   108/28.1*93.1,
            },
        pressure = 1.01325, # bar
        temp = 25, 
        RR = 0.5,
    ): 
    # db = PhreeqcDatabase("thermoddem-v1.10.dat")
    # db = PhreeqcDatabase("pitzer.dat")
    # db = PhreeqcDatabase("minteq.v4.dat")
    # db = PhreeqcDatabase("phreeqc.dat")
    db = SupcrtDatabase("supcrtbl-organics")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S Si"))
    solution.set(ActivityModelPitzer())
    # solution.set(ActivityModelDebyeHuckel())

    # Create a gaseous phase
    gaseousphase = GaseousPhase("CO2(g)")
    gaseousphase.set(ActivityModelPengRobinson())

    Magnesite = MineralPhase("Magnesite")
    
    Calcite = MineralPhase("Calcite")
    Halite = MineralPhase("Halite")
    # Gypsum = MineralPhase("Gypsum")
    Anhydrite = MineralPhase("Anhydrite")
    # Epsomite = MineralPhase("Epsomite")
    Brucite = MineralPhase("Brucite")
    Portlandite = MineralPhase("Portlandite")
    Chalcedony = MineralPhase("Chalcedony")
    

    system = ChemicalSystem(db, 
                            gaseousphase,
                            solution,
                            Magnesite,
                            Calcite,
                            Halite,
                            Anhydrite,
                            # Gypsum,
                            # Epsomite,
                            Brucite,
                            # Portlandite,
                            Chalcedony,
                            )

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    # specs.phaseAmount("GaseousPhase")
    specs.pH()
    specs.charge()
    specs.openTo("Cl-")
    specs.volume()
    specs.openTo("H2O(aq)")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)
    

    state = ChemicalState(system)
    state.temperature(temp, "celsius")

    state.set("H2O(aq)", 1.0, "kg")
    for ion, conc in feed_props.items():
        state.add(ion, conc/(1-RR), "mg")

    conditions = EquilibriumConditions(specs)
    conditions.temperature(temp, "celsius")
    conditions.pressure(pressure, "bar")

    # conditions.phaseAmount("GaseousPhase", 0.0001, "umol")
    # conditions.setLowerBoundPressure(0.01, "bar")
    # conditions.setUpperBoundPressure(1000.0, "bar")
    conditions.setLowerBoundTemperature(25, "celsius")
    conditions.setUpperBoundTemperature(200, "celsius")
    conditions.pH(6.6)
    conditions.charge(0.0)
    conditions.volume(0.001, "m3")

    result = solver.solve(state, conditions)


    # print(state)
    
    props = ChemicalProps(state)
    aprops = AqueousProps(system)
    aprops.update(state)
    props = ChemicalProps(state)

    output = {'pH': float(aprops.pH()),
              'density': float(props.phaseProps("AqueousPhase").density()),
              'temperature': float(aprops.temperature()) - 273.15,
              'pressure': float(aprops.pressure()),
              }
    
    outflow_species = ["H2O(aq)", "Ca+2", "Cl-", "SO4-2","MgSO4(aq)","CaSO4(aq)"
                       "Mg+2", "Na+", "OH-", "K+", "HCO3-", "CO3-2"]
    outflow_mass = dict(zip(outflow_species,
                            [float(props.speciesAmount(i)) * molar_mass[i] for i in outflow_species])) # mg

    # print('Original feed pH:', float(aprops.pH()))

    return output, state, system, props

ori_output, ori_state, ori_system, ori_props = permian_brine(RR=0.5)


# %% Step 1: add NaOH
def permian_step_1(
        state, system,
        pressure = 1.01325, # bar
        add_NaOH_vol = 0.0365, # L
        add_NaOH_conc = 5, # mol/L
        add_NaOH_mass = 5, # gram
        temp = 25, 
    ): 
    # solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S"))
    # solution.set(ActivityModelPitzer())

    # Create a gaseous phase
    # gaseousphase = GaseousPhase("CO2(g)")
    # gaseousphase.set(ActivityModelPengRobinson())

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    specs.charge()
    specs.openTo("Cl-")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)


    # Additional chemicals (5 M NaOH)
    if add_NaOH_conc == 'solid':
        Na_mass = add_NaOH_mass * 23/40 # gram
        OH_mass = add_NaOH_mass * 17/40 # gram
        H2O_mass = 0

        state.add("Na+"  ,  Na_mass  / 1000, "kg")
        state.add("OH-"  ,  OH_mass  / 1000, "kg")
        state.add("H2O(aq)"  ,  H2O_mass  / 1000, "kg")
    else:
        NaOH_props, masses = NaOH(conc = add_NaOH_conc)
        H2O_mass, Na_mass, OH_mass = masses # g/L

        state.add("Na+"  ,  Na_mass * add_NaOH_vol / 1000, "kg")
        state.add("OH-"  ,  OH_mass * add_NaOH_vol / 1000, "kg")
        state.add("H2O(aq)"  ,  H2O_mass * add_NaOH_vol / 1000, "kg")

    conditions = EquilibriumConditions(specs)
    conditions.temperature(temp, "celsius")
    conditions.pressure(pressure, "bar")

    conditions.setLowerBoundTemperature(25, "celsius")
    conditions.setUpperBoundTemperature(200, "celsius")
    conditions.charge(0.0)

    result = solver.solve(state, conditions)

    props = ChemicalProps(state)
    aprops = AqueousProps(system)
    aprops.update(state)

    props = ChemicalProps(state)
    # print(props)
    output = {'pH': float(aprops.pH()),
              'density': props.phaseProps("AqueousPhase").density(),
              'temperature': float(aprops.temperature()) - 273.15,
              'pressure': float(aprops.pressure()),
              'Amount': state.speciesAmounts().asarray()
              }
    out_brucite = float(props.speciesAmount("Brucite")) * 58.3197 * 1000 # mg
    outflow_species = ["H2O(aq)", "Ca+2", "Cl-", "SO4-2", "CaOH+", "MgOH+", "NaSO4-",
                       "Mg+2", "Na+", "OH-", "K+", "HCO3-", "CO3-2"]
    outflow_mass = dict(zip(outflow_species,
                            [float(props.speciesAmount(i)) * molar_mass[i] for i in outflow_species])) # mg

    return output, props, outflow_mass, out_brucite, state, system


permian_brine_output, permian_brine_state, permian_brine_system, permian_brine_props = permian_brine(RR=0.5)
output_1_PB, props_1_PB, outflow_mass_1_PB, out_brucite_1_PB, state_1_PB, system_1_PB = permian_step_1(state = permian_brine_state,
                                                                                     system=permian_brine_system,
                                                                                     add_NaOH_conc=5,
                                                                                     add_NaOH_vol=0.025,
                                                                                     add_NaOH_mass=3,
                                                                             )
print('Step One simulation')
print('pH=', output_1_PB['pH'])
print('Ca in CaCO3 (mg)', float(props_1_PB.speciesAmount('Calcite')) * molar_mass['Ca+2']* 1000)
print('Ca+2 (mg)', float(props_1_PB.speciesAmount('Ca+2')) * molar_mass['Ca+2']* 1000)
print('Mg in MgCO3 (mg)', (float(props_1_PB.speciesAmount('Magnesite')))* molar_mass['Mg+2']* 1000 )
print('Mg in MgOH2 (mg)', (float(props_1_PB.speciesAmount('Brucite')))* molar_mass['Mg+2'] * 1000 )
print('Mg+2 (mg)', (float(props_1_PB.speciesAmount('Mg+2')))* molar_mass['Mg+2'] * 1000 )
print('Ca rec rate', float(props_1_PB.speciesAmount('Anhydrite')) / float(props_1_PB.elementAmount('Ca')))
print('Mg rec rate', float(props_1_PB.speciesAmount('Brucite')) / float(props_1_PB.elementAmount('Mg')))

#%% Step 1: Plot (add solid)
#%% Step 1: Plot (add solution)
def plot_step_1_permian(       
        maxNaOH_vol = 0.5, # L
        NaOH_conc = 1, # mol/L
        RR = 0.5, 
        runs = 50,
        ):
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
    Mg_rr = []
    vols = []

    NaOHs = [i * maxNaOH_vol / runs for i in range(0,runs+1,1)]
    for v in NaOHs:
        ori_output, ori_state, ori_system, ori_props = permian_brine(RR=RR)

        output, props, outflow_mass, out_brucite, state, system = permian_step_1(
            state=ori_state, 
            system=ori_system,
            pressure = 1.01325, # bar
            add_NaOH_vol = v, # L
            add_NaOH_conc = NaOH_conc, # mol/L
            temp = 25, 
            )

        Calcites.append(float(props.speciesAmount("Calcite"))*1000)
        Brucites.append(float(props.speciesAmount("Brucite"))*1000)
        Magnesites.append(float(props.speciesAmount("Magnesite"))*1000)
        Anhydrites.append(float(props.speciesAmount("Anhydrite"))*1000)
        # Gypsums.append(float(props.speciesAmount("Gypsum"))*1000)

        temperatures.append(output['temperature'])
        pressures.append(output['pressure'])
        pHs.append(output['pH'])

        Cas.append(float(props.speciesAmount("Ca+2"))
                #    +float(props.speciesAmount("CaSO4"))+
                #    float(props.speciesAmount("CaCl+"))+
                #    float(props.speciesAmount("CaCl2"))
                   )
        Nas.append(float(props.speciesAmount("Na+"))
                #    +float(props.speciesAmount("NaSO4-"))
                   )
        Mgs.append(float(props.speciesAmount("Mg+2"))
                #    +float(props.speciesAmount("MgSO4"))
                #    +float(props.speciesAmount("MgCl+"))
                #    +float(props.speciesAmount("MgOH+"))
                   )

        vols.append(float(props.phaseProps("AqueousPhase").volume())*1000)

        Mg_rr.append(props.speciesAmount("Brucite")/float(props.elementAmount("Mg"))*100)

        # print(masses)
    fig, ax1 = plt.subplots()

    Ca_concs = [ Cas[i]*40.087/vols[i] for i in range(len(NaOHs))]
    Na_concs = [ Nas[i]*22.99/vols[i] for i in range(len(NaOHs))]
    Mg_concs = [ Mgs[i]*24.305/vols[i] for i in range(len(NaOHs))]

    calcites_conc = [ Calcites[i] / vols[i] for i in range(len(vols))]
    brucites_conc = [ Brucites[i] / vols[i] for i in range(len(vols))]
    magnesites_conc = [ Magnesites[i] / vols[i] for i in range(len(vols))]
    Anhydrites_conc = [ Anhydrites[i] / vols[i] for i in range(len(vols)) ]

    # ax2 = ax1.twinx()
    ax1.plot(pHs, Anhydrites, 'b-', label='Ca(OH)2')
    ax1.plot(pHs, Brucites, 'r-', label='Mg(OH)2')
    # ax1.plot(pHs, Gypsums, 'g-', label='CaSO4')


    ax1.set_xlabel('pH')
    ax1.set_ylabel('Solid (mmol/L)', color='k')
    # ax2.set_ylabel('Vol (L)', color='k')

    ax1.legend(loc='upper center')
    # ax2.legend(loc = 'best')
    plt.show()

    print('Max vol: ', maxNaOH_vol, "Conc: ", NaOH_conc)
    return Brucites, Portlandites, Mg_rr, pHs, vols, NaOHs


# Run function
maxNaOH_vol = 0.1 # L
NaOH_conc = 5 # mol/L
RR = 0.5
Brucites, Anhydrites, Mg_rr, pHs, vols, NaOHs = plot_step_1_permian(
                                                    maxNaOH_vol=maxNaOH_vol,
                                                    NaOH_conc = NaOH_conc,
                                                    RR = RR,
                                                    runs = 100,
                                                    )
# %%
