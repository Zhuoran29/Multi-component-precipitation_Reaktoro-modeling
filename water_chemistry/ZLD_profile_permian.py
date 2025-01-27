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
    "CaSO4": 136.14,
    "CaOH+": 57,
    "MgOH+": 41.305,
    "MgSO4": 120.366,
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
    db = PhreeqcDatabase("pitzer.dat")

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

    state.set("H2O", 1.0, "kg")
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
    H2O_mass = float(props.speciesAmount("H2O")) * molar_mass["H2O"] # g/L
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
    db = PhreeqcDatabase("minteq.v4.dat")
    # db = PhreeqcDatabase("phreeqc.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S Si"))
    solution.set(ActivityModelPitzer())
    # solution.set(ActivityModelDebyeHuckel())

    # Create a gaseous phase
    gaseousphase = GaseousPhase("CO2(g)")
    gaseousphase.set(ActivityModelPengRobinson())

    Magnesite = MineralPhase("Magnesite")
    
    Calcite = MineralPhase("Calcite")
    Halite = MineralPhase("Halite")
    Gypsum = MineralPhase("Gypsum")
    Anhydrite = MineralPhase("Anhydrite")
    Epsomite = MineralPhase("Epsomite")
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
                            Gypsum,
                            Epsomite,
                            Brucite,
                            Portlandite,
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
    specs.openTo("H2O")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)
    

    state = ChemicalState(system)
    state.temperature(temp, "celsius")

    state.set("H2O", 1.0, "kg")
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
    
    outflow_species = ["H2O", "Ca+2", "Cl-", "SO4-2",
                       "Mg+2", "Na+", "OH-", "K+", "HCO3-", "CO3-2"]
    outflow_mass = dict(zip(outflow_species,
                            [float(props.speciesAmount(i)) * molar_mass[i] for i in outflow_species])) # mg

    # print('Original feed pH:', float(aprops.pH()))

    return output, state, system, props

ori_output, ori_state, ori_system, ori_props = permian_brine(RR=0.5)
#%% Step 1: Add NaOH to separate Mg
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
        state.add("H2O"  ,  H2O_mass  / 1000, "kg")
    else:
        NaOH_props, masses = NaOH(conc = add_NaOH_conc)
        H2O_mass, Na_mass, OH_mass = masses # g/L

        state.add("Na+"  ,  Na_mass * add_NaOH_vol / 1000, "kg")
        state.add("OH-"  ,  OH_mass * add_NaOH_vol / 1000, "kg")
        state.add("H2O"  ,  H2O_mass * add_NaOH_vol / 1000, "kg")

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
    outflow_species = ["H2O", "Ca+2", "Cl-", "SO4-2", "CaOH+", "CaSO4", "MgOH+", "MgSO4", "NaSO4-",
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
print('Ca rec rate', float(props_1_PB.speciesAmount('Portlandite')) / float(props_1_PB.elementAmount('Ca')))
print('Mg rec rate', float(props_1_PB.speciesAmount('Brucite')) / float(props_1_PB.elementAmount('Mg')))

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
        Portlandites.append(float(props.speciesAmount("Portlandite"))*1000)
        Gypsums.append(float(props.speciesAmount("Gypsum"))*1000)

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
    Portlandites_conc = [ Portlandites[i] / vols[i] for i in range(len(vols)) ]

    # ax2 = ax1.twinx()
    ax1.plot(pHs, Portlandites, 'b-', label='Ca(OH)2')
    ax1.plot(pHs, Brucites, 'r-', label='Mg(OH)2')
    ax1.plot(pHs, Gypsums, 'g-', label='CaSO4')


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
Brucites, Portlandites, Mg_rr, pHs, vols, NaOHs = plot_step_1_permian(
                                                    maxNaOH_vol=maxNaOH_vol,
                                                    NaOH_conc = NaOH_conc,
                                                    RR = RR,
                                                    runs = 100,
                                                    )
#%% Step 1: Plot (add solid)
def plot_step_1_permian_solid(     
        maxNaOH_mass = 5, # g
        NaOH_conc = 'solid', # mol/L
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

    NaOHs = [i * maxNaOH_mass / runs for i in range(0,runs+1,1)]
    for v in NaOHs:
        ori_output, ori_state, ori_system, ori_props = permian_brine(RR=RR)

        output, props, outflow_mass, out_brucite, state, system = permian_step_1(
            state=ori_state, 
            system=ori_system,
            pressure = 1.01325, # bar
            add_NaOH_mass = v, # L
            add_NaOH_conc = 'solid', # mol/L
            temp = 25, 
            )

        Calcites.append(float(props.speciesAmount("Calcite"))*1000)
        Brucites.append(float(props.speciesAmount("Brucite"))*1000)
        Magnesites.append(float(props.speciesAmount("Magnesite"))*1000)
        Portlandites.append(float(props.speciesAmount("Portlandite"))*1000)
        Gypsums.append(float(props.speciesAmount("Gypsum"))*1000)

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
    Portlandites_conc = [ Portlandites[i] / vols[i] for i in range(len(vols)) ]

    # ax2 = ax1.twinx()
    ax1.plot(pHs, Portlandites, 'b-', label='Ca(OH)2')
    ax1.plot(pHs, Brucites, 'r-', label='Mg(OH)2')
    ax1.plot(pHs, Gypsums, 'g-', label='CaSO4')


    ax1.set_xlabel('pH')
    ax1.set_ylabel('Solid (mmol/L)', color='k')
    # ax2.set_ylabel('Vol (L)', color='k')

    ax1.legend(loc='upper center')
    # ax2.legend(loc = 'best')
    plt.show()

    print('Max vol: ', maxNaOH_vol, "Conc: ", NaOH_conc)
    return Brucites, Portlandites, Mg_rr, pHs, vols, NaOHs


# Run function
maxNaOH_mass = 10 # g
NaOH_conc = 'solid' # mol/L
RR = 0.5
Brucites, Portlandites, Mg_rr, pHs, vols, NaOHs = plot_step_1_permian_solid(
                                                    maxNaOH_mass=maxNaOH_mass,
                                                    NaOH_conc = NaOH_conc,
                                                    RR = RR,
                                                    runs = 50,
                                                    )

#%% Step 1 Export csv
results = {"added_NaOH (L)": NaOHs,
        "Mg(OH)2 (mmol/L)": Brucites,
        "Ca(OH)2 (mmol/L)":Portlandites, 
        "Mg recovery":Mg_rr, 
        "pH":pHs, 
        "Volume (L)":vols,
        }
     
df = pd.DataFrame(results)
outpath = '/Users/zhuoranzhang/Documents/Crystallization_paper/ZLD_profile_results/'
filename = f'Permian_S1_RR{RR}_NaOH{NaOH_conc}M.csv'

df.to_csv(outpath+filename)
print(f'csv created for RR={RR} and NaOH conc of {NaOH_conc}')
#%% Step 2: Add CO2 and NaOH to precipitate Ca
def permian_step_2(
        state, system,
        pressure = 1.01325, # bar
        add_NaOH_vol = 0.0365, # L
        add_NaOH_conc = 5, # mol/L
        CO2_fug = 0, #atm
        temp = 25, 
    ): 
    # db = PhreeqcDatabase("thermoddem-v1.10.dat")
    db = PhreeqcDatabase("pitzer.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg S"))
    solution.set(ActivityModelPitzer())

    # Create a gaseous phase
    gaseousphase = GaseousPhase("CO2(g)")
    gaseousphase.set(ActivityModelPengRobinson())

    # Hydromagnesite = MineralPhase("Hydromagnesite")
    # Magnesite = MineralPhase("Magnesite"

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    # specs.phaseAmount("GaseousPhase")
    # specs.pH()
    specs.charge()
    specs.openTo("Cl-")
    # specs.openTo("OH-")
    specs.fugacity("CO2(g)")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)

    # Additional chemicals (5 M NaOH)
    NaOH_props, masses = NaOH(conc = add_NaOH_conc)

    H2O_mass, Na_mass, OH_mass = masses # g/L

    state.add("Na+"  ,  Na_mass * add_NaOH_vol / 1000, "kg")
    state.add("OH-"  ,  OH_mass * add_NaOH_vol / 1000, "kg")
    state.add("H2O"  ,  H2O_mass * add_NaOH_vol / 1000, "kg")
    # state.add("CO2(g)"  ,  add_CO2, "mg")

    # equilibrate(state)

    conditions = EquilibriumConditions(specs)
    conditions.temperature(temp, "celsius")
    # conditions.pressure(pressure, "bar")

    # conditions.phaseAmount("GaseousPhase", 0.0001, "umol")
    # conditions.setLowerBoundPressure(0.01, "bar")
    # conditions.setUpperBoundPressure(1000.0, "bar")
    conditions.setLowerBoundTemperature(25, "celsius")
    conditions.setUpperBoundTemperature(200, "celsius")
    # conditions.pH(8.2)
    conditions.charge(0.0)
    conditions.fugacity("CO2(g)", CO2_fug * 1.01325, "bar")

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
    out_brucite = float(props.speciesAmount("Brucite")) * 58.3197 * 1000 # mg
    outflow_species = ["H2O", "Ca+2", "Cl-", "SO4-2",
                       "Mg+2", "Na+", "OH-", "HCO3-", "CO3-2"]
    outflow_mass = dict(zip(outflow_species,
                            [float(props.speciesAmount(i)) * molar_mass[i] for i in outflow_species])) # mg

    return output, props, out_brucite, state, system


# Run for case
permian_brine_output, permian_brine_state, permian_brine_system, permian_brine_props = permian_brine(RR=0.5)
output_1_PB, props_1_PB, outflow_mass_1_PB, out_brucite_1_PB, state_1_PB, system_1_PB = permian_step_1(state = permian_brine_state,
                                                                                     system=permian_brine_system,
                                                                                     add_NaOH_conc=5,
                                                                                     add_NaOH_vol=0.025,
                                                                             )
output_2_PB, props_2_PB, out_brucite_2, state_2_PB, system_2_PB = permian_step_2(
                                                                             state= state_1_PB,
                                                                             system=system_1_PB,
                                                                             add_NaOH_conc=5,
                                                                             add_NaOH_vol=0.052,
                                                                             CO2_fug=0.0000000001,
                                                                             )
print('Step 2 simulation')
print('pH=', output_2_PB['pH'])
print('Ca precipitated (g)', float(props_2_PB.speciesAmount('Calcite')) * molar_mass['Ca+2'])
print('Mg precipitated (g)', (float(props_2_PB.speciesAmount('Magnesite')) + float(props_2_PB.speciesAmount('Brucite')))* molar_mass['Mg+2'])
print('Ca rec rate', float(props_2_PB.speciesAmount('Calcite')) / float(props_2_PB.elementAmount('Ca')))

#%% Step 2: Plot
def plot_step_2(         
    maxNaOH_vol = 0.1, # L
    NaOH_conc = 5, # mol/L
    max_CO2 = 0.001, # bar 
    runs_NaOH = 20,
    runs_CO2 = 20):

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

    Ca_rr = []

    CO2s = [ i * max_CO2 / runs_CO2 for i in range(0,runs_CO2+1,1)]
    NaOHs = [i * maxNaOH_vol / runs_NaOH for i in range(0,runs_NaOH+1,1)]

    for i in NaOHs:
        brucites = []
        magnesites = []
        calcites = []

        temps = []
        pres = []
        phs = []

        nas = []
        mgs = []
        cas = []
        ca_rr =[]
        for j in CO2s:
            permian_brine_output, permian_brine_state, permian_brine_system, permian_brine_props = permian_brine(RR=0.5)
            output_1_PB, props_1_PB, outflow_mass_1_PB, out_brucite_1_PB, state_1_PB, system_1_PB = permian_step_1(state = permian_brine_state,
                                                                                                system=permian_brine_system,
                                                                                                add_NaOH_conc=5,
                                                                                                add_NaOH_vol=0.025,
                                                                                        )
            output_2_PB, props_2_PB, out_brucite_2, state_2_PB, system_2_PB = permian_step_2(
                                                                                        state= state_1_PB,
                                                                                        system=system_1_PB,
                                                                                        add_NaOH_conc=NaOH_conc,
                                                                                        add_NaOH_vol=i,
                                                                                        CO2_fug= j
                                                                                        )
            
            
            brucites.append(float(props_2_PB.speciesAmount("Brucite")))
            calcites.append(float(props_2_PB.speciesAmount("Calcite")))
            magnesites.append(float(props_2_PB.speciesAmount("Magnesite")))

            temps.append(output_2_PB['temperature'])
            pres.append(output_2_PB['pressure'])
            phs.append(output_2_PB['pH'])

            cas.append(float(props_2_PB.speciesAmount("Ca+2")))
            nas.append(float(props_2_PB.speciesAmount("Na+")))
            mgs.append(float(props_2_PB.speciesAmount("Mg+2")))

            ca_rr.append(float(props_2_PB.speciesAmount("Calcite"))/float(props_2_PB.elementAmount("Ca"))*100)

        Brucites.append(brucites)
        Calcites.append(calcites)
        Magnesites.append(magnesites)

        temperatures.append(temps)
        pressures.append(pres)
        pHs.append(phs)

        Cas.append(cas)
        Nas.append(nas)
        Mgs.append(mgs)

        Ca_rr.append(ca_rr)

    return CO2s, NaOHs, pHs, Ca_rr, Magnesites, Brucites
    
# make plot
CO2s, NaOHs, pHs, Ca_rr, Magnesites, Brucites = plot_step_2(      
                                        maxNaOH_vol = 0.1, # L
                                        NaOH_conc = 5, # mol/L
                                        max_CO2 = 0.001, # bar 
                                        runs_NaOH = 20,
                                        runs_CO2 = 20
                                        )

# Plot pH
fig,  axes= plt.subplots(1,1)
ax1 = axes

im1 = ax1.imshow([pHs[i-1] for i in range(len(pHs),0,-1)], cmap='YlGnBu', interpolation='none', vmax=9.5)
# ax1.set_aspect(1)
ax1.set(xticks=np.arange(0, len(CO2s)+1, 5), xticklabels=np.arange(0, np.max(CO2s), np.max(CO2s)/5))
ax1.set(yticks=np.arange(0, len(NaOHs)+1, 5), yticklabels=[0.05,0.04,0.03,0.02,0.01])
ax1.set_xlabel('CO2 pressure (bar)')
ax1.set_ylabel('5M NaOH addition (L)')
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im1,  ticks = [7, 7.5, 8, 8.5, 9, 9.5])
cbar.set_label('pH')

# Plot Ca_rr
fig,  axes= plt.subplots(1,1)
ax1 = axes

im1 = ax1.imshow([Ca_rr[i-1] for i in range(len(pHs),0,-1)], cmap=plt.cm.Reds, interpolation='none')
# ax1.set_aspect(1)
# ax1.set(xticks=np.arange(0, len(CO2s)+1, 5), xticklabels=np.arange(0, np.max(CO2s)+1, np.max(CO2s)/5))
# ax1.set(yticks=np.arange(0, len(NaOHs)+1, 5), yticklabels=[0.05,0.04,0.03,0.02,0.01,0])
ax1.set_xlabel('CO2 addition (mg)')
ax1.set_ylabel('5M NaOH addition (L)')
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im1,  ticks = [0,10,20,30,40,50,60,70,80,90,100])
cbar.set_label('Ca recovery ratio')

# Plot MgCO3
fig,  axes= plt.subplots(1,1)
ax1 = axes

im1 = ax1.imshow([Magnesites[i-1] for i in range(len(pHs),0,-1)], cmap="YlOrBr", interpolation='none')
# ax1.set_aspect(1)
# ax1.set(xticks=np.arange(0, len(CO2s)+1, 5), xticklabels=np.arange(0, np.max(CO2s)+1, np.max(CO2s)/5))
# ax1.set(yticks=np.arange(0, len(NaOHs)+1, 5), yticklabels=[0.05,0.04,0.03,0.02,0.01,0])
ax1.set_xlabel('CO2 addition (mg)')
ax1.set_ylabel('5M NaOH addition (L)')
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im1,  ticks = [0, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006 ])
cbar.set_label('MgCO3 (mol)')

# Plot Mg(OH)2
fig,  axes= plt.subplots(1,1)
ax1 = axes

im1 = ax1.imshow([Brucites[i-1] for i in range(len(pHs),0,-1)], cmap="YlOrBr", interpolation='none')
# ax1.set_aspect(1)
# ax1.set(xticks=np.arange(0, len(CO2s)+1, 5), xticklabels=np.arange(0, np.max(CO2s)+1, np.max(CO2s)/5))
# ax1.set(yticks=np.arange(0, len(NaOHs)+1, 5), yticklabels=[0.05,0.04,0.03,0.02,0.01,0])
ax1.set_xlabel('CO2 addition (mg)')
ax1.set_ylabel('5M NaOH addition (L)')
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im1,  ticks = [0, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007 ])
cbar.set_label('Mg(OH)2 (mol)')


#%% Step 2: Use H2CO3
def permian_step_2_alt(
        state, system,
        pressure = 1.01325, # bar
        add_NaOH_vol = 0.0365, # L
        add_NaOH_conc = 5, # mol/L
        add_CO2 = 0, #atm
        temp = 25, 
    ): 
    # db = PhreeqcDatabase("thermoddem-v1.10.dat")
    db = PhreeqcDatabase("pitzer.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg S"))
    solution.set(ActivityModelPitzer())

    # Create a gaseous phase
    gaseousphase = GaseousPhase("CO2(g)")
    gaseousphase.set(ActivityModelPengRobinson())

    # Hydromagnesite = MineralPhase("Hydromagnesite")
    # Magnesite = MineralPhase("Magnesite"

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    # specs.phaseAmount("GaseousPhase")
    # specs.pH()
    specs.charge()
    specs.openTo("Cl-")
    # specs.openTo("OH-")
    # specs.fugacity("CO2(g)")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)

    # Additional chemicals (5 M NaOH)
    NaOH_props, masses = NaOH(conc = add_NaOH_conc)

    H2O_mass, Na_mass, OH_mass = masses # g/L

    # state.add("CO2(g)"  ,  add_CO2, "mg")

    add_H2CO3 = add_CO2 / 44.009 * 62.03 # mg
    remove_H2O = add_CO2 /44.009 * 18.021 # mg
    
    state.add("H2CO3", add_H2CO3, "mg")

    state.add("Na+"  ,  Na_mass * add_NaOH_vol / 1000, "kg")
    state.add("OH-"  ,  OH_mass * add_NaOH_vol / 1000, "kg")
    state.add("H2O"  ,  H2O_mass * add_NaOH_vol / 1000 - remove_H2O / 1e6, "kg")

    state.set("Brucite", 0, "kg")

    # equilibrate(state)

    conditions = EquilibriumConditions(specs)
    conditions.temperature(temp, "celsius")
    # conditions.pressure(pressure, "bar")

    # conditions.phaseAmount("GaseousPhase", 0.0001, "umol")
    # conditions.setLowerBoundPressure(0.01, "bar")
    # conditions.setUpperBoundPressure(1000.0, "bar")
    conditions.setLowerBoundTemperature(25, "celsius")
    conditions.setUpperBoundTemperature(200, "celsius")
    # conditions.pH(8.2)
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
    outflow_species = ["H2O", "Ca+2", "Cl-", "SO4-2",
                       "Mg+2", "Na+", "OH-", "HCO3-", "CO3-2"]
    outflow_mass = dict(zip(outflow_species,
                            [float(props.speciesAmount(i)) * molar_mass[i] for i in outflow_species])) # mg

    return output, props, state, system


# Run for case
permian_brine_output, permian_brine_state, permian_brine_system, permian_brine_props = permian_brine(RR=0.5)
output_1_PB, props_1_PB, outflow_mass_1_PB, out_brucite_1_PB, state_1_PB, system_1_PB = permian_step_1(state = permian_brine_state,
                                                                                     system=permian_brine_system,
                                                                                     add_NaOH_conc=5,
                                                                                     add_NaOH_vol=0.025,
                                                                             )
output_2_PB, props_2_PB, state_2_PB, system_2_PB = permian_step_2_alt(
                                                                    state= state_1_PB,
                                                                    system=system_1_PB,
                                                                    add_NaOH_conc=5,
                                                                    add_NaOH_vol=0.08,
                                                                    add_CO2 = 9000,
                                                                #  CO2_fug=0.0000000001,
                                                                    )
print('Step 2 simulation')
print('pH=', output_2_PB['pH'])
print('Ca precipitated (g)', float(props_2_PB.speciesAmount('Calcite')) * molar_mass['Ca+2'])
print('Mg precipitated (g)', (float(props_2_PB.speciesAmount('Magnesite')) + float(props_2_PB.speciesAmount('Brucite')))* molar_mass['Mg+2'])
print('Ca rec rate', float(props_2_PB.speciesAmount('Calcite')) / float(props_2_PB.elementAmount('Ca')))
total_solid = (float(props_2_PB.speciesAmount("Calcite"))*molar_mass["Calcite"]
                + float(props_2_PB.speciesAmount("Brucite"))*molar_mass["Brucite"]
                + float(props_2_PB.speciesAmount("Magnesite"))* molar_mass["Magnesite"]
)
print("Ca purity (%)", float(props_2_PB.speciesAmount("Calcite"))*molar_mass["Calcite"]/total_solid*100)

#%% Step 2: Plot
def plot_step_2_alt(         
    maxNaOH_vol = 0.1, # L
    NaOH_conc = 5, # mol/L
    max_CO2 = 0.001, # bar 
    runs_NaOH = 20,
    runs_CO2 = 20):

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

    Ca_rr = []
    Ca_Purity = []
    Vols =[]

    CO2s = [ i * max_CO2 / runs_CO2 for i in range(0,runs_CO2+1,1)]
    NaOHs = [i * maxNaOH_vol / runs_NaOH for i in range(0,runs_NaOH+1,1)]

    for i in NaOHs:
        brucites = []
        magnesites = []
        calcites = []
        portlandites = []

        temps = []
        pres = []
        phs = []
        vols = []

        nas = []
        mgs = []
        cas = []
        ca_rr =[]
        ca_purity=[]
        for j in CO2s:
            permian_brine_output, permian_brine_state, permian_brine_system, permian_brine_props = permian_brine(RR=0.5)
            output_1_PB, props_1_PB, outflow_mass_1_PB, out_brucite_1_PB, state_1_PB, system_1_PB = permian_step_1(state = permian_brine_state,
                                                                                                system=permian_brine_system,
                                                                                                add_NaOH_conc=5,
                                                                                                add_NaOH_vol=0.025,
                                                                                        )
            output_2_PB, props_2_PB, state_2_PB, system_2_PB = permian_step_2_alt(
                                                                                        state= state_1_PB,
                                                                                        system=system_1_PB,
                                                                                        add_NaOH_conc=NaOH_conc,
                                                                                        add_NaOH_vol=i,
                                                                                        add_CO2 = j
                                                                                        )
            
            
            brucites.append(float(props_2_PB.speciesAmount("Brucite"))*1000) # mol convert to mmol
            calcites.append(float(props_2_PB.speciesAmount("Calcite")))
            magnesites.append(float(props_2_PB.speciesAmount("Magnesite")))
            portlandites.append(float(props_2_PB.speciesAmount("Portlandite")))

            vols.append(float(props_2_PB.phaseProps("AqueousPhase").volume())*1000)

            temps.append(output_2_PB['temperature'])
            pres.append(output_2_PB['pressure'])
            phs.append(output_2_PB['pH'])

            cas.append(float(props_2_PB.speciesAmount("Ca+2")))
            nas.append(float(props_2_PB.speciesAmount("Na+")))
            mgs.append(float(props_2_PB.speciesAmount("Mg+2")))

            ca_rr.append(float(props_2_PB.speciesAmount("Calcite"))/float(props_2_PB.elementAmount("Ca"))*100)

            total_solid = (float(props_2_PB.speciesAmount("Calcite"))*molar_mass["Calcite"]
                           + float(props_2_PB.speciesAmount("Brucite"))*molar_mass["Brucite"]
                           + float(props_2_PB.speciesAmount("Magnesite"))* molar_mass["Magnesite"]
            )
            ca_purity.append(float(props_2_PB.speciesAmount("Calcite"))*molar_mass["Calcite"]/total_solid*100)


        Brucites.append(brucites)
        Calcites.append(calcites)
        Magnesites.append(magnesites)
        Portlandites.append(portlandites)

        temperatures.append(temps)
        pressures.append(pres)
        pHs.append(phs)
        Vols.append(vols)

        Cas.append(cas)
        Nas.append(nas)
        Mgs.append(mgs)

        Ca_rr.append(ca_rr)
        Ca_Purity.append(ca_purity)

    return CO2s, NaOHs, pHs, Ca_rr, Magnesites, Brucites, Calcites, Portlandites, Ca_Purity, Vols
    
# make plot
CO2s, NaOHs, pHs, Ca_rr, Magnesites, Brucites, Calcites, Portlandites, Ca_Purity, Vols = plot_step_2_alt(      
                                        maxNaOH_vol = 0.1, # L
                                        NaOH_conc = 5, # mol/L
                                        max_CO2 =  15000, # mg 
                                        runs_NaOH = 20,
                                        runs_CO2 = 20
                                        )

CO2s = [i/1000 for i in CO2s]

# Plot pH
fig,  axes= plt.subplots(1,1)
ax1 = axes


im1 = ax1.imshow([pHs[i-1] for i in range(len(pHs),0,-1)], cmap='YlGnBu', interpolation='none', vmax =13)
# ax1.set_aspect(1)
ax1.set(xticks=np.arange(0, len(CO2s)+1, 4), xticklabels=np.arange(0, int(np.max(CO2s))+1, int(np.max(CO2s)/5)))
ax1.set(yticks=np.arange(0, len(NaOHs)+1, 4), yticklabels=[0.1, 0.08, 0.06, 0.04, 0.02, 0])

ax1.set_xlabel('CO2 addition (g)')
ax1.set_ylabel('5M NaOH addition (L)')
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im1,  ticks = [4,5,6,7,8,9,10,11,12,13])
cbar.set_label('pH')
plt.rcParams['figure.dpi']=300
plt.show()

# Plot Ca_rr
fig,  axes= plt.subplots(1,1)
ax1 = axes

im1 = ax1.imshow([Ca_rr[i-1] for i in range(len(pHs),0,-1)], cmap=plt.cm.Reds, interpolation='none')
# ax1.set_aspect(1)
ax1.set(xticks=np.arange(0, len(CO2s)+1, 4), xticklabels=np.arange(0, int(np.max(CO2s))+1, int(np.max(CO2s)/5)))
ax1.set(yticks=np.arange(0, len(NaOHs)+1, 4), yticklabels=[0.1, 0.08, 0.06, 0.04, 0.02, 0])
ax1.set_xlabel('CO2 addition (g)')
ax1.set_ylabel('5M NaOH addition (L)')
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im1,  ticks = [0,10,20,30,40,50,60,70,80,90,100])
cbar.set_label('Ca recovery ratio (%)')
plt.rcParams['figure.dpi']=300
plt.show()

# # Plot MgCO3
# fig,  axes= plt.subplots(1,1)
# ax1 = axes

# im1 = ax1.imshow([Magnesites[i-1] for i in range(len(pHs),0,-1)], cmap="YlOrBr", interpolation='none')
# # ax1.set_aspect(1)
# ax1.set(xticks=np.arange(0, len(CO2s)+1, 4), xticklabels=np.arange(0, int(np.max(CO2s))+1, int(np.max(CO2s)/5)))
# ax1.set(yticks=np.arange(0, len(NaOHs)+1, 4), yticklabels=[0.1, 0.08, 0.06, 0.04, 0.02, 0])
# ax1.set_xlabel('CO2 addition (g)')
# ax1.set_ylabel('5M NaOH addition (L)')
# # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
# cbar = fig.colorbar(im1)
# cbar.set_label('MgCO3 (mol)')
# plt.rcParams['figure.dpi']=300
# plt.show()

# Plot Mg(OH)2
fig,  axes= plt.subplots(1,1)
ax1 = axes

im1 = ax1.imshow([Brucites[i-1] for i in range(len(pHs),0,-1)], cmap="YlOrBr", interpolation='none')
# ax1.set_aspect(1)
ax1.set(xticks=np.arange(0, len(CO2s)+1, 4), xticklabels=np.arange(0, np.max(CO2s)*1.01, int(np.max(CO2s)/5)))
ax1.set(yticks=np.arange(0, len(NaOHs)+1, 4), yticklabels=[0.1, 0.08, 0.06, 0.04, 0.02, 0])
ax1.set_xlabel('CO2 addition (mg)')
ax1.set_ylabel('5M NaOH addition (L)')
# cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im1)
cbar.set_label('Mg(OH)2 (mmol)')
plt.rcParams['figure.dpi']=300
plt.show()

# # Plot purity of CaCO3
# fig,  axes= plt.subplots(1,1)
# ax1 = axes

# im1 = ax1.imshow([Ca_Purity[i-1] for i in range(len(pHs),0,-1)], cmap="YlOrBr", interpolation='none')
# # ax1.set_aspect(1)
# ax1.set(xticks=np.arange(0, len(CO2s)+1, 4), xticklabels=np.arange(0, np.max(CO2s)*1.01, int(np.max(CO2s)/5)))
# ax1.set(yticks=np.arange(0, len(NaOHs)+1, 4), yticklabels=[0.1, 0.08, 0.06, 0.04, 0.02, 0])
# ax1.set_xlabel('CO2 addition (mg)')
# ax1.set_ylabel('5M NaOH addition (L)')
# # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
# cbar = fig.colorbar(im1)
# cbar.set_label('Calcite purity (%)')
# plt.rcParams['figure.dpi']=300
# plt.show()

# #%%  best rr=99.88% (9000mg CO2 + 0.08L NaOH)
# i, j = 12, 16
# print('CO2 (mg):  ', CO2s[i])
# print('NaOH(L):   ', round(NaOHs[j],3))
# print('pH:        ', round(pHs[j][i],2))
# print('Ca rr:     ', round(Ca_rr[j][i],2), '%')
# print('MgCO3:     ', round(Magnesites[j][i],4), 'mol')
# print('Mg(OH)2:   ', round(Brucites[j][i],4), 'mol')
#%% Step 2 Export csv
pH_df = pd.DataFrame(pHs, columns=CO2s, index=NaOHs)
Ca_df = pd.DataFrame(Ca_rr, columns=CO2s, index=NaOHs)
Mag_df = pd.DataFrame(Magnesites, columns=CO2s, index=NaOHs)
Bru_df = pd.DataFrame(Brucites, columns=CO2s, index=NaOHs)
Cal_df = pd.DataFrame(Calcites, columns=CO2s, index=NaOHs)
Por_df = pd.DataFrame(Portlandites, columns=CO2s, index=NaOHs)
Ca_Purity_df = pd.DataFrame(Ca_Purity, columns=CO2s, index=NaOHs)
Vol_df = pd.DataFrame(Vols, columns=CO2s, index=NaOHs)
    
outpath = '/Users/zhuoranzhang/Documents/Crystallization_paper/ZLD_profile_results/'
pH_filename = f'Permian_S2_pH.csv'
Ca_filename = f'Permian_S2_Ca_rr.csv'
Mag_filename = f'Permian_S2_MgCO3.csv'
Bru_filename = f'Permian_S2_MgOH2.csv'
Cal_filename = f'Permian_S2_CaCO3.csv'
Por_filename = f'Permian_S2_CaOH2.csv'
Ca_purity_filename =  f'Permian_S2_purity.csv'
Vol_filename = f'Permian_S2_volume.csv'

pH_df.to_csv(outpath+pH_filename)
Ca_df.to_csv(outpath+Ca_filename)
Mag_df.to_csv(outpath+Mag_filename)
Bru_df.to_csv(outpath+Bru_filename)
Cal_df.to_csv(outpath+Cal_filename)
Por_df.to_csv(outpath+Por_filename)
Ca_Purity_df.to_csv(outpath+Ca_purity_filename)
Vol_df.to_csv(outpath+Vol_filename)
print(f'csv created')
#%% Step 3: Add CO2 to carbonize Mg
def permian_step_3_alt(
        state, system,
        pressure = 1.01325, # bar
        add_NaOH_vol = 0.0365, # L
        add_NaOH_conc = 5, # mol/L
        add_CO2 = 0, #atm
        add_brucite =0,
        temp = 25, 
    ): 
    # db = PhreeqcDatabase("thermoddem-v1.10.dat")
    # db = PhreeqcDatabase("pitzer.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg S"))
    solution.set(ActivityModelPitzer())

    # Create a gaseous phase
    gaseousphase = GaseousPhase("CO2(g)")
    gaseousphase.set(ActivityModelPengRobinson())

    # Hydromagnesite = MineralPhase("Hydromagnesite")
    # Magnesite = MineralPhase("Magnesite"

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    # specs.phaseAmount("GaseousPhase")
    # specs.pH()
    specs.charge()
    specs.openTo("Cl-")
    # specs.openTo("OH-")
    # specs.fugacity("CO2(g)")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)

    # Additional chemicals (5 M NaOH)
    NaOH_props, masses = NaOH(conc = add_NaOH_conc)

    H2O_mass, Na_mass, OH_mass = masses # g/L

    # state.add("CO2(g)"  ,  add_CO2, "mg")

    add_H2CO3 = add_CO2 / 44.009 * 62.03 # mg
    remove_H2O = add_CO2 /44.009 * 18.021 # mg
    
    state.add("H2CO3", add_H2CO3, "mg")

    state.add("Na+"  ,  Na_mass * add_NaOH_vol / 1000, "kg")
    state.add("OH-"  ,  OH_mass * add_NaOH_vol / 1000, "kg")
    state.add("H2O"  ,  H2O_mass * add_NaOH_vol / 1000 - remove_H2O / 1e6, "kg")

    state.add("Brucite", add_brucite, "mg")
    state.set("Calcite", 0, "mg")

    # equilibrate(state)

    conditions = EquilibriumConditions(specs)
    conditions.temperature(temp, "celsius")
    # conditions.pressure(pressure, "bar")

    # conditions.phaseAmount("GaseousPhase", 0.0001, "umol")
    # conditions.setLowerBoundPressure(0.01, "bar")
    # conditions.setUpperBoundPressure(1000.0, "bar")
    conditions.setLowerBoundTemperature(25, "celsius")
    conditions.setUpperBoundTemperature(200, "celsius")
    # conditions.pH(8.2)
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
              'Amount': state.speciesAmounts().asarray(),
              'Volume': float(props.phaseProps("AqueousPhase").volume())
              }
    out_brucite = float(props.speciesAmount("Brucite")) * 58.3197 * 1000 # mg
    outflow_species = ["H2O", "Ca+2", "Cl-", "SO4-2",
                       "Mg+2", "Na+", "OH-", "HCO3-", "CO3-2"]
    outflow_mass = dict(zip(outflow_species,
                            [float(props.speciesAmount(i)) * molar_mass[i] for i in outflow_species])) # mg

    return output, props, outflow_mass, state, system


# Run for case
permian_brine_output, permian_brine_state, permian_brine_system, permian_brine_props = permian_brine(RR=0.5)
output_1_PB, props_1_PB, outflow_mass_1_PB, out_brucite_1_PB, state_1_PB, system_1_PB = permian_step_1(state = permian_brine_state,
                                                                                     system=permian_brine_system,
                                                                                     add_NaOH_conc=5,
                                                                                     add_NaOH_vol=0.025,
                                                                             )
output_2_PB, props_2_PB, state_2_PB, system_2_PB = permian_step_2_alt(
                                                                             state= state_1_PB,
                                                                             system=system_1_PB,
                                                                             add_NaOH_conc=5,
                                                                             add_NaOH_vol=0.08,
                                                                             add_CO2 = 9000,
                                                                            #  CO2_fug=0.0000000001,
                                                                             )
output_3_PB, props_3_PB, outflow_mass_3, state_3_PB, system_3_PB = permian_step_3_alt(
                                                                             state= state_2_PB,
                                                                             system=system_2_PB,
                                                                             add_CO2 = 2600,
                                                                             add_NaOH_conc=5,
                                                                             add_NaOH_vol=0,
                                                                             add_brucite=out_brucite_1_PB,
                                                                            #  CO2_fug=0.0000000001,
                                                                             )
print('Step 3 simulation')
print('pH: ', output_3_PB['pH'])
print('Mg in MgCO3:', float(props_3_PB.speciesAmount("Magnesite")) * molar_mass["Mg+2"])
print('Mg in MgOH2:', float(props_3_PB.speciesAmount("Brucite")) * molar_mass["Mg+2"])
print('MgCO3(%)', float(props_3_PB.speciesAmount("Magnesite"))/(float(props_3_PB.elementAmount("Mg")))*100)
print('total Mg:', (float(props_3_PB.elementAmount("Mg"))) * molar_mass["Mg+2"])
print('Ca in CaCO3:', float(props_3_PB.speciesAmount("Calcite")) * molar_mass["Ca+2"])
print('total Ca: ', float(props_3_PB.elementAmount("Ca")) * molar_mass["Ca+2"])
print('TDS (mg/L): ', (outflow_mass_3['Na+'] + outflow_mass_3['Cl-'] )/ output_3_PB['Volume'] / 1000)

#%% Step 3: Plot
def plot_step_3_alt(
    max_CO2 = 10000,
    runs_CO2 = 20,
    NaOH_vol = 0,
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
    Hydromagnesites = []
    Periclases = []
    vols = []
    Nas = []
    Mgs = []
    Cas = []
    CO3s = []
    HCO3s = []
    H2CO3s = []

    Mg_rr = []

    CO2s = [ i * max_CO2 / runs_CO2 for i in range(0,runs_CO2+1,1)]
    # NaOHs = [i * 0.001 for i in range(1,50,1)]

    # for i in CO2s:
    for j in CO2s:
        permian_brine_output, permian_brine_state, permian_brine_system, permian_brine_props = permian_brine(RR=0.5)
        output_1_PB, props_1_PB, outflow_mass_1_PB, out_brucite_1_PB, state_1_PB, system_1_PB = permian_step_1(state = permian_brine_state,
                                                                                            system=permian_brine_system,
                                                                                            add_NaOH_conc=5,
                                                                                            add_NaOH_vol=0.025,
                                                                                    )
        output_2_PB, props_2_PB, state_2_PB, system_2_PB = permian_step_2_alt(
                                                                                    state= state_1_PB,
                                                                                    system=system_1_PB,
                                                                                    add_NaOH_conc=5,
                                                                                    add_NaOH_vol=0.08,
                                                                                    add_CO2 = 9000,
                                                                                    #  CO2_fug=0.0000000001,
                                                                                    )

        output_3_PB, props_3_PB, out_brucite_3, state_3_PB, system_3_PB = permian_step_3_alt(
                                                                                    state= state_2_PB,
                                                                                    system=system_2_PB,
                                                                                    add_CO2 = j,
                                                                                    add_NaOH_conc=5,
                                                                                    add_NaOH_vol=NaOH_vol,
                                                                                    add_brucite=out_brucite_1_PB,
                                                                                    #  CO2_fug=0.0000000001,
                                                                                    )
        # rrs.append(i)

        vol = float(props_3_PB.phaseProps("AqueousPhase").volume())

        Brucites.append(float(props_3_PB.speciesAmount("Brucite")) / vol)
        # Hydromagnesites.append(float(props.speciesAmount("Hydromagnesite")))
        Magnesites.append(float(props_3_PB.speciesAmount("Magnesite")) / vol)
        Calcites.append(float(props_3_PB.speciesAmount("Calcite")) / vol)
        CO3s.append(float(props_3_PB.speciesAmount("CO3-2")) / vol)
        HCO3s.append(float(props_3_PB.speciesAmount("HCO3-")) / vol)
        H2CO3s.append(float(props_3_PB.speciesAmount("H2CO3")) / vol)
        # Periclases.append(float(props.speciesAmount("Periclase")))

        temperatures.append(output_3_PB['temperature'])
        pressures.append(output_3_PB['pressure'])
        pHs.append(output_3_PB['pH'])

        vols.append(float(props_3_PB.phaseProps("AqueousPhase").volume()))

        Cas.append(float(props_3_PB.speciesAmount("Ca+2")))
        Nas.append(float(props_3_PB.speciesAmount("Na+")))
        Mgs.append(float(props_3_PB.speciesAmount("Mg+2")))

        Mg_rr.append(props_3_PB.speciesAmount("Magnesite")/float(props_3_PB.elementAmount("Mg"))*100)

    fig, ax1 = plt.subplots()

    # ax2 = ax1.twinx()
    ax1.plot(CO2s, Magnesites, 'b-', label='MgCO3')
    # ax1.plot(CO2s, Calcites, 'r-', label='CaCO3')
    ax1.plot(CO2s, Brucites, 'g-', label='Mg(OH)2')
    ax1.plot(CO2s, CO3s, 'r-', label="CO3")
    ax1.plot(CO2s, HCO3s, 'y-', label="HCO3")
    ax1.plot(CO2s, H2CO3s, 'p-', label="H2CO3")

    # ax2.plot(CO2s, pHs, 'y-', label='pH')

    ax1.set_xlabel('Added CO2 (mg)')
    ax1.set_ylabel('Solids (mol)', color='k')
    # ax2.set_ylabel('pH', color='k')

    ax1.legend(loc='upper right')
    # ax2.legend(loc = 'best')
    plt.rcParams['figure.dpi']=300
    plt.show()

    # Plot Mg recovery ratio
    fig, ax1 = plt.subplots()
    ax1.plot(CO2s, Mg_rr, label = "Mg carbonization rate (%)")

    ax1.set_xlabel('Added CO2 (mg)')
    ax1.set_ylabel('Mg precipitation rate (%)', color='k')
    plt.rcParams['figure.dpi']=300
    plt.show()

    return Mgs, Magnesites, Brucites, Calcites, Mg_rr, pHs, vols,  CO3s, HCO3s, H2CO3s


NaOH_vol = 0.02
max_CO2 = 10000
runs_CO2 = 50
CO2s = [ i * max_CO2 / runs_CO2 for i in range(0,runs_CO2+1,1)]

Mgs_3, Magnesites_3, Brucites_3, Calcites_3, Mg_rr_3, pHs_3, vols_3,CO3s_3, HCO3s_3, H2CO3s_3 = plot_step_3_alt(
                                                                max_CO2=max_CO2,
                                                                runs_CO2=runs_CO2,
                                                                NaOH_vol = NaOH_vol,
                                                                )
#%% Step 3 Export csv
results = {"added_CO2 (mg)": CO2s,
        "MgCO3 (mol/L)": Magnesites_3,
        "CaCO3 (mol/L)":Calcites_3,
        "Mg(OH)2 (mol/L)":Brucites_3,  
        "Mg carbonization rate (%)":Mg_rr_3, 
        "CO3": CO3s_3,
        "HCO3": HCO3s_3,
        "H2CO3": H2CO3s_3,
        "pH":pHs_3, 
        "Volume (L)":vols_3,
        }
     
df = pd.DataFrame(results)
outpath = '/Users/zhuoranzhang/Documents/Crystallization_paper/ZLD_profile_results/'
filename = f'Permian_S3_NaOH{NaOH_vol}L.csv'

df.to_csv(outpath+filename)
print(f'csv created for step 3 at {filename}')


# %%
