#%%
from reaktoro import *
from reaktplot import *
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


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
}
#%% Run for CO2 addition
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
state.temperature(25, "celsius")

# Seawater composition (mg/kg)
Na  = 18738
Cl  = 27695
SO4 = 4096
Mg = 1913
Ca = 614
K  = 738
CO2 = 2500


# state.pressure(1.0, "atm")
state.set("H2O", 1.0, "kg")
state.add("Ca+2"   ,   Ca , "mg")
state.add("Mg+2"   ,  Mg , "mg")
state.add("Na+"    ,  Na , "mg")
state.add("K+"     ,   K , "mg")
state.add("Cl-"    , Cl , "mg")
state.add("SO4-2"  ,  SO4 , "mg")
state.add("CO2"  ,  CO2, "mg")


# equilibrate(state)

conditions = EquilibriumConditions(specs)
conditions.temperature(25, "celsius")
conditions.pressure(1.01325, "bar")

# conditions.phaseAmount("GaseousPhase", 0.0001, "umol")
# conditions.setLowerBoundPressure(0.01, "bar")
# conditions.setUpperBoundPressure(1000.0, "bar")
# conditions.setLowerBoundTemperature(25, "celsius")
# conditions.setUpperBoundTemperature(200, "celsius")
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

#%% Original solution (pH=8.2)
def original_state(
        feed_props = {    
            "Na+"  : 17987,
            "Cl-"  : 32220,
            "SO4-2" : 4766,
            "Mg+2" : 2226,
            "Ca+2" : 714,
            "K+"  : 858,
            "HCO3-": 0,
            # "CO3-2": 0.6,
            "OH-": 0,
            },
        pressure = 1.01325, # bar
        add_NaOH_vol = 0, # L
        add_NaOH_conc = 5, # mol/L
        temp = 25, 
    ): 
    # db = PhreeqcDatabase("thermoddem-v1.10.dat")
    db = PhreeqcDatabase("pitzer.dat")
    # db = PhreeqcDatabase("minteq.v4.dat")
    # db = PhreeqcDatabase("phreeqc.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S"))
    solution.set(ActivityModelPitzer())
    # solution.set(ActivityModelDebyeHuckel())

    # Create a gaseous phase
    gaseousphase = GaseousPhase("CO2(g)")
    gaseousphase.set(ActivityModelPengRobinson())

    # Hydromagnesite = MineralPhase("Hydromagnesite")
    Magnesite = MineralPhase("Magnesite")
    
    Calcite = MineralPhase("Calcite")
    # Polyhalite = MineralPhase("Polyhalite")
    # Glauberite = MineralPhase("Glauberite")
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
                            # Polyhalite,
                            # Glauberite,
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
    specs.pH()
    specs.charge()
    specs.openTo("Cl-")
    # specs.openTo("OH-")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)
    

    state = ChemicalState(system)
    state.temperature(temp, "celsius")

    # Seawater composition (mg/kg)
    Na  = 17987
    Cl  = 32220
    SO4 = 4766
    Mg = 2226
    Ca = 714
    K  = 858
    # HCO3= 300

    state.set("H2O", 1.0, "kg")
    for ion, conc in feed_props.items():
        state.add(ion, conc, "mg")
    # state.pressure(1.0, "atm")
    # state.add("Ca+2"   ,   Ca , "mg")
    # state.add("Mg+2"   ,  Mg , "mg")
    # state.add("Na+"    ,  Na , "mg")
    # state.add("K+"     ,   K , "mg")
    # state.add("Cl-"    , Cl , "mg")
    # state.add("SO4-2"  ,  SO4 , "mg")
    # state.add("HCO3-"  ,  HCO3 , "mg")
    # equilibrate(state)

    conditions = EquilibriumConditions(specs)
    conditions.temperature(temp, "celsius")
    conditions.pressure(pressure, "bar")

    # conditions.phaseAmount("GaseousPhase", 0.0001, "umol")
    # conditions.setLowerBoundPressure(0.01, "bar")
    # conditions.setUpperBoundPressure(1000.0, "bar")
    conditions.setLowerBoundTemperature(25, "celsius")
    conditions.setUpperBoundTemperature(200, "celsius")
    conditions.pH(8.2)
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
                       "Mg+2", "Na+", "OH-", "K+", "HCO3-", "CO3-2"]
    outflow_mass = dict(zip(outflow_species,
                            [float(props.speciesAmount(i)) * molar_mass[i] for i in outflow_species])) # mg

    # print('Original feed pH:', float(aprops.pH()))

    return output, state, system, props

ori_output, ori_state, ori_system, ori_props = original_state()
print(ori_output['pH'])
#%% Step 1 using original solution
def step1_add_NaOH_to_original(
        state, system,
        pressure = 1.01325, # bar
        add_NaOH_vol = 0.0365, # L
        add_NaOH_conc = 5, # mol/L
        add_CO2 = 0,
        temp = 25, 
    ): 
    # db = PhreeqcDatabase("thermoddem-v1.10.dat")
    db = PhreeqcDatabase("pitzer.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S"))
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

    solver = EquilibriumSolver(specs)
    # solver = SmartEquilibriumSolver(specs)

    # Additional chemicals (5 M NaOH)
    NaOH_props, masses = NaOH(conc = add_NaOH_conc)

    H2O_mass, Na_mass, OH_mass = masses # g/L

    state.add("Na+"  ,  Na_mass * add_NaOH_vol / 1000, "kg")
    state.add("OH-"  ,  OH_mass * add_NaOH_vol / 1000, "kg")
    state.add("H2O"  ,  H2O_mass * add_NaOH_vol / 1000, "kg")
    # state.add("CO2"  ,  add_CO2, "mg")

    # equilibrate(state)

    conditions = EquilibriumConditions(specs)
    conditions.temperature(temp, "celsius")
    conditions.pressure(pressure, "bar")

    # conditions.phaseAmount("GaseousPhase", 0.0001, "umol")
    # conditions.setLowerBoundPressure(0.01, "bar")
    # conditions.setUpperBoundPressure(1000.0, "bar")
    conditions.setLowerBoundTemperature(25, "celsius")
    conditions.setUpperBoundTemperature(200, "celsius")
    # conditions.pH(8.2)
    conditions.charge(0.0)
    # conditions.fugacity("CO2(g)", 0.00038, "bar")

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
                       "Mg+2", "Na+", "OH-", "K+", "HCO3-", "CO3-2"]
    outflow_mass = dict(zip(outflow_species,
                            [float(props.speciesAmount(i)) * molar_mass[i] for i in outflow_species])) # mg

    return output, props, outflow_mass, out_brucite, state, system


ori_output, ori_state, ori_system, ori_props = original_state()
output_1, props_1, outflow_mass_1, out_brucite_1, state_1, system_1 = step1_add_NaOH_to_original(state= ori_state, system=ori_system,
                                                                             add_NaOH_conc=5,
                                                                             add_NaOH_vol=0.0365,
                                                                             add_CO2 =0,
                                                                             )
print('Step One simulation')
print('pH=', output_1['pH'])
print('Ca in CaCO3 (mg)', float(props_1.speciesAmount('Calcite')) * molar_mass['Ca+2']* 1000)
print('Ca+2 (mg)', float(props_1.speciesAmount('Ca+2')) * molar_mass['Ca+2']* 1000)
print('Mg in MgCO3 (mg)', (float(props_1.speciesAmount('Magnesite')))* molar_mass['Mg+2']* 1000 )
print('Mg in MgOH2 (mg)', (float(props_1.speciesAmount('Brucite')))* molar_mass['Mg+2'] * 1000 )
print('Mg+2 (mg)', (float(props_1.speciesAmount('Mg+2')))* molar_mass['Mg+2'] * 1000 )
print('Ca rec rate', float(props_1.speciesAmount('Calcite')) * molar_mass['Ca+2'] / 0.714)


#%% kinetic
from reaktoro import *
db = PhreeqcDatabase("pitzer.dat")
Abar = 6.0

def ratefn(props: ChemicalProps):
    aprops = AqueousProps(props)  # we need an AqueousProps object to compute the saturation ratio of the mineral
    k0 = pow(10.0, -0.21)  # the reaction rate constant at 25 °C from Palandri and Kharaka (2004)
    q = props.phaseProps("Brucite").volume()  # the current volume of the mineral (in m3)
    Omega = aprops.saturationRatio("Calcite")  # the current saturation ratio of the mineral Ω = IAP/K
    return q * Abar * k0 * (1 - Omega) 

params = Params.embedded("PalandriKharaka.yaml")

Brucite = MineralPhase("Brucite")
solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S"))
solution.set(ActivityModelPitzer())

system = ChemicalSystem(db,
    solution,
    Brucite,
    # GeneralReaction("Brucite = Ca+2 + 2*OH-").setRateModel(ratefn),  # define the reaction and set its reaction rate model!
    MineralReaction("Brucite").setRateModel(ReactionRateModelPalandriKharaka(params)),

)
state = ChemicalState(system)
state.temperature(25.0, "C")
state.pressure(1.0, "bar")
state.set("H2O", 1.0, "kg")

Na  = 17987
Cl  = 32220
SO4 = 4766
Mg = 2226
Ca = 714
K  = 858

state.add("Ca+2"   ,   Ca , "mg")
state.add("Mg+2"   ,  Mg , "mg")
state.add("Na+"    ,  Na , "mg")
state.add("K+"     ,   K , "mg")
state.add("Cl-"    , Cl , "mg")
state.add("SO4-2"  ,  SO4 , "mg")
# state.add("OH-",      3114, "mg")
NaOH_props, masses = NaOH(conc = 5)

H2O_mass, Na_mass, OH_mass = masses # g/L
add_NaOH_vol = 1 # L

state.add("Na+"  ,  Na_mass * add_NaOH_vol / 1000, "kg")
state.add("OH-"  ,  OH_mass * add_NaOH_vol / 1000, "kg")
state.add("H2O"  ,  H2O_mass * add_NaOH_vol / 1000, "kg")

state.scalePhaseVolume("Brucite", 0, "cm3") 
solver = KineticsSolver(system)  # the chemical kinetics solver
# solver = SmartKineticsSolver(system)
table = Table()  # used to create table of data for later output and plotting

dt = 60  # time step (in seconds)

# Initiate the time stepping for the kinetics modeling
for i in range(60):
    result = solver.solve(state, dt)  # compute the chemical state of the system after it reacted for given time length

    assert result.succeeded(), f"Calculation did not succeed at time step #{i}."

    props = state.props()  # get the current thermodynamic and chemical properties of the system

    table.column("Time")   << i*dt / 60  # from seconds to minutes
    table.column("Brucite") << props.speciesAmount("Brucite") 
    table.column("Mg+2")    << props.speciesAmount("Mg+2")  # in mol
    table.column("OH-")    << props.speciesAmount("OH-")  # in mol
print(table)
#%% Step 1 Plot using original feed
def plot_step_1_original(       
        maxNaOH_vol = 0.5, # L
        NaOH_conc = 1, # mol/L
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

    CO2s = [i * 250 for i in range(100)]

    NaOHs = [i * maxNaOH_vol / runs for i in range(0,runs+1,1)]
    # for i in CO2s:
    for v in NaOHs:
        ori_output, ori_state, ori_system, ori_props = original_state()

        output, props, outflow_mass, out_brucite, state, system = step1_add_NaOH_to_original(
            state=ori_state, 
            system=ori_system,
            pressure = 1.01325, # bar
            add_NaOH_vol = v, # L
            add_NaOH_conc = NaOH_conc, # mol/L
            add_CO2= 0, #mg
            temp = 25, 
                        )

        Calcites.append(float(props.speciesAmount("Calcite")))
        Brucites.append(float(props.speciesAmount("Brucite")))
        Magnesites.append(float(props.speciesAmount("Magnesite")))

        temperatures.append(output['temperature'])
        pressures.append(output['pressure'])
        pHs.append(output['pH'])

        # Nas.append(output['Amount'][-12] * 23 / (output['Amount'][1] * 18 / 1000) * 1000)
        # Cas.append(output['Amount'][4] * 40 / (output['Amount'][1] * 18 / 1000) * 1000)
        # Mgs.append(output['Amount'][9] * 24 / (output['Amount'][1] * 18 / 1000) * 1000)
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

        vols.append(float(props.phaseProps("AqueousPhase").volume()))

        Mg_rr.append(props.speciesAmount("Brucite")/float(props.elementAmount("Mg"))*100)

        # print(masses)
    fig, ax1 = plt.subplots()

    Ca_concs = [ Cas[i]*40.087/vols[i] for i in range(len(NaOHs))]
    Na_concs = [ Nas[i]*22.99/vols[i] for i in range(len(NaOHs))]
    Mg_concs = [ Mgs[i]*24.305/vols[i] for i in range(len(NaOHs))]

    calcites_conc = [ Calcites[i] / vols[i] for i in range(len(vols))]
    brucites_conc = [ Brucites[i] / vols[i] for i in range(len(vols))]
    magnesites_conc = [ Magnesites[i] / vols[i] for i in range(len(vols))]


    ax2 = ax1.twinx()
    ax1.plot(pHs, Ca_concs, 'b-', label='Ca+2')
    ax1.plot(pHs, Mg_concs, 'r-', label='Mg+2')
    ax2.plot(pHs, Na_concs, 'g-', label='Na+')

    # ax2.plot(NaOHs, pHs, 'k-', label='pH')

    ax1.set_xlabel('pH')
    ax1.set_ylabel('Mg,Ca (mg/L)', color='k')
    ax2.set_ylabel('Na (mg/L)', color='k')

    # ax1.set_ylim(-100, 2300)
    #ax2.set_ylim(17000, 20000)
    # ax1.set_xlim(0,maxNaOH_vol)
    # ax1.set_yticks([i * 200 for i in range(13)])

    ax1.legend(loc='upper center')
    ax2.legend(loc = 'best')
    plt.show()

    # Plot Mg recovery ratio
    # fig, ax1 = plt.subplots()
    # ax1.plot(pHs, Mg_rr, label = "Mg precipitation rate")

    # ax1.set_xlabel('pH')
    # ax1.set_ylabel('Mg precipitation rate (%)', color='k')
    # plt.show()
    print('Max vol: ', maxNaOH_vol, "Conc: ", NaOH_conc)
    return Ca_concs, Na_concs, Mg_concs, pHs, vols, Nas

Ca_concs, Na_concs, Mg_concs, pHs, vols, Nas = plot_step_1_original(maxNaOH_vol=0.5,
                                                                    NaOH_conc = 1,
                                                                    )
#%% Step 1 design
def step1_add_NaOH(
        feed_props = {    
            "Na+"  : 17987,
            "Cl-"  : 32220,
            "SO4-2" : 4766,
            "Mg+2" : 2226,
            "Ca+2" : 714,
            "K+"  : 858,
            "HCO3-": 0,
            "CO3-2": 0,
            },
        pressure = 1.01325, # bar
        add_NaOH_vol = 1, # L
        add_NaOH_conc = 5, # mol/L
        temp = 25, 
    ): 
    # db = PhreeqcDatabase("thermoddem-v1.10.dat")
    db = PhreeqcDatabase("pitzer.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S Si Li Sr B Br"))
    solution.set(ActivityModelPitzer())

    # Create a gaseous phase
    # gaseousphase = GaseousPhase("H2O(g)")
    # gaseousphase.set(ActivityModelPengRobinson())

    # Hydromagnesite = MineralPhase("Hydromagnesite")
    # Magnesite = MineralPhase("Magnesite"
    
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
                            # Magnesite,
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
    # specs.openTo("OH-")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)
    

    state = ChemicalState(system)
    state.temperature(temp, "celsius")

    # Seawater composition (mg/kg)
    Na  = 17987
    Cl  = 32220
    SO4 = 4766
    Mg = 2226
    Ca = 714
    K  = 858
    # HCO3= 300

    state.set("H2O", 1.0, "kg")
    for ion, conc in feed_props.items():
        state.add(ion, conc, "mg")
    # state.pressure(1.0, "atm")
    # state.add("Ca+2"   ,   Ca , "mg")
    # state.add("Mg+2"   ,  Mg , "mg")
    # state.add("Na+"    ,  Na , "mg")
    # state.add("K+"     ,   K , "mg")
    # state.add("Cl-"    , Cl , "mg")
    # state.add("SO4-2"  ,  SO4 , "mg")
    # state.add("HCO3-"  ,  HCO3 , "mg")

    state.add("Li+", 0.6, "mg")
    state.add("Sr+2", 12.6, "mg")
    # state.add("B+3", 5.87, "mg")
    state.add("Br-", 110, "mg")


    # Additional chemicals (5 M NaOH)
    NaOH_props, masses = NaOH(conc = add_NaOH_conc)

    H2O_mass, Na_mass, OH_mass = masses # g/L

    state.add("Na+"  ,  Na_mass * add_NaOH_vol / 1000, "kg")
    state.add("OH-"  ,  OH_mass * add_NaOH_vol / 1000, "kg")
    state.add("H2O"  ,  H2O_mass * add_NaOH_vol / 1000, "kg")

    # equilibrate(state)

    conditions = EquilibriumConditions(specs)
    conditions.temperature(temp, "celsius")
    conditions.pressure(pressure, "bar")

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
                       "Mg+2", "Na+", "OH-", "K+", "HCO3-", "CO3-2"]
    outflow_mass = dict(zip(outflow_species,
                            [float(props.speciesAmount(i)) * molar_mass[i] for i in outflow_species])) # mg

    return output, props, outflow_mass, masses, state

props_experiment_2 = {    
            "Na+"  : 20200,
            "Cl-"  : 32000,
            "SO4-2" : 5400,
            "Mg+2" : 2420,
            "Ca+2" : 770,
            "K+"  : 711,
            "HCO3-": 0,
            "CO3-2": 0,
            }
output1, props1, outflow_mass1, masses1, state1 = step1_add_NaOH(feed_props=props_experiment_2,
                                                            add_NaOH_conc=5,
                                                            add_NaOH_vol=0.0397)
print(output1['pH'])

#%% Step 1 Plot
def plot_step_1(        
        feed_props = {    
            "Na+"  : 17987,
            "Cl-"  : 32220,
            "SO4-2" : 4766,
            "Mg+2" : 2226,
            "Ca+2" : 714,
            "K+"  : 858,
            "HCO3-": 0,
            "CO3-2": 0,
            },
        maxNaOH_vol = 0.6, # L
        runs = 30,
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

    # CO2s = [i * 500 for i in range(31)]

    NaOHs = [i * maxNaOH_vol / runs for i in range(0,runs,1)]

    # for i in CO2s:
    for j in NaOHs:
        output,props,outflow, masses, state = step1_add_NaOH(
            feed_props= feed_props,
            pressure = 1.01325, # bar
            add_NaOH_vol = j, # L
            add_NaOH_conc = 1, # mol/L
            temp = 25, 
                        )

        Calcites.append(float(props.speciesAmount("Calcite")))
        Brucites.append(float(props.speciesAmount("Brucite")))

        temperatures.append(output['temperature'])
        pressures.append(output['pressure'])
        pHs.append(output['pH'])

        # Nas.append(output['Amount'][-12] * 23 / (output['Amount'][1] * 18 / 1000) * 1000)
        # Cas.append(output['Amount'][4] * 40 / (output['Amount'][1] * 18 / 1000) * 1000)
        # Mgs.append(output['Amount'][9] * 24 / (output['Amount'][1] * 18 / 1000) * 1000)
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

        vols.append(float(props.phaseProps("AqueousPhase").volume()))

        Mg_rr.append(props.speciesAmount("Brucite")/float(props.elementAmount("Mg"))*100)

        # print(masses)
    fig, ax1 = plt.subplots()

    Ca_concs = [ Cas[i]*40.087/vols[i] for i in range(len(NaOHs))]
    Na_concs = [ Nas[i]*23.99/vols[i] for i in range(len(NaOHs))]
    Mg_concs = [ Mgs[i]*24.305/vols[i] for i in range(len(NaOHs))]


    ax2 = ax1.twinx()
    ax1.plot(pHs, Mg_concs, 'b-', label='Mg+2')
    ax1.plot(pHs, Ca_concs, 'r-', label='Ca+2')
    ax2.plot(pHs, Na_concs, 'g-', label='Na+')

    # ax2.plot(NaOHs, pHs, 'k-', label='pH')

    ax1.set_xlabel('pH')
    ax1.set_ylabel('Mg,Ca (mg/L)', color='k')
    ax2.set_ylabel('Na (mg/L)', color='k')

    ax1.set_ylim(-100, 2300)
    # ax2.set_ylim(17000, 20000)
    # ax1.set_xlim(0,maxNaOH_vol)
    # ax1.set_yticks([i * 200 for i in range(13)])

    ax1.legend(loc='upper center')
    ax2.legend(loc = 'best')
    plt.show()

    # Plot Mg recovery ratio
    fig, ax1 = plt.subplots()
    ax1.plot(pHs, Mg_rr, label = "Mg precipitation rate")

    ax1.set_xlabel('pH')
    ax1.set_ylabel('Mg precipitation rate (%)', color='k')
    plt.show()

    return props

props = plot_step_1()

#%% Validation paper 2
def original_state_2(
        feed_props = {    
            "Na+"  : 20000,
            "Cl-"  : 35000,
            "SO4-2" : 5000,
            "Mg+2" : 2950,
            "Ca+2" : 900,
            # "Fe+2"  : 1,
            # "NH4+": 0.4,
            # "CO3-2": 0.6,
            "OH-": 0,
            },
        pressure = 1.01325, # bar
        add_NaOH_vol = 0, # L
        add_NaOH_conc = 5, # mol/L
        temp = 25, 
    ): 
    # db = PhreeqcDatabase("thermoddem-v1.10.dat")
    # db = PhreeqcDatabase("pitzer.dat")
    db = PhreeqcDatabase("minteq.v4.dat")
    # db = PhreeqcDatabase("phreeqc.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg S "))
    solution.set(ActivityModelPitzer())
    # solution.set(ActivityModelDebyeHuckel())

    # Create a gaseous phase
    gaseousphase = GaseousPhase("CO2(g)")
    gaseousphase.set(ActivityModelPengRobinson())

    # Hydromagnesite = MineralPhase("Hydromagnesite")
    Magnesite = MineralPhase("Magnesite")
    
    Calcite = MineralPhase("Calcite")
    # Polyhalite = MineralPhase("Polyhalite")
    # Glauberite = MineralPhase("Glauberite")
    Halite = MineralPhase("Halite")
    Gypsum = MineralPhase("Gypsum")
    Anhydrite = MineralPhase("Anhydrite")
    Epsomite = MineralPhase("Epsomite")
    Brucite = MineralPhase("Brucite")
    Portlandite = MineralPhase("Portlandite")
    

    system = ChemicalSystem(db, 
                            gaseousphase,
                            solution,
                            Magnesite,
                            Calcite,
                            # Polyhalite,
                            # Glauberite,
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
    # specs.charge()
    specs.volume()
    # specs.openTo("Cl-")
    specs.openTo("H2O")
    # specs.openTo("OH-")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)
    

    state = ChemicalState(system)
    state.temperature(temp, "celsius")

    # Seawater composition (mg/kg)
    Na  = 17987
    Cl  = 32220
    SO4 = 4766
    Mg = 2226
    Ca = 714
    K  = 858
    # HCO3= 300

    state.set("H2O", 1.0, "kg")
    for ion, conc in feed_props.items():
        state.add(ion, conc, "mg")
    # state.pressure(1.0, "atm")
    # state.add("Ca+2"   ,   Ca , "mg")
    # state.add("Mg+2"   ,  Mg , "mg")
    # state.add("Na+"    ,  Na , "mg")
    # state.add("K+"     ,   K , "mg")
    # state.add("Cl-"    , Cl , "mg")
    # state.add("SO4-2"  ,  SO4 , "mg")
    # state.add("HCO3-"  ,  HCO3 , "mg")
    # equilibrate(state)

    conditions = EquilibriumConditions(specs)
    conditions.temperature(temp, "celsius")
    conditions.pressure(pressure, "bar")

    # conditions.phaseAmount("GaseousPhase", 0.0001, "umol")
    # conditions.setLowerBoundPressure(0.01, "bar")
    # conditions.setUpperBoundPressure(1000.0, "bar")
    conditions.setLowerBoundTemperature(25, "celsius")
    conditions.setUpperBoundTemperature(200, "celsius")
    # conditions.pH(8.2)
    # conditions.charge(0.0)
    conditions.volume(0.001, "m3")

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

    # print('Original feed pH:', float(aprops.pH()))

    return output, state, system, props

ori_output, ori_state, ori_system, ori_props = original_state_2()
print(ori_output['pH'])

#%% Validation case 2 Add NaOH
def validate_add_NaOH_to_original(
        state, system,
        pressure = 1.01325, # bar
        add_NaOH_vol = 0.0365, # L
        add_NaOH_conc = 5, # mol/L
        add_CO2 = 0,
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
    conditions.fugacity("CO2(g)", 0.00076 * 1.01325, "bar")

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

    return output, props, outflow_mass, out_brucite, state, system


val_output, val_state, val_system, val_props = original_state_2()
output_val, props_val, outflow_mass_val, out_brucite_val, state_val, system_val = validate_add_NaOH_to_original(state= val_state, system=val_system,
                                                                             add_NaOH_conc=1,
                                                                             add_NaOH_vol=0,
                                                                             )
print('Step One simulation')
print('pH=', output_val['pH'])
print('Ca in CaCO3 (mg)', float(props_val.speciesAmount('Calcite')) * molar_mass['Ca+2']* 1000)
print('Ca+2 (mg)', float(props_val.speciesAmount('Ca+2')) * molar_mass['Ca+2']* 1000)
print('Mg in MgCO3 (mg)', (float(props_val.speciesAmount('Magnesite')))* molar_mass['Mg+2']* 1000 )
print('Mg in MgOH2 (mg)', (float(props_val.speciesAmount('Brucite')))* molar_mass['Mg+2'] * 1000 )
print('Mg+2 (mg)', (float(props_val.speciesAmount('Mg+2')))* molar_mass['Mg+2'] * 1000 )
print('Ca rec rate', float(props_val.speciesAmount('Calcite')) * molar_mass['Ca+2'] / 0.714)

#%% Validation plot
def plot_val_original(       
        maxNaOH_vol = 0.5, # L
        NaOH_conc = 1, # mol/L
        runs = 200,
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

    CO2s = [i * 250 for i in range(100)]

    NaOHs =  [0.005 * i for i in range(40)]+  [i * maxNaOH_vol / runs for i in range(1,runs+1,1)]

    # for i in CO2s:
    for v in NaOHs:
        ori_output, ori_state, ori_system, ori_props = original_state_2()

        output, props, outflow_mass, out_brucite, state, system = validate_add_NaOH_to_original(
            state=ori_state, 
            system=ori_system,
            pressure = 1.01325, # bar
            add_NaOH_vol = v, # L
            add_NaOH_conc = NaOH_conc, # mol/L
            add_CO2= 0, #mg
            temp = 25, 
                        )

        Calcites.append(float(props.speciesAmount("Calcite")))
        Brucites.append(float(props.speciesAmount("Brucite")))
        Magnesites.append(float(props.speciesAmount("Magnesite")))

        temperatures.append(output['temperature'])
        pressures.append(output['pressure'])
        pHs.append(output['pH'])

        # Nas.append(output['Amount'][-12] * 23 / (output['Amount'][1] * 18 / 1000) * 1000)
        # Cas.append(output['Amount'][4] * 40 / (output['Amount'][1] * 18 / 1000) * 1000)
        # Mgs.append(output['Amount'][9] * 24 / (output['Amount'][1] * 18 / 1000) * 1000)
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

        vols.append(float(props.phaseProps("AqueousPhase").volume()))

        Mg_rr.append(props.speciesAmount("Brucite")/float(props.elementAmount("Mg"))*100)

        # print(masses)
    fig, ax1 = plt.subplots()

    Ca_concs = [ Cas[i]*40.087/vols[i] for i in range(len(NaOHs))]
    Na_concs = [ Nas[i]*22.99/vols[i] for i in range(len(NaOHs))]
    Mg_concs = [ Mgs[i]*24.305/vols[i] for i in range(len(NaOHs))]

    calcites_conc = [ Calcites[i] / vols[i] for i in range(len(vols))]
    brucites_conc = [ Brucites[i] / vols[i] for i in range(len(vols))]
    magnesites_conc = [ Magnesites[i] / vols[i] for i in range(len(vols))]


    # ax2 = ax1.twinx()
    ax1.plot(pHs, calcites_conc, 'r-', label='CaCO3')
    ax1.plot(pHs, brucites_conc, 'g-', label='MgOH2')
    ax1.plot(pHs, magnesites_conc, 'b-', label='MgCO3')

    # ax2.plot(NaOHs, pHs, 'k-', label='pH')

    ax1.set_xlabel('pH')
    ax1.set_ylabel('Concentration (mmol/L)', color='k')
    # ax2.set_ylabel('Na (mg/L)', color='k')

    # ax1.set_ylim(-0, 32)
    #ax2.set_ylim(17000, 20000)
    # ax1.set_xlim(0,maxNaOH_vol)
    # ax1.set_yticks([i * 200 for i in range(13)])

    ax1.legend(loc='upper center')
    # ax2.legend(loc = 'best')
    plt.show()

    # Plot Mg recovery ratio
    # fig, ax1 = plt.subplots()
    # ax1.plot(pHs, Mg_rr, label = "Mg precipitation rate")

    # ax1.set_xlabel('pH')
    # ax1.set_ylabel('Mg precipitation rate (%)', color='k')
    # plt.show()
    print('Max vol: ', maxNaOH_vol, "Conc: ", NaOH_conc)
    return calcites_conc, brucites_conc, magnesites_conc, pHs, vols

calcites_conc, brucites_conc, magnesites_conc, pHs, vols = plot_val_original(maxNaOH_vol=40,
                                                                NaOH_conc = 0.01,
                                                                runs = 200,
                                                                )
#%% Step 2 design
def step2_add_NaOH_CO2(
        # state, system,
        feed_props = {    
            "Na+"  : 22.18313226220165,
            "Cl-"  : 32.756795904174304,
            "SO4-2" : 4.765747117899852,
            "Mg+2" : 0.015354423011699113,
            "Ca+2" : 0.7139839160295243,
            "K+"  : 0.8580120385839307,
            "OH-" : 0.009934676623053135,
            "H2O" : 1037.3791483025027,
            # "MgOH+": 0.0002884,
            "HCO3-": 0,
            "CO3-2": 0,
            },
        pressure = 1.01325, # bar
        add_NaOH_vol = 1, # L
        add_NaOH_conc = 5, # mol/L
        add_CO2 = 0, # mg
        temp = 25, 
    ): 
    db = PhreeqcDatabase("pitzer.dat")
    # db = PhreeqcDatabase("minteq.v4.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S"))
    solution.set(ActivityModelPitzer())

    # Create a gaseous phase
    gaseousphase = GaseousPhase("CO2(g)")
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
    Brucite = MineralPhase("Brucite")
    Portlandite = MineralPhase("Portlandite")
    

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
    state.temperature(temp, "celsius")

    state.set("H2O", feed_props["H2O"] / 1000, "kg")
    for ion, conc in feed_props.items():
        if ion == "H2O":
            continue
        state.add(ion, conc * 1000, "mg")


    # Additional chemicals (5 M NaOH)
    NaOH_props, masses = NaOH(conc = add_NaOH_conc)

    H2O_mass, Na_mass, OH_mass = masses # g/L

    state.add("Na+"  ,  Na_mass * add_NaOH_vol / 1000, "kg")
    state.add("OH-"  ,  OH_mass * add_NaOH_vol / 1000, "kg")
    state.add("H2O"  ,  H2O_mass * add_NaOH_vol / 1000, "kg")
    state.add("CO2"  ,  add_CO2, "mg")

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
    
    outflow_species = ["H2O", "Ca+2", "Cl-", "SO4-2",
                       "Mg+2", "Na+", "OH-", "K+", "HCO3-", "CO3-2"]
    outflow_mass = dict(zip(outflow_species,
                            [float(props.speciesAmount(i)) * molar_mass[i] for i in outflow_species]))

    return output, props, outflow_mass

output2, props2, outflow_mass2 = step2_add_NaOH_CO2(feed_props=outflow_mass_1,
                                                    add_NaOH_vol=0.012,
                                                    add_NaOH_conc=5,
                                                    add_CO2=1000)


# Test on repeatedly adding NaOH and CO2
# ori_output, ori_state, ori_system, ori_props = original_state()
# outflow_species = ["H2O", "Ca+2", "Cl-", "SO4-2",
#                     "Mg+2", "Na+", "OH-", "K+", "HCO3-", "CO3-2"]
# ori_outflow  = dict(zip(outflow_species,
#                         [float(props.speciesAmount(i)) * molar_mass[i] for i in outflow_species])) # mg

# output2, props2, outflow_mass2 = step2_add_NaOH_CO2(feed_props=ori_outflow,
#                                                     add_NaOH_conc=1,
#                                                     add_NaOH_vol=0, add_CO2=0)
print('Step 2 simulation')
print('pH=', output2['pH'])
print('Ca precipitated (g)', float(props2.speciesAmount('Calcite')) * molar_mass['Ca+2'])
print('Mg precipitated (g)', (float(props2.speciesAmount('Magnesite')) + float(props2.speciesAmount('Brucite')))* molar_mass['Mg+2'])
print('Ca rec rate', float(props2.speciesAmount('Calcite')) * molar_mass['Ca+2'] / 0.714)

#%% Plot
def plot_step_2(        
        feed_props = {    
            "Na+"  : 22.18313226220165,
            "Cl-"  : 32.756795904174304,
            "SO4-2" : 4.765747117899852,
            "Mg+2" : 0.015354423011699113,
            "Ca+2" : 0.7139839160295243,
            "K+"  : 0.8580120385839307,
            "OH-" : 0.009934676623053135,
            "H2O" : 1037.3791483025027,
            # "MgOH+": 0.0002884,
            "HCO3-": 0,
            "CO3-2": 0,
            },):

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

    CO2s = [i * 200 for i in range(26)]
    NaOHs = [i * 0.002 for i in range(26)]

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
            output,props, outflow = step2_add_NaOH_CO2(
                feed_props = feed_props,
                pressure = 1.01325, # bar
                add_NaOH_vol = i, # L
                add_NaOH_conc = 5, # mol/L
                add_CO2 = j, # mg
                temp = 25, 
                            )
            
            brucites.append(float(props.speciesAmount("Brucite")))
            calcites.append(float(props.speciesAmount("Calcite")))
            magnesites.append(float(props.speciesAmount("Magnesite")))

            temps.append(output['temperature'])
            pres.append(output['pressure'])
            phs.append(output['pH'])

            cas.append(float(props.speciesAmount("Ca+2")))
            nas.append(float(props.speciesAmount("Na+")))
            mgs.append(float(props.speciesAmount("Mg+2")))

            ca_rr.append(float(props.speciesAmount("Calcite"))/float(props.elementAmount("Ca"))*100)

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

    print('')
    print('Calcites (mol)')
    for i in Calcites:
        print([round(j,3) for j in i])

    print('')
    print('Brucites (mol)')
    for i in Brucites:
        print([round(j,3) for j in i])

    print('')
    print('Magnesites (mol)')
    for i in Magnesites:
        print([round(j,3) for j in i])


    print('')
    print('pH')
    for i in pHs:
        print([round(j,3) for j in i])
# Plot pH
    fig,  axes= plt.subplots(1,1)
    ax1 = axes

    im1 = ax1.imshow([pHs[i-1] for i in range(len(pHs),0,-1)], cmap='YlGnBu', interpolation='none')
    # ax1.set_aspect(1)
    ax1.set(xticks=np.arange(0, len(CO2s)+1, 5), xticklabels=np.arange(0, np.max(CO2s)+1, np.max(CO2s)/5))
    ax1.set(yticks=np.arange(0, len(NaOHs)+1, 5), yticklabels=[0.05,0.04,0.03,0.02,0.01,0])
    ax1.set_xlabel('CO2 addition (mg)')
    ax1.set_ylabel('5M NaOH addition (L)')
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(im1,  ticks = [4,5,6,7,8,9,10,11,12,13])
    cbar.set_label('pH')
    
# Plot Ca_rr
    fig,  axes= plt.subplots(1,1)
    ax1 = axes

    im1 = ax1.imshow([Ca_rr[i-1] for i in range(len(pHs),0,-1)], cmap=plt.cm.Reds, interpolation='none')
    # ax1.set_aspect(1)
    ax1.set(xticks=np.arange(0, len(CO2s)+1, 5), xticklabels=np.arange(0, np.max(CO2s)+1, np.max(CO2s)/5))
    ax1.set(yticks=np.arange(0, len(NaOHs)+1, 5), yticklabels=[0.05,0.04,0.03,0.02,0.01,0])
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
    ax1.set(xticks=np.arange(0, len(CO2s)+1, 5), xticklabels=np.arange(0, np.max(CO2s)+1, np.max(CO2s)/5))
    ax1.set(yticks=np.arange(0, len(NaOHs)+1, 5), yticklabels=[0.05,0.04,0.03,0.02,0.01,0])
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
    ax1.set(xticks=np.arange(0, len(CO2s)+1, 5), xticklabels=np.arange(0, np.max(CO2s)+1, np.max(CO2s)/5))
    ax1.set(yticks=np.arange(0, len(NaOHs)+1, 5), yticklabels=[0.05,0.04,0.03,0.02,0.01,0])
    ax1.set_xlabel('CO2 addition (mg)')
    ax1.set_ylabel('5M NaOH addition (L)')
    # cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cbar = fig.colorbar(im1,  ticks = [0, 0.0001, 0.0002, 0.0003, 0.0004, 0.0005, 0.0006, 0.0007 ])
    cbar.set_label('Mg(OH)2 (mol)')

plot_step_2()

#%% Step 3 design
def step3_add_CO2(
        feed_props = {    
            "Na+"  : 23.56304888864093,
            "Cl-"  : 32.757606881477145,
            "SO4-2" : 4.765494250623454,
            "Mg+2" : 2.736647107392137e-5,
            "Ca+2" : 0.0020614526693185464,
            "K+"  : 0.8580240773367747,
            "OH-" : 0.23640685634626776,
            "H2O" : 1050.048811961014,
            # "MgOH+": 0.0002884,
            "HCO3-": 0.0013107795201177737,
            "CO3-2": 0.29625660487721833,
            },
        brucite = 5304.26, # mg Mg(OH)2
        pressure = 1.01325, # bar
        add_NaOH_vol = 1, # L
        add_NaOH_conc = 5, # mol/L
        add_CO2 = 0, # mg
        CO2_fug = 0, # bar
        temp = 25, 
    ): 
    # db = PhreeqcDatabase("pitzer.dat")
    db = PhreeqcDatabase("minteq.v4.dat")
    # db = PhreeqcDatabase("thermoddem-v1.10.dat")

    solution = AqueousPhase(speciate("H O C Na Cl Ca Mg K S"))
    solution.set(ActivityModelPitzer())

    # Create a gaseous phase
    gaseousphase = GaseousPhase("CO2(g)")
    gaseousphase.set(ActivityModelPengRobinson())

    # Hydromagnesite = MineralPhase("Hydromagnesite")
    Magnesite = MineralPhase("Magnesite")
    Calcite = MineralPhase("Calcite")
    Halite = MineralPhase("Halite")
    Gypsum = MineralPhase("Gypsum")
    Anhydrite = MineralPhase("Anhydrite")
    Epsomite = MineralPhase("Epsomite")
    Brucite = MineralPhase("Brucite")
    # Hydromagnesite =  MineralPhase("Hydromagnesite")
    # Periclase = MineralPhase("Periclase")
    

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
                            # Hydromagnesite,
                            # Periclase
                            )

    specs = EquilibriumSpecs(system)
    specs.temperature()
    specs.pressure()
    # specs.phaseAmount("GaseousPhase")
    # specs.pH()
    specs.charge()
    specs.openTo("Cl-")
    # specs.fugacity("CO2(g)")

    # solver = EquilibriumSolver(specs)
    solver = SmartEquilibriumSolver(specs)
    

    state = ChemicalState(system)
    state.temperature(temp, "celsius")

    add_H2CO3 = add_CO2 / 44.009 * 62.03 # mg
    remove_H2O = add_CO2 /44.009 * 18.021 # mg

    state.set("H2O", feed_props["H2O"] / 1000 - remove_H2O / 1e6, "kg")
    for ion, conc in feed_props.items():
        if ion == "H2O":
            continue
        state.add(ion, conc * 1000, "mg")
    
    state.add("H2CO3", add_H2CO3, "mg")


    # Add back Mg(OH)2

    state.add("Brucite", brucite, "mg")
    # state.add("CO2"  ,  add_CO2, "mg")

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
    # conditions.fugacity("CO2(g)", CO2_fug, "bar")

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
                       "Mg+2", "Na+", "OH-", "K+", "HCO3-", "CO3-2"]
    outflow_mass = dict(zip(outflow_species,
                            [float(props.speciesAmount(i)) * molar_mass[i] for i in outflow_species]))

    return output, props, outflow_mass

# validation

output3, props3, outflow_mass3 = step3_add_CO2(feed_props=outflow_mass2,
                                               add_CO2=8000,
                                               brucite = 5304.26
                                            #    CO2_fug=CO2_fug,
                                               )
# CO2=4066
print('pH: ', output3['pH'])
print('Mg in MgCO3:', float(props3.speciesAmount("Magnesite")) * molar_mass["Mg+2"])
print('MgCO3 %:', float(props3.speciesAmount("Magnesite")) / (float(props3.elementAmount("Mg"))))
print('Mg in MgOH2:', float(props3.speciesAmount("Brucite")) * molar_mass["Mg+2"])
print('Mg+2:', (float(props3.speciesAmount("Mg+2"))) * molar_mass["Mg+2"])
print('Ca in CaCO3:', float(props3.speciesAmount("Calcite")) * molar_mass["Ca+2"])
print('Ca+2: ', float(props3.speciesAmount("Ca+2")) * molar_mass["Ca+2"])

#%% Plot
def plot_step_3(
        feed_props = {    
            "Na+"  : 23.56304888864093,
            "Cl-"  : 32.757606881477145,
            "SO4-2" : 4.765494250623454,
            "Mg+2" : 2.736647107392137e-5,
            "Ca+2" : 0.0020614526693185464,
            "K+"  : 0.8580240773367747,
            "OH-" : 0.23640685634626776,
            "H2O" : 1050.048811961014,
            # "MgOH+": 0.0002884,
            "HCO3-": 0.0013107795201177737,
            "CO3-2": 0.29625660487721833,
            },
        brucite = 5304.26,
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

    Mg_rr = []

    CO2s = [i * 100 for i in range(81)]
    # NaOHs = [i * 0.001 for i in range(1,50,1)]

    # for i in CO2s:
    for j in CO2s:
        output,props,outflow = step3_add_CO2(
            feed_props= outflow_mass2,
            brucite= brucite,
            pressure = 1.01325, # bar
            add_CO2= j, # mg
            temp = 25, 
                        )
        # rrs.append(i)

        

        Brucites.append(float(props.speciesAmount("Brucite")))
        # Hydromagnesites.append(float(props.speciesAmount("Hydromagnesite")))
        Magnesites.append(float(props.speciesAmount("Magnesite")))
        Calcites.append(float(props.speciesAmount("Calcite")))
        # Periclases.append(float(props.speciesAmount("Periclase")))

        temperatures.append(output['temperature'])
        pressures.append(output['pressure'])
        pHs.append(output['pH'])

        vols.append(float(props.phaseProps("AqueousPhase").volume()))

        Cas.append(float(props.speciesAmount("Ca+2")))
        Nas.append(float(props.speciesAmount("Na+")))
        Mgs.append(float(props.speciesAmount("Mg+2")))

        Mg_rr.append(props.speciesAmount("Brucite")/float(props.elementAmount("Mg"))*100)

    fig, ax1 = plt.subplots()

    ax2 = ax1.twinx()
    ax1.plot(CO2s, Magnesites, 'b-', label='MgCO3')
    ax1.plot(CO2s, Calcites, 'r-', label='CaCO3')
    ax1.plot(CO2s, Brucites, 'g-', label='Mg(OH)2')
    # ax1.plot(CO2s, Mgs, 'y-', label="Mg+2")

    ax2.plot(CO2s, pHs, 'y-', label='pH')

    ax1.set_xlabel('Added CO2 (mg)')
    ax1.set_ylabel('Solids (mol)', color='k')
    ax2.set_ylabel('pH', color='k')

    # ax1.set_ylim(-100, 2200)
    # ax2.set_ylim(17000, 20000)
    # ax1.set_xlim(0,0.05)
    # ax1.set_yticks([i * 200 for i in range(13)])

    ax1.legend(loc='upper center')
    ax2.legend(loc = 'best')
    plt.show()

    # Plot Mg recovery ratio
    fig, ax1 = plt.subplots()
    ax1.plot(pHs, Mg_rr, label = "Mg precipitation rate")

    ax1.set_xlabel('pH')
    ax1.set_ylabel('Mg precipitation rate (%)', color='k')
    plt.show()

    return Mgs, Magnesites, Brucites, Calcites, pHs, vols

feed = {'H2O': 1050.049,
 'Ca+2': 0.00206,
 'Cl-': 32.756,
 'SO4-2': 4.766,
 'Mg+2': 0.0154,
 'Na+': 23.5624,
 'OH-': 0.2364,
 'K+': 0.858,
 'HCO3-': 0.00131,
 'CO3-2': 0.2963}

Mgs_3, Magnesites_3, Brucites_3, Calcites_3, pHs_3, vols_3 = plot_step_3(feed_props=feed)
# %% Applied to NF brine
plot_step_1(
    feed_props = {    
    "Na+"  : 5384,
    "Cl-"  : 16586.5,
    "SO4-2" : 11797,
    "Mg+2" : 6063,
    "Ca+2" : 1628.5,
    "K+"  : 199.55,
    "HCO3-": 531.5,
    "CO3-2": 0,
    },
    maxNaOH_vol=0.15,
    runs = 50,
)
# %%
