from reaktoro import * 

db = PhreeqcDatabase("pitzer.dat")

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

# Decalre the equilibrium state properties
specs = EquilibriumSpecs(system)
specs.pressure()
specs.phaseAmount("GaseousPhase")
# specs.pH()
specs.charge()
specs.openTo("Cl-")

solver = EquilibriumSolver(specs)


state = ChemicalState(system)
state.temperature(90, "celsius") # temperature here doesn't matter
state.set("H2O", 1.0, "kg")

state.add("Ca+2"   ,   21100, "mg")
state.add("Mg+2"   ,   65940, "mg")
state.add("Na+"    ,  130680, "mg")
state.add("K+"     ,   20240, "mg")
state.add("Cl-"    ,  351138, "mg")
# state.add("HCO3-"  ,   7490, "mg")
state.add("SO4-2"  ,  137250, "mg")


# Speicify the equilibrum condition
conditions = EquilibriumConditions(specs)

conditions.pressure(0.99, "bar")
# Small amount of gas phase to indicate it starts bubbling
conditions.phaseAmount("GaseousPhase", 0.01, "umol") 
conditions.setLowerBoundTemperature(25, "celsius")
conditions.setUpperBoundTemperature(200, "celsius")
conditions.charge(0.0)

result = solver.solve(state, conditions)

# Check calculation succeeded
assert result.succeeded()
 
props = ChemicalProps(state)
aprops = AqueousProps(system)
aprops.update(state)

# print(props)
print(state)
print('density = ', props.phaseProps("AqueousPhase").density())
print('pH = ', float(aprops.pH()))