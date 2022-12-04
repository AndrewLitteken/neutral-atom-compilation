# Neutral Atom Compilation

## Introduction
This is the compiler from the paper _Exploiting Long-Distance Intractions and Tolerating Atom Loss in Neutral Atom Architectures_ from ISCA 2021: [arxiv]().

It is also the compiler from the paper _Exploiting Long-Distance Intractions and Tolerating Atom Loss in Neutral Atom Architectures_ from QCE 2022: [arxiv]().

## Installation
This can be installed using `pip install .` from the main directory.

## Examples

### Hardware Model
The Neutral Atom Architecture is modelled as an 2 or 3 dimensional grid of equally spaced qubits in each dimension.  It is instantiated as such:
```
nac.Hardware(num_dimensions=2, dimensions_length=(10, 10), dimensions_spacing=(1,1))
```

This data structure contains several methods to check physical and manhattan distances between qubits, as well as whether two qubits are within specific radius.

### Interaction Model
We model the possible long range interactions using a `nx.Graph` in the `InteractionGraph` data structure.  We specify the underlying hardware, the function of the physical distance to restriction zone radius, and the maximum interaction distance. It is instantiated as:
```
nac.InteractionModel(hardware=hw, d_to_r=lambda x: x / 2, max_int_dist=3)
```
where our restriction zone radius is half of the physical distance and a maximum interaction distance of 3.

### Compiler
The compiler used for neutral atom architectures is instatiated by specifying the `Hardware` and `InteractionModel`:
```
comp = nac.LookaheadCompiler(interaction_model=im, hardware=hw)
```
We compile a circuit using this compiler, including swaps as such:
```
c = comp.compile(circuit, lookahead_distance=float('inf'), weighting_function=lambda x: np.e ** (-x))
```

### Atom Loss Handlers
For Atom Loss handlers, we provide an interaction model and hardware model:
```
hh = NAC.HoleHandler(hw, im)
```
To adjust a mapped and routed circuit to lost atoms, we provide the circuit, the lost list of atoms in terms of physical locations, a strategy for shifting the mapped atoms, and rerouting the shifted atoms, if at all.
```
hh.reroute_with_holes(compiled_circuit, lost_list,
                      shift_strategy=NAC.ShiftStrategy.Strategy,
                      route_strategy=NAC.ReRouteStrategy.Strategy)
```

To reset the atom loss use:
```
hh.reset()
hw.reset()
im.reset_graph()
```

The shift strategies are:
- `NAC.ShiftStrategy.NaiveMinMovement`: Move atoms in the direction where there are the most unused atoms along the Hardware Graph.
- `NAC.ShiftStrategy.InteractionGraph`: Move atoms along the shortest Interaction Graph path to a spare qubit.

The routing strategies are:
- `NAC.Reroute.Fail`: Do not attempt to reroute the circuit when qubits are out of range.
- `NAC.RerouteShiftStrategy.Swap`: Attempt to reroute the circuit when qubits are out of range.

### Atom Loss Handlers: Relocation
To use relocation, we need to use either the `LookaheadCompilerConstrained` or `ParallelLookaheadCompilerConstrained` compilers.  These make sure that the circuit is constrianed to a specific section of the circuit.  The `LookaheadCompilerConstrained` takes one extra argument, "wide" or "tight". "wide" uses a bounding box that is the square root of the cirucit by the square root of the size of the circuit.  "tight" is the square root of the size of the circuit by the number of qubits in the circuit divided by the square root.  The `ParallelLookaheadCompilerConstrained` has one additional argument, "percentage" which is a decimal value from 0 to 1, and denotes a maximum of how much the architecture is allowed to be used. If the size of the circuit exceeds the percentage usage, only one instance of the cirucit is run.

The constrained compiler has an extra attribute `.viable_node_sets`, which is a list of lists of equal numbers of qubits.  To perform a relocation of the circuit on the architecture, pick one sections from this list and use:
```
new_set = compiler.viable_node_sets[i]
hh.readjust_starting_loc(circuit, compiler.viable_node_sets[0], new_set, shift_strategy=ShiftStrategy, route_strategy=RerouteStrategy)
```
to get the new circuit.

The parallel compiles has an extra attribute `.included_constrained_sets`, which is a list of lists of equal numbers of qubits.  To perform a relocation of the circuit on the architecture for multiple parallel instances, pick one sections from this list and use:
```
new_set = compiler.included_constrained_sets[i][0]
hh.readjust_starting_loc(circuit, compiler.included_constrained_sets[0][0], new_set, shift_strategy=ShiftStrategy, route_strategy=RerouteStrategy)
```
to get the new circuit.

### Swap Gates Decomposition
To insert barriers into a circuit such that follows the specified restrictions for the architecture we do:
`new_circuit = NAC.utilities.decompose_swaps.decompose_swap(circuit, hardware, interaction_model, mapping=fixed_mapping)`
where the mapping is the remapping from physical qubits to the new physical qubit locations if needed.