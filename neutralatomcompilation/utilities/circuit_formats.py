import qiskit
import networkx as nx
from collections import defaultdict
from typing import List, Tuple, DefaultDict, Dict, Set

def create_circuit_digraph(circuit: qiskit.circuit.QuantumCircuit) -> Tuple[nx.DiGraph,
                            Dict[int, qiskit.circuit.Instruction]]:

    prev_instructions = defaultdict(lambda: -1)
    id_to_instruction = {}

    graph = nx.DiGraph()
    for index, instruction in enumerate(circuit):
        id_to_instruction[index] = instruction

        gate = instruction[0]
        qubits = instruction[1]
        cbits = instruction[2]

        graph.add_node(index)
        for q in qubits:
            prev = prev_instructions[q]
            prev_instructions[q] = index
            graph.add_edge(prev, index)

    return graph, id_to_instruction