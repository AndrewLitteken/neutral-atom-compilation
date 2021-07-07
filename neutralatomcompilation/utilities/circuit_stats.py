import qiskit
import networkx as nx
from collections import defaultdict
from typing import List, Tuple, DefaultDict, Dict, Set

from .circuit_formats import create_circuit_digraph
from ..hardware.hardware import Hardware
from ..compiler.lookahead_compiler_new import circuit_to_dag, LabeledOp

def swap_num(circuit: qiskit.QuantumCircuit) -> int:
  '''
    Get the number of swaps in a given quantum circuit
  '''
  num_swaps = 0
  for instruction in circuit:
    if isinstance(instruction[0], qiskit.circuit.library.standard_gates.swap.SwapGate):
      num_swaps += 1
  return num_swaps

def swap_lengths(circuit: qiskit.QuantumCircuit,
                 hardware: Hardware,
                 physical_to_hardware):
  '''
    Get the number of each length of swaps
  '''
  mapping = {}
  for p in physical_to_hardware:
    mapping[physical_to_hardware[p]] = p

  swap_counts = defaultdict(int)

  for instruction in circuit:
    if isinstance(instruction[0], qiskit.circuit.library.standard_gates.swap.SwapGate):
      qubits = instruction[1]
      distance = hardware.integer_distance(mapping[qubits[0].register], mapping[qubits[1].register])
      swap_counts[distance] += 1

  return swap_counts

def digraph_to_time_steps(digraph: nx.DiGraph,
                          id_to_instruction: Dict[int, qiskit.circuit.Instruction]) \
                           -> List[List[qiskit.circuit.Instruction]]:
  time_steps = []

  visited = set()
  if id_to_instruction:
    frontier = [-1]
  else:
    frontier = [LabeledOp([None, []], -1)]

  while len(frontier) > 0:

    new_frontier = []
    for f in frontier:
      visited.add(f)
      out_edges = digraph.out_edges(f)
      for o in out_edges:
        in_edges = digraph.in_edges(o[1])
        cleared = True
        for i in in_edges:
          if i[0] not in visited:
            cleared = False;
        if cleared:
          new_frontier.append(o[1])
    if id_to_instruction != None:
      time_steps.append([id_to_instruction[f] for f in new_frontier])
    else:
      time_steps.append([f for f in new_frontier])
    frontier = new_frontier

  return time_steps


