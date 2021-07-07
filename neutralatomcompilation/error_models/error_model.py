import qiskit
import numpy as np

from ..utilities.circuit_formats import create_circuit_digraph
from ..utilities.circuit_stats import digraph_to_time_steps

class ErrorModel:

  def __init__(self,
               one_qubit_success: float,
               two_qubit_success: float,
               three_qubit_success: float,
               one_qubit_time: float,
               two_qubit_time: float,
               three_qubit_time: float,
               t1_time: float,
               transition_time: float):
    self.one_qubit_success = one_qubit_success
    self.two_qubit_success = two_qubit_success
    self.three_qubit_success = three_qubit_success
    self.one_qubit_time = one_qubit_time
    self.two_qubit_time = two_qubit_time
    self.three_qubit_time = three_qubit_time
    self.success_dict = {1: self.one_qubit_success, 2: self.two_qubit_success, 3: self.three_qubit_success}
    self.time_dict = {1: self.one_qubit_time, 2: self.two_qubit_time, 3: self.three_qubit_time}
    self.t1_time = t1_time
    self.transition_time = transition_time

  def calculate_gate_error(self, circuit: qiskit.circuit.QuantumCircuit) -> float:
    cumulative_error = 1
    for instruction in circuit:
      if isinstance(instruction[0], qiskit.circuit.library.Barrier):
        continue
      num_qubits = len(instruction[1])
      cumulative_error *= self.success_dict[num_qubits]
    return cumulative_error

  def calculate_coherence_error(self, circuit: qiskit.circuit.QuantumCircuit):
    digraph, id_to_instruction = create_circuit_digraph(circuit)

    time_steps = digraph_to_time_steps(digraph, id_to_instruction)
    overall_time = 0
    for step in time_steps:
      max_size = 0
      for inst in step:
        if isinstance(inst[0], qiskit.circuit.library.Barrier):
          continue
        num_qubits = len(inst[1])
        if num_qubits > max_size:
          max_size = num_qubits
      if max_size != 0:
        overall_time += self.time_dict[max_size] + self.transition_time

    error = np.exp(-overall_time / self.t1_time)

    return error


