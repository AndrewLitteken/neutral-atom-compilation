import qiskit
import networkx as nx
import time

class LabeledOp:
    def __init__(self, instruction, index):
        self.instruction = instruction
        self.op = instruction[0]
        self.qubits = instruction[1]
        self.index = index
        
    def __eq__(self, other):
        return self.index == other.index and self.instruction == other.instruction
    
    def __hash__(self):
        return self.index.__hash__()
    
    def __repr__(self):
        return f'{self.op} ({self.index})'

def circuit_to_dag(c : qiskit.circuit.QuantumCircuit) -> nx.DiGraph:
    
    dag = nx.DiGraph()
    
    qubit_last_use = {}
    
    null_op = LabeledOp([None, []], -1)
    
    for index, instruction in enumerate(c):
        op = LabeledOp(instruction, index)
        
        for qubit in op.qubits:
            if qubit not in qubit_last_use:
                dag.add_edge(null_op, op)
            else:
                dag.add_edge(qubit_last_use[qubit], op)
                
            qubit_last_use[qubit] = op
            
    final_op = LabeledOp([None, []], float('inf'))
    for q in qubit_last_use:
        dag.add_edge(qubit_last_use[qubit], final_op)
            
    return nx.transitive_reduction(dag)

def decompose_swap(circuit, hw, im, record_largest_gate, mapping = None):

    def _update_dag_(dag, op_to_remove):
        for target in dag[op_to_remove]:
            dag.add_edge(start_node, target)
        dag.remove_node(op_to_remove)

    new_c = circuit.copy()
    new_c.data = []
    #new_c = qiskit.circuit.QuantumCircuit(*circuit.qregs)
    
    swap = qiskit.circuit.library.standard_gates.swap.SwapGate
    cx = qiskit.circuit.library.standard_gates.x.CXGate
    
    # 1. Decompose swaps
    for gate in circuit:
        if isinstance(gate[0], qiskit.circuit.library.standard_gates.swap.SwapGate):
            new_c.append(cx(), gate[1])
            new_c.append(cx(), [gate[1][1], gate[1][0]])
            new_c.append(cx(), gate[1])
        elif isinstance(gate[0], qiskit.circuit.barrier.Barrier):
            pass
        else:
            new_c.append(gate[0], gate[1])
            
    # 2. Schedule
    dag = circuit_to_dag(new_c)
    
    start_node = LabeledOp([None, []], -1)
    final_op = LabeledOp([None, []], float('inf'))
    
    frontier = []
    for op in dag[start_node]:
        if op != final_op:
            frontier.append(op)
            
    final_circuit = new_c.copy()
    final_circuit.data = []
    
    # Want largest non overlapping set, but that's hard :(
    largest_per_time_step = []
    while len(frontier) > 0:
        to_execute = []
        largest = float('-inf')
        
        for op in frontier:
            in_edges = dag.in_edges(op)
            runnable = True
            for edge in in_edges:
                if edge[0] != start_node:
                    runnable = False
                    break
            if not runnable:
                continue
        
            if mapping is None:
                if im.valid_interaction_constrained(
                    [hw.qiskit_to_hardware[qq] for qq in op.qubits], 
                    [[hw.qiskit_to_hardware[qq] for qq in o.qubits] for o in to_execute]):
                    to_execute.append(op)
            else:
                if im.valid_interaction_constrained(
                    [mapping[hw.qiskit_to_hardware[qq]] for qq in op.qubits], 
                    [[mapping[hw.qiskit_to_hardware[qq]] for qq in o.qubits] for o in to_execute]):
                    to_execute.append(op)
            
        for op in to_execute:
            if len(op.qubits) > largest:
                largest = len(op.qubits)
            final_circuit.append(op.op, op.qubits)
            
            _update_dag_(dag, op)
        
            # Update the frontier
        largest_per_time_step.append(largest)

        frontier = []
        for op in dag[start_node]:
            if op != final_op:
                frontier.append(op)
        
        qubits = list(final_circuit.qubits)
        final_circuit.append(qiskit.circuit.barrier.Barrier(1), [qubits[0]], [])
    
    if record_largest_gate:
        return final_circuit, largest_per_time_step        
    return final_circuit
            
            
            