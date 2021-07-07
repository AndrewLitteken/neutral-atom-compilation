import qiskit
from ..compiler.lookahead_compiler_new import circuit_to_dag
from ..utilities.circuit_stats import digraph_to_time_steps
import networkx as nx
'''
    Metrics to Compute:
        1. Number of gates
            a. # Logic
            b. # from SWAPS / Communication
        2. Circuit depth
        3. Parallelism
            a. Given as either a raw score (1 <= score <= num_qubits) computed as
            #ops / len(cp) or ratio of score / max_score. The upper limit is saturated 
            when every qubit has gates equal to the length of the critical path
        4. Parallelism Preservation
            How well does a circuit post compilation preserve the parallelsim of an input circuit
            Computed 2 ways: 
                a. Compute the above paralleism score, but consider the CP as one _including_ barriers
                b. Check to see if ops which could be done in parallel in the original circuit
                    are done in parallel in the final circuit
        5. More?

'''

def gate_count(circuit, include=None, exclude=None):
    '''
        Counts the number of gates
            1. of a given set of types, i.e. include
            2. given as any gate not in exclude
        Does this exclusively, and prioritizes include first
    '''
    count = 0
    
    if include is not None:
        for gate in circuit:
            if type(gate) in include:
                count += 1
    elif exclude is not None:
        for gate in circuit:
            if type(gate) not in exclude:
                count += 1
    else:
        count = len(circuit)
        
    return count

class LabeledOp:
    def __init__(self, instruction, index):
        self.instruction = instruction
        self.op = instruction[0]
        self.qubits = instruction[1]
        self.index = index
        
    def __eq__(self, other):
        if self is None or other is None:
            return False
        return self.index == other.index and self.instruction == other.instruction
    
    def __hash__(self):
        return self.index.__hash__()
    
    def __repr__(self):
        return f'{self.op} ({self.index})'

def gen_dag(circuit):
    '''
        Special Dag generator which accounts for Barriers without inserting them
        as gates
    '''
    dag = nx.DiGraph()
    
    qubit_last_use = {}
    
    null_op = LabeledOp([None, []], -1)
    
    for index, instruction in enumerate(circuit):        
        op = LabeledOp(instruction, index)
        
        for qubit in op.qubits:
            if qubit not in qubit_last_use:
                dag.add_edge(null_op, op)
            else:
                dag.add_edge(qubit_last_use[qubit], op)
                
            qubit_last_use[qubit] = op
            
    final_op = LabeledOp([None, []], float('inf'))
    for q in qubit_last_use:
        dag.add_edge(qubit_last_use[q], final_op)
        
    # Finally prune out the barriers
    for node in list(dag.nodes):
        if isinstance(op.instruction, qiskit.circuit.barrier.Barrier):
            # All things in adjacent to all things out
            for inedge in nx.in_edges(dag, node):
                for outedge in nx.out_edges(dag, node):
                    dag.add_edge(inedge[0], outedge[1])
            dag.remove_node(node)
      
    return dag
    # return nx.transitive_reduction(dag)
def gate_count(qc : qiskit.QuantumCircuit) -> int:
    total_count = 0
    for key in qc.count_ops():
        if key != 'barrier':
            total_count += qc.count_ops()[key]
    
    return total_count        
        
def circuit_depth(circuit):
    '''
        Computes the depth of an input circuit
        Equivalent to the length of the critical path in the circuit dag
    '''
    start_node = LabeledOp([None, []], -1)
    final_node = LabeledOp([None, []], float('inf'))
    
    dag = gen_dag(circuit)
    return nx.algorithms.dag.dag_longest_path_length(dag) - 1

def parallelism_score(circuit, raw=True):
    start_node = LabeledOp([None, []], -1)
    final_node = LabeledOp([None, []], float('inf'))
    
    dag = gen_dag(circuit)

    if raw:
        return (len(dag) - 2) / (nx.algorithms.dag.dag_longest_path_length(dag) - 1)
    else:
        P = (len(dag) - 2) / (nx.algorithms.dag.dag_longest_path_length(dag) - 1)
        mx = len(circuit.qubits)
        return P / mx

def fred_parallelism_score(circuit):
    return gate_count(circuit) / circuit_depth(circuit)

def parallel_lists(circuit: qiskit.QuantumCircuit):
    digraph = circuit_to_dag(circuit)
  
    parallel_ops = {}
    nodes = list(digraph.nodes())
    for node in nodes:
        parallel_ops[node] = set()

    for i in range(len(nodes)):
        for j in range(i):
            if not nx.has_path(digraph, nodes[i], nodes[j]):
                if not nx.has_path(digraph, nodes[j], nodes[i]):
                    parallel_ops[nodes[i]].add(nodes[j])
                    parallel_ops[nodes[j]].add(nodes[i])

    return parallel_ops
        
def preserved_parallelism(original_circuit, new_circuit):
  original_lists = parallel_lists(original_circuit)
  count = 0
  original_digraph = circuit_to_dag(original_circuit)
  int_to_op = {}
  for node in original_digraph.nodes():
    int_to_op[node.index] = node
  digraph = circuit_to_dag(new_circuit)
  time_steps = digraph_to_time_steps(digraph, None)
  for time_step in time_steps:
    for op_1 in time_step:
      if op_1.op is None:
        continue

      if op_1.op.name == "barrier":
        continue
      elif op_1.op.label == None:
        continue
      op_1_int = int(op_1.op.label)
      for op_2 in time_step:
        if op_2.op is None:
          continue

        if op_2.op.name == "barrier":
          continue
        elif op_2.op.label == None:
          continue
        op_2_int = int(op_2.op.label)
        if op_1 == op_2:
          continue
        if int_to_op[op_2_int] in original_lists[int_to_op[op_1_int]]:
          count += 1
          break
  list_length = sum([1 for i in original_lists if len(original_lists[i]) > 0])
  return count / list_length
    