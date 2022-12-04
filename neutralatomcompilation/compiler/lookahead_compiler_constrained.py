import qiskit
import networkx as nx
import numpy as np
from collections import defaultdict

from neutralatomcompilation.compiler.compiler import Compiler

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
    
class SwapOp:
    def __init__(self, source, dest):
        self.qubits = [source, dest]
        
        
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
        dag.add_edge(qubit_last_use[q], final_op)
    
    return nx.algorithms.transitive_reduction(dag)

import qiskit
import networkx as nx
from collections import defaultdict
from typing import List, Tuple, DefaultDict, Dict, Set, Callable

import neutralatomcompilation as NAC
from neutralatomcompilation.compiler import Compiler

class LookaheadCompilerConstrained(Compiler):
    start_node = LabeledOp([None, []], -1)
    final_op = LabeledOp([None, []], float('inf'))

    def __init__(self,
                 interaction_model,
                 hardware,
                 constrained_method="wide") -> None:
        self.constrained_method = constrained_method
        super().__init__(interaction_model, hardware)
    
    def build_interaction_edges(self, dist_dict, frontier, weighting_function):
        edges = defaultdict(float)
        counted_los = {}
        
        for lo1 in frontier:
            for lo2 in dist_dict[lo1]:
                if len(lo2.qubits) > 1:
                    if lo2 not in counted_los:
                        for q1 in lo2.qubits:
                            for q2 in lo2.qubits:
                                if q1 != q2:
                                    edges[frozenset({q1, q2})] += weighting_function(dist_dict[lo1][lo2])
                    elif dist_dict[lo1][lo2] < counted_los[lo2]:
                        for q1 in lo2.qubits:
                            for q2 in lo2.qubits:
                                if q1 != q2:
                                    edges[frozenset({q1, q2})] += weighting_function(dist_dict[lo1][lo2])
                                    edges[frozenset({q1, q2})] -= weighting_function(counted_los[lo2])
                                    
                    counted_los[lo2] = dist_dict[lo1][lo2]
        return edges

    def _dist_norm(self, A, B):
        if A == B:
            return 0
        if self.hardware.distance(A, B) <= self.interaction_model.max_interaction_distance:
            return 0
        
        # This will probably break
        # if self.interaction_model.max_interaction_distance < 2:
        #    return self.hardware.integer_distance(A, B) - 1
        return np.floor(self.hardware.distance(A, B) / (self.interaction_model.max_interaction_distance ))

    def _weight_to_(self, wedges):
        weight_to = {}
        for edge in wedges:
            ledge = list(edge)
            if ledge[0] not in weight_to:
                weight_to[ledge[0]] = {}
            if ledge[1] not in weight_to:
                weight_to[ledge[1]] = {}
            
            weight_to[ledge[1]][ledge[0]] = wedges[edge]
            weight_to[ledge[0]][ledge[1]] = wedges[edge]

        return weight_to
    
    def _initial_mapping(self,
                         circuit,
                         dist_dict,
                         frontier,
                         weighting_function):
        
        logical_to_physical = {}
        physical_to_logical = {}
        placed_qubits = set()
        
        edges = self.build_interaction_edges(dist_dict, frontier, weighting_function)
        sorted_edges = sorted(edges, key=lambda t : -edges[t])

        self.viable_nodes_sets = self.hardware.get_constrained_qubit_sets(len(circuit.qubits), self.constrained_method)
        first_set = self.viable_nodes_sets[0]

        constrained_graph = self.interaction_model.interaction_graph.copy()
        for node in list(constrained_graph.nodes()):
            if node not in first_set:
                constrained_graph.remove_node(node)
        
        try:
            center_qubit = nx.center(constrained_graph)[0]
        except:
            return None, None
        
        if len(sorted_edges) == 0:
            return {}, {}
        edge_to_place = list(sorted_edges[0])
        
        logical_to_physical[edge_to_place[0]] = center_qubit
        physical_to_logical[center_qubit] = edge_to_place[0]
        
        placed_qubits.add(edge_to_place[0])
        
        weight_to = self._weight_to_(edges) 

        if len(sorted_edges) == 0:
            return logical_to_physical, physical_to_logical

        placed_qubits.add(edge_to_place[0])
        
        def _place_qubit_(q):
            best_location, best_score = None, float('inf')
            for node in first_set:
                if node not in physical_to_logical:
                    # If this node is unoccuped
                    score = 0
                    for target in weight_to[q]:
                        if target in placed_qubits:
                            score += weight_to[q][target] * self._dist_norm(logical_to_physical[target], node)
                    if score < best_score:
                        best_location, best_score = node, score
            
            logical_to_physical[q] = best_location
            physical_to_logical[best_location] = q
            placed_qubits.add(q)
        
        def _find_next_to_place_():
            '''
                Choose qubit with largest weight to existing placed qubits
            '''
            best_weight, which = float('-inf'), None

            for i in weight_to:
                if i not in placed_qubits:
                    weight = 0
                    for j in weight_to[i]:
                        if j in placed_qubits:
                            weight += weight_to[i][j]
                    if weight > best_weight:
                        best_weight, which = weight, i
            return which
        
        
        while len(placed_qubits) < len(weight_to):
            which = _find_next_to_place_()
            _place_qubit_(which)
        
        return logical_to_physical, physical_to_logical
    
    def _map_qubit(self, q, weight_to, l2p, p2l):
        best_location, best_score = None, float('inf')
        for node in self.hardware.connectivity:
            if node not in p2l:
                # If this node is unoccuped
                score = 0
                if q in weight_to:
                    for target in wedges[q]:
                        if target in placed_qubits:
                            score += wedges[q][target] * self._dist_norm(l2p[target], node)
                if score < best_score:
                    best_location, best_score = node, score
            
            l2p[q] = best_location
            p2l[best_location] = q
    
    def compile(self,
                 circuit,
                 lookahead_distance,
                 weighting_function,
                 extra_mapping_params=None
                ):
        dag = circuit_to_dag(circuit)
        frontier = []
        for op in dag[self.start_node]:
            if op != self.final_op:
                frontier.append(op)
                
        # Based on this circuit dag compute as much as we can about future interactions
        all_paths = nx.algorithms.shortest_paths.all_pairs_shortest_path_length(dag, 
                                                                                cutoff=lookahead_distance)
        
        dist_dict = {}
        for item in all_paths:
            to_del = []
            for i2 in item[1].keys():
                if len(i2.qubits) < 2:
                    to_del.append(i2)
            for i2 in to_del:
                item[1].pop(i2)
                
            dist_dict[item[0]] = item[1]
        if len(circuit.qubits) > len(self.hardware.connectivity.nodes()):
            #print(len(circuit.qubits), len(self.hardware.connectivity.nodes()))
            #print("size")
            return None
            
        l2p, p2l = self._initial_mapping(circuit, dist_dict, frontier, weighting_function)
        self.original_l2p = l2p.copy()
        self.original_p2l = p2l.copy()
        if l2p is None:
            print("no mapping")
            return None
        return self._schedule_and_route(dag, l2p, p2l, lookahead_distance, weighting_function, dist_dict, frontier)
    

    def _schedule_and_route(self,
                            dag,
                            l2p,
                            p2l,
                            lookahead_distance,
                            weighting_function,
                            dist_dict,
                            frontier):

        ccc = 0
        
        def _update_dag_(op_to_remove):
            for target in dag[op_to_remove]:
                dag.add_edge(self.start_node, target)
            dag.remove_node(op_to_remove)
        
        def _score_swaps_(qubits,
                          wig,
                          l2p, 
                          p2l,
                          to_execute):
            cur_dist = {}
            literal_dist = {}
            psuedo_dists = {}
            for q1 in qubits:
                cur_dist[q1] = 0
                literal_dist[q1] = 0
                psuedo_dists[q1] = 0
                
                for q2 in qubits:
                    if q1 == q2:
                        continue
                    cur_dist[q1] += self._dist_norm(l2p[q1], l2p[q2])
                    literal_dist[q1] += self.hardware.distance(l2p[q1], l2p[q2])
                    psuedo_dists[q1] += self.hardware.integer_distance(l2p[q1], l2p[q2])
            
            def _get_valid_swaps_(q):
                # q is a logical qubit
                # locs is a list of places q can go to
                locs = []
                pq = l2p[q]
                pqs = [l2p[qq] for qq in qubits if qq != q]
                for node in self.hardware.qubits_in_radius(pq, self.interaction_model.max_d()):
                    dist = 0
                    psuedo_dist = 0
                    if node in pqs:
                        continue
                    if not node in self.viable_nodes_sets[0]:
                        continue
                    for pq2 in pqs:
                        if pq2 == pq:
                            continue
                        if node != pq2:
                            psuedo_dist += self.hardware.integer_distance(pq2, node)
                            dist += self.hardware.distance(pq2, node)
                    if dist < literal_dist[q]:
                        if self.interaction_model.valid_interaction_constrained(
                            [node, pq],
                            [[l2p[q] for q in o.qubits] for o in to_execute]):
                            locs.append(node)
                return locs
                                
            def _score_(q, p):
                s = 0
                
                d_score = cur_dist[q]
                for q2 in qubits:
                    if q2 != q:
                        d_score -= self._dist_norm(l2p[q2], p)
                        
                future = 0
                for q in wig:
                    for targ in wig[q]:
                        if targ in l2p:
                            p_to_targ = self._dist_norm(l2p[targ], p)
                            q_to_targ = self._dist_norm(l2p[targ], l2p[q])
                            
                            future += wig[q][targ] * (q_to_targ - p_to_targ)
                            
                if p in p2l:
                    if p2l[p] in wig:
                        for targ in wig[p2l[p]]:
                            p_to_targ = self._dist_norm(l2p[targ], p)
                            q_to_targ = self._dist_norm(l2p[targ], l2p[q])
                            
                            future += wig[p2l[p]][targ] * (p_to_targ - q_to_targ)
                            
                # -- Flagged for needing possible weighting --
                return d_score + future
            
            output_dict = {}
            for q in qubits:
                best_score, best_loc = -float('inf'), None
                for p in _get_valid_swaps_(q):
                    score = _score_(q, p)
                    if score > best_score:
                        best_score, best_loc = score, p
                        
                output_dict[q] = (best_score, best_loc)
            return output_dict

            
        output_circuit = qiskit.QuantumCircuit(self.hardware.qiskit_qubits)
        swapped_previously = defaultdict(set)

        while len(frontier) > 0:
            # ---- FLAGGED FOR POTENTIAL SLOWNESS ----
            wedges = None
            weight_to = None
            # ----
            
            to_execute = []
            currently_executable = []
            requires_swaps = []
            for op in frontier:
                in_edges = dag.in_edges(op)
                runnable = True
                for edge in in_edges:
                    if edge[0] != self.start_node:
                        runnable = False
                        break
                if not runnable:
                    continue
                    
                # Map qubits which are unmapped
                for q in op.qubits:
                    if q not in l2p:
                        if weight_to is None:
                            wedges = self.build_interaction_edges(dist_dict, frontier, weighting_function)
                            weight_to = self._weight_to_(wedges)
                        self._map_qubit(q, weight_to, l2p, p2l)

                if len(op.qubits) == 1:
                    to_execute.append(op)
                else:
                    # Divide up what can be done and what can't be done right now
                    if self.interaction_model.valid_interaction_distance([l2p[qq] for qq in op.qubits]):
                        currently_executable.append(op)
                    else:
                        requires_swaps.append(op)
                    
                    
            # Choose a set of ops to do, 
            # -- Flagged for Naive / Could be Improved --
            for op in currently_executable:
                if self.interaction_model.valid_interaction_constrained(
                    [l2p[qq] for qq in op.qubits],
                    [[l2p[qq] for qq in o.qubits] for o in to_execute]):
                    to_execute.append(op)

            if len(requires_swaps) > 0:
                frontier_no_exec = set()
                for op in frontier:
                    if op not in to_execute:
                        frontier_no_exec.add(op)
                    else:
                        for o2 in dag[op]:
                            frontier_no_exec.add(o2)
                        
                # Don't care about what the circuit is already going to execute
                wedges_no_exec = self.build_interaction_edges(dist_dict,
                                                              frontier_no_exec,
                                                              weighting_function)
                
                weight_to = self._weight_to_(wedges_no_exec)
                
                # Find all the possible SWAP scores
                all_scores = {}
                for op in requires_swaps:
                    scores = _score_swaps_(op.qubits, weight_to, l2p, p2l, to_execute)
                    for q in scores:
                        if scores[q][1] in swapped_previously[q]:
                            continue
                        if scores[q][1] in p2l:
                            if l2p[q] in swapped_previously[p2l[scores[q][1]]]:
                                continue
                        all_scores[q] = scores[q]
                        
                sorted_scores_keys = sorted(all_scores.keys(), key=lambda x: all_scores[x][0], reverse=True)
                        
                def _pick_best_swap_(sorted_scores_keys, all_scores):
                    for key in sorted_scores_keys:
                        source = l2p[key]
                        dest = all_scores[key][1]
                        
                        if dest is not None:
                            constraints = []
                            for o in to_execute:
                                if isinstance(o, SwapOp):
                                    constraints.append(op.qubits)
                                else:
                                    constraints.append([l2p[qq] for qq in o.qubits])
                            if self.interaction_model.valid_interaction_constrained([source, dest], constraints):
                                to_execute.append(SwapOp(source, dest))
                                return
                prev_len = len(to_execute)   
                _pick_best_swap_(sorted_scores_keys, all_scores)
                if prev_len == 0:
                    if prev_len >= len(to_execute):
                        print("no swaps")
                        return None
                    assert prev_len < len(to_execute), 'no swap found.'

            for op in to_execute:
                if isinstance(op, SwapOp):
                    swapped_previously[p2l[op.qubits[0]]].add(op.qubits[0])
                    qiskit_qubits = [q.qiskit_qubit for q in op.qubits]
                    output_circuit.swap(qiskit_qubits[0], qiskit_qubits[1])
                    
                    l1 = p2l[op.qubits[0]] if op.qubits[0] in p2l else None
                    l2 = p2l[op.qubits[1]] if op.qubits[1] in p2l else None
                    
                    if l1 is not None and l2 is not None:
                        # Atomic swap
                        p2l[op.qubits[1]], p2l[op.qubits[0]] = p2l[op.qubits[0]], p2l[op.qubits[1]]
                        l2p[l1], l2p[l2] = l2p[l2], l2p[l1]
                    elif l1 is not None:
                        p2l.pop(op.qubits[0])
                        p2l[op.qubits[1]] = l1
                        l2p[l1] = op.qubits[1]
                    elif l2 is not None:
                        p2l.pop(op.qubits[1])
                        p2l[op.qubits[0]] = l2
                        l2p[l2] = op.qubits[0]
                    else:
                        print("ERRRRRR")
                else:
                    for q in op.qubits: swapped_previously[q].clear()
                    qiskit_qubits = [l2p[q].qiskit_qubit for q in op.qubits]
                    op.op.label = str(op.index)
                    output_circuit.append(op.op, qiskit_qubits)
                    new_op = LabeledOp(output_circuit[-1], op.index)
                    
                    _update_dag_(op)
                    
            frontier.clear()
            for op in dag[self.start_node]:
                if len(dag.in_edges(op)) == 1 and op != self.final_op:
                    frontier.append(op)
                    
            ccc += 1
            
            #output_circuit.barrier()
            
        return output_circuit
