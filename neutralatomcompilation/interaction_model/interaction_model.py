import qiskit
import networkx as nx
import copy
from typing import Callable, List, Tuple, Union, Iterable

from ..hardware import Hardware

class InteractionModel:

    def __init__(self,
                 hardware: Hardware,
                 d_to_r : Callable[[float], float],
                 max_int_dist: float = float('inf'),
                 ) -> None:
        '''
            Gives detail on how interactions are done on the hardware

            Currently queries what interactions can be done given a piece of
            hardware

            TODO:
                Can query error costs of performing certain operations

            Inputs:
                hardware: specifies the dimensions of the underlying hardware

                d_to_r: A function taking a interaction distance d and
                producing the restricted zone of radius r around each qubit

                max_int_dist: Fixes the maximum size of interactions in the
                hardware;

        '''
        self.hardware = hardware
        self.d_to_r = d_to_r
        self.max_interaction_distance = max_int_dist
        self.interaction_graph = nx.Graph()
        for hardware_qubit in self.hardware.hardware_qubits:
            for qubit in self.hardware.qubits_in_radius(hardware_qubit, self.max_d()):
                if not self.interaction_graph.has_edge(hardware_qubit, qubit):
                    self.interaction_graph.add_edge(hardware_qubit, qubit)
        self.original_graph = self.interaction_graph.copy()

    def reset_graph(self):
        self.interaction_graph = self.original_graph.copy()

    def remove_qubit(self, qubit):
        self.interaction_graph.remove_node(qubit)

    def max_d(self) -> Union[float, None]:
        '''
            Returns max interaction distance
        '''
        return self.max_interaction_distance

    def path_using_max_d(self, q1, q2):
        return nx.dijkstra_path(self.interaction_graph, q1, q2)

    def valid_interaction_constrained(self,
                                      qs : List[Tuple[int, int]],
                                      constraints : Iterable[List[Tuple[int, int]]]
                                      ) -> bool:
        '''
            Given a set of qubits (qs) designated to interact
            And a set of constraining interactions (constraints)

            Decide if the interaction on qs can be done in parallel

            i.e. The zone of use for the new operation cannot intersect

        '''

        illegal_qubits_constraints = set()
        
        new_constraints = []
        for c in constraints:
            if len(c) != 0:
                new_constraints.append(c)
                
        constraints = new_constraints
        
        if len(constraints) == 0:
            return self.valid_interaction_distance(qs)


        for c in constraints:

            if len(c) > 1:
                max_dist_pairwise = max(self.hardware.distance(q1, q2) for q1 in c for q2 in c if q1 != q2)
            else:
                max_dist_pairwise = 0

            # Verify the qubit sets satisfy the max_interaction_distance constraint
            assert max_dist_pairwise <= self.max_d(), 'Illegal qubit set in constraints.'

            r = self.d_to_r(max_dist_pairwise)
            for q in c:
                for q_in_r in self.hardware.qubits_in_radius(q, r):
                    illegal_qubits_constraints.add(q_in_r)

        if len(qs) == 1:
            return not qs[0] in illegal_qubits_constraints

        max_dist_qs = max(self.hardware.distance(q1, q2) for q1 in qs for q2 in qs if q1 != q2)

        # Verify the qubit sets satisfy the max_interaction_distance constraint
        assert max_dist_qs <= self.max_d(), 'Illegal qubit set in qs.'

        r = self.d_to_r(max_dist_qs)
        illegal_qubits_qs = set()
        for q in qs:
            for item in self.hardware.qubits_in_radius(q, r):
                illegal_qubits_qs.add(item)

        # If the intersection is empty => valid parallel ops
        return len(illegal_qubits_constraints.intersection(illegal_qubits_qs)) == 0

    def valid_interaction_distance(self,
                                   qs: List[Tuple[int, int]]
                                   ) -> bool:
        '''
            Verifies the qubits of an operation are within the max distance
        '''
        if len(qs) == 1:
            return True
        max_dist_qs = max(self.hardware.distance(q1, q2) for q1 in qs for q2 in qs if q1 != q2)
        return max_dist_qs <= self.max_d()