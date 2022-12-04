import qiskit
from typing import Callable, Iterable, Tuple, List

import networkx as nx
import numpy as np
import copy

from collections import defaultdict


class DimensionData:

    def __init__(self,
                 size: int,
                 spacing: float):
        self.size = size
        self.spacing = spacing

class HardwareQubit:

    def __init__(self,
                 hardware_tuple: Tuple[int],
                 qiskit_qubit: qiskit.circuit.QuantumRegister):

        self.dims = len(hardware_tuple)
        self.hardware_tuple = hardware_tuple
        self.qiskit_qubit = qiskit_qubit

    def dim_loc(self, dim: int):
        assert dim < self.dims
        return self.hardware_tuple[dim]

    def __hash__(self):
        return hash(self.hardware_tuple)

    def __eq__(self, other):
        if self.dims != other.dims:
            return False

        for i in range(self.dims):
            if self.hardware_tuple[i] != other.hardware_tuple[i]:
                return False

        return True
    
    def __repr__(self):
        return f'{self.hardware_tuple}'


class Hardware:
    def __init__(self,
                 num_dimensions: int,
                 dimensions_length: List[int],
                 dimensions_spacing: List[float],
                 **kwargs) -> None:

        assert len(dimensions_length) == num_dimensions
        assert len(dimensions_spacing) == num_dimensions

        hardware_qubits = [[]]
        for length in dimensions_length:
            new_hardware_qubits = []
            for i in range(length):
                for qubit in hardware_qubits:
                    new_qubit = list(qubit)
                    new_qubit.append(i)
                    new_hardware_qubits.append(new_qubit)
            hardware_qubits = new_hardware_qubits
            
        self.qiskit_to_hardware = {}
        self.hardware_to_qiskit = {}
        self.hardware_loc_to_obj = {}
        self.active_holes = set()

        self.num_dimensions = num_dimensions
        self.qiskit_qubits = qiskit.circuit.QuantumRegister(len(hardware_qubits))
        self.hardware_qubits = []
        for index, qubit_tuple in enumerate(hardware_qubits):
            hwq = HardwareQubit(tuple(qubit_tuple), self.qiskit_qubits[index])
            self.qiskit_to_hardware[self.qiskit_qubits[index]] = hwq
            self.hardware_to_qiskit[hwq] = self.qiskit_qubits[index]
            self.hardware_loc_to_obj[tuple(qubit_tuple)] = hwq
            self.hardware_qubits.append(hwq)

        self.qiskit_to_hardware_orig = copy.copy(self.qiskit_to_hardware)
        self.hardware_to_qiskit_orig = copy.copy(self.hardware_to_qiskit)

        self.dimension_info = {}
        for dim in range(num_dimensions):
            self.dimension_info[dim] = DimensionData(dimensions_length[dim],
                                                     dimensions_spacing[dim])

        self.distances = defaultdict(dict)
        self.manhattan_distances = defaultdict(dict)
        for index, qubit1 in enumerate(self.hardware_qubits):
            for qubit2 in self.hardware_qubits[index:]:
                if qubit1 == qubit2:
                    self.distances[qubit1][qubit2] = 0
                    self.manhattan_distances[qubit1][qubit2] = 0
                else:
                    squared_sum = sum([(qubit1.dim_loc(i) - qubit2.dim_loc(i)) ** 2 for i in range(self.num_dimensions)])
                    distance = np.sqrt(squared_sum)
                    manhattan_distance = sum([abs(qubit1.dim_loc(i) - qubit2.dim_loc(i)) for i in range(self.num_dimensions)])
                    self.distances[qubit1][qubit2] = distance
                    self.distances[qubit2][qubit1] = distance
                    self.manhattan_distances[qubit1][qubit2] = manhattan_distance
                    self.manhattan_distances[qubit2][qubit1] = manhattan_distance
        self.orig_manhattan_distances = copy.copy(self.manhattan_distances)
        self.orig_distances_distances = copy.copy(self.distances)

        self.connectivity = nx.Graph()
        for qubit in self.hardware_qubits:
            self.connectivity.add_node(qubit)
            for other_qubit in self.manhattan_distances[qubit]:
                if self.manhattan_distances[qubit][other_qubit] == 1:
                    if not self.connectivity.has_edge(qubit, other_qubit):
                        spacing = None
                        for i in range(self.num_dimensions):
                            val = abs(qubit.dim_loc(i) - other_qubit.dim_loc(i))
                            if val == 1:
                                spacing = self.dimension_info[i].spacing
                        assert spacing is not None
                        self.connectivity.add_edge(qubit, other_qubit, weight=spacing)
        self.orig_connectivity = self.connectivity.copy()

    def get_constrained_qubit_sets(self, number_qubits, method="wide"):
        square_dimension = np.sqrt(number_qubits)
        dims = [None, None]
        if square_dimension.is_integer():
            dims = [int(square_dimension), int(square_dimension)]
        else:
            dims = [int(np.ceil(square_dimension)), int(np.ceil(square_dimension))]

        if method == "tight":
            dims[-1] = int(np.ceil(number_qubits / dims[0]))

        qubit_sets = []
        curr_first = 0
        break_out_outer = False
        while curr_first + dims[0] <= self.dimension_info[0].size:

            curr_second = 0
            break_out_inner = False
            while curr_second + dims[1] <= self.dimension_info[1].size:

                curr_set = []
                for i in range(curr_first, curr_first + dims[0]):
                    for j in range(curr_second, curr_second + dims[1]):
                        curr_set.append(self.hardware_loc_to_obj[(i, j)])
                qubit_sets.append(curr_set)
                curr_second += dims[1]

                if curr_second + dims[1] >= self.dimension_info[1].size and not break_out_inner:
                    curr_second = self.dimension_info[1].size - dims[1]
                    break_out_inner = True

            curr_first += dims[0]
            if curr_first + dims[0] > self.dimension_info[0].size and not break_out_outer:
                curr_first = self.dimension_info[0].size - dims[0]
                break_out_outer = True
        return qubit_sets

    def get_non_overlapping_constrained_qubit_sets(self, number_qubits, method="wide", spacing=0):
        num_sections = len(self.qiskit_qubits) // number_qubits
        satisified = None
        if num_sections < 3:
            satisified = (2, 1)
        elif num_sections < 4:
            satisified = (3, 1)
        curr_num = num_sections
        while satisified is None:
            factor_list = []
            for i in range(1, int(curr_num/2) + 1):
                if curr_num % i == 0:
                    factor_list.append((i, curr_num // i))
            factor_pair = factor_list[-1]
            if factor_pair[0] / factor_pair[1] >= 0.5:
                satisified = factor_pair
            else:
                curr_num = curr_num - 1
        dims = [None, None]
        if self.dimension_info[0].size > self.dimension_info[1].size:
            if np.floor(self.dimension_info[0].size / satisified[1]) * np.floor(self.dimension_info[1].size / satisified[0]) >= number_qubits:
                dims[0] = int(np.floor(self.dimension_info[0].size / satisified[1]))
                dims[1] = int(np.floor(self.dimension_info[1].size / satisified[0]))
            else:
                dims[0] = int(np.ceil(self.dimension_info[0].size / satisified[1]))
                dims[1] = int(np.ceil(self.dimension_info[1].size / satisified[0]))
        else:
            if np.floor(self.dimension_info[1].size / satisified[1]) * np.floor(self.dimension_info[0].size / satisified[0]) >= number_qubits:
                dims[1] = int(np.floor(self.dimension_info[1].size / satisified[1]))
                dims[0] = int(np.floor(self.dimension_info[0].size / satisified[0]))
            else:
                dims[1] = int(np.ceil(self.dimension_info[1].size / satisified[1]))
                dims[0] = int(np.ceil(self.dimension_info[0].size / satisified[0]))
        qubit_sets = []
        overlapping_sets = []
        already_captured = set()
        curr_first = 0
        break_out_outer = False
        while curr_first + dims[0] <= self.dimension_info[0].size:

            curr_second = 0
            break_out_inner = False
            set_add_outer = False
            while curr_second + dims[1] <= self.dimension_info[1].size:

                curr_set = []
                for i in range(curr_first, curr_first + dims[0]):
                    for j in range(curr_second, curr_second + dims[1]):
                        curr_set.append(self.hardware_loc_to_obj[(i, j)])
                qubit_sets.append(curr_set)
                set_add = True
                for i in curr_set:
                    if i in already_captured:
                        qubit_sets = qubit_sets[:-1]
                        set_add = False
                        break
                set_add_outer = set_add or set_add_outer
                curr_second += dims[1] + (spacing if set_add else 0)

            curr_first += dims[0] + (spacing if set_add_outer else 0)

        return qubit_sets

    def reset(self):
        self.active_holes = set()
        self.connectivity = self.orig_connectivity.copy()
    
    def add_hole(self, hole):
        self.active_holes.add(hole)
        self.connectivity.remove_node(hole)    

    def qubits_in_radius(self, q : HardwareQubit, r : float) -> List[HardwareQubit]:
        if q in self.active_holes:
            return []

        self.validate_qubits(q)

        in_radius = []

        for qubit in self.distances[q]:
            if self.distances[q][qubit] <= r:
                if qubit not in self.active_holes: in_radius.append(qubit)

        return in_radius

    def validate_qubits(self, *qs : Iterable[HardwareQubit]) -> None:
        '''
            Takes a list of possible qubits (i, j) and ensures they are in the hardware map
            i.e. self.connectivity
        '''
        for q in qs:
            assert q in self.connectivity, f'{q} is not in hardware.'

    def distance(self, q1, q2):
        '''
            Returns the float distance between qubits 1 and 2 if in the hardware
        '''
        self.validate_qubits(q1, q2)

        return self.distances[q1][q2]

    def integer_distance(self, q1, q2):
        '''
            Returns the integer distance (L2 norm in the underlying 2D mesh)
            between two qubits
        '''
        self.validate_qubits(q1, q2)

        return self.manhattan_distances[q1][q2]
