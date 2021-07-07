import qiskit
import networkx as nx
import enum
import copy
from collections import defaultdict
import time
import numpy as np

from typing import List

from ..hardware.hardware_general import Hardware2, HardwareQubit

class ReRouteStrategy(enum.Enum):
  Fail = 1
  Swap = 2

  def __eq__(self, other):
    return self.value == other.value

class ShiftStrategy(enum.Enum):
  NaiveMaxSpace = 1
  NaiveMinMovement = 2
  Clever = 3
  InteractionGraph = 4

  def __eq__(self, other):
    return self.value == other.value

class Remapper:

  def __init__(self):
    self.data = {}
    self.access_counter = 0
    self.read_count = 0
    self.write_count = 0
    self.access_time = 0

  def __getitem__(self, key):
    s = time.time()
    self.access_counter += 1
    self.read_count += 1
    v = self.data[key]
    e = time.time()
    self.access_time += e - s
    return v

  def __setitem__(self, key, value):
    s = time.time()
    self.write_count += 1
    self.access_counter += 1
    self.data[key] = value
    e = time.time()
    self.access_time += e - s

  def reset(self):
    self.data = {}
    self.access_counter = 0
    self.read_count = 0
    self.write_count = 0
    self.access_time = 0

class HoleHandler:
  def __init__(self, hardware, interaction_model):

    self.active_holes = set()
    self.hardware = hardware
    self.interaction_model = interaction_model
    self.interaction_graph = self.interaction_model.interaction_graph.copy()
    
    self.qh_mapping = copy.copy(hardware.qiskit_to_hardware)
    self.hq_mapping = copy.copy(hardware.hardware_to_qiskit)

    # Where is my qubit now?
    self.hh_mapping = Remapper()
    # What hardware qubit is here now?
    self.hh_mapping_reverse = Remapper()
    for h in self.hq_mapping:
      self.hh_mapping[h] = h
      self.hh_mapping_reverse[h] = h
    self.hh_mapping_reverse.read_count = 0
    self.hh_mapping_reverse.write_count = 0
    self.hh_mapping_reverse.access_counter = 0
    self.added_swaps = 0

  def reset(self):
    self.active_holes = set()
    self.qh_mapping = copy.copy(self.hardware.qiskit_to_hardware)
    self.hq_mapping = copy.copy(self.hardware.hardware_to_qiskit)
    self.hh_mapping.reset()
    self.hh_mapping_reverse.reset()
    for h in self.hq_mapping:
      self.hh_mapping[h] = h
      self.hh_mapping_reverse[h] = h
    self.hh_mapping_reverse.read_count = 0
    self.hh_mapping_reverse.write_count = 0
    self.hh_mapping_reverse.access_counter = 0
    self.interaction_graph = self.interaction_model.interaction_graph.copy()



  def shift(self, location, used_qubits, dim, direct, strat="short"):
    to_remove = list(location)
    shift_loc = list(location)
    to_shift_from = list(location)
    dimension = dim
    dim_end = self.hardware.dimension_info[dim].size - 1
    start_loc = 0 if direct == -1 else dim_end
    end_loc = shift_loc[dimension]
    first_qubit = True
    for i in range(start_loc, end_loc, -1*direct):
      shift_loc[dimension] = i
      hw_qubit = self.hardware.hardware_loc_to_obj[tuple(shift_loc)]
      #print(hw_qubit)
      if hw_qubit in self.active_holes:
        #print("hole")
        continue
      curr_loc = i + -1 * direct
      to_shift_from[dimension] = curr_loc
      while self.hardware.hardware_loc_to_obj[tuple(to_shift_from)] in self.active_holes:
        curr_loc += -1 * direct
        to_shift_from[dimension] = curr_loc
      to_shift_from_hw_qubit = self.hh_mapping_reverse[self.hardware.hardware_loc_to_obj[tuple(to_shift_from)]]
      if strat == "short":
        if to_shift_from_hw_qubit not in used_qubits and first_qubit:
          continue
      if first_qubit:
        first_qubit = False
        #print("first", hw_qubit)
        self.hh_mapping_reverse[hw_qubit] = None
      #print(to_shift_from_hw_qubit, ":", self.hh_mapping[to_shift_from_hw_qubit], "-->", hw_qubit)
      self.hh_mapping[to_shift_from_hw_qubit] = hw_qubit
      self.hh_mapping_reverse[hw_qubit] = to_shift_from_hw_qubit
    #print("step done")
    
    hole = self.hardware.hardware_loc_to_obj[tuple(to_remove)]
    self.hh_mapping_reverse[hole] = None

  def im_shift(self, qiskit_qubits, hole, stategy="space"):
    if self.hh_mapping_reverse[hole] not in qiskit_qubits:
      #print("not used", hole, self.hh_mapping_reverse[hole])
      self.hh_mapping_reverse[hole] = None
      return True

    frontier = [hole]
    found_node = None
    visited = set()
    paths = {hole: None}
    while found_node is None and len(frontier) > 0:
      new_frontier = []
      for f in frontier:
        visited.add(f)
        ns = sorted(list(self.interaction_graph.neighbors(f)), key= lambda x: self.hardware.distance(f, x))
        self.hh_mapping.access_counter += np.log2(len(ns))
        for n in ns:
          if n in visited:
            continue
          paths[n] = f
          currently_here = self.hh_mapping_reverse[n]
          self.hh_mapping.access_counter += np.log2(len(qiskit_qubits))
          if currently_here not in qiskit_qubits and currently_here not in self.active_holes:
            found_node = n
            break
          new_frontier.append(n)
        if not found_node is None:
          break
      frontier = new_frontier
    if found_node is None:
      return False

    path = [found_node]
    curr = found_node
    while curr is not None:
      next_n = paths[curr]
      path.append(next_n)
      curr = next_n

    current_node = path[0]
    self.hh_mapping_reverse[current_node] = currently_here
    self.hh_mapping[currently_here] = None
    for index, h_qubit in enumerate(path[1:]):
      pulling_from = self.hh_mapping_reverse[h_qubit]
      self.hh_mapping_reverse[path[index - 1]] = pulling_from
      self.hh_mapping[pulling_from] = path[index - 1]
    self.hh_mapping_reverse[hole] = None

    return True


  # This shifts in whichever direction has the most open spots
  def naive_shift(self, qiskit_qubits, hole, stategy="space"):
    if self.hh_mapping_reverse[hole] not in qiskit_qubits:
      #print("not used", hole, self.hh_mapping_reverse[hole])
      self.hh_mapping_reverse[hole] = None
      return True
    location = hole.hardware_tuple

    if stategy == "space":
      best_dist = 0
      best_shift = None
    else:
      best_dist = float("inf")
      best_shift = None

    for dimension, l in enumerate(location):
      dim_start = 0
      dim_end = self.hardware.dimension_info[dimension].size - 1

      all_holes = True
      for i in range(dim_start, location[dimension] + 1):
        loc = list(location)
        loc[dimension] = i
        hw_qubit = self.hardware.hardware_loc_to_obj[tuple(loc)]
        if hw_qubit not in self.active_holes and self.hh_mapping_reverse[hw_qubit] not in qiskit_qubits:
          all_holes = False
        if self.hh_mapping_reverse[hw_qubit] in qiskit_qubits:
          if stategy == "space":
            if i > best_dist and not all_holes:
              best_shift = (dimension, -1)
              best_dist = i
          else:
            if location[dimension] - i < best_dist and not all_holes:
              best_shift = (dimension, -1)
              best_dist = location[dimension] - i
          break

      all_holes = True
      for i in range(dim_end - 1, location[dimension] - 1, -1):
        loc = list(location)
        loc[dimension] = i
        hw_qubit = self.hardware.hardware_loc_to_obj[tuple(loc)]
        if hw_qubit not in self.active_holes and self.hh_mapping_reverse[hw_qubit] not in qiskit_qubits:
          all_holes = False
        if self.hh_mapping_reverse[hw_qubit] in qiskit_qubits:
          if stategy == "space":
            if dim_end - i > best_dist and not all_holes:
              best_shift = (dimension, 1)
              best_dist = dim_end - (i + 1)
          else:
            if i - location[dimension] < best_dist and not all_holes:
              best_shift = (dimension, 1)
              best_dist = i - location[dimension]
          break

    if best_shift is None:
      return False

    self.shift(location, qiskit_qubits, best_shift[0], best_shift[1], "short")

    return True

  def fail_reroute(self, circuit) -> bool:
    for instruction in circuit:
      if len(instruction[1]) < 2:
        continue
      qubits = [self.hh_mapping[self.qh_mapping[i]] for i in instruction[1]]

      try:
        if not self.interaction_model.valid_interaction_distance(qubits):
          return False
      except Exception as e:
        print(instruction[1], [self.qh_mapping[i] for i in instruction[1]], qubits)
        raise e
    return True

  def insert_swaps(self, circuit):
    c = qiskit.QuantumCircuit(*circuit.qregs)
    for instruction in circuit:
      if len(instruction[1]) < 2:
        c.append(*instruction)
        continue
      qubits = [self.hh_mapping[self.qh_mapping[i]] for i in instruction[1]]
      try:
        self.interaction_model.valid_interaction_distance(qubits)
      except Exception as e:
        print([self.qh_mapping[i] for i in instruction[1]], qubits)
        raise e
      if self.interaction_model.valid_interaction_distance(qubits):
        c.append(*instruction)
        continue
      else:
        try:
          path = nx.dijkstra_path(self.interaction_graph, qubits[0], qubits[1])
        except nx.NetworkXNoPath:
          return None
        except Exception as e:
          raise e
        for index, p1 in enumerate(path[:-2]):
          p2 = path[index + 1]
          qq1 = self.hh_mapping_reverse[p1].qiskit_qubit
          qq2 = self.hh_mapping_reverse[p2].qiskit_qubit
          try:
            self.added_swaps += 1
            c.swap(qq1, qq2)
          except Exception as e:
            print(path)
            print(self.hh_mapping_reverse[p1], self.hh_mapping_reverse[p2])
            print(qq1, qq2)
            raise e
        qq1 = self.hh_mapping_reverse[path[-2]].qiskit_qubit
        qq2 = self.hh_mapping_reverse[path[-1]].qiskit_qubit
        c.append(instruction[0], (qq1, qq2))
        for index in range(len(path) - 2, 0, -1):
          p1 = path[index]
          p2 = path[index - 1]
          qq1 = self.hh_mapping_reverse[p1].qiskit_qubit
          qq2 = self.hh_mapping_reverse[p2].qiskit_qubit
          c.swap(qq1, qq2)
    return c

  def reroute_with_holes(self,
                         circuit,
                         holes,
                         qiskit_qubits=None,
                         shift_strategy=ShiftStrategy.NaiveMinMovement,
                         route_strategy=ReRouteStrategy.Fail):

    if qiskit_qubits is None:
      qiskit_qubits = set()
      for instruction in circuit:
        for qubit in instruction[1]:
          qiskit_qubits.add(self.qh_mapping[qubit])
    #print(qiskit_qubits)
    if len(qiskit_qubits) > len(self.hardware.hardware_qubits) - len(self.active_holes) - len(holes):
      return None, "size"

    if route_strategy == ReRouteStrategy.Fail:
      for h in holes:
        if not self.naive_shift(qiskit_qubits, h, ""):
          self.active_holes.add(h)
          return None, "shift"
        self.active_holes.add(h)
      works = self.fail_reroute(circuit)
      if works: return circuit, "good"
      else: return None, "routing"
    elif route_strategy == ReRouteStrategy.Swap:
      curr_circuit = circuit
      for h in holes:
        if shift_strategy == ShiftStrategy.NaiveMinMovement:
          if not self.naive_shift(qiskit_qubits, h, ""):
            self.active_holes.add(h)
            return None, "shift"
        elif shift_strategy == ShiftStrategy.InteractionGraph:
          if not self.im_shift(qiskit_qubits, h, ""):
            self.active_holes.add(h)
            return None, "shift"
        self.active_holes.add(h)
        self.interaction_graph.remove_node(h)
      curr_circuit = self.insert_swaps(curr_circuit)
      if curr_circuit is None:
        return None, "routing"
      return curr_circuit, "good"
    else:
      print("nothing")
