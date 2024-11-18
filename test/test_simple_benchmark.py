import pickle
import numpy as np
import inspect as ins
from numba import jit, njit
from openalea.mtg.traversal import post_order2, pre_order2
import timeit


# FAKE PARAMETERS
Km_Nm_root_LATS = 1e1
Km_Nm_root_HATS = 1e-3
begin_N_regulation = 1
span_N_regulation = 6e-5
vmax_Nm_root = 1e-6
transport_C_regulation = 7e-3
soil_Nm = 1e2


# FAKE FUNCTIONS FOR BELLOW TESTS

# Testing with the simplest function
def phloem_exchange_surface(radius, length):
    return 2 * np.pi * radius * length * 0.1


# Testing with a complex function
def import_Nm(Nm, root_exchange_surface, C_hexose_root=1e-4):
        # We define mineral nitrogen active uptake from soil
        precision = 0.99
        Km_Nm_root = (Km_Nm_root_LATS - Km_Nm_root_HATS) / (
                1 + (precision / ((1 - precision) * np.exp(-begin_N_regulation))
                     * np.exp(-Nm / span_N_regulation))
        ) + Km_Nm_root_HATS
        # (Michaelis-Menten kinetic, surface dependency, active transport C requirements)
        return ((soil_Nm * vmax_Nm_root / (soil_Nm + Km_Nm_root)) * root_exchange_surface * (
            C_hexose_root / (C_hexose_root + transport_C_regulation)))

# Testing with an iterating funtion
def fake_complex_iterator():
    root_gen = g.component_roots_at_scale_iter(g.root, scale=1)
    root = next(root_gen)
    for vid in post_order2(g, root):
        if vid in vertices and vid != 5528:
            n = g.node(vid)
            n.import_Nm = import_Nm(n.Nm, n.root_exchange_surface, n.C_hexose_root)
            p = n.parent()
            if p:
                p.import_Nm = import_Nm(p.Nm, p.root_exchange_surface, p.C_hexose_root)
    
    for vid in pre_order2(g, root):
        if vid in vertices and vid != 5528:
            n = g.node(vid)
            n.import_Nm = import_Nm(n.Nm, n.root_exchange_surface, n.C_hexose_root)
            p = n.parent()
            if p:
                p.import_Nm = import_Nm(p.Nm, p.root_exchange_surface, p.C_hexose_root)

@njit 
def phloem_exchange_surface_numba(radius, length):
    return 2 * np.pi * radius * length * 0.1

@njit
def import_Nm_numba(Nm, root_exchange_surface, C_hexose_root=1e-4):
        # We define mineral nitrogen active uptake from soil
        precision = 0.99
        Km_Nm_root = (Km_Nm_root_LATS - Km_Nm_root_HATS) / (
                1 + (precision / ((1 - precision) * np.exp(-begin_N_regulation))
                     * np.exp(-Nm / span_N_regulation))
        ) + Km_Nm_root_HATS
        # (Michaelis-Menten kinetic, surface dependency, active transport C requirements)
        return ((soil_Nm * vmax_Nm_root / (soil_Nm + Km_Nm_root)) * root_exchange_surface * (
            C_hexose_root / (C_hexose_root + transport_C_regulation)))

class InputClass:
     def __init__(self, g):
          self.props = g.properties()
          for name in self.props:
               setattr(self, name, self.props[name]) 

class ReceiverClass:
     def __init__(self, g, input_instance):
          self.inputs = input_instance
          self.props = g.properties()
          self.vertices = list(self.props["struct_mass"].keys())
          for name in self.props:
               setattr(self.__class__, name, property(eval(f"""lambda self: dict(zip(self.vertices, [self.inputs.{name}[vid] for vid in self.vertices]))"""))) 

# STRATEGIES TO TEST
def strategy_0():
    """
    Like strategy 1 but without self iterating function to assess its impact.
    """
    for _ in range(n_processes):
        for func in parallel_functions:
            props[func.__name__].update({vid:func(*(props[name][vid] for name in input_names[func.__name__])) for vid in vertices})

def strategy_1():
    for _ in range(n_processes):
        for func in parallel_functions:
            props[func.__name__].update({vid:func(*(props[name][vid] for name in input_names[func.__name__])) for vid in vertices})
    fake_complex_iterator()

def strategy_2():
    numpy_props = {name: np.array(list(props[name].values())) for name in props if isinstance(props[name], dict)}
    for _ in range(n_processes):
        for func in parallel_functions:
            props[func.__name__] = func(*(numpy_props[name] for name in input_names[func.__name__]))
    fake_complex_iterator()
    for name in numpy_props:
            props[name] = dict(zip(vertices, numpy_props[name]))

def strategy_3():
    numpy_props = {name: np.array(list(props[name].values())) for name in props if isinstance(props[name], dict)}
    for _ in range(n_processes):
        props["phloem_exchange_surface"] = phloem_exchange_surface_numba(*(numpy_props[name] for name in input_names["phloem_exchange_surface"]))
        props["import_Nm"] = import_Nm_numba(*(numpy_props[name] for name in input_names["import_Nm"]))
    fake_complex_iterator()
    for name in numpy_props:
            props[name] = dict(zip(vertices, numpy_props[name]))

def strategy_4():
    """
    Twist of strategy 1 but with properties retreived from a class to see if the bottleneck might be here.
    """
    for _ in range(n_processes):
        for func in parallel_functions:
            props[func.__name__].update({vid:func(*(getattr(class_holding_properties, name)[vid] for name in input_names[func.__name__])) for vid in vertices})
    fake_complex_iterator()

def strategy_5():
    """
    Twist of strategy 4 but with properties retreived from an "Input class" whose values are retrieved from the other class as properties at each call.
    """
    for _ in range(n_processes):
        for func in parallel_functions:
            props[func.__name__].update({vid:func(*(getattr(receiver_instance, name)[vid] for name in input_names[func.__name__])) for vid in vertices})
    fake_complex_iterator()


if __name__ == "__main__":
    mtg_path = "inputs/root_1080.pckl"

    with open(mtg_path, "rb") as f:
        g = pickle.load(f)

    props = g.properties()
    class_holding_properties = InputClass(g=g)
    receiver_instance = ReceiverClass(g=g, input_instance=class_holding_properties)
    vertices = list(props["struct_mass"].keys())
    
    n_processes = 50
    parallel_functions = [phloem_exchange_surface, import_Nm]
    input_names = {func.__name__: ins.getfullargspec(func)[0] for func in parallel_functions}

    execution_times = dict(
        strat_0=timeit.timeit(strategy_0, number=10),
        strat_1=timeit.timeit(strategy_1, number=10),
        strat_2=timeit.timeit(strategy_2, number=10),
        strat_3=timeit.timeit(strategy_3, number=10),
        strat_4=timeit.timeit(strategy_4, number=10),
        strat_5=timeit.timeit(strategy_5, number=10)
    )
    
    print(execution_times)
