import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from file_manager import load_system, select_atoms_by_Z, json_to_multiple_txt
from information_measures import *
import json

def compute_entropy(name, cls, kwargs, start_time):
        print(f"{name} iniciada en t+{time.time() - start_time:.3f} s", flush = True)
        result = cls(**kwargs)
        return name, result

if __name__ == "__main__":
    BASE_PATH = "COEFS"

    start_time = time.time()

    neutral_atom = load_system(1, BASE_PATH)
    atoms_selected = select_atoms_by_Z(neutral_atom, "all")
    atoms_selected = atoms_selected[3:]
   

    alpha_values = [num / 4 for num in range(2, 13)]

    tasks = []
    for a in alpha_values:
        tasks += [
            (f"JTD - posicion alpha={a:.3f}-cat", JenTsallisDivergence, {"atoms": atoms_selected, "alpha": a, "space": "position"}),
            (f"JTD - momento alpha={a:.3f}-cat", JenTsallisDivergence, {"atoms": atoms_selected, "alpha": a, "space": "momentum"})
        ]

    resultsJTD = {}

    with ProcessPoolExecutor(max_workers=15) as executor:
        future_to_name = {
            executor.submit(compute_entropy, name, cls, kwargs, start_time): name
            for name, cls, kwargs in tasks
        }

        for future in as_completed(future_to_name):
            name = future_to_name[future]
            task_name, result = future.result()
            print(f"{task_name} completada en t+{time.time() - start_time:.3f} s", flush = True)
            resultsJTD[task_name] = result

    end_time = time.time()
    print(f"\nTiempo total de cómputo: {end_time - start_time:.2f} segundos", flush = True)
    print(resultsJTD, flush = True)

    with open("resultsJTD-cat.json", "w") as json_file:
        json.dump(resultsJTD, json_file, indent=4)

    json_to_multiple_txt("resultsJTD-cat.json")
