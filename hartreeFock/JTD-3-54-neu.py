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

    neutral_atom = load_system(0, BASE_PATH)
    atoms_selected = select_atoms_by_Z(neutral_atom, range(3, 55))

   

    alpha_values = [num / 4 for num in [2, 4]]

    tasks = []
    for a in alpha_values:
        tasks += [
            (f"JTD - posicion alpha={a:.3f}-neu", JenTsallisDivergence, {"atoms": atoms_selected[0:3], "alpha": a, "space": "position"}),
            (f"JTD - momento alpha={a:.3f}-neu", JenTsallisDivergence, {"atoms": atoms_selected[0:3], "alpha": a, "space": "momentum"})
        ]

    resultsJTD = {}

    with ProcessPoolExecutor(max_workers=2) as executor:
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

    with open("resultsJTD-neu.json", "w") as json_file:
        json.dump(resultsJTD, json_file, indent=4)

    json_to_multiple_txt("resultsJTD-neu.json")
