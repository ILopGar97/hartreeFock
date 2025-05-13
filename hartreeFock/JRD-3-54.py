if __name__ == "__main__":
    """
    Para que la ejecución sea correcta este archivo debe estar en el mismo directorio que los siguientes:
        - COEFS (la carpeta de coeficientes, no es necesario modificar nada)  
        - file_manager.py (encargado de leer coeficientes y exportar resultados)
        - wave_function_constructor.py (encargado de construir las funciones de onda, densidades, etc en ambos espacios, con armonicos esfericos o sin ellos etc)
        - information_measures.py (encargado de calcular las medidas de información (de los sitemas atomicos)que necesitemos)
        - data_structures.py (es un fichero que contiene la clase AtomicData y Orbitals)


    """
    import time 
    from concurrent.futures import ThreadPoolExecutor, as_completed
    from file_manager import load_system, select_atoms_by_Z, json_to_multiple_txt
    from information_measures import *
    import json
    BASE_PATH = "COEFS"

    start_time = time.time()

    neutral_atom = load_system(0, BASE_PATH) #<-------- Aquí se puede cambiar los sitemas que queremos estudiar si queremos estudiar atomos neutros 0, -1 aniones, +1 cationes (con eso hacemos referencia a la carga). Si introducimos un número entre 2 y 9 iremos a las series isoelectronicas (estamos indicando el número de electrones).
    atoms_selected = select_atoms_by_Z(neutral_atom, range(3, 55)) # <------ Indicamos (con Z) el rango de átomos/iones/cationes... que queremos estudiar (en este caso está indicando atomos con Z entre 3 y 54, dado que range genera una lista que es [3, 4, ..., n-1]. Si queremos estudiar todos los iones por ejemplo podemos poner "all" (con las comillas incluidas) i.e select_atoms_by_Z(neutral_atom, "all").

    def compute_entropy(name, cls, kwargs): #no tocar
        print(f"{name} iniciada en t+{time.time() - start_time:.3f} s", flush = True)
        result = cls(**kwargs)
        return name, result

    alpha_values = [0.5] # <---------- Aquí incluimos el/los valor/es de alpha/q que queremos estudiar para JRDq

    tasks = []
    for a in alpha_values:
        for atom in atoms_selected:
            tasks += [ #<----- Aqui estamos haciendo una lista con las tareas que vamos a realizar para poder paralelizarlo, en este caso el separamos el calculo de JRD en espacios y en valores de a y atomos
                (f"JRD posicion - alpha={a:.3f}", JenRenyiDivergence, {"atoms": atom, "alpha": a, "space": "position"}), #<---- saldrá un warning indicando que por defecto se ha cogido el espacio indicado (ignorar)
                (f"JRD momento - alpha={a:.3f}", JenRenyiDivergence, {"atoms": atom, "alpha": a, "space": "momentum"})
            ]

    resultsJRD = {} #<---------Aquí se iran guardando los datos obtenidos y obtendremos al final un formato json

    with ThreadPoolExecutor(max_workers=10) as executor: #no tocar <------- el max_workers=10 dice que como MAXIMO cogera 10 CPUS para hacer el calculo paralelo. Hay que tener cuidado de indicar -c 10 en proteus, sino no hará nada en paralelo...
        future_to_name = { 
            executor.submit(compute_entropy, name, cls, kwargs): name
            for name, cls, kwargs in tasks
        }

        for future in as_completed(future_to_name):
            name = future_to_name[future]
            task_name, result = future.result()
            print(f"{task_name} completada en t+{time.time() - start_time:.3f} s", flush = True)
            if task_name in resultsJRD:
                resultsJRD[task_name].append(result[0])
            else:
                resultsJRD[task_name] = result


    end_time = time.time()
    print(f"\nTiempo total de cómputo: {end_time - start_time:.2f} segundos", flush = True) 
    print(resultsJRD, flush = True) #<---------- por seguridad se imprimira el resultado tambien por pantalla, en el caso de proteus irá al archivo .out 

    with open("resultsJTD.json", "w") as json_file: #<------------- guardamos el json
        json.dump(resultsJRD, json_file, indent=4)

    json_to_multiple_txt("resultsJTD.json")  #<--------- lo ponemos bonito y ordenados en tablas (se genera una carpeta nueva, pueden salir desordenadas, se arreglará pronto...)
