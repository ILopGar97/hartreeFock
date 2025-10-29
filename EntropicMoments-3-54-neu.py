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
    """
    Para que la ejecución sea correcta este archivo debe estar en el mismo directorio que los siguientes:
        - COEFS (la carpeta de coeficientes, no es necesario modificar nada)  
        - file_manager.py (encargado de leer coeficientes y exportar resultados)
        - wave_function_constructor.py (encargado de construir las funciones de onda, densidades, etc en ambos espacios, con armonicos esfericos o sin ellos etc)
        - information_measures.py (encargado de calcular las medidas de información (de los sitemas atomicos)que necesitemos)
        - data_structures.py (es un fichero que contiene la clase AtomicData y Orbitals)

    CUIDADO: No se si probarlo en tu local es buena idea, tira de mucha CPU, puede que se te congele el PC (no me ha pasado pero el mío tiene bastante). Si quieres probarlo te recomiendo tocar las variables rg (reduces cuantos atomos quieres), y el max_workers a 3-4 no más. 
    """

    BASE_PATH = "COEFS"

    start_time = time.time()

    Q = 0 #<---- aqui elige que quieres estudiar neutros = 0, aniones = -1, cationes = 1. No te salgas de esos valores, si no dará error o calculará isoelectronicas. 
    if Q == 0:
        ext = "neu" #<------- es la extension que añadiremos para diferenciar los diferentes archivos se selecciona automaticamente al darle un valor a Q
    elif Q==1:
        ext = "cat"
    else:
        ext = "ani"
    q_values = [num / 4 for num in range(2, 13)] # <---------- Aquí incluimos el/los valor/es de alpha/q que queremos estudiar para EntropicMomentq (si vas a incluir valores enteros como el 1, pon el punto i.e 1.0). También puedes definirlo con [num/den for num in range(ini, final-1)]
    if Q == 0: 
        rg = range(3, 55)
    else:
        rg = "all"
    atoms_list = load_system(Q, BASE_PATH) #<-------- Aquí se puede cambiar los sitemas que queremos estudiar si queremos estudiar atomos neutros 0, -1 aniones, +1 cationes (con eso hacemos referencia a la carga). Si introducimos un número entre 2 y 9 iremos a las series isoelectronicas (estamos indicando el número de electrones).
    atoms_selected = select_atoms_by_Z(atoms_list, rg) # <------ Indicamos (con Z) el rango de átomos/iones/cationes... que queremos estudiar (en este caso está indicando atomos con Z entre 3 y 54, dado que range genera una lista que es [3, 4, ..., n-1]. Si queremos estudiar todos los iones por ejemplo podemos poner "all" (con las comillas incluidas) i.e select_atoms_by_Z(neutral_atom, "all").
    if Q != 0:
        atoms_selected = atoms_selected[1:] #<---- para los aniones descartamos el hidrogeno, el helio y el litio

    tasks = []
    for q in q_values:
        for atom in atoms_selected:
            tasks += [ #<----- Aqui estamos haciendo una lista con las tareas que vamos a realizar para poder paralelizarlo, en este caso el separamos el calculo de EntropicMoment en espacios y en valores de a y atomos
                (f"EntropicMoment posicion - q={q:.3f} -"+ext, EntropicMoment, {"atoms": atom, "qs": q, "space": "position", "two_electron_density": True, "product_density":True, "average_density":True}), #<---- saldrá un warning indicando que por defecto se ha cogido el espacio indicado (ignorar)
                (f"EntropicMoment momento - q={q:.3f} -"+ext, EntropicMoment, {"atoms": atom, "qs": q, "space": "momentum", "two_electron_density": True, "product_density":True, "average_density":True})
            ]

    resultsEntropicMoment = {} #<---------Aquí se iran guardando los datos obtenidos y obtendremos al final un formato json

    with ProcessPoolExecutor(max_workers=15) as executor: #no tocar <------- el max_workers=15 dice que como MAXIMO cogera 15 CPUS para hacer el calculo paralelo. Hay que tener cuidado de indicar -c 15 en proteus, sino no hará nada en paralelo... Y por supuesto no cojas más de los que indiques aquí, no debería pasar nada... hasta que pasa, mejor no lo hagas. Con 15 va a ir muy rápido
        future_to_name = { 
            executor.submit(compute_entropy, name, cls, kwargs, start_time): name
            for name, cls, kwargs in tasks
        }

        for future in as_completed(future_to_name):
            name = future_to_name[future]
            task_name, result = future.result()
            print(f"{task_name} completada en t+{time.time() - start_time:.3f} s", flush = True)
            if task_name in resultsEntropicMoment:
                resultsEntropicMoment[task_name].append(result[0])
            else:
                resultsEntropicMoment[task_name] = result


    end_time = time.time()
    print(f"\nTiempo total de cómputo: {end_time - start_time:.2f} segundos", flush = True) 
    print(resultsEntropicMoment, flush = True) #<---------- por seguridad se imprimira el resultado tambien por pantalla, en el caso de proteus irá al archivo .out 

    with open("resultsEntropicMoment-"+ext+".json", "w") as json_file: #<------------- guardamos el json
        json.dump(resultsEntropicMoment, json_file, indent=4)

    json_to_multiple_txt("resultsEntropicMoment-"+ext+".json")  #<--------- lo ponemos bonito y ordenados en tablas (se genera una carpeta nueva, pueden salir desordenadas, se arreglará pronto... Mientras tengamos el json todo ok)
