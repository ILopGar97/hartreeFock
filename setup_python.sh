#!/bin/bash

# Cargar el modulo de Python 3.10.4
module load Python/3.10.4-GCCcore-11.3.0-bare

# Verificar versión de Python
echo "Python cargado:"
python3 --version

# Instalar numpy y scipy localmente
echo "Instalando numpy y scipy..."
python3 -m pip install --user --upgrade pip
python3 -m pip install --user numpy scipy

# Verificar instalacion
echo "Verificando paquetes instalados:"
python3 -c "import numpy, scipy; print('NumPy:', numpy.__version__, '| SciPy:', scipy.__version__)"

echo "Entorno Python configurado correctamente."
