## Instructions

### Deps

 * `python` 3 with `numpy`, `scipy`, `matplotlib`, `pandas` and `cffi` modules installed
 * Jupyter
 * SFML

### Lennard-Jones gas simulation

 * Compile the simulation core and the CFFI module with `make brownian-motion-pysimul` (with SFML) or `make brownian-motion-pysimul-headless` (without display)
 * Run simulation with the Jupyter notebook `main.ipynb`

### Langevin simulations for computing survival probabilities and MFPT

 * Compile the simulation core and the CFFI module with `make langevin-survival-pysimul` after adjusting the various flags in `langevin-survival.cpp`
 * For computing the survival probabilities, run the simulation with one of the Jupyter notebook in `langevin-survival/`
 * For computing the MFPT, run the simulation with `langevin-mfpt/langevin-ft-automated.ipynb`
