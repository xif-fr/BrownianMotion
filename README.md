## Instructions

### Deps

 * `python` 3 with `numpy`, `matplotlib`, `pandas` and `cffi` modules installed
 * Jupyter
 * SFML

### Build

 * Compile the simulation core and the CFFI module with `make brownian-motion-pysimul`

### Usage for studying diffusion

 * Open the Jupyter notebook `main.ipynb`
 * Adjust `rea_number` (typ. `0`), `target_T`, `t_beg` and `t_end`
 * Set `savepath` and create the corresponding directory
 * Run cells 0, 1, 2, 3, 5, 6, 7, 8
 * Let the cycle run and spaghettis form
 * To stop, run `simul.end()`