# monte_carlo_diagonalization
This package uses a Monte Carlo sampling of the Fock states and computes the ground state and the Green function for this truncated space, using the Hubbard model.


# Table of contents

- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Authors](#authors)
- [License](#license)


# Dependencies
python3, openmp, lapack-blas are required. 
### Linus (Ubuntu)
```shell
sudo apt update
sudo apt install git liblapack-dev libblas-dev libopenmpt-dev python3 
sudo apt install lablapacke-dev
```
The python3 packages numpy, matplotlib and scipy are needed.
```shell
pip3 install numpy matplotlib scipy
```

# Installation
1. Clone this repository
```shell
git clone https://github.com/ben65314/monte_carlo_diagonalization.git
```
2. Make the .exe 
```shell
cd monte-carlo-diagonalization
make
```
3. You may want to put the `build` directory in your `$PATH`

# Usage
The `mcd_solver` and `subspace_finder` are the two main executables, while `qMatrixGreenCompute.py`, `graphGreen.py` and `av_combine.py` are used to plot Green functions.
1. `mcd_solver` samples states, finds the ground state and computes the Green functions of a given $N_e$- $S_z$ bloc. Produces a `qMatrices.txt` file countaining the Q-matrices needed to plot the Green functions. Needs to be given a parameter file :
```shell
mcd_solver parameters.txt
```
2. `subspace_finder` computes the ground state searching through all $N_e$- $S_z$ sectors with ED. Needs to be given a parameter file : 
```shell
mcd_solver parameters.txt
```
The parameters needed in the parameter file are presented in `/examples/PARAMETERS_DOC.md`

3. `qMatrixGreenCompute.py` computes the Green function with the given `qMatrices.txt` file. Produces a `greenFunctionValueQM.txt`. There are multiple way to call it,
- Will search in the `pwd` for a `qMatrices.txt` file : 
```shell
qMatrixGreenCompute.py
```
- Will use the given file (`qMatrices.txt` format) :
```shell
qMatrixGreenCompute.py qMatrices.txt
```
- Will use the given file (`qMatrices.txt` format) and use the given plot limits :
```shell
qMatrixGreenCompute.py qMatrices.txt -10 10
```

4. `graphGreen.py` uses matplotlib to plot the computed Green function : 
```shell
graphGreen.py greenFunctionValueQM.txt
```

5. `av_combine.py` combines the multiple Green functions into one (DOS). The produced file can be plot with `graphGreen.py`: 
```shell
av_combine.py greenFunctionValueQM.txt
```
# Examples
In `/examples/` there are preset file parameters ready to be execute. There is documentation about the parameters in `/examples/PARAMETERS_DOC.md` so you can make your own systems.
# Authors
This package was made and published by **Benjamin Bernard** and **Maxime Charelbois**. If you have any question about the code, contact us.
- <benjamin.bernard@uqtr.ca>
- <maxime.charlebois@uqtr.ca>

# License
MIT

