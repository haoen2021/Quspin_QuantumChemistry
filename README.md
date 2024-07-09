# Quspin_QuantumChemistry
A very simple python script to construct molecular Hamiltonian matrix using quspin.This project demonstrates how to construct the Hamiltonian under the fermionic basis using the QuSpin library, along with molecular integrals computed using PySCF. 

## Installation

Clone the repository:

  ```
    git clone https://github.com/haoen2021/Quspin_QuantumChemistry.git
    cd molecular_hamiltonian
  ```
Install the required Python packages:

  ```
    pip install -r requirements.txt
  ```

## Usage

To run the script, simply execute:

  ```
  python main.py --atom "H 0 0 0; Li 0 0 1.546" --basis "sto3g" --mol_dub "LiH"
  ```
By default, the Hamiltonian matrix will be saved as `Hmat_{mol_dub}.npy` in the current directory. You can specify a different filename using the `--save_file` option.

## Example Output

  ```
  converged SCF energy = -7.86313368869442
  GS energy from the exact diagonalization of the matrix we constructed: -7.882761848745534
  FCI from pyscf = -7.882761848745533
  Hamiltonian matrix saved to Hmat_LiH.npy
  ```
