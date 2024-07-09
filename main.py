"""
trying to construct the Hamiltonian under the fermionic basis using quspin
"""
import argparse
from quspin.operators import hamiltonian # Hamiltonians and operators
from quspin.basis import spinful_fermion_basis_1d, spinless_fermion_basis_1d
import numpy as np 
 
from pyscf import gto, scf, fci, ao2mo
import sys
import os
from utils.spatial2spin import spinorb_from_spatial

def main(atom, basis, mol_dub, save_file):
    # define molecular Hamiltonians
    mol = gto.M(
        atom = atom,
        basis = basis,
        #verbose = 4,
    )
    dub = mol_dub


    # read 1- and 2- electron integrals
    e1 = mol.intor('int1e_nuc')+mol.intor('int1e_kin')
    mf = scf.RHF(mol)
    mf.kernel()
    mo1e = mf.mo_coeff.T @ e1 @ mf.mo_coeff
    mo2e = ao2mo.get_mo_eri(mol, mf.mo_coeff, aosym = 1)
    mo2e = mo2e.reshape(mol.nao, mol.nao, mol.nao, mol.nao)
    mo2e = np.transpose(mo2e, (0,2,3,1)) # need to be checked later!!!

    mo1e_spin, mo2e_spin = spinorb_from_spatial(mo1e, mo2e)

    # construct the fermionic operators using the data structure provided by quspin
    ferm_SingExcit, ferm_DoubExcit = [], []
    for i in range(mo1e_spin.shape[0]):
        for j in range(mo1e_spin.shape[0]):
                if abs(mo1e_spin[i][j]) > 1e-15:
                    ferm_SingExcit.append([mo1e_spin[i][j], i, j]) # a_i^+ a_j, '+-'

    for i in range(mo2e_spin.shape[0]):
        for j in range(mo2e_spin.shape[0]):
            for k in range(mo2e_spin.shape[0]):
                for l in range(mo2e_spin.shape[0]):
                    if abs(mo2e_spin[i][j][k][l]) > 1e-15:
                        ferm_DoubExcit.append([0.5*mo2e_spin[i][j][k][l], i, j, k, l]) # a_i^+ a_j^+ a_k a_l, '++--'

    # spinless fermionic basis, since we have already converted the mo1e and mo2e tensors into spinful ones
    fermion_basis_mol = spinless_fermion_basis_1d(mo1e_spin.shape[0], Nf=mol.nelectron)

    # hamiltonian constructed using quspin
    no_checks = dict(check_pcon=False,check_symm=False,check_herm=False)
    H_spinf = hamiltonian([['+-', ferm_SingExcit]]+[['++--', ferm_DoubExcit]],[], basis=fermion_basis_mol, dtype=np.float64,**no_checks)
    H_matf = np.array(H_spinf.todense())
    H_matf += mol.energy_nuc()*np.eye(H_matf.shape[0])
    e, evc = np.linalg.eigh(H_matf)
    print(f"GS energy from the exact diagonalization of the matrix we constructed: {e[0]}")
    print("FCI from pyscf = {}".format(fci.FCI(mf).kernel()[0]))

    # Save the Hamiltonian matrix to a .npy file
    np.save(save_file, H_matf)
    print(f"Hamiltonian matrix saved to {save_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Construct molecular Hamiltonian using QuSpin and PySCF.')
    parser.add_argument('--atom', type=str, required=True, help='Atomic configuration, e.g., "H 0 0 0; Li 0 0 1.546"')
    parser.add_argument('--basis', type=str, required=True, help='Basis set, e.g., "sto3g"')
    parser.add_argument('--mol_dub', type=str, required=True, help='Molecule description, e.g., "LiH"')

    # Default save file name
    args = parser.parse_args()
    default_save_file = f"Hmat_{args.mol_dub}.npy"
    parser.add_argument('--save_file', type=str, default=default_save_file, help='Filename to save the Hamiltonian matrix, e.g., "hamiltonian.npy"')
    args = parser.parse_args()

    main(args.atom, args.basis, args.mol_dub, args.save_file)