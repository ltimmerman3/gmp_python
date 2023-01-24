from ase.io import read
from gmp import *
import os
from tqdm import tqdm

sigmas = np.logspace(-1,0.6,3).round(2)
max_order = 4
rs_scale = 1.0
cutoff = 20
images = read('../water_2d/water_2d.traj',index=':')
pbc_bools = [1,1,1]
# first entry is max msch order
params_i = [[i,0,1] for i in range(max_order+1) for x in range(len(sigmas))]
params_d = [[sigma, 1.0, (1.0 / (sigma * np.sqrt(2.0 * np.pi))) ** 3, 1.0 / (2.0 * sigma * sigma), cutoff, (1.0 / (rs_scale * sigma))] for x in range(max_order+1) for sigma in sigmas]
nmcsh = len(sigmas)*(max_order+1)
O_psd = np.genfromtxt('O_pseudodensity_4.g')
H_psd = np.genfromtxt('H_pseudodensity_2.g')
atom_gaussians = np.zeros((2,8),dtype=float)
ngaussians = [4,2]
element_index_to_order = np.zeros(120,dtype=int)
element_index_to_order[8] = 0
element_index_to_order[1] = 1
atom_gaussians[0,::2] = O_psd[:,0]
atom_gaussians[0,1::2] = O_psd[:,1]
atom_gaussians[1,::2] = np.concatenate((H_psd[:,0],np.zeros(2)))
atom_gaussians[1,1::2] = np.concatenate((H_psd[:,1],np.zeros(2)))
# print('Input params: cell {}\ncart {}\nscale {}\npbc_bools {}\natom_i {}\nnatoms {}\ncal_atoms {}\ncal_num {}\nparams_i {}\nparams_d {}\nnmcsh {}\natom_gaussians {}\nngaussians {}\nelement_index_to_order {}'.format(cell, cart, scale, pbc_bools,
#                                         atom_i, natoms, cal_atoms, cal_num,
#                                         params_i, params_d, nmcsh, atom_gaussians, ngaussians, element_index_to_order))
for image in tqdm(images, desc="Images"):
    cell = np.copy(image.cell)
    cart = image.get_positions(wrap=True)
    scale = image.get_scaled_positions()
    atom_i = np.array(image.get_atomic_numbers())
    natoms = len(image)
    cal_atoms = np.linspace(0,len(image)-1,len(image),dtype=int)
    cal_num = len(image)
    mcsh, dmcsh = calculate_solid_gmpordernorm(cell, cart, scale, pbc_bools,
                                        atom_i, natoms, cal_atoms, cal_num,
                                        params_i, params_d, nmcsh, atom_gaussians, ngaussians, element_index_to_order)
    for idx,i in enumerate(['O','H','H']):
        if not os.path.exists('./all_py_water.txt'):
            with open('all_py_water.txt','w') as fp:
                print(':{}:\n{}'.format(i,mcsh[idx,:]),file=fp)
        else:
            with open('all_py_water.txt','a') as fp:
                print(':{}:\n{}'.format(i,mcsh[idx,:]),file=fp)
