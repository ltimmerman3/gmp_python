import numpy as np
from math import ceil
from mcsh import *

M_PI = 3.14159265358979323846

NUM_IMPLEMENTED_TYPE = 5
IMPLEMENTED_MCSH_TYPE = [0,1,2,3,4]

def calculate_solid_gmpordernorm(cell, cart, scale, pbc_bools,
                                        atom_i, natoms, cal_atoms, cal_num,
                                        params_i, params_d, nmcsh, atom_gaussian, ngaussians, element_index_to_order):
    mcsh = np.zeros((cal_num,nmcsh),dtype=float)
    dmcsh = np.zeros((cal_num*nmcsh,natoms*3),dtype=float)
    bin_range = np.zeros(3,dtype=int)
    nbins = np.zeros(3,dtype=int)
    cell_shift = np.zeros(3,dtype=int)
    max_bin = np.zeros(3,dtype=int)
    min_bin = np.zeros(3,dtype=int)
    pbc_bin = np.zeros(3,dtype=int)
    plane_d = np.zeros(3,dtype=float)
    total_shift = np.zeros(3,dtype=float)
    cross = np.zeros((3,3),dtype=float)
    reci = np.zeros((3,3),dtype=float)
    inv = np.zeros((3,3),dtype=float)

    # Check for not implemented mcsh type.
    for m in range(nmcsh):
        implemented = False
        for i in range(NUM_IMPLEMENTED_TYPE):
            # I think there's an issue here - doesn't actually check what it should check - removing for simplification
            if (params_i[m][0] == IMPLEMENTED_MCSH_TYPE[i]):
                implemented = True
                break
        if (not implemented): 
            return 1
    
    bin_i = np.zeros((natoms,4),dtype=int)

    cutoff = 0.0
    #let cutoff equal to the maximum of Rc
    for m in range(nmcsh):
        if (cutoff < params_d[m][4]):
            cutoff = params_d[m][4]
    cutoff_sqr = cutoff * cutoff

    total_bins = 1

    #calculate the inverse matrix of cell and the distance between cell plane
    cross[0,0] = cell[1,1]*cell[2,2] - cell[1,2]*cell[2,1]
    cross[0,1] = cell[1,2]*cell[2,0] - cell[1,0]*cell[2,2]
    cross[0,2] = cell[1,0]*cell[2,1] - cell[1,1]*cell[2,0]
    cross[1,0] = cell[2,1]*cell[0,2] - cell[2,2]*cell[0,1]
    cross[1,1] = cell[2,2]*cell[0,0] - cell[2,0]*cell[0,2]
    cross[1,2] = cell[2,0]*cell[0,1] - cell[2,1]*cell[0,0]
    cross[2,0] = cell[0,1]*cell[1,2] - cell[0,2]*cell[1,1]
    cross[2,1] = cell[0,2]*cell[1,0] - cell[0,0]*cell[1,2]
    cross[2,2] = cell[0,0]*cell[1,1] - cell[0,1]*cell[1,0]

    vol = cross[0,0]*cell[0,0] + cross[0,1]*cell[0,1] + cross[0,2]*cell[0,2]

    inv[0,0] = cross[0,0]/vol
    inv[0,1] = cross[1,0]/vol
    inv[0,2] = cross[2,0]/vol
    inv[1,0] = cross[0,1]/vol
    inv[1,1] = cross[1,1]/vol
    inv[1,2] = cross[2,1]/vol
    inv[2,0] = cross[0,2]/vol
    inv[2,1] = cross[1,2]/vol
    inv[2,2] = cross[2,2]/vol

    # bin: number of repetitive cells?
    for i in range(3):
        tmp = 0
        for j in range(3):
            reci[i,j] = cross[i,j]/vol
            tmp += reci[i,j]*reci[i,j]
        plane_d[i] = 1/np.sqrt(tmp)
        nbins[i] = ceil(plane_d[i]/cutoff)
        total_bins *= nbins[i]
    atoms_bin = np.zeros(total_bins,dtype=int)

    # assign the bin index to each atom
    for i in range(natoms):
        for j in range(3):
            bin_i[i,j] = scale[i,j] * float(nbins[j])
        bin_i[i,3] = bin_i[i,0] + nbins[0]*bin_i[i,1] + nbins[0]*nbins[1]*bin_i[i,2]
        atoms_bin[bin_i[i,3]] += 1
    max_atoms_bin = 0
    for i in range(total_bins):
        if (atoms_bin[i] > max_atoms_bin):
            max_atoms_bin = atoms_bin[i]
    # number of bins in each direction
    neigh_check_bins = 1
    for i in range(3):
        bin_range[i] = ceil(cutoff * nbins[i] / plane_d[i])
        neigh_check_bins *= 2*bin_range[i]
    for ii in range(cal_num):
        i=cal_atoms[ii]
        # calculate neighbor atoms
        nei_list_d = np.zeros(max_atoms_bin * 4 * neigh_check_bins,dtype=float)
        nei_list_i = np.zeros(max_atoms_bin * 2 * neigh_check_bins, dtype=int)
        nneigh = 0

        for j in range(3):
            max_bin[j] = bin_i[i,j] + bin_range[j]
            min_bin[j] = bin_i[i,j] - bin_range[j]
        dx = min_bin[0]
        for dummy_x in range(max_bin[0]-min_bin[0]+1):
            dy = min_bin[1]
            for dummy_j in range(max_bin[1]-min_bin[1]+1):
                dz = min_bin[2]
                for dummy_k in range(max_bin[2]-min_bin[2]+1):
                    pbc_bin[0] = (dx%nbins[0] + nbins[0]) % nbins[0]
                    pbc_bin[1] = (dy%nbins[1] + nbins[1]) % nbins[1]
                    pbc_bin[2] = (dz%nbins[2] + nbins[2]) % nbins[2]
                    cell_shift[0] = (dx-pbc_bin[0]) / nbins[0]
                    cell_shift[1] = (dy-pbc_bin[1]) / nbins[1]
                    cell_shift[2] = (dz-pbc_bin[2]) / nbins[2]
                    bin_num = pbc_bin[0] + nbins[0]*pbc_bin[1] + nbins[0]*nbins[1]*pbc_bin[2]

                    for j in range(natoms):
                        if (bin_i[j,3] != bin_num):
                            continue

                        # take care of pbc
                        if (not pbc_bools[0] and cell_shift[0] != 0):
                            continue

                        if (not pbc_bools[1] and cell_shift[1] != 0):
                            continue

                        if (not pbc_bools[2] and cell_shift[2] != 0):
                            continue

                        for a in range(3):
                            total_shift[a] = cell_shift[0]*cell[0,a] + cell_shift[1]*cell[1,a] + cell_shift[2]*cell[2,a] + cart[j,a] - cart[i,a]

                        # tmp = sqrt(total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2])
                        tmp_r2 = total_shift[0]*total_shift[0] + total_shift[1]*total_shift[1] + total_shift[2]*total_shift[2]
                        # if (tmp < cutoff) 
                        if (tmp_r2 < cutoff_sqr):
                            for a in range(3):
                                nei_list_d[nneigh*4 + a] = total_shift[a]
                            nei_list_d[nneigh*4 + 3] = tmp_r2
                            nei_list_i[nneigh*2]    = atom_i[j]
                            nei_list_i[nneigh*2 + 1] = j
                            nneigh += 1
                    dz += 1
                dy += 1
            dx += 1
        for m in range(nmcsh):
            mcsh_order = params_i[m][0]
            square = params_i[m][1]
            num_groups = get_num_groups(mcsh_order)
            # params_d: sigma, weight, A, alpha, cutoff, inv_rs
            A = params_d[m][2]
            alpha = params_d[m][3]
            #weight = 1.0
            sum_square = 0.0
            group_index = 1
            for dummy_gi in range(num_groups):
                mcsh_function = get_solid_mcsh_function(mcsh_order, group_index)
                group_coefficient = get_group_coefficients(mcsh_order, group_index)
                mcsh_type = get_mcsh_type(mcsh_order, group_index)
                if (mcsh_type == 1):
                    sum_miu = 0.0
                    sum_dmiu_dxj = np.zeros(nneigh,dtype=float)
                    sum_dmiu_dyj = np.zeros(nneigh,dtype=float)
                    sum_dmiu_dzj = np.zeros(nneigh,dtype=float)

                    m_desc = np.zeros(1,dtype=float)
                    deriv = np.zeros(3,dtype=float)

                    for j in range(nneigh):

                        neigh_atom_element_index = nei_list_i[j*2]
                        neigh_atom_element_order = element_index_to_order[neigh_atom_element_index]
                        x0 = nei_list_d[j*4]
                        y0 = nei_list_d[j*4+1]
                        z0 = nei_list_d[j*4+2]
                        r0_sqr = nei_list_d[j*4+3]
                        for g in range(ngaussians[neigh_atom_element_order]):
                            B = atom_gaussian[neigh_atom_element_order,g*2]
                            beta = atom_gaussian[neigh_atom_element_order,g*2+1]
                            m_desc, deriv = mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta)
                            sum_miu += m_desc
                            sum_dmiu_dxj[j] += deriv[0]
                            sum_dmiu_dyj[j] += deriv[1]
                            sum_dmiu_dzj[j] += deriv[2]
                    sum_square += group_coefficient * sum_miu * sum_miu

                    dmdx = np.zeros(1,dtype=float)
                    dmdy = np.zeros(1,dtype=float)
                    dmdz = np.zeros(1,dtype=float)

                    for j in range(nneigh):
                        dmdx = (sum_miu * sum_dmiu_dxj[j]) * group_coefficient * 2.0
                        dmdy = (sum_miu * sum_dmiu_dyj[j]) * group_coefficient * 2.0
                        dmdz = (sum_miu * sum_dmiu_dzj[j]) * group_coefficient * 2.0

                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3] += dmdx
                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3 + 1] += dmdy
                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3 + 2] += dmdz

                        dmcsh[ii*nmcsh + m,i*3]     -= dmdx
                        dmcsh[ii*nmcsh + m,i*3 + 1] -= dmdy
                        dmcsh[ii*nmcsh + m,i*3 + 2] -= dmdz

                if (mcsh_type == 2):
                    sum_miu1 = 0.0
                    sum_miu2 = 0.0
                    sum_miu3 = 0.0

                    sum_dmiu1_dxj = np.zeros(nneigh,dtype=float)
                    sum_dmiu2_dxj = np.zeros(nneigh,dtype=float)
                    sum_dmiu3_dxj = np.zeros(nneigh,dtype=float)
                    sum_dmiu1_dyj = np.zeros(nneigh,dtype=float)
                    sum_dmiu2_dyj = np.zeros(nneigh,dtype=float)
                    sum_dmiu3_dyj = np.zeros(nneigh,dtype=float)
                    sum_dmiu1_dzj = np.zeros(nneigh,dtype=float)
                    sum_dmiu2_dzj = np.zeros(nneigh,dtype=float)
                    sum_dmiu3_dzj = np.zeros(nneigh,dtype=float)

                    miu = np.zeros(3,dtype=float)
                    deriv = np.zeros(9,dtype=float)
                    for j in range(nneigh):
                        neigh_atom_element_index = nei_list_i[j*2]
                        neigh_atom_element_order = element_index_to_order[neigh_atom_element_index]
                        x0 = nei_list_d[j*4]
                        y0 = nei_list_d[j*4+1]
                        z0 = nei_list_d[j*4+2]
                        r0_sqr = nei_list_d[j*4+3]
                        for g in range(ngaussians[neigh_atom_element_order]):
                            B = atom_gaussian[neigh_atom_element_order,g*2]
                            beta = atom_gaussian[neigh_atom_element_order,g*2+1]
                            miu, deriv = mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta)
                            # miu: miu_1, miu_2, miu_3
                            # deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]
                            sum_miu2 += miu[1]
                            sum_miu3 += miu[2]
                            sum_dmiu1_dxj[j] += deriv[0]
                            sum_dmiu1_dyj[j] += deriv[1]
                            sum_dmiu1_dzj[j] += deriv[2]
                            sum_dmiu2_dxj[j] += deriv[3]
                            sum_dmiu2_dyj[j] += deriv[4]
                            sum_dmiu2_dzj[j] += deriv[5]
                            sum_dmiu3_dxj[j] += deriv[6]
                            sum_dmiu3_dyj[j] += deriv[7]
                            sum_dmiu3_dzj[j] += deriv[8]
                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3)

                    dmdx = dmdy = dmdz = 0.0
                    for j in range(nneigh):
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] + sum_miu3 * sum_dmiu3_dxj[j]) * group_coefficient * 2.0
                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] + sum_miu3 * sum_dmiu3_dyj[j]) * group_coefficient * 2.0
                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] + sum_miu3 * sum_dmiu3_dzj[j]) * group_coefficient * 2.0

                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3] += dmdx
                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3 + 1] += dmdy
                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3 + 2] += dmdz

                        dmcsh[ii*nmcsh + m,i*3]     -= dmdx
                        dmcsh[ii*nmcsh + m,i*3 + 1] -= dmdy
                        dmcsh[ii*nmcsh + m,i*3 + 2] -= dmdz

                if (mcsh_type == 3):
                    sum_miu1 = sum_miu2 = sum_miu3 = sum_miu4 = sum_miu5 = sum_miu6 = 0.0

                    sum_dmiu1_dxj = np.zeros(nneigh,dtype=float)
                    sum_dmiu2_dxj = np.zeros(nneigh,dtype=float)
                    sum_dmiu3_dxj = np.zeros(nneigh,dtype=float)
                    sum_dmiu4_dxj = np.zeros(nneigh,dtype=float)
                    sum_dmiu5_dxj = np.zeros(nneigh,dtype=float)
                    sum_dmiu6_dxj = np.zeros(nneigh,dtype=float)
                    sum_dmiu1_dyj = np.zeros(nneigh,dtype=float)
                    sum_dmiu2_dyj = np.zeros(nneigh,dtype=float)
                    sum_dmiu3_dyj = np.zeros(nneigh,dtype=float)
                    sum_dmiu4_dyj = np.zeros(nneigh,dtype=float)
                    sum_dmiu5_dyj = np.zeros(nneigh,dtype=float)
                    sum_dmiu6_dyj = np.zeros(nneigh,dtype=float)
                    sum_dmiu1_dzj = np.zeros(nneigh,dtype=float)
                    sum_dmiu2_dzj = np.zeros(nneigh,dtype=float)
                    sum_dmiu3_dzj = np.zeros(nneigh,dtype=float)
                    sum_dmiu4_dzj = np.zeros(nneigh,dtype=float)
                    sum_dmiu5_dzj = np.zeros(nneigh,dtype=float)
                    sum_dmiu6_dzj = np.zeros(nneigh,dtype=float)

                    miu = np.zeros(6,dtype=float)
                    deriv = np.zeros(18,dtype=float)
                    for j in range(nneigh):
                        neigh_atom_element_index = nei_list_i[j*2]
                        neigh_atom_element_order = element_index_to_order[neigh_atom_element_index]
                        x0 = nei_list_d[j*4]
                        y0 = nei_list_d[j*4+1]
                        z0 = nei_list_d[j*4+2]
                        r0_sqr = nei_list_d[j*4+3]
                        for g in range(ngaussians[neigh_atom_element_order]):
                            B = atom_gaussian[neigh_atom_element_order,g*2]
                            beta = atom_gaussian[neigh_atom_element_order,g*2+1]
                            miu, deriv = mcsh_function(x0, y0, z0, r0_sqr, A, B, alpha, beta)
                            # miu: miu_1, miu_2, miu_3
                            # deriv: dmiu1_dxj, dmiu1_dyj, dmiu1_dzj, dmiu2_dxj, dmiu2_dyj, dmiu2_dzj, dmiu3_dxj, dmiu3_dyj, dmiu3_dzj
                            sum_miu1 += miu[0]
                            sum_miu2 += miu[1]
                            sum_miu3 += miu[2]
                            sum_miu4 += miu[3]
                            sum_miu5 += miu[4]
                            sum_miu6 += miu[5]
                            sum_dmiu1_dxj[j] += deriv[0]
                            sum_dmiu1_dyj[j] += deriv[1]
                            sum_dmiu1_dzj[j] += deriv[2]
                            sum_dmiu2_dxj[j] += deriv[3]
                            sum_dmiu2_dyj[j] += deriv[4]
                            sum_dmiu2_dzj[j] += deriv[5]
                            sum_dmiu3_dxj[j] += deriv[6]
                            sum_dmiu3_dyj[j] += deriv[7]
                            sum_dmiu3_dzj[j] += deriv[8]
                            sum_dmiu4_dxj[j] += deriv[9]
                            sum_dmiu4_dyj[j] += deriv[10]
                            sum_dmiu4_dzj[j] += deriv[11]
                            sum_dmiu5_dxj[j] += deriv[12]
                            sum_dmiu5_dyj[j] += deriv[13]
                            sum_dmiu5_dzj[j] += deriv[14]
                            sum_dmiu6_dxj[j] += deriv[15]
                            sum_dmiu6_dyj[j] += deriv[16]
                            sum_dmiu6_dzj[j] += deriv[17]

                    sum_square += group_coefficient * (sum_miu1*sum_miu1 + sum_miu2*sum_miu2 + sum_miu3*sum_miu3 +
                                                       sum_miu4*sum_miu4 + sum_miu5*sum_miu5 + sum_miu6*sum_miu6)
                    dmdx = dmdy = dmdz = 0.0

                    for j in range(nneigh):
                        dmdx = (sum_miu1 * sum_dmiu1_dxj[j] + sum_miu2 * sum_dmiu2_dxj[j] +
                                sum_miu3 * sum_dmiu3_dxj[j] + sum_miu4 * sum_dmiu4_dxj[j] +
                                sum_miu5 * sum_dmiu5_dxj[j] + sum_miu6 * sum_dmiu6_dxj[j]) * group_coefficient * 2.0

                        dmdy = (sum_miu1 * sum_dmiu1_dyj[j] + sum_miu2 * sum_dmiu2_dyj[j] +
                                sum_miu3 * sum_dmiu3_dyj[j] + sum_miu4 * sum_dmiu4_dyj[j] +
                                sum_miu5 * sum_dmiu5_dyj[j] + sum_miu6 * sum_dmiu6_dyj[j]) * group_coefficient * 2.0

                        dmdz = (sum_miu1 * sum_dmiu1_dzj[j] + sum_miu2 * sum_dmiu2_dzj[j] +
                                sum_miu3 * sum_dmiu3_dzj[j] + sum_miu4 * sum_dmiu4_dzj[j] +
                                sum_miu5 * sum_dmiu5_dzj[j] + sum_miu6 * sum_dmiu6_dzj[j]) * group_coefficient * 2.0

                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3] += dmdx
                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3 + 1] += dmdy
                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3 + 2] += dmdz

                        dmcsh[ii*nmcsh + m,i*3]     -= dmdx
                        dmcsh[ii*nmcsh + m,i*3 + 1] -= dmdy
                        dmcsh[ii*nmcsh + m,i*3 + 2] -= dmdz

                group_index += 1
            # sum_square = sum_square * weight
            if (square != 0):
                mcsh[ii,m] = sum_square
            else:
                temp = np.sqrt(sum_square)
                if (abs(temp) < 1e-2):
                    mcsh[ii,m] = 0.0
                    for j in range(nneigh):
                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3] = 0.0
                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3 + 1] = 0.0
                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3 + 2] = 0.0

                else:
                    mcsh[ii,m] = temp
                    for j in range(nneigh):
                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3] *= (0.5 / temp)
                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3 + 1] *= (0.5 / temp)
                        dmcsh[ii*nmcsh + m,nei_list_i[j*2 + 1]*3 + 2] *= (0.5 / temp)

    return mcsh, dmcsh 