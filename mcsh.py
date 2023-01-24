import numpy as np
from math import ceil
from helper import *

M_PI = 3.14159265358979323846

NUM_IMPLEMENTED_TYPE = 5
IMPLEMENTED_MCSH_TYPE = [0,1,2,3,4]

def get_solid_mcsh_function(mcsh_order, group_num):
    result = 0

    if mcsh_order == 0:
        if group_num == 1:
            result = calc_solid_MCSH_0_1
    elif mcsh_order == 1:
        if group_num == 1:
            result = calc_solid_MCSH_1_1
    elif mcsh_order == 2:
        if group_num == 1:
            result = calc_solid_MCSH_2_1
        elif group_num == 2:
            result = calc_solid_MCSH_2_2
    elif mcsh_order == 3:
        if group_num == 1:
            result = calc_solid_MCSH_3_1
        elif group_num == 2:
            result = calc_solid_MCSH_3_2
        elif group_num == 3:
            result = calc_solid_MCSH_3_3
    elif mcsh_order == 4:
        if group_num == 1:
            result = calc_solid_MCSH_4_1
        elif group_num == 2:
            result = calc_solid_MCSH_4_2
        elif group_num == 3:
            result = calc_solid_MCSH_4_3
        elif group_num == 4:
            result = calc_solid_MCSH_4_4
    # elif mcsh_order == 5:
    #     if group_num == 1:
    #         result = calc_solid_MCSH_5_1
    #     elif group_num == 2:
    #         result = calc_solid_MCSH_5_2
    #     elif group_num == 3:
    #         result = calc_solid_MCSH_5_3
    #     elif group_num == 4:
    #         result = calc_solid_MCSH_5_4
    #     elif group_num == 5:
    #         result = calc_solid_MCSH_5_5
    # elif mcsh_order == 6:
    #     if group_num == 1:
    #         result = calc_solid_MCSH_6_1
    #     elif group_num == 2:
    #         result = calc_solid_MCSH_6_2
    #     elif group_num == 3:
    #         result = calc_solid_MCSH_6_3
    #     elif group_num == 4:
    #         result = calc_solid_MCSH_6_4
    #     elif group_num == 5:
    #         result = calc_solid_MCSH_6_5
    #     elif group_num == 6:
    #         result = calc_solid_MCSH_6_6
    #     elif group_num == 7:
    #         result = calc_solid_MCSH_6_7
    # elif mcsh_order == 7:
    #     if group_num == 1:
    #         result = calc_solid_MCSH_7_1
    #     elif group_num == 2:
    #         result = calc_solid_MCSH_7_2
    #     elif group_num == 3:
    #         result = calc_solid_MCSH_7_3
    #     elif group_num == 4:
    #         result = calc_solid_MCSH_7_4
    #     elif group_num == 5:
    #         result = calc_solid_MCSH_7_5
    #     elif group_num == 6:
    #         result = calc_solid_MCSH_7_6
    #     elif group_num == 7:
    #         result = calc_solid_MCSH_7_7
    #     elif group_num == 8:
    #         result = calc_solid_MCSH_7_8
    # elif mcsh_order == 8:
    #     if group_num == 1:
    #         result = calc_solid_MCSH_8_1
    #     elif group_num == 2:
    #         result = calc_solid_MCSH_8_2
    #     elif group_num == 3:
    #         result = calc_solid_MCSH_8_3
    #     elif group_num == 4:
    #         result = calc_solid_MCSH_8_4
    #     elif group_num == 5:
    #         result = calc_solid_MCSH_8_5
    #     elif group_num == 6:
    #         result = calc_solid_MCSH_8_6
    #     elif group_num == 7:
    #         result = calc_solid_MCSH_8_7
    #     elif group_num == 8:
    #         result = calc_solid_MCSH_8_8
    #     elif group_num == 9:
    #         result = calc_solid_MCSH_8_9
    #     elif group_num == 10:
    #         result = calc_solid_MCSH_8_10
    # elif mcsh_order == 9:
    #     if group_num == 1:
    #         result = calc_solid_MCSH_9_1
    #     elif group_num == 2:
    #         result = calc_solid_MCSH_9_2
    #     elif group_num == 3:
    #         result = calc_solid_MCSH_9_3
    #     elif group_num == 4:
    #         result = calc_solid_MCSH_9_4
    #     elif group_num == 5:
    #         result = calc_solid_MCSH_9_5
    #     elif group_num == 6:
    #         result = calc_solid_MCSH_9_6
    #     elif group_num == 7:
    #         result = calc_solid_MCSH_9_7
    #     elif group_num == 8:
    #         result = calc_solid_MCSH_9_8
    #     elif group_num == 9:
    #         result = calc_solid_MCSH_9_9
    #     elif group_num == 10:
    #         result = calc_solid_MCSH_9_10
    #     elif group_num == 11:
    #         result = calc_solid_MCSH_9_11
    #     elif group_num == 12:
    #         result = calc_solid_MCSH_9_12
    
    return result




def calc_solid_MCSH_0_1(x0, y0, z0, r0_sqr, A, B, alpha, beta):
    C1 = calc_C1(A,B,alpha,beta)
    C2 = calc_C2(alpha, beta)
    m_0_1 = C1 * np.exp(C2 * r0_sqr)

    deriv = [m_0_1 * (2.0 * C2 * x0), m_0_1 * (2.0 * C2 * y0), m_0_1 * (2.0 * C2 * z0)]
    value = m_0_1

    return value, deriv

def calc_solid_MCSH_1_1(x0, y0, z0, r0_sqr, A, B, alpha, beta):
    C1 = calc_C1(A, B, alpha, beta)
    C2 = calc_C2(alpha, beta)
    temp = C1 * np.exp(C2 * r0_sqr)

    lamda = calc_lambda(alpha, beta)
    gamma = calc_gamma(alpha, beta)

    P1x = P1(lamda, x0, gamma)
    P1y = P1(lamda, y0, gamma)
    P1z = P1(lamda, z0, gamma)

    term_1 = (1.0 * P1x)
    term_2 = (1.0 * P1y)
    term_3 = (1.0 * P1z)

    miu_1 = temp * term_1
    miu_2 = temp * term_2
    miu_3 = temp * term_3

    value = [miu_1, miu_2, miu_3]

    dP0x_exp = 2.0 * C2 * x0
    dP0y_exp = 2.0 * C2 * y0
    dP0z_exp = 2.0 * C2 * z0
    
    dP1x_exp = dP1_exp(P1x, C2, lamda, x0, gamma)
    dP1y_exp = dP1_exp(P1y, C2, lamda, y0, gamma)
    dP1z_exp = dP1_exp(P1z, C2, lamda, z0, gamma)

    dterm_1_dx = (1.0 * dP1x_exp)
    dterm_1_dy = (1.0 * P1x * dP0y_exp)
    dterm_1_dz = (1.0 * P1x * dP0z_exp)

    dterm_2_dx = (1.0 * dP0x_exp * P1y)
    dterm_2_dy = (1.0 * dP1y_exp)
    dterm_2_dz = (1.0 * P1y * dP0z_exp)

    dterm_3_dx = (1.0 * dP0x_exp * P1z)
    dterm_3_dy = (1.0 * dP0y_exp * P1z)
    dterm_3_dz = (1.0 * dP1z_exp)

    deriv = [temp * dterm_1_dx, temp * dterm_1_dy, temp * dterm_1_dz, temp * dterm_2_dx, temp * dterm_2_dy, temp * dterm_2_dz, temp * dterm_3_dx, temp * dterm_3_dy, temp * dterm_3_dz]

    return value, deriv

def calc_solid_MCSH_2_1(x0, y0, z0, r0_sqr, A, B, alpha, beta):
    C1 = calc_C1(A,B,alpha,beta)
    C2 = calc_C2(alpha, beta)
    temp = C1 * np.exp( C2 * r0_sqr)

    lamda = calc_lambda(alpha, beta)
    gamma = calc_gamma(alpha, beta)

    P2x = P2(lamda, x0, gamma)
    P2y = P2(lamda, y0, gamma)
    P2z = P2(lamda, z0, gamma)

    term_1 = (2.0 * P2x) - (1.0 * P2y) - (1.0 * P2z)
    term_2 = (2.0 * P2y) - (1.0 * P2x) - (1.0 * P2z)
    term_3 = (2.0 * P2z) - (1.0 * P2x) - (1.0 * P2y)

    miu_1 = temp * term_1
    miu_2 = temp * term_2
    miu_3 = temp * term_3

    value = [miu_1, miu_2, miu_3]

    dP0x_exp = 2.0 * C2 * x0
    dP0y_exp = 2.0 * C2 * y0
    dP0z_exp = 2.0 * C2 * z0
    
    dP2x_exp = dP2_exp(P2x, C2, lamda, x0, gamma)
    dP2y_exp = dP2_exp(P2y, C2, lamda, y0, gamma)
    dP2z_exp = dP2_exp(P2z, C2, lamda, z0, gamma)

    dterm_1_dx = (2.0 * dP2x_exp) - (1.0 * dP0x_exp * P2y) - (1.0 * dP0x_exp * P2z)
    dterm_1_dy = (2.0 * P2x * dP0y_exp) - (1.0 * dP2y_exp) - (1.0 * dP0y_exp * P2z)
    dterm_1_dz = (2.0 * P2x * dP0z_exp) - (1.0 * P2y * dP0z_exp) - (1.0 * dP2z_exp)

    dterm_2_dx = (2.0 * dP0x_exp * P2y) - (1.0 * dP2x_exp) - (1.0 * dP0x_exp * P2z)
    dterm_2_dy = (2.0 * dP2y_exp) - (1.0 * P2x * dP0y_exp) - (1.0 * dP0y_exp * P2z)
    dterm_2_dz = (2.0 * P2y * dP0z_exp) - (1.0 * P2x * dP0z_exp) - (1.0 * dP2z_exp)

    dterm_3_dx = (2.0 * dP0x_exp * P2z) - (1.0 * dP2x_exp) - (1.0 * dP0x_exp * P2y)
    dterm_3_dy = (2.0 * dP0y_exp * P2z) - (1.0 * P2x * dP0y_exp) - (1.0 * dP2y_exp)
    dterm_3_dz = (2.0 * dP2z_exp) - (1.0 * P2x * dP0z_exp) - (1.0 * P2y * dP0z_exp)

    deriv = [temp * dterm_1_dx, temp * dterm_1_dy, temp * dterm_1_dz,
            temp * dterm_2_dx, temp * dterm_2_dy, temp * dterm_2_dz,
            temp * dterm_3_dx, temp * dterm_3_dy, temp * dterm_3_dz]

    return value, deriv

def calc_solid_MCSH_2_2(x0, y0, z0, r0_sqr, A, B, alpha, beta):
    C1 = calc_C1(A, B, alpha, beta)
    C2 = calc_C2(alpha, beta)
    temp = C1 * np.exp(C2 * r0_sqr)

    lamda = calc_lambda(alpha, beta)
    gamma = calc_gamma(alpha, beta)

    P1x = P1(lamda, x0, gamma)
    P1y = P1(lamda, y0, gamma)
    P1z = P1(lamda, z0, gamma)

    term_1 = (3.0 * P1x * P1y)
    term_2 = (3.0 * P1x * P1z)
    term_3 = (3.0 * P1y * P1z)

    miu_1 = temp * term_1
    miu_2 = temp * term_2
    miu_3 = temp * term_3

    value = [miu_1,miu_2,miu_3]


    dP0x_exp = 2.0 * C2 * x0
    dP0y_exp = 2.0 * C2 * y0
    dP0z_exp = 2.0 * C2 * z0
    
    dP1x_exp = dP1_exp(P1x, C2, lamda, x0, gamma)
    dP1y_exp = dP1_exp(P1y, C2, lamda, y0, gamma)
    dP1z_exp = dP1_exp(P1z, C2, lamda, z0, gamma)

    dterm_1_dx = (3.0 * dP1x_exp * P1y)
    dterm_1_dy = (3.0 * P1x * dP1y_exp)
    dterm_1_dz = (3.0 * P1x * P1y * dP0z_exp)

    dterm_2_dx = (3.0 * dP1x_exp * P1z)
    dterm_2_dy = (3.0 * P1x * dP0y_exp * P1z)
    dterm_2_dz = (3.0 * P1x * dP1z_exp)

    dterm_3_dx = (3.0 * dP0x_exp * P1y * P1z)
    dterm_3_dy = (3.0 * dP1y_exp * P1z)
    dterm_3_dz = (3.0 * P1y * dP1z_exp)

    deriv = [temp * dterm_1_dx,temp * dterm_1_dy,temp * dterm_1_dz,temp * dterm_2_dx,temp * dterm_2_dy,temp * dterm_2_dz,temp * dterm_3_dx,temp * dterm_3_dy,temp * dterm_3_dz]

    return value, deriv

def calc_solid_MCSH_3_1(x0, y0, z0, r0_sqr, A, B, alpha, beta):
    C1 = calc_C1(A,B,alpha,beta)
    C2 = calc_C2(alpha, beta)
    temp = C1 * np.exp( C2 * r0_sqr)

    lamda = calc_lambda(alpha, beta)
    gamma = calc_gamma(alpha, beta)

    P1x = P1(lamda, x0, gamma)
    P1y = P1(lamda, y0, gamma)
    P1z = P1(lamda, z0, gamma)

    P2x = P2(lamda, x0, gamma)
    P2y = P2(lamda, y0, gamma)
    P2z = P2(lamda, z0, gamma)

    P3x = P3(lamda, x0, gamma)
    P3y = P3(lamda, y0, gamma)
    P3z = P3(lamda, z0, gamma)

    term_1 = (6.0 * P3x) - (9.0 * P1x * P2y) - (9.0 * P1x * P2z)
    term_2 = (6.0 * P3y) - (9.0 * P2x * P1y) - (9.0 * P1y * P2z)
    term_3 = (6.0 * P3z) - (9.0 * P2x * P1z) - (9.0 * P2y * P1z)

    miu_1 = temp * term_1
    miu_2 = temp * term_2
    miu_3 = temp * term_3

    value = [miu_1,miu_2,miu_3]


    dP0x_exp = 2.0 * C2 * x0
    dP0y_exp = 2.0 * C2 * y0
    dP0z_exp = 2.0 * C2 * z0
    
    dP1x_exp = dP1_exp(P1x, C2, lamda, x0, gamma)
    dP1y_exp = dP1_exp(P1y, C2, lamda, y0, gamma)
    dP1z_exp = dP1_exp(P1z, C2, lamda, z0, gamma)

    dP2x_exp = dP2_exp(P2x, C2, lamda, x0, gamma)
    dP2y_exp = dP2_exp(P2y, C2, lamda, y0, gamma)
    dP2z_exp = dP2_exp(P2z, C2, lamda, z0, gamma)

    dP3x_exp = dP3_exp(P3x, C2, lamda, x0, gamma)
    dP3y_exp = dP3_exp(P3y, C2, lamda, y0, gamma)
    dP3z_exp = dP3_exp(P3z, C2, lamda, z0, gamma)

    dterm_1_dx = (6.0 * dP3x_exp) - (9.0 * dP1x_exp * P2y) - (9.0 * dP1x_exp * P2z)
    dterm_1_dy = (6.0 * P3x * dP0y_exp) - (9.0 * P1x * dP2y_exp) - (9.0 * P1x * dP0y_exp * P2z)
    dterm_1_dz = (6.0 * P3x * dP0z_exp) - (9.0 * P1x * P2y * dP0z_exp) - (9.0 * P1x * dP2z_exp)

    dterm_2_dx = (6.0 * dP0x_exp * P3y) - (9.0 * dP2x_exp * P1y) - (9.0 * dP0x_exp * P1y * P2z)
    dterm_2_dy = (6.0 * dP3y_exp) - (9.0 * P2x * dP1y_exp) - (9.0 * dP1y_exp * P2z)
    dterm_2_dz = (6.0 * P3y * dP0z_exp) - (9.0 * P2x * P1y * dP0z_exp) - (9.0 * P1y * dP2z_exp)

    dterm_3_dx = (6.0 * dP0x_exp * P3z) - (9.0 * dP2x_exp * P1z) - (9.0 * dP0x_exp * P2y * P1z)
    dterm_3_dy = (6.0 * dP0y_exp * P3z) - (9.0 * P2x * dP0y_exp * P1z) - (9.0 * dP2y_exp * P1z)
    dterm_3_dz = (6.0 * dP3z_exp) - (9.0 * P2x * dP1z_exp) - (9.0 * P2y * dP1z_exp)

    deriv = [temp * dterm_1_dx,temp * dterm_1_dy,temp * dterm_1_dz,temp * dterm_2_dx,temp * dterm_2_dy,temp * dterm_2_dz,temp * dterm_3_dx,temp * dterm_3_dy,temp * dterm_3_dz]

    return value, deriv

def calc_solid_MCSH_3_2(x0, y0, z0, r0_sqr, A, B, alpha, beta):
    C1 = calc_C1(A,B,alpha,beta)
    C2 = calc_C2(alpha, beta)
    temp = C1 * np.exp( C2 * r0_sqr)
    
    lamda = calc_lambda(alpha, beta)
    gamma = calc_gamma(alpha, beta)
    
    P1x = P1(lamda, x0, gamma)
    P1y = P1(lamda, y0, gamma)
    P1z = P1(lamda, z0, gamma)
    
    P2x = P2(lamda, x0, gamma)
    P2y = P2(lamda, y0, gamma)
    P2z = P2(lamda, z0, gamma)
    
    P3x = P3(lamda, x0, gamma)
    P3y = P3(lamda, y0, gamma)
    P3z = P3(lamda, z0, gamma)
    
    term_1 = (12.0 * P2x * P1y) - (3.0 * P3y) - (3.0 * P1y * P2z)
    term_2 = (12.0 * P1x * P2y) - (3.0 * P3x) - (3.0 * P1x * P2z)
    term_3 = (12.0 * P2x * P1z) - (3.0 * P3z) - (3.0 * P2y * P1z)
    term_4 = (12.0 * P1x * P2z) - (3.0 * P3x) - (3.0 * P1x * P2y)
    term_5 = (12.0 * P2y * P1z) - (3.0 * P3z) - (3.0 * P2x * P1z)
    term_6 = (12.0 * P1y * P2z) - (3.0 * P3y) - (3.0 * P2x * P1y)
    
    miu_1 = temp * term_1
    miu_2 = temp * term_2
    miu_3 = temp * term_3
    miu_4 = temp * term_4
    miu_5 = temp * term_5
    miu_6 = temp * term_6
    
    value = [miu_1, miu_2, miu_3, miu_4, miu_5, miu_6]
    
    dP0x_exp = 2.0 * C2 * x0
    dP0y_exp = 2.0 * C2 * y0
    dP0z_exp = 2.0 * C2 * z0
    
    dP1x_exp = dP1_exp(P1x, C2, lamda, x0, gamma)
    dP1y_exp = dP1_exp(P1y, C2, lamda, y0, gamma)
    dP1z_exp = dP1_exp(P1z, C2, lamda, z0, gamma)
    
    dP2x_exp = dP2_exp(P2x, C2, lamda, x0, gamma)
    dP2y_exp = dP2_exp(P2y, C2, lamda, y0, gamma)
    dP2z_exp = dP2_exp(P2z, C2, lamda, z0, gamma)
    
    dP3x_exp = dP3_exp(P3x, C2, lamda, x0, gamma)
    dP3y_exp = dP3_exp(P3y, C2, lamda, y0, gamma)
    dP3z_exp = dP3_exp(P3z, C2, lamda, z0, gamma)
    
    dterm_1_dx = (12.0 * dP2x_exp * P1y) - (3.0 * dP0x_exp * P3y) - (3.0 * dP0x_exp * P1y * P2z)
    dterm_1_dy = (12.0 * P2x * dP1y_exp) - (3.0 * dP3y_exp) - (3.0 * dP1y_exp * P2z)
    dterm_1_dz = (12.0 * P2x * P1y * dP0z_exp) - (3.0 * P3y * dP0z_exp) - (3.0 * P1y * dP2z_exp)
    
    dterm_2_dx = (12.0 * dP1x_exp * P2y) - (3.0 * dP3x_exp) - (3.0 * dP1x_exp * P2z)
    dterm_2_dy = (12.0 * P1x * dP2y_exp) - (3.0 * P3x * dP0y_exp) - (3.0 * P1x * dP0y_exp * P2z)
    dterm_2_dz = (12.0 * P1x * P2y * dP0z_exp) - (3.0 * P3x * dP0z_exp) - (3.0 * P1x * dP2z_exp)
    
    dterm_3_dx = (12.0 * dP2x_exp * P1z) - (3.0 * dP0x_exp * P3z) - (3.0 * dP0x_exp * P2y * P1z)
    dterm_3_dy = (12.0 * P2x * dP0y_exp * P1z) - (3.0 * dP0y_exp * P3z) - (3.0 * dP2y_exp * P1z)
    dterm_3_dz = (12.0 * P2x * dP1z_exp) - (3.0 * dP3z_exp) - (3.0 * P2y * dP1z_exp)

    dterm_4_dx = (12.0 * dP1x_exp * P2z) - (3.0 * dP3x_exp) - (3.0 * dP1x_exp * P2y)
    dterm_4_dy = (12.0 * P1x * dP0y_exp * P2z) - (3.0 * P3x * dP0y_exp) - (3.0 * P1x * dP2y_exp)
    dterm_4_dz = (12.0 * P1x * dP2z_exp) - (3.0 * P3x * dP0z_exp) - (3.0 * P1x * P2y * dP0z_exp)

    dterm_5_dx = (12.0 * dP0x_exp * P2y * P1z) - (3.0 * dP0x_exp * P3z) - (3.0 * dP2x_exp * P1z)
    dterm_5_dy = (12.0 * dP2y_exp * P1z) - (3.0 * dP0y_exp * P3z) - (3.0 * P2x * dP0y_exp * P1z)
    dterm_5_dz = (12.0 * P2y * dP1z_exp) - (3.0 * dP3z_exp) - (3.0 * P2x * dP1z_exp)

    dterm_6_dx = (12.0 * dP0x_exp * P1y * P2z) - (3.0 * dP0x_exp * P3y) - (3.0 * dP2x_exp * P1y)
    dterm_6_dy = (12.0 * dP1y_exp * P2z) - (3.0 * dP3y_exp) - (3.0 * P2x * dP1y_exp)
    dterm_6_dz = (12.0 * P1y * dP2z_exp) - (3.0 * P3y * dP0z_exp) - (3.0 * P2x * P1y * dP0z_exp)


    deriv = [temp * dterm_1_dx,temp * dterm_1_dy,temp * dterm_1_dz,temp * dterm_2_dx,temp * dterm_2_dy,temp * dterm_2_dz,temp * dterm_3_dx,
    temp * dterm_3_dy,temp * dterm_3_dz,temp * dterm_4_dx,temp * dterm_4_dy,temp * dterm_4_dz,temp * dterm_5_dx,temp * dterm_5_dy,temp * dterm_5_dz,
    temp * dterm_6_dx,temp * dterm_6_dy,temp * dterm_6_dz]

    return value, deriv

def calc_solid_MCSH_3_3(x0, y0, z0, r0_sqr, A, B, alpha, beta):
    C1 = calc_C1(A,B,alpha,beta)
    C2 = calc_C2(alpha, beta)
    temp = C1 * np.exp( C2 * r0_sqr)

    lamda = calc_lambda(alpha, beta)
    gamma = calc_gamma(alpha, beta)

    P1x = P1(lamda, x0, gamma)
    P1y = P1(lamda, y0, gamma)
    P1z = P1(lamda, z0, gamma)

    term_1 = (15.0 * P1x * P1y * P1z)

    m = temp * term_1

    value = m


    # dP0x_exp = 2.0 * C2 * x0
    # dP0y_exp = 2.0 * C2 * y0
    # dP0z_exp = 2.0 * C2 * z0
    
    dP1x_exp = dP1_exp(P1x, C2, lamda, x0, gamma)
    dP1y_exp = dP1_exp(P1y, C2, lamda, y0, gamma)
    dP1z_exp = dP1_exp(P1z, C2, lamda, z0, gamma)

    dterm_1_dx = (15.0 * dP1x_exp * P1y * P1z)
    dterm_1_dy = (15.0 * P1x * dP1y_exp * P1z)
    dterm_1_dz = (15.0 * P1x * P1y * dP1z_exp)

    deriv = [temp * dterm_1_dx, temp * dterm_1_dy, temp * dterm_1_dz]

    return value, deriv


def calc_solid_MCSH_4_1(x0, y0, z0, r0_sqr, A, B, alpha, beta):
    C1 = calc_C1(A,B,alpha,beta)
    C2 = calc_C2(alpha, beta)
    temp = C1 * np.exp( C2 * r0_sqr)

    lamda = calc_lambda(alpha, beta)
    gamma = calc_gamma(alpha, beta)

    P2x = P2(lamda, x0, gamma)
    P2y = P2(lamda, y0, gamma)
    P2z = P2(lamda, z0, gamma)

    P4x = P4(lamda, x0, gamma)
    P4y = P4(lamda, y0, gamma)
    P4z = P4(lamda, z0, gamma)

    term_1 = (24.0 * P4x) - (72.0 * P2x * P2y) - (72.0 * P2x * P2z) + (9.0 * P4y) + (18.0 * P2y * P2z) + (9.0 * P4z)
    term_2 = (24.0 * P4y) - (72.0 * P2x * P2y) - (72.0 * P2y * P2z) + (9.0 * P4x) + (18.0 * P2x * P2z) + (9.0 * P4z)
    term_3 = (24.0 * P4z) - (72.0 * P2x * P2z) - (72.0 * P2y * P2z) + (9.0 * P4x) + (18.0 * P2x * P2y) + (9.0 * P4y)

    miu_1 = temp * term_1
    miu_2 = temp * term_2
    miu_3 = temp * term_3

    value = [miu_1, miu_2, miu_3]


    dP0x_exp = 2.0 * C2 * x0
    dP0y_exp = 2.0 * C2 * y0
    dP0z_exp = 2.0 * C2 * z0
    
    dP2x_exp = dP2_exp(P2x, C2, lamda, x0, gamma)
    dP2y_exp = dP2_exp(P2y, C2, lamda, y0, gamma)
    dP2z_exp = dP2_exp(P2z, C2, lamda, z0, gamma)

    dP4x_exp = dP4_exp(P4x, C2, lamda, x0, gamma)
    dP4y_exp = dP4_exp(P4y, C2, lamda, y0, gamma)
    dP4z_exp = dP4_exp(P4z, C2, lamda, z0, gamma)

    dterm_1_dx = (24.0 * dP4x_exp) - (72.0 * dP2x_exp * P2y) - (72.0 * dP2x_exp * P2z) + (9.0 * dP0x_exp * P4y) + (18.0 * dP0x_exp * P2y * P2z) + (9.0 * dP0x_exp * P4z)
    dterm_1_dy = (24.0 * P4x * dP0y_exp) - (72.0 * P2x * dP2y_exp) - (72.0 * P2x * dP0y_exp * P2z) + (9.0 * dP4y_exp) + (18.0 * dP2y_exp * P2z) + (9.0 * dP0y_exp * P4z)
    dterm_1_dz = (24.0 * P4x * dP0z_exp) - (72.0 * P2x * P2y * dP0z_exp) - (72.0 * P2x * dP2z_exp) + (9.0 * P4y * dP0z_exp) + (18.0 * P2y * dP2z_exp) + (9.0 * dP4z_exp)

    dterm_2_dx = (24.0 * dP0x_exp * P4y) - (72.0 * dP2x_exp * P2y) - (72.0 * dP0x_exp * P2y * P2z) + (9.0 * dP4x_exp) + (18.0 * dP2x_exp * P2z) + (9.0 * dP0x_exp * P4z)
    dterm_2_dy = (24.0 * dP4y_exp) - (72.0 * P2x * dP2y_exp) - (72.0 * dP2y_exp * P2z) + (9.0 * P4x * dP0y_exp) + (18.0 * P2x * dP0y_exp * P2z) + (9.0 * dP0y_exp * P4z)
    dterm_2_dz = (24.0 * P4y * dP0z_exp) - (72.0 * P2x * P2y * dP0z_exp) - (72.0 * P2y * dP2z_exp) + (9.0 * P4x * dP0z_exp) + (18.0 * P2x * dP2z_exp) + (9.0 * dP4z_exp)

    dterm_3_dx = (24.0 * dP0x_exp * P4z) - (72.0 * dP2x_exp * P2z) - (72.0 * dP0x_exp * P2y * P2z) + (9.0 * dP4x_exp) + (18.0 * dP2x_exp * P2y) + (9.0 * dP0x_exp * P4y)
    dterm_3_dy = (24.0 * dP0y_exp * P4z) - (72.0 * P2x * dP0y_exp * P2z) - (72.0 * dP2y_exp * P2z) + (9.0 * P4x * dP0y_exp) + (18.0 * P2x * dP2y_exp) + (9.0 * dP4y_exp)
    dterm_3_dz = (24.0 * dP4z_exp) - (72.0 * P2x * dP2z_exp) - (72.0 * P2y * dP2z_exp) + (9.0 * P4x * dP0z_exp) + (18.0 * P2x * P2y * dP0z_exp) + (9.0 * P4y * dP0z_exp)

    # dmiu1 dx/dy/dz
    deriv = [temp * dterm_1_dx, temp * dterm_1_dy, temp * dterm_1_dz, # dmiu1 dx/dy/dz
             temp * dterm_2_dx, temp * dterm_2_dy, temp * dterm_2_dz, # dmiu2 dx/dy/dz
             temp * dterm_3_dx, temp * dterm_3_dy, temp * dterm_3_dz] # dmiu3 dx/dy/dz

    return value, deriv

def calc_solid_MCSH_4_2(x0, y0, z0, r0_sqr, A, B, alpha, beta):
    C1 = calc_C1(A,B,alpha,beta)
    C2 = calc_C2(alpha, beta)
    temp = C1 * np.exp( C2 * r0_sqr)

    lamda = calc_lambda(alpha, beta)
    gamma = calc_gamma(alpha, beta)

    P1x = P1(lamda, x0, gamma)
    P1y = P1(lamda, y0, gamma)
    P1z = P1(lamda, z0, gamma)

    P2x = P2(lamda, x0, gamma)
    P2y = P2(lamda, y0, gamma)
    P2z = P2(lamda, z0, gamma)

    P3x = P3(lamda, x0, gamma)
    P3y = P3(lamda, y0, gamma)
    P3z = P3(lamda, z0, gamma)

    term_1 = (60.0 * P3x * P1y) - (45.0 * P1x * P3y) - (45.0 * P1x * P1y * P2z)
    term_2 = (60.0 * P1x * P3y) - (45.0 * P3x * P1y) - (45.0 * P1x * P1y * P2z)
    term_3 = (60.0 * P3x * P1z) - (45.0 * P1x * P3z) - (45.0 * P1x * P2y * P1z)
    term_4 = (60.0 * P1x * P3z) - (45.0 * P3x * P1z) - (45.0 * P1x * P2y * P1z)
    term_5 = (60.0 * P3y * P1z) - (45.0 * P1y * P3z) - (45.0 * P2x * P1y * P1z)
    term_6 = (60.0 * P1y * P3z) - (45.0 * P3y * P1z) - (45.0 * P2x * P1y * P1z)

    miu_1 = temp * term_1
    miu_2 = temp * term_2
    miu_3 = temp * term_3
    miu_4 = temp * term_4
    miu_5 = temp * term_5
    miu_6 = temp * term_6

    value = [miu_1, miu_2, miu_3, miu_4, miu_5, miu_6]

    dP0x_exp = 2.0 * C2 * x0
    dP0y_exp = 2.0 * C2 * y0
    dP0z_exp = 2.0 * C2 * z0
    
    dP1x_exp = dP1_exp(P1x, C2, lamda, x0, gamma)
    dP1y_exp = dP1_exp(P1y, C2, lamda, y0, gamma)
    dP1z_exp = dP1_exp(P1z, C2, lamda, z0, gamma)

    dP2x_exp = dP2_exp(P2x, C2, lamda, x0, gamma)
    dP2y_exp = dP2_exp(P2y, C2, lamda, y0, gamma)
    dP2z_exp = dP2_exp(P2z, C2, lamda, z0, gamma)

    dP3x_exp = dP3_exp(P3x, C2, lamda, x0, gamma)
    dP3y_exp = dP3_exp(P3y, C2, lamda, y0, gamma)
    dP3z_exp = dP3_exp(P3z, C2, lamda, z0, gamma)

    dterm_1_dx = (60.0 * dP3x_exp * P1y) - (45.0 * dP1x_exp * P3y) - (45.0 * dP1x_exp * P1y * P2z)
    dterm_1_dy = (60.0 * P3x * dP1y_exp) - (45.0 * P1x * dP3y_exp) - (45.0 * P1x * dP1y_exp * P2z)
    dterm_1_dz = (60.0 * P3x * P1y * dP0z_exp) - (45.0 * P1x * P3y * dP0z_exp) - (45.0 * P1x * P1y * dP2z_exp)

    dterm_2_dx = (60.0 * dP1x_exp * P3y) - (45.0 * dP3x_exp * P1y) - (45.0 * dP1x_exp * P1y * P2z)
    dterm_2_dy = (60.0 * P1x * dP3y_exp) - (45.0 * P3x * dP1y_exp) - (45.0 * P1x * dP1y_exp * P2z)
    dterm_2_dz = (60.0 * P1x * P3y * dP0z_exp) - (45.0 * P3x * P1y * dP0z_exp) - (45.0 * P1x * P1y * dP2z_exp)

    dterm_3_dx = (60.0 * dP3x_exp * P1z) - (45.0 * dP1x_exp * P3z) - (45.0 * dP1x_exp * P2y * P1z)
    dterm_3_dy = (60.0 * P3x * dP0y_exp * P1z) - (45.0 * P1x * dP0y_exp * P3z) - (45.0 * P1x * dP2y_exp * P1z)
    dterm_3_dz = (60.0 * P3x * dP1z_exp) - (45.0 * P1x * dP3z_exp) - (45.0 * P1x * P2y * dP1z_exp)

    dterm_4_dx = (60.0 * dP1x_exp * P3z) - (45.0 * dP3x_exp * P1z) - (45.0 * dP1x_exp * P2y * P1z)
    dterm_4_dy = (60.0 * P1x * dP0y_exp * P3z) - (45.0 * P3x * dP0y_exp * P1z) - (45.0 * P1x * dP2y_exp * P1z)
    dterm_4_dz = (60.0 * P1x * dP3z_exp) - (45.0 * P3x * dP1z_exp) - (45.0 * P1x * P2y * dP1z_exp)

    dterm_5_dx = (60.0 * dP0x_exp * P3y * P1z) - (45.0 * dP0x_exp * P1y * P3z) - (45.0 * dP2x_exp * P1y * P1z)
    dterm_5_dy = (60.0 * dP3y_exp * P1z) - (45.0 * dP1y_exp * P3z) - (45.0 * P2x * dP1y_exp * P1z)
    dterm_5_dz = (60.0 * P3y * dP1z_exp) - (45.0 * P1y * dP3z_exp) - (45.0 * P2x * P1y * dP1z_exp)

    dterm_6_dx = (60.0 * dP0x_exp * P1y * P3z) - (45.0 * dP0x_exp * P3y * P1z) - (45.0 * dP2x_exp * P1y * P1z)
    dterm_6_dy = (60.0 * dP1y_exp * P3z) - (45.0 * dP3y_exp * P1z) - (45.0 * P2x * dP1y_exp * P1z)
    dterm_6_dz = (60.0 * P1y * dP3z_exp) - (45.0 * P3y * dP1z_exp) - (45.0 * P2x * P1y * dP1z_exp)

    # dmiu1 dx/dy/dz
    deriv = [temp * dterm_1_dx,temp * dterm_1_dy,temp * dterm_1_dz,temp * dterm_2_dx,temp * dterm_2_dy,temp * dterm_2_dz,temp * dterm_3_dx,
    temp * dterm_3_dy,temp * dterm_3_dz,temp * dterm_4_dx,temp * dterm_4_dy,temp * dterm_4_dz,temp * dterm_5_dx,temp * dterm_5_dy,temp * dterm_5_dz,
    temp * dterm_6_dx,temp * dterm_6_dy,temp * dterm_6_dz]

    return value, deriv

def calc_solid_MCSH_4_3(x0, y0, z0, r0_sqr, A, B, alpha, beta):
    C1 = calc_C1(A,B,alpha,beta)
    C2 = calc_C2(alpha,beta)
    temp = C1 * np.exp( C2 * r0_sqr)

    lamda = calc_lambda(alpha, beta)
    gamma = calc_gamma(alpha, beta)

    P2x = P2(lamda, x0, gamma)
    P2y = P2(lamda, y0, gamma)
    P2z = P2(lamda, z0, gamma)

    P4x = P4(lamda, x0, gamma)
    P4y = P4(lamda, y0, gamma)
    P4z = P4(lamda, z0, gamma)

    term_1 = (81.0 * P2x * P2y) - (12.0 * P4x) - (9.0 * P2x * P2z) - (12.0 * P4y) - (9.0 * P2y * P2z) + (3.0 * P4z)
    term_2 = (81.0 * P2x * P2z) - (12.0 * P4x) - (9.0 * P2x * P2y) - (12.0 * P4z) - (9.0 * P2y * P2z) + (3.0 * P4y)
    term_3 = (81.0 * P2y * P2z) - (12.0 * P4y) - (9.0 * P2x * P2y) - (12.0 * P4z) - (9.0 * P2x * P2z) + (3.0 * P4x)

    miu_1 = temp * term_1
    miu_2 = temp * term_2
    miu_3 = temp * term_3

    value = [miu_1, miu_2, miu_3]

    dP0x_exp = 2.0 * C2 * x0
    dP0y_exp = 2.0 * C2 * y0
    dP0z_exp = 2.0 * C2 * z0
    
    dP2x_exp = dP2_exp(P2x, C2, lamda, x0, gamma)
    dP2y_exp = dP2_exp(P2y, C2, lamda, y0, gamma)
    dP2z_exp = dP2_exp(P2z, C2, lamda, z0, gamma)

    dP4x_exp = dP4_exp(P4x, C2, lamda, x0, gamma)
    dP4y_exp = dP4_exp(P4y, C2, lamda, y0, gamma)
    dP4z_exp = dP4_exp(P4z, C2, lamda, z0, gamma)

    dterm_1_dx = (81.0 * dP2x_exp * P2y) - (12.0 * dP4x_exp) - (9.0 * dP2x_exp * P2z) - (12.0 * dP0x_exp * P4y) - (9.0 * dP0x_exp * P2y * P2z) + (3.0 * dP0x_exp * P4z)
    dterm_1_dy = (81.0 * P2x * dP2y_exp) - (12.0 * P4x * dP0y_exp) - (9.0 * P2x * dP0y_exp * P2z) - (12.0 * dP4y_exp) - (9.0 * dP2y_exp * P2z) + (3.0 * dP0y_exp * P4z)
    dterm_1_dz = (81.0 * P2x * P2y * dP0z_exp) - (12.0 * P4x * dP0z_exp) - (9.0 * P2x * dP2z_exp) - (12.0 * P4y * dP0z_exp) - (9.0 * P2y * dP2z_exp) + (3.0 * dP4z_exp)

    dterm_2_dx = (81.0 * dP2x_exp * P2z) - (12.0 * dP4x_exp) - (9.0 * dP2x_exp * P2y) - (12.0 * dP0x_exp * P4z) - (9.0 * dP0x_exp * P2y * P2z) + (3.0 * dP0x_exp * P4y)
    dterm_2_dy = (81.0 * P2x * dP0y_exp * P2z) - (12.0 * P4x * dP0y_exp) - (9.0 * P2x * dP2y_exp) - (12.0 * dP0y_exp * P4z) - (9.0 * dP2y_exp * P2z) + (3.0 * dP4y_exp)
    dterm_2_dz = (81.0 * P2x * dP2z_exp) - (12.0 * P4x * dP0z_exp) - (9.0 * P2x * P2y * dP0z_exp) - (12.0 * dP4z_exp) - (9.0 * P2y * dP2z_exp) + (3.0 * P4y * dP0z_exp)

    dterm_3_dx = (81.0 * dP0x_exp * P2y * P2z) - (12.0 * dP0x_exp * P4y) - (9.0 * dP2x_exp * P2y) - (12.0 * dP0x_exp * P4z) - (9.0 * dP2x_exp * P2z) + (3.0 * dP4x_exp)
    dterm_3_dy = (81.0 * dP2y_exp * P2z) - (12.0 * dP4y_exp) - (9.0 * P2x * dP2y_exp) - (12.0 * dP0y_exp * P4z) - (9.0 * P2x * dP0y_exp * P2z) + (3.0 * P4x * dP0y_exp)
    dterm_3_dz = (81.0 * P2y * dP2z_exp) - (12.0 * P4y * dP0z_exp) - (9.0 * P2x * P2y * dP0z_exp) - (12.0 * dP4z_exp) - (9.0 * P2x * dP2z_exp) + (3.0 * P4x * dP0z_exp)

    deriv = [temp * dterm_1_dx, temp * dterm_1_dy, temp * dterm_1_dz, temp * dterm_2_dx, temp * dterm_2_dy, temp * dterm_2_dz, temp * dterm_3_dx, temp * dterm_3_dy, temp * dterm_3_dz]

    return value, deriv


def calc_solid_MCSH_4_4(x0, y0, z0, r0_sqr, A, B, alpha, beta):
    C1 = calc_C1(A,B,alpha,beta)
    C2 = calc_C2(alpha, beta)
    temp = C1 * np.exp( C2 * r0_sqr)

    lamda = calc_lambda(alpha, beta)
    gamma = calc_gamma(alpha, beta)

    P1x = P1(lamda, x0, gamma)
    P1y = P1(lamda, y0, gamma)
    P1z = P1(lamda, z0, gamma)

    P2x = P2(lamda, x0, gamma)
    P2y = P2(lamda, y0, gamma)
    P2z = P2(lamda, z0, gamma)

    P3x = P3(lamda, x0, gamma)
    P3y = P3(lamda, y0, gamma)
    P3z = P3(lamda, z0, gamma)

    term_1 = (90.0 * P2x * P1y * P1z) - (15.0 * P3y * P1z) - (15.0 * P1y * P3z)
    term_2 = (90.0 * P1x * P2y * P1z) - (15.0 * P3x * P1z) - (15.0 * P1x * P3z)
    term_3 = (90.0 * P1x * P1y * P2z) - (15.0 * P3x * P1y) - (15.0 * P1x * P3y)

    miu_1 = temp * term_1
    miu_2 = temp * term_2
    miu_3 = temp * term_3

    values = [miu_1, miu_2, miu_3]

    dP0x_exp = 2.0 * C2 * x0
    dP0y_exp = 2.0 * C2 * y0
    dP0z_exp = 2.0 * C2 * z0

    dP1x_exp = dP1_exp(P1x, C2, lamda, x0, gamma)
    dP1y_exp = dP1_exp(P1y, C2, lamda, y0, gamma)
    dP1z_exp = dP1_exp(P1z, C2, lamda, z0, gamma)

    dP2x_exp = dP2_exp(P2x, C2, lamda, x0, gamma)
    dP2y_exp = dP2_exp(P2y, C2, lamda, y0, gamma)
    dP2z_exp = dP2_exp(P2z, C2, lamda, z0, gamma)

    dP3x_exp = dP3_exp(P3x, C2, lamda, x0, gamma)
    dP3y_exp = dP3_exp(P3y, C2, lamda, y0, gamma)
    dP3z_exp = dP3_exp(P3z, C2, lamda, z0, gamma)

    dterm_1_dx = (90.0 * dP2x_exp * P1y * P1z) - (15.0 * dP0x_exp * P3y * P1z) - (15.0 * dP0x_exp * P1y * P3z)
    dterm_1_dy = (90.0 * P2x * dP1y_exp * P1z) - (15.0 * dP3y_exp * P1z) - (15.0 * dP1y_exp * P3z)
    dterm_1_dz = (90.0 * P2x * P1y * dP1z_exp) - (15.0 * P3y * dP1z_exp) - (15.0 * P1y * dP3z_exp)

    dterm_2_dx = (90.0 * dP1x_exp * P2y * P1z) - (15.0 * dP3x_exp * P1z) - (15.0 * dP1x_exp * P3z)
    dterm_2_dy = (90.0 * P1x * dP2y_exp * P1z) - (15.0 * P3x * dP0y_exp * P1z) - (15.0 * P1x * dP0y_exp * P3z)
    dterm_2_dz = (90.0 * P1x * P2y * dP1z_exp) - (15.0 * P3x * dP1z_exp) - (15.0 * P1x * dP3z_exp)

    dterm_3_dx = (90.0 * dP1x_exp * P1y * P2z) - (15.0 * dP3x_exp * P1y) - (15.0 * dP1x_exp * P3y)
    dterm_3_dy = (90.0 * P1x * dP1y_exp * P2z) - (15.0 * P3x * dP1y_exp) - (15.0 * P1x * dP3y_exp)
    dterm_3_dz = (90.0 * P1x * P1y * dP2z_exp) - (15.0 * P3x * P1y * dP0z_exp) - (15.0 * P1x * P3y * dP0z_exp)

    # dmiu1 dx/dy/dz
    deriv = [temp * dterm_1_dx, temp * dterm_1_dy, temp * dterm_1_dz, 
             temp * dterm_2_dx, temp * dterm_2_dy, temp * dterm_2_dz, 
             temp * dterm_3_dx, temp * dterm_3_dy, temp * dterm_3_dz]

    return values, deriv