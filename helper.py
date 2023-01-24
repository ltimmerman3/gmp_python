import numpy as np
import math
M_PI = 3.14159265358979323846

NUM_IMPLEMENTED_TYPE = 5
IMPLEMENTED_MCSH_TYPE = [0,1,2,3,4]

def calc_C1(A, B,  alpha, beta):
    temp = np.sqrt(M_PI / (alpha + beta))
    return A * B * temp * temp * temp


def calc_C2(alpha, beta):
    return -1.0 * (alpha * beta / (alpha + beta))



def calc_lambda( alpha,  beta):
    return beta / (alpha + beta)


def calc_gamma( alpha,  beta):
    return alpha + beta



def P1( lamda,  x0,  gamma):
    return lamda * x0


def P2( lamda,  x0,  gamma):
    lambda_x0_2 = lamda * lamda * x0 * x0
    return (0.5 / gamma) + lambda_x0_2


def P3( lamda,  x0,  gamma):
    lambda_x0 = lamda * x0
    lambda_x0_3 = lambda_x0 * lambda_x0 * lambda_x0
    return (1.5 * lambda_x0 / gamma) + lambda_x0_3


def P4( lamda,  x0,  gamma):
    lambda_x0_2 = lamda * lamda * x0 * x0
    lambda_x0_4 = lambda_x0_2 * lambda_x0_2
    return (0.75 / (gamma * gamma)) + (3.0 * lambda_x0_2 / gamma) + lambda_x0_4


def P5( lamda,  x0,  gamma):
    lambda_x0 = lamda * x0
    lambda_x0_2 = lambda_x0 * lambda_x0
    lambda_x0_3 = lambda_x0 * lambda_x0_2
    lambda_x0_5 = lambda_x0_3 * lambda_x0_2
    return ((15.0 * lambda_x0) / (4.0 * gamma * gamma)) + (5.0 * lambda_x0_3 / gamma) + lambda_x0_5


def P6(lamda, x0, gamma):
    lambda_x0_2 = lamda * lamda * x0 * x0
    lambda_x0_4 = lambda_x0_2 * lambda_x0_2
    lambda_x0_6 = lambda_x0_4 * lambda_x0_2
    return (15.0 / (8.0 * gamma * gamma * gamma)) + ((11.25 * lambda_x0_2) / (gamma * gamma)) + (7.5 * lambda_x0_4 / gamma) + lambda_x0_6


def P7(lamda, x0, gamma):
    lambda_x0 = lamda * x0
    lambda_x0_2 = lambda_x0 * lambda_x0
    lambda_x0_3 = lambda_x0 * lambda_x0_2
    lambda_x0_5 = lambda_x0_3 * lambda_x0_2
    lambda_x0_7 = lambda_x0_5 * lambda_x0_2
    term1 = 105.0 * lambda_x0 / (8.0 * gamma * gamma * gamma)
    term2 = 105.0 * lambda_x0_3 / (4.0 * gamma * gamma)
    term3 = 10.5 * lambda_x0_5 / gamma
    return term1 + term2 + term3 + lambda_x0_7


def P8(lamda, x0, gamma):
    lambda_x0_2 = lamda * lamda * x0 * x0
    lambda_x0_4 = lambda_x0_2 * lambda_x0_2
    lambda_x0_6 = lambda_x0_4 * lambda_x0_2
    lambda_x0_8 = lambda_x0_6 * lambda_x0_2
    term1 = 105.0 / (16.0 * gamma * gamma * gamma * gamma)
    term2 = 105.0 * lambda_x0_2 / (2.0 * gamma * gamma * gamma) 
    term3 = 105.0 * lambda_x0_4 / (2.0 * gamma * gamma)
    term4 = 14.0 * lambda_x0_6 / gamma
    return term1 + term2 + term3 + term4 + lambda_x0_8


def P9(lamda, x0, gamma):
    lambda_x0 = lamda * x0
    lambda_x0_2 = lambda_x0 * lambda_x0
    lambda_x0_3 = lambda_x0 * lambda_x0_2
    lambda_x0_5 = lambda_x0_3 * lambda_x0_2
    lambda_x0_7 = lambda_x0_5 * lambda_x0_2
    lambda_x0_9 = lambda_x0_7 * lambda_x0_2
    term1 = 945.0 * lambda_x0 / (16.0 * gamma * gamma * gamma * gamma)
    term2 = 315.0 * lambda_x0_3 / (2.0 * gamma * gamma * gamma)
    term3 = 189.0 * lambda_x0_5 / (2.0 * gamma * gamma)
    term4 = 18.0 * lambda_x0_7 / gamma
    return term1 + term2 + term3 + term4 + lambda_x0_9


def dP1(lamda, x0, gamma):
    return lamda


def dP2(lamda, x0, gamma):
    return 2.0 * lamda * lamda * x0


def dP3(lamda, x0, gamma):
    lambda_x0 = lamda * x0
    lambda_x0_2 = lambda_x0 * lambda_x0
    return lamda * (3.0 * lambda_x0_2 + (3.0 / (2.0 * gamma)))


def dP4(lamda, x0, gamma):
    lambda_x0 = lamda * x0
    lambda_x0_3 = lambda_x0 * lambda_x0 * lambda_x0
    return lamda * (4.0 * lambda_x0_3 + ((6.0 * lambda_x0) / gamma))


def dP5(lamda, x0, gamma):
    lambda_x0 = lamda * x0
    lambda_x0_2 = lambda_x0 * lambda_x0
    lambda_x0_4 = lambda_x0_2 * lambda_x0_2
    return lamda * (5.0 * lambda_x0_4 + ((15.0 * lambda_x0_2) / gamma) + (15.0 / (4.0 * gamma * gamma)))


def dP6(lamda, x0, gamma):
    lambda_x0 = lamda * x0
    lambda_x0_2 = lambda_x0 * lambda_x0
    lambda_x0_3 = lambda_x0 * lambda_x0_2
    lambda_x0_5 = lambda_x0_3 * lambda_x0_2
    return lamda * (6.0 * lambda_x0_5 + ((30.0 * lambda_x0_3) / gamma) + (45.0 * lambda_x0 / (2.0 * gamma * gamma)))


def dP7(lamda, x0, gamma):
    lambda_x0 = lamda * x0
    lambda_x0_2 = lambda_x0 * lambda_x0
    lambda_x0_4 = lambda_x0_2 * lambda_x0_2
    lambda_x0_6 = lambda_x0_4 * lambda_x0_2
    return lamda * (7.0 * lambda_x0_6 + ((105.0 * lambda_x0_4) / (2.0 * gamma)) + ((315.0 * lambda_x0_2) / (4.0 * gamma * gamma)) + (105.0 / (8.0 * gamma * gamma * gamma)))


def dP8(lamda, x0, gamma):
    lambda_x0 = lamda * x0
    lambda_x0_2 = lambda_x0 * lambda_x0
    lambda_x0_3 = lambda_x0 * lambda_x0_2
    lambda_x0_5 = lambda_x0_3 * lambda_x0_2
    lambda_x0_7 = lambda_x0_5 * lambda_x0_2
    return lamda * (8.0 * lambda_x0_7 + ((84.0 * lambda_x0_5) / gamma) + (210.0 * lambda_x0_3 / (gamma * gamma)) + ((105.0 * lambda_x0) / (gamma * gamma * gamma)))


def dP9(lamda, x0, gamma):
    lambda_x0 = lamda * x0
    lambda_x0_2 = lambda_x0 * lambda_x0
    lambda_x0_4 = lambda_x0_2 * lambda_x0_2
    lambda_x0_6 = lambda_x0_4 * lambda_x0_2
    lambda_x0_8 = lambda_x0_6 * lambda_x0_2
    return lamda * (9.0 * lambda_x0_8 + ((126.0 * lambda_x0_6) / gamma) + ((945.0 * lambda_x0_4) / (2.0 * gamma * gamma)) + ((945.0 * lambda_x0_2) / (2.0 * gamma * gamma * gamma)) + (945.0 / (16.0 * gamma * gamma * gamma * gamma)))

 

def dP1_exp(P, C2, lamda, x0, gamma):
    temp = 2.0 * C2 * x0
    return P * temp + dP1(lamda, x0, gamma)


def dP2_exp(P, C2, lamda, x0, gamma):
    temp = 2.0 * C2 * x0
    return P * temp + dP2(lamda, x0, gamma)


def dP3_exp(P, C2, lamda, x0, gamma):
    temp = 2.0 * C2 * x0
    return P * temp + dP3(lamda, x0, gamma)


def dP4_exp(P, C2, lamda, x0, gamma):
    temp = 2.0 * C2 * x0
    return P * temp + dP4(lamda, x0, gamma)


def dP5_exp(P, C2, lamda, x0, gamma):
    temp = 2.0 * C2 * x0
    return P * temp + dP5(lamda, x0, gamma)


def dP6_exp(P, C2, lamda, x0, gamma):
    temp = 2.0 * C2 * x0
    return P * temp + dP6(lamda, x0, gamma)


def dP7_exp(P, C2, lamda, x0, gamma):
    temp = 2.0 * C2 * x0
    return P * temp + dP7(lamda, x0, gamma)


def dP8_exp(P, C2, lamda, x0, gamma):
    temp = 2.0 * C2 * x0
    return P * temp + dP8(lamda, x0, gamma)


def dP9_exp(P, C2, lamda, x0, gamma):
    temp = 2.0 * C2 * x0
    return P * temp + dP9(lamda, x0, gamma)



def get_num_groups(mcsh_order):
    if (mcsh_order == 0): 
        return 1
    elif (mcsh_order == 1): 
        return 1
    elif (mcsh_order == 2): 
        return 2
    elif (mcsh_order == 3): 
        return 3
    elif (mcsh_order == 4): 
        return 4
    elif (mcsh_order == 5): 
        return 5
    elif (mcsh_order == 6): 
        return 7
    elif (mcsh_order == 7): 
        return 8
    elif (mcsh_order == 8): 
        return 10
    elif (mcsh_order == 9): 
        return 12
    else:
        return 0


def get_group_coefficients(mcsh_order, group_num):
    if (mcsh_order == 0):
        if (group_num == 1):
            return 1.0
        else: 
            return 0.0
        
    elif (mcsh_order == 1):
        if (group_num == 1):
            return 1.0
        else: 
            return 0.0
        
    elif (mcsh_order == 2):
        if (group_num == 1):
            return 1.0
        elif (group_num == 2):
            return 2.0
        else: 
            return 0.0
        
    elif (mcsh_order == 3):
        if (group_num == 1):
            return 1.0
        elif (group_num == 2):
            return 3.0
        elif (group_num == 3):
            return 6.0
        else: 
            return 0.0
        
    elif (mcsh_order == 4):
        if (group_num == 1):
            return 1.0
        elif (group_num == 2):
            return 4.0
        elif (group_num == 3):
            return 6.0
        elif (group_num == 4):
            return 12.0
        else: 
            return 0.0
        
    elif (mcsh_order == 5):
        if (group_num == 1):
            return 1.0
        elif (group_num == 2):
            return 5.0
        elif (group_num == 3):
            return 10.0
        elif (group_num == 4):
            return 20.0
        elif (group_num == 5):
            return 30.0
        else:
            return 0.0
        
    elif (mcsh_order == 6):
        if (group_num == 1):
            return 1.0
        elif (group_num == 2):
            return 6.0
        elif (group_num == 3):
            return 15.0
        elif (group_num == 4):
            return 30.0
        elif (group_num == 5):
            return 20.0
        elif (group_num == 6):
            return 60.0
        elif (group_num == 7):
            return 90.0
        else:
            return 0.0
        
    elif (mcsh_order == 7):
        if (group_num == 1):
            return 1.0
        elif (group_num == 2):
            return 7.0
        elif (group_num == 3):
            return 21.0
        elif (group_num == 4):
            return 42.0
        elif (group_num == 5):
            return 35.0
        elif (group_num == 6):
            return 105.0
        elif (group_num == 7):
            return 140.0
        elif (group_num == 8):
            return 210.0
        else:
            return 0.0
        
    elif (mcsh_order == 8):
        if (group_num == 1):
            return 1.0
        elif (group_num == 2):
            return 8.0
        elif (group_num == 3):
            return 28.0
        elif (group_num == 4):
            return 56.0
        elif (group_num == 5):
            return 56.0
        elif (group_num == 6):
            return 168.0
        elif (group_num == 7):
            return 70.0
        elif (group_num == 8):
            return 280.0
        elif (group_num == 9):
            return 420.0
        elif (group_num == 10):
            return 560.0
        else:
            return 0.0
        
    elif (mcsh_order == 9):
        if (group_num == 1):
            return 1.0
        elif (group_num == 2):
            return 9.0
        elif (group_num == 3):
            return 36.0
        elif (group_num == 4):
            return 72.0
        elif (group_num == 5):
            return 84.0
        elif (group_num == 6):
            return 252.0
        elif (group_num == 7):
            return 126.0
        elif (group_num == 8):
            return 504.0
        elif (group_num == 9):
            return 756.0
        elif (group_num == 10):
            return 630.0
        elif (group_num == 11):
            return 1260.0
        elif (group_num == 12):
            return 1680.0
        else:
            return 0.0
        
    else:
        return 0.0
    



def get_mcsh_type(mcsh_order, group_num):

    if (mcsh_order == 0):
        if (group_num == 1):
            return 1
        else:
            return 0
        
    elif (mcsh_order == 1):
        if (group_num == 1):
            return 2
        else:
            return 0
        
    elif (mcsh_order == 2):
        if (group_num == 1):
            return 2
        elif (group_num == 2):
            return 2
        else:
            return 0
        
    elif (mcsh_order == 3):
        if (group_num == 1):
            return 2
        elif (group_num == 2):
            return 3
        elif (group_num == 3):
            return 1
        else:
            return 0
        
    elif (mcsh_order == 4):
        if (group_num == 1):
            return 2
        elif (group_num == 2):
            return 3
        elif (group_num == 3):
            return 2
        elif (group_num == 4):
            return 2
        else:
            return 0
        
    elif (mcsh_order == 5):
        if (group_num == 1):
            return 2
        elif (group_num == 2):
            return 3
        elif (group_num == 3):
            return 3
        elif (group_num == 4):
            return 2
        elif (group_num == 5):
            return 2
        else:
            return 0
        
    elif (mcsh_order == 6):
        if (group_num == 1):
            return 2
        elif (group_num == 2):
            return 3
        elif (group_num == 3):
            return 3
        elif (group_num == 4):
            return 2
        elif (group_num == 5):
            return 2
        elif (group_num == 6):
            return 3
        elif (group_num == 7):
            return 1
        else:
            return 0
        
    elif (mcsh_order == 7):
        if (group_num == 1):
            return 2
        elif (group_num == 2):
            return 3
        elif (group_num == 3):
            return 3
        elif (group_num == 4):
            return 2
        elif (group_num == 5):
            return 3
        elif (group_num == 6):
            return 3
        elif (group_num == 7):
            return 2
        elif (group_num == 8):
            return 2
        else:
            return 0
        
    elif (mcsh_order == 8):
        if (group_num == 1):
            return 2
        elif (group_num == 2):
            return 3
        elif (group_num == 3):
            return 3
        elif (group_num == 4):
            return 2
        elif (group_num == 5):
            return 3
        elif (group_num == 6):
            return 3
        elif (group_num == 7):
            return 2
        elif (group_num == 8):
            return 3
        elif (group_num == 9):
            return 2
        elif (group_num == 10):
            return 2
        else:
            return 0
        
    elif (mcsh_order == 9):
        if (group_num == 1):
            return 2
        elif (group_num == 2):
            return 3
        elif (group_num == 3):
            return 3
        elif (group_num == 4):
            return 2
        elif (group_num == 5):
            return 3
        elif (group_num == 6):
            return 3
        elif (group_num == 7):
            return 3
        elif (group_num == 8):
            return 3
        elif (group_num == 9):
            return 2
        elif (group_num == 10):
            return 2
        elif (group_num == 11):
            return 3
        elif (group_num == 12):
            return 1
        else:
            return 0
        
    else:
        return 0