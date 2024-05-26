import numpy as np
import matplotlib.pyplot as plt

# P parametrs Definition
p1 = 8.8
p2 = 440
p3 = 100

# d paramters definition
d1 = 1.375e-14
d2 = 1.375e-4
d3 = 3e-5

# d paramters definition
k1 = 1.925e-5
k2 = 1e5
k3 = 1.5e5

init_parameters = {
    "p53" : 100,
    "NDMm" : 150,
    "NDMst" : 200,
    "PTEN" : 100,
    "hop" : 20
}

def calculate_dev_p53(p53_initial_value, NDMm_value):
    
    second_part = d1 * p53_initial_value * NDMm_value

    return p1 - second_part

def calculate_dev_NDMst(p53_value, PTEN_value, NDMst_initial_value):
    
    first_part = p2 * (pow(p53_value, 4) / (pow(p53_value, 4) + pow(k2, 4)))
    second_part = k1 * (pow(k3, 2) /(pow(k3, 2) + pow(PTEN_value, 2))) * NDMst_initial_value
    third_pard = d2 * NDMst_initial_value

    return first_part - second_part - third_pard

def calculate_dev_NDMm(NDMst_value, PTEN_value, NDMm_initial_value):
    
    first_part = k1 * (pow(k3, 2) /(pow(k3, 2) + pow(PTEN_value, 2))) * NDMst_value
    second_part = d2 * NDMm_initial_value
    
    return first_part - second_part

def calculate_dev_PTEN(p53_value, PTEN_initial_value):

    first_part = p3 * (pow(p53_value, 4) / (pow(p53_value, 4) + pow(k2, 4)))
    second_part = d3 * PTEN_initial_value

    return first_part - second_part

def A_value_update(h_value, A_initial_value, kx_value, k_nnumber):
    if k_nnumber in (2, 3):
        h_value = h_value / 2
    
    return A_initial_value + (h_value * kx_value)

def hop_result_p53(p53, NDMm, hop):

    print("Calculating P53 Value after the hop = ", hop)
    k1 = calculate_dev_p53(p53, NDMm)
    print("k1: ", k1)
    new_p53 = A_value_update(hop, p53, k1, 2)
    k2 = calculate_dev_p53(new_p53, NDMm)
    print("k2: ", k2)
    new_p53 = A_value_update(hop, p53, k2, 3)
    k3 = calculate_dev_p53(new_p53, NDMm)
    print("k3: ", k3)    
    new_p53 = A_value_update(hop, p53, k3, 4)
    k4 = calculate_dev_p53(new_p53, NDMm)
    print("k4: ", k4)

    return p53 + (hop / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

def hop_result_NDMm(NDMm, NDMst, PTEN, hop):

    print("Calculating NDMm Value after the hop = ", hop)
    k1 = calculate_dev_NDMm(NDMst, PTEN, NDMm)
    print("k1: ", k1)
    new_NDMm = A_value_update(hop, NDMm, k1, 2)
    k2 = calculate_dev_NDMm(NDMst, PTEN, new_NDMm)
    print("k2: ", k2)
    new_NDMm = A_value_update(hop, NDMm, k2, 3)
    k3 = calculate_dev_NDMm(NDMst, PTEN, new_NDMm)
    print("k3: ", k3)    
    new_NDMm = A_value_update(hop, NDMm, k3, 4)
    k4 = calculate_dev_NDMm(NDMst, PTEN, new_NDMm)
    print("k4: ", k4)
    
    return NDMm + (hop / 6) * (k1 + 2 * k2 + 2 * k3 + k4)



def hop_result_NDMst(p53, NDMst, PTEN, hop):

    print("Calculating NDMst Value after the hop = ", hop)
    k1 = calculate_dev_NDMst(p53, PTEN, NDMst)
    print("k1: ", k1)
    new_NDMst = A_value_update(hop, NDMst, k1, 2)
    k2 = calculate_dev_NDMst(p53, PTEN, new_NDMst)
    print("k2: ", k2)
    new_NDMst = A_value_update(hop, NDMst, k2, 3)
    k3 = calculate_dev_NDMst(p53, PTEN, new_NDMst)
    print("k3: ", k3)    
    new_NDMst = A_value_update(hop, NDMst, k3, 4)
    k4 = calculate_dev_NDMst(p53, PTEN, new_NDMst)
    print("k4: ", k4)
    
    return NDMst + (hop / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def hop_result_PTEN(p53, PTEN, hop):

    print("Calculating PTEN Value after the hop = ", hop)
    k1 = calculate_dev_PTEN(p53, PTEN)
    print("k1: ", k1)
    new_PTEN = A_value_update(hop, PTEN, k1, 2)
    k2 = calculate_dev_PTEN(p53, new_PTEN)
    print("k2: ", k2)
    new_PTEN = A_value_update(hop, PTEN, k2, 3)
    k3 = calculate_dev_PTEN(p53, new_PTEN)
    print("k3: ", k3)    
    new_PTEN = A_value_update(hop, PTEN, k3, 4)
    k4 = calculate_dev_PTEN(p53, new_PTEN)
    print("k4: ", k4)
    
    return PTEN + (hop / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def calculate_symulation(p53, NDMm, NDMst, PTEN, hop, time):

    p53_array = np.array([p53])
    NDMm_array = np.array([NDMm])
    NDMst_array = np.array([NDMst])
    PTEN_array = np.array([PTEN])

    for i in range(time):
        i += hop

        p53 = hop_result_p53(p53, NDMm, hop)
        p53_array = np.append(p53_array, p53)

        NDMm = hop_result_NDMm(NDMm, NDMst, PTEN, hop)
        NDMm_array = np.append(NDMm_array, NDMm)

        NDMst = hop_result_NDMst(p53, NDMst, PTEN, hop)
        NDMst_array = np.append(NDMst_array, NDMst)

        PTEN = hop_result_PTEN(p53, PTEN, hop)
        PTEN_array = np.append(PTEN_array, PTEN)

    
    print("P53: ", p53_array)
    print("NDMm: ", NDMm_array)
    print("NDMst: ", NDMst_array)
    print("PTEN: ", PTEN_array)
    return p53_array, NDMm_array, NDMst_array, PTEN_array


p53_array, NDMm_array, NDMst_array, PTEN_array = calculate_symulation(**init_parameters, time=17280) 

plt.plot(p53_array)
plt.plot(NDMm_array)
plt.plot(NDMst_array)
plt.plot(PTEN_array)
plt.ylabel("P53 value")
plt.xlabel("Time")
plt.show()
    