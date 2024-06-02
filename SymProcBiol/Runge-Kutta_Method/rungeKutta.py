import numpy as np
import matplotlib.pyplot as plt

    
def calculate_dev_p53(p53_initial_value, NDMm_value, parameters):
    d1 = parameters["d1"]
    p1 = parameters["p1"]
    
    second_part = d1 * p53_initial_value * NDMm_value

    return p1 - second_part

def calculate_dev_NDMst(p53_value, PTEN_value, NDMst_initial_value, parameters):
    p2 = parameters["p2"]
    k1 = parameters["k1"]
    k2 = parameters["k2"]
    k3 = parameters["k3"]
    d2 = parameters["d2"]
    
    first_part = p2 * (p53_value**4 / (p53_value**4 + pow(k2, 4)))
    second_part = k1 * (pow(k3, 2) /(pow(k3, 2) + pow(PTEN_value, 2))) * NDMst_initial_value
    third_pard = d2 * NDMst_initial_value

    return first_part - second_part - third_pard

def calculate_dev_NDMm(NDMst_value, PTEN_value, NDMm_initial_value, parameters):
    
    k1 = parameters["k1"]
    k3 = parameters["k3"]
    d2 = parameters["d2"]

    first_part = k1 * (pow(k3, 2) /(pow(k3, 2) + pow(PTEN_value, 2))) * NDMst_value
    second_part = d2 * NDMm_initial_value
    
    return first_part - second_part

def calculate_dev_PTEN(p53_value, PTEN_initial_value, parameters):

    p3 = parameters["p3"]
    k2 = parameters["k2"]
    d3 = parameters["d3"]

    first_part = p3 * (pow(p53_value, 4) / (pow(p53_value, 4) + pow(k2, 4)))
    second_part = d3 * PTEN_initial_value

    return first_part - second_part



def A_value_update(h_value, A_initial_value, kx_value, k_nnumber):
    if k_nnumber in (2, 3):
        h_value = h_value / 2
    
    return A_initial_value + (h_value * kx_value)

def hop_result_p53(p53, NDMm, hop, parameters):

    print("Calculating P53 Value after the hop = ", hop)
    k1 = calculate_dev_p53(p53, NDMm, parameters)
    print("k1: ", k1)
    new_p53 = A_value_update(hop, p53, k1, 2)
    k2 = calculate_dev_p53(new_p53, NDMm, parameters)
    print("k2: ", k2)
    new_p53 = A_value_update(hop, p53, k2, 3)
    k3 = calculate_dev_p53(new_p53, NDMm, parameters)
    print("k3: ", k3)    
    new_p53 = A_value_update(hop, p53, k3, 4)
    k4 = calculate_dev_p53(new_p53, NDMm, parameters)
    print("k4: ", k4)

    return p53 + (hop / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

def hop_result_NDMm(NDMm, NDMst, PTEN, hop, parameters):

    print("Calculating NDMm Value after the hop = ", hop)
    k1 = calculate_dev_NDMm(NDMst, PTEN, NDMm, parameters)
    print("k1: ", k1)
    new_NDMm = A_value_update(hop, NDMm, k1, 2)
    k2 = calculate_dev_NDMm(NDMst, PTEN, new_NDMm, parameters)
    print("k2: ", k2)
    new_NDMm = A_value_update(hop, NDMm, k2, 3)
    k3 = calculate_dev_NDMm(NDMst, PTEN, new_NDMm, parameters)
    print("k3: ", k3)    
    new_NDMm = A_value_update(hop, NDMm, k3, 4)
    k4 = calculate_dev_NDMm(NDMst, PTEN, new_NDMm, parameters)
    print("k4: ", k4)
    
    return NDMm + (hop / 6) * (k1 + 2 * k2 + 2 * k3 + k4)



def hop_result_NDMst(p53, NDMst, PTEN, hop, parameters):

    print("Calculating NDMst Value after the hop = ", hop)
    k1 = calculate_dev_NDMst(p53, PTEN, NDMst, parameters)
    print("k1: ", k1)
    new_NDMst = A_value_update(hop, NDMst, k1, 2)
    k2 = calculate_dev_NDMst(p53, PTEN, new_NDMst, parameters)
    print("k2: ", k2)
    new_NDMst = A_value_update(hop, NDMst, k2, 3)
    k3 = calculate_dev_NDMst(p53, PTEN, new_NDMst, parameters)
    print("k3: ", k3)    
    new_NDMst = A_value_update(hop, NDMst, k3, 4)
    k4 = calculate_dev_NDMst(p53, PTEN, new_NDMst, parameters)
    print("k4: ", k4)
    
    return NDMst + (hop / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def hop_result_PTEN(p53, PTEN, hop, parameters):

    print("Calculating PTEN Value after the hop = ", hop)
    k1 = calculate_dev_PTEN(p53, PTEN, parameters)
    print("k1: ", k1)
    new_PTEN = A_value_update(hop, PTEN, k1, 2)
    k2 = calculate_dev_PTEN(p53, new_PTEN, parameters)
    print("k2: ", k2)
    new_PTEN = A_value_update(hop, PTEN, k2, 3)
    k3 = calculate_dev_PTEN(p53, new_PTEN, parameters)
    print("k3: ", k3)    
    new_PTEN = A_value_update(hop, PTEN, k3, 4)
    k4 = calculate_dev_PTEN(p53, new_PTEN, parameters)
    print("k4: ", k4)
    
    return PTEN + (hop / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


def calculate_symulation(p53, NDMm, NDMst, PTEN, hop, time, parameters):
    print("d2 value", parameters["d2"])
    p53_array = np.array([p53])
    NDMm_array = np.array([NDMm])
    NDMst_array = np.array([NDMst])
    PTEN_array = np.array([PTEN])

    for i in range(time):
        i += hop

        p53 = hop_result_p53(p53, NDMm, hop, parameters)
        p53_array = np.append(p53_array, p53)

        NDMm = hop_result_NDMm(NDMm, NDMst, PTEN, hop, parameters)
        NDMm_array = np.append(NDMm_array, NDMm)

        NDMst = hop_result_NDMst(p53, NDMst, PTEN, hop, parameters)
        NDMst_array = np.append(NDMst_array, NDMst)

        PTEN = hop_result_PTEN(p53, PTEN, hop, parameters)
        PTEN_array = np.append(PTEN_array, PTEN)

    
    # print("P53: ", p53_array)
    # print("NDMm: ", NDMm_array)
    # print("NDMst: ", NDMst_array)
    # print("PTEN: ", PTEN_array)
    return p53_array, NDMm_array, NDMst_array, PTEN_array

# p53_array, NDMm_array, NDMst_array, PTEN_array = calculate_symulation(100, 100, 100, 100, 100, 17800)
# plt.plot(p53_array)
# plt.plot(NDMm_array)
# plt.plot(NDMst_array)
# plt.plot(PTEN_array)
# plt.ylabel("P53 value")
# plt.xlabel("Time")
# plt.show()


def RK_method(p53=100, NDMm=100, NDMst=100, PTEN=100, hop=100, time=17280, PTEN_off=False, is_siRNA=False, DNA_damage=False):

    # utworzenie słownika (HashMap) dla każdego z parametrów
    parameters = {"p1": 8.8,
                  "d1": 1.375e-14,
                  "d3": 3e-5,
                  "k1": 1.925e-5,
                  "k2": 1e5,
                  "k3": 1.5e5,
                  }

    # warunek, gdy PTEN ne działa zmienia wartość parametru p3
    if PTEN_off:
        parameters["p3"] = 0.0
    else:
        parameters["p3"] = 100

    # warunek na obecność siRNA, zmienia wartość paramteru p2
    if is_siRNA:
        parameters["p2"] = 0.02
    else:
        parameters["p2"] = 440
  
    # warunek sprawdza czy w mamy uszkodzenia DNA i odpowiednio zmienia wartość parametru d2
    if DNA_damage:
        parameters["d2"] = 1.375e-4
    else:
        parameters["d2"] = 0.1

    # Wywołanie głównej funkcji wykonującej całą symulację
    # przekazanie wszytkich danych wejściowych białek oraz parametrów
    p53_array, NDMm_array, NDMst_array, PTEN_array = calculate_symulation(p53=p53, 
                                                                          NDMm=NDMm, 
                                                                          NDMst=NDMst, 
                                                                          PTEN=PTEN, 
                                                                          hop=hop, 
                                                                          time=time,
                                                                          parameters=parameters)
    
    # rysowanie funkcji ilości białek w czasie
    plt.plot(p53_array)
    plt.plot(NDMm_array)
    plt.plot(NDMst_array)
    plt.plot(PTEN_array)
    plt.ylabel("P53 value")
    plt.xlabel("Time")
    plt.show()


# wywołanie funkcji
RK_method(DNA_damage=True)