import random

# generuj randomowe parametry
def create_parameters_set_left():
    parameters = {}

    parameters["p1"] = random.uniform(0, 20)
    parameters["p2"] = random.uniform(0, 1000)
    parameters["p3"] = random.uniform(0, 1000)
    print(parameters["p1"], parameters["p2"], parameters["p3"])

    parameters["d1"] = random.uniform(1e-14, 1e-13)
    parameters["d2"] = random.uniform(1e-2, 1e-1)
    parameters["d3"] = random.uniform(1e-5, 1e-4)
    print(parameters["d1"], parameters["d2"], parameters["d3"])

    parameters["k1"] = random.uniform(1e-5, 1e-4)
    parameters["k2"] = random.uniform(1e5, 5e5)
    parameters["k3"] = random.uniform(1e5, 5e5)
    print(parameters["k1"], parameters["k2"], parameters["k3"])
    
    return parameters


# generuj randomowe parametry
def create_parameters_set_right():
    parameters = {}

    parameters["p1"] = random.uniform(0, 20)
    parameters["p2"] = random.uniform(0, 1000)
    parameters["p3"] = 0
    print(parameters["p1"], parameters["p2"], parameters["p3"])

    parameters["d1"] = random.uniform(1e-14, 1e-13)
    parameters["d2"] = random.uniform(1e-4, 1e-3)
    parameters["d3"] = random.uniform(1e-5, 1e-4)
    print(parameters["d1"], parameters["d2"], parameters["d3"])

    parameters["k1"] = random.uniform(1e-5, 1e-4)
    parameters["k2"] = random.uniform(1e5, 5e5)
    parameters["k3"] = random.uniform(1e5, 5e5)
    print(parameters["k1"], parameters["k2"], parameters["k3"])
    
    return parameters






# Obliczanie rkIV ze stałym krokiem 
# wykonane tak samo jak w I zadaniu

import numpy as np
import matplotlib.pyplot as plt

# Funkcja obliczająca pochodną dla p53
def calculate_dev_p53(p53_initial_value, NDMm_value, parameters):
    
    # przypisanie wartości parametrów za pomocą słownika  
    d1 = parameters["d1"]
    p1 = parameters["p1"]

    # obliczenie drugiej części pochodnej
    second_part = d1 * p53_initial_value * pow(NDMm_value, 2)

    return p1 - second_part

# Funkcja obliczająca pochodną dla NDMst
def calculate_dev_NDMst(p53_value, PTEN_value, NDMst_initial_value, parameters):
    
    # przypisanie wartości parametrów za pomocą słownika  
    p2 = parameters["p2"]
    k1 = parameters["k1"]
    k2 = parameters["k2"]
    k3 = parameters["k3"]
    d2 = parameters["d2"]
    
    # obliczanie kolejnych części pochodnej
    first_part = p2 * (pow(p53_value, 4) / (pow(p53_value,4) + pow(k2, 4)))
    second_part = k1 * (pow(k3, 2) /(pow(k3, 2) + pow(PTEN_value, 2))) * NDMst_initial_value
    third_pard = d2 * NDMst_initial_value

    return first_part - second_part - third_pard


# Funkcja obliczająca pochodną dla NDMm
def calculate_dev_NDMm(NDMst_value, PTEN_value, NDMm_initial_value, parameters):
    
    # przypisanie wartości parametrów za pomocą słownika  
    k1 = parameters["k1"]
    k3 = parameters["k3"]
    d2 = parameters["d2"]

    # obliczanie kolejnych części pochodnej
    first_part = k1 * (pow(k3, 2) /(pow(k3, 2) + pow(PTEN_value, 2))) * NDMst_value
    second_part = d2 * NDMm_initial_value
    
    return first_part - second_part

# Funkcja obliczająca pochodną dla PTEN
def calculate_dev_PTEN(p53_value, PTEN_initial_value, parameters):

    # przypisanie wartości parametrów za pomocą słownika  
    p3 = parameters["p3"]
    k2 = parameters["k2"]
    d3 = parameters["d3"]

    # obliczanie kolejnych części pochodnej

    first_part = p3 * (pow(p53_value, 4) / (pow(p53_value, 4) + pow(k2, 4)))
    second_part = d3 * PTEN_initial_value

    return first_part - second_part


# Funkcja obliczająca wartości białek dla kolejnych k -> k2, k3, k4
# h_value -> hop value / skok
# A_initial_value -> wartość białka A (PTEN, p53, NDMst, NDMm)
# kx_value -> wartości odpowiednich liczb k np. k2
# k_nnumber -> numer liczby k, która jest obliczana
def A_value_update(h_value, A_initial_value, kx_value, k_nnumber):
    
    # dla k2 i k3 -> wartości h dzielimy o połowę
    if k_nnumber in (2, 3):
        h_value = h_value / 2
    
    # zwrócenie rówaninia np. k3 = A + (h * k2)
    return A_initial_value + (h_value * kx_value)


# funkcja wykonująca skok i obliczająca ilość białka p53 po skoku
# parameters -> słownik parametrów k, d, p
def hop_result_p53(p53, NDMm, hop, parameters):

    # obliczanie kolejnych wartości k oraz nowym A
    # k1
    k1 = calculate_dev_p53(p53, NDMm, parameters)

    # k2
    new_p53 = A_value_update(hop, p53, k1, 2)
    k2 = calculate_dev_p53(new_p53, NDMm, parameters)
    
    # k3
    new_p53 = A_value_update(hop, p53, k2, 3)
    k3 = calculate_dev_p53(new_p53, NDMm, parameters)
    
    # k4
    new_p53 = A_value_update(hop, p53, k3, 4)
    k4 = calculate_dev_p53(new_p53, NDMm, parameters)

    # zwrócenie wartości równania na nową ilość p53 po skoku 
    return p53 + (hop / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


# funkcja wykonująca skok i obliczająca ilość białka NDMm po skoku
# parameters -> słownik parametrów k, d, p
def hop_result_NDMm(NDMm, NDMst, PTEN, hop, parameters):

    k1 = calculate_dev_NDMm(NDMst, PTEN, NDMm, parameters)

    new_NDMm = A_value_update(hop, NDMm, k1, 2)
    k2 = calculate_dev_NDMm(NDMst, PTEN, new_NDMm, parameters)

    new_NDMm = A_value_update(hop, NDMm, k2, 3)
    k3 = calculate_dev_NDMm(NDMst, PTEN, new_NDMm, parameters)
  
    new_NDMm = A_value_update(hop, NDMm, k3, 4)
    k4 = calculate_dev_NDMm(NDMst, PTEN, new_NDMm, parameters)

    # zwrócenie wartości równania na nową ilość NDMm po skoku 
    return NDMm + (hop / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


# funkcja wykonująca skok i obliczająca ilość białka NDMst po skoku
# parameters -> słownik parametrów k, d, p
def hop_result_NDMst(p53, NDMst, PTEN, hop, parameters):

    k1 = calculate_dev_NDMst(p53, PTEN, NDMst, parameters)

    new_NDMst = A_value_update(hop, NDMst, k1, 2)
    k2 = calculate_dev_NDMst(p53, PTEN, new_NDMst, parameters)

    new_NDMst = A_value_update(hop, NDMst, k2, 3)
    k3 = calculate_dev_NDMst(p53, PTEN, new_NDMst, parameters)
  
    new_NDMst = A_value_update(hop, NDMst, k3, 4)
    k4 = calculate_dev_NDMst(p53, PTEN, new_NDMst, parameters)

    
    # zwrócenie wartości równania na nową ilość NDMst po skoku
    return NDMst + (hop / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

# funkcja wykonująca skok i obliczająca ilość białka PTEN po skoku
# parameters -> słownik parametrów k, d, p
def hop_result_PTEN(p53, PTEN, hop, parameters):

    k1 = calculate_dev_PTEN(p53, PTEN, parameters)

    new_PTEN = A_value_update(hop, PTEN, k1, 2)
    k2 = calculate_dev_PTEN(p53, new_PTEN, parameters)

    new_PTEN = A_value_update(hop, PTEN, k2, 3)
    k3 = calculate_dev_PTEN(p53, new_PTEN, parameters)
  
    new_PTEN = A_value_update(hop, PTEN, k3, 4)
    k4 = calculate_dev_PTEN(p53, new_PTEN, parameters)

    # zwrócenie wartości równania na nową ilość PTEN po skoku
    return PTEN + (hop / 6) * (k1 + 2 * k2 + 2 * k3 + k4)


# funkcja odpowiadająca za obliczanie symulacji
# p53, NDMm, NDMst, PTEN -> białka
# hop -> odległość skoku
# time -> czas symulacji
# parameters -> parametry p, k, d
def calculate_symulation(p53, NDMm, NDMst, PTEN, hop, time, parameters):
    print("starting")
    # inicjacja macierzy z ilościami początkowymi białek
    p53_array = np.array([p53])
    NDMm_array = np.array([NDMm])
    NDMst_array = np.array([NDMst])
    PTEN_array = np.array([PTEN])

    # pętla przy stałym skoku
    # zaczynaj od wartości hop do czasu - time, skok o wartość hop
    for i in range(hop, time + 1, hop):

        # obliczenie wartości p53 po skoku i wpisanie jej do macierzy
        p53 = hop_result_p53(p53, NDMm, hop, parameters)
        p53_array = np.append(p53_array, p53)

        # obliczenie wartości NDMm po skoku i wpisanie jej do macierzy
        NDMm = hop_result_NDMm(NDMm, NDMst, PTEN, hop, parameters)
        NDMm_array = np.append(NDMm_array, NDMm)

        # obliczenie wartości NDMst po skoku i wpisanie jej do macierzy
        NDMst = hop_result_NDMst(p53, NDMst, PTEN, hop, parameters)
        NDMst_array = np.append(NDMst_array, NDMst)

        # obliczenie wartości PTEN po skoku i wpisanie jej do macierzy
        PTEN = hop_result_PTEN(p53, PTEN, hop, parameters)
        PTEN_array = np.append(PTEN_array, PTEN)

    # zwrócenie macierzy każej z białek
    return p53_array, NDMm_array, NDMst_array, PTEN_array


# Główna funkcja Symulacji Runge-Kutty IV
# podane parametry są wartościami domyślnymi
# można je zmienić podając nowe odpowiednio NDMm=50
# czas to 48 h = 60 s * 60 min * 24 h * 2 dni
# domyślenie PTEN jest aktywny, brak siRNA oraz nie ma uszkodzeń DNA
def Global_analysis(p53=1, NDMm=1, NDMst=1, PTEN=1, hop=10, time=172800):
    print("Global analysis starting!")
    first_simulations_set = np.array([])
    second_simulation_set = np.array([])
    third_simulation_set = np.array([])

    for i in range(50):
        print("starting set: ", i)
        # utworzenie słownika (HashMap) dla każdego z parametrów
        left_parameters = create_parameters_set_left()
        right_parameters = create_parameters_set_right()

        first_parameters_set = left_parameters

        second_parameters_set = {}
        second_parameters_set["p1"] = right_parameters["p1"]
        second_parameters_set["p2"] = left_parameters["p2"]
        second_parameters_set["p3"] = right_parameters["p3"]

        second_parameters_set["d1"] = right_parameters["d1"]
        second_parameters_set["d2"] = left_parameters["d2"]
        second_parameters_set["d3"] = right_parameters["d3"]

        second_parameters_set["k1"] = right_parameters["k1"]
        second_parameters_set["k2"] = right_parameters["k2"]
        second_parameters_set["k3"] = right_parameters["k3"]

        third_parameters_set = {}
        third_parameters_set["p1"] = left_parameters["p1"]
        third_parameters_set["p2"] = right_parameters["p2"]
        third_parameters_set["p3"] = left_parameters["p3"]

        third_parameters_set["d1"] = left_parameters["d1"]
        third_parameters_set["d2"] = right_parameters["d2"]
        third_parameters_set["d3"] = left_parameters["d3"]

        third_parameters_set["k1"] = left_parameters["k1"]
        third_parameters_set["k2"] = left_parameters["k2"]
        third_parameters_set["k3"] = left_parameters["k3"]


        simulation_first = {}
        simulation_first["p53_array"], simulation_first["NDMm_array"], simulation_first["NDMst_array"], simulation_first["PTEN_array"] = calculate_symulation(p53=p53, NDMm=NDMm, NDMst=NDMst, PTEN=PTEN, 
                                                                            hop=hop, time=time, parameters=first_parameters_set)
        first_simulations_set = np.append(first_simulations_set, simulation_first)


        simulation_second = {}
        simulation_second["p53_array"], simulation_second["NDMm_array"], simulation_second["NDMst_array"], simulation_second["PTEN_array"] = calculate_symulation(p53=p53, NDMm=NDMm, NDMst=NDMst, PTEN=PTEN, 
                                                                            hop=hop, time=time, parameters=second_parameters_set)
        second_simulation_set = np.append(second_simulation_set, simulation_second)


        simulation_third = {}
        simulation_third["p53_array"], simulation_third["NDMm_array"], simulation_third["NDMst_array"], simulation_third["PTEN_array"] = calculate_symulation(p53=p53, NDMm=NDMm, NDMst=NDMst, PTEN=PTEN, 
                                                                            hop=hop, time=time, parameters=third_parameters_set)
        third_simulation_set = np.append(third_simulation_set, simulation_third)

    print("Global analysis step1")

    total_lengh = len(first_simulations_set[0]["p53_array"])

    p53_sum_array = np.zeros(total_lengh)
    NDMm_sum_array = np.zeros(total_lengh)
    NDMst_sum_array = np.zeros(total_lengh)
    PTEN_sum_array = np.zeros(total_lengh)

    for i in range(len(first_simulations_set)):

        p53_array = first_simulations_set[i]["p53_array"]
        NDMm_array = first_simulations_set[i]["NDMm_array"]
        NDMst_array = first_simulations_set[i]["NDMst_array"]
        PTEN_array = first_simulations_set[i]["PTEN_array"]

        for i in range(total_lengh):
            p53_sum_array[i] += p53_array[i]
            NDMm_sum_array[i] += NDMm_array[i]
            NDMst_sum_array[i] += NDMst_array[i]
            PTEN_sum_array[i] += PTEN_array[i]

    # y
    p53_in_time = p53_sum_array/ len(first_simulations_set)
    NDMm_in_time = NDMm_sum_array/ len(first_simulations_set)
    NDMst_in_time = NDMst_sum_array/ len(first_simulations_set)
    PTEN_in_time = PTEN_sum_array/ len(first_simulations_set)
    print("Global analysis step1 done!")


    print("Global analysis step2")

    p53_sum_array = np.zeros(total_lengh)
    NDMm_sum_array = np.zeros(total_lengh)
    NDMst_sum_array = np.zeros(total_lengh)
    PTEN_sum_array = np.zeros(total_lengh)

    for i in range(len(first_simulations_set)):

        p53_array = first_simulations_set[i]["p53_array"]
        NDMm_array = first_simulations_set[i]["NDMm_array"]
        NDMst_array = first_simulations_set[i]["NDMst_array"]
        PTEN_array = first_simulations_set[i]["PTEN_array"]

        for i in range(total_lengh):
            p53_sum_array[i] = p53_sum_array[i] + (p53_array[i]**2)
            NDMm_sum_array[i] = NDMm_sum_array[i] + (NDMm_array[i]**2)
            NDMst_sum_array[i] = NDMst_sum_array[i] + (NDMst_array[i]**2)
            PTEN_sum_array[i] = PTEN_sum_array[i] + (PTEN_array[i]**2)

    # D + y^2
    p53_squeared_in_time = p53_sum_array/ len(first_simulations_set)
    NDMm_squeared_in_time = NDMm_sum_array/ len(first_simulations_set)
    NDMst_squeared_in_time = NDMst_sum_array/ len(first_simulations_set)
    PTEN_squeared_in_time = PTEN_sum_array/ len(first_simulations_set)

    the_D_p53 = p53_squeared_in_time - (p53_in_time**2)
    the_D_NDMm = NDMm_squeared_in_time - (NDMm_in_time**2)
    the_D_NDMst = NDMst_squeared_in_time - (NDMst_in_time**2)
    the_D_PTEN = PTEN_squeared_in_time - (PTEN_in_time**2)

    print("Global analysis step2 done!")

    print("Global analysis step3")

    p53_sum_array = np.zeros(total_lengh)
    NDMm_sum_array = np.zeros(total_lengh)
    NDMst_sum_array = np.zeros(total_lengh)
    PTEN_sum_array = np.zeros(total_lengh)

    for i in range(len(first_simulations_set)):

        p53_array = first_simulations_set[i]["p53_array"]
        NDMm_array = first_simulations_set[i]["NDMm_array"]
        NDMst_array = first_simulations_set[i]["NDMst_array"]
        PTEN_array = first_simulations_set[i]["PTEN_array"]

        prim2_p53_array = second_simulation_set[i]["p53_array"]
        prim2_NDMm_array = second_simulation_set[i]["NDMm_array"]
        prim2_NDMst_array = second_simulation_set[i]["NDMst_array"]
        prim2_PTEN_array = second_simulation_set[i]["PTEN_array"]



        for i in range(total_lengh):
            p53_sum_array[i] = p53_sum_array[i] + (p53_array[i] * prim2_p53_array[i])
            NDMm_sum_array[i] = NDMm_sum_array[i] + (NDMm_array[i] * prim2_NDMm_array[i])
            NDMst_sum_array[i] = NDMst_sum_array[i] + (NDMst_array[i] * prim2_NDMst_array[i])
            PTEN_sum_array[i] = PTEN_sum_array[i] + (PTEN_array[i] * prim2_PTEN_array[i])
        
    p53_prim2_in_time = p53_sum_array/ len(first_simulations_set)
    NDMm_prim2_in_time = NDMm_sum_array/ len(first_simulations_set)
    NDMst_prim2_in_time = NDMst_sum_array/ len(first_simulations_set)
    PTEN_prim2_in_time = PTEN_sum_array/ len(first_simulations_set)

    the_Dy_p53 = p53_prim2_in_time - (p53_in_time**2)
    the_Dy_NDMm = NDMm_prim2_in_time - (NDMm_in_time**2)
    the_Dy_NDMst = NDMst_prim2_in_time - (NDMst_in_time**2)
    the_Dy_PTEN = PTEN_prim2_in_time - (PTEN_in_time**2)
    print("Global analysis step3 done")

    print("Global analysis step4")
    p53_sum_array = np.zeros(total_lengh)
    NDMm_sum_array = np.zeros(total_lengh)
    NDMst_sum_array = np.zeros(total_lengh)
    PTEN_sum_array = np.zeros(total_lengh)

    for i in range(len(first_simulations_set)):

        p53_array = first_simulations_set[i]["p53_array"]
        NDMm_array = first_simulations_set[i]["NDMm_array"]
        NDMst_array = first_simulations_set[i]["NDMst_array"]
        PTEN_array = first_simulations_set[i]["PTEN_array"]

        prim1_p53_array = third_simulation_set[i]["p53_array"]
        prim1_NDMm_array = third_simulation_set[i]["NDMm_array"]
        prim1_NDMst_array = third_simulation_set[i]["NDMst_array"]
        prim1_PTEN_array = third_simulation_set[i]["PTEN_array"]



        for i in range(total_lengh):
            p53_sum_array[i] = p53_sum_array[i] + (p53_array[i] * prim1_p53_array[i])
            NDMm_sum_array[i] = NDMm_sum_array[i] + (NDMm_array[i] * prim1_NDMm_array[i])
            NDMst_sum_array[i] = NDMst_sum_array[i] + (NDMst_array[i] * prim1_NDMst_array[i])
            PTEN_sum_array[i] = PTEN_sum_array[i] + (PTEN_array[i] * prim1_PTEN_array[i])
        
    p53_prim1_in_time = p53_sum_array/ len(first_simulations_set)
    NDMm_prim1_in_time = NDMm_sum_array/ len(first_simulations_set)
    NDMst_prim1_in_time = NDMst_sum_array/ len(first_simulations_set)
    PTEN_prim1_in_time = PTEN_sum_array/ len(first_simulations_set)

    the_Dz_p53 = p53_prim1_in_time - (p53_in_time**2)
    the_Dz_NDMm = NDMm_prim1_in_time - (NDMm_in_time**2)
    the_Dz_NDMst = NDMst_prim1_in_time - (NDMst_in_time**2)
    the_Dz_PTEN = PTEN_prim1_in_time - (PTEN_in_time**2)
    print("Global analysis step4 done!")

    print("Final calculations!")
    Dtot_p53_p1 = the_D_p53 - the_Dz_p53
    S_p53_p1 = the_Dy_p53 / the_D_p53
    Stot_p53_p1 = Dtot_p53_p1 / the_D_p53

    Dtot_NDMm_p1 = the_D_NDMm - the_Dz_NDMm
    S_NDMm_p1 = the_Dy_NDMm / the_D_NDMm
    Stot_NDMm_p1 = Dtot_NDMm_p1 / the_D_NDMm

    Dtot_NDMst_p1 = the_D_NDMst - the_Dz_NDMst
    S_NDMst_p1 = the_Dy_NDMst / the_D_NDMst
    Stot_NDMst_p1 = Dtot_NDMst_p1 / the_D_NDMst

    Dtot_PTEN_p1 = the_D_PTEN - the_Dz_PTEN
    S_PTEN_p1 = the_Dy_PTEN / the_D_PTEN
    Stot_PTEN_p1 = Dtot_PTEN_p1 / the_D_PTEN
    print("Calculations Done!")
    print(S_p53_p1, S_NDMm_p1, S_NDMst_p1, S_PTEN_p1)

    plt.plot(S_p53_p1, label="p53", color="#cc6699")
    plt.ylabel("Protein value")
    plt.xlabel("Time [s]")
    plt.legend(loc="upper left")
    plt.show()

    plt.plot(S_NDMm_p1, label="p53", color="#cc6699")
    plt.ylabel("Protein value")
    plt.xlabel("Time [s]")
    plt.legend(loc="upper left")
    plt.show()

    plt.plot(S_NDMst_p1, label="p53", color="#cc6699")
    plt.ylabel("Protein value")
    plt.xlabel("Time [s]")
    plt.legend(loc="upper left")
    plt.show()

    plt.plot(S_PTEN_p1, label="p53", color="#cc6699")
    plt.ylabel("Protein value")
    plt.xlabel("Time [s]")
    plt.legend(loc="upper left")
    plt.show()

Global_analysis()

    

