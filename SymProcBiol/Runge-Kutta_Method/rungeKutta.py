import numpy as np
import matplotlib.pyplot as plt

# Funkcja obliczająca pochodną dla p53
def calculate_dev_p53(p53_initial_value, NDMm_value, parameters):
    
    # przypisanie wartości parametrów za pomocą słownika  
    d1 = parameters["d1"]
    p1 = parameters["p1"]

    # obliczenie drugiej części pochodnej
    second_part = d1 * p53_initial_value * NDMm_value

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

    # inicjacja macierzy z ilościami początkowymi białek
    p53_array = np.array([p53])
    NDMm_array = np.array([NDMm])
    NDMst_array = np.array([NDMst])
    PTEN_array = np.array([PTEN])

    # pętla przy stałym skoku
    # zaczynaj od wartości hop do czasu - time, skok o wartość hop
    for i in range(hop, time, hop):

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
def RK_method(p53=100, NDMm=100, NDMst=100, PTEN=100, hop=10, time=17280, PTEN_off=False, is_siRNA=False, DNA_damage=False):

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
        print("Parameters: PTEN is off")
        parameters["p3"] = 0.0
    else:
        print("Parameters: PTEN is on")
        parameters["p3"] = 100

    # warunek na obecność siRNA, zmienia wartość paramteru p2
    if is_siRNA:
        print("Parameters: siRNA")
        parameters["p2"] = 0.02
    else:
        print("Parameters: no siRNA")
        parameters["p2"] = 440
  
    # warunek sprawdza czy w mamy uszkodzenia DNA i odpowiednio zmienia wartość parametru d2
    if DNA_damage:
        print("Parameters: DNA damage")
        parameters["d2"] = 1.375e-4
    else:
        print("Parameters: no DNA damage")
        parameters["d2"] = 0.1

    # Wywołanie funkcji wykonującej całą symulację
    # przekazanie wszytkich danych wejściowych białek oraz parametrów
    p53_array, NDMm_array, NDMst_array, PTEN_array = calculate_symulation(p53=p53, NDMm=NDMm, NDMst=NDMst, PTEN=PTEN, 
                                                                          hop=hop, time=time, parameters=parameters)
    # utworzenie macierzy czasu
    time = np.array([x for x in range(0, time, hop)])

    # rysowanie funkcji ilości białek w czasie
    plt.plot(p53_array, label="p53", color="#0033cc")
    plt.plot(NDMm_array, label="NDMm", color="#ffff00")
    plt.plot(NDMst_array, label="NDMst", color="#003300")
    plt.plot(PTEN_array, label="PTEN", color="#cc6699")
    plt.ylabel("P53 value")
    plt.xlabel("Time")
    plt.legend(loc="upper left")
    plt.show()

    # Zapis wyników do pliku
    with open("./university/SymProcBiol/Runge-Kutta_Method/result.txt", "w") as file:
        line = "p53 \t NDMm \t NDMst \t PTEN \n"
        file.write(line)
        for i in range(len(PTEN_array)):
            line = (str(p53_array[i]) + "\t" + str(NDMm_array[i]) + "\t" + str(NDMst_array[i]) + "\t" + str(PTEN_array[i]) + "\n")
            file.write(line)


# funkcja okre                 
def calculate_hop_length_p53(max_percent_change, hop, p53, NDMm, parameters):
    after_hop_p53 = hop_result_p53(p53, NDMm, hop, parameters)
    p53_diff = (abs(after_hop_p53 - p53) * 100) / p53
    
    if (p53 != after_hop_p53) and (p53_diff > max_percent_change):
        
        if hop/2 >= 0.75:
            return calculate_hop_length_p53(max_percent_change=max_percent_change, hop=hop/2, p53=p53, NDMm=NDMm, parameters=parameters)
        else:
            return hop

    else:
        return hop

def calculate_hop_length_p53(max_percent_change, hop, p53, NDMm, parameters):
    after_hop_p53 = hop_result_p53(p53, NDMm, hop, parameters)
    p53_diff = (abs(after_hop_p53 - p53) * 100) / p53
    
    if (p53 != after_hop_p53) and (p53_diff > max_percent_change):
        
        if hop/2 >= 0.75:
            return calculate_hop_length_p53(max_percent_change=max_percent_change, hop=hop/2, p53=p53, NDMm=NDMm, parameters=parameters)
        else:
            return hop

    else:
        return hop



def RK_method_V2(p53=100, NDMm=100, NDMst=100, PTEN=100, hop=10, time=17280, PTEN_off=False, is_siRNA=False, DNA_damage=False):

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
        print("Parameters: PTEN is off")
        parameters["p3"] = 0.0
    else:
        print("Parameters: PTEN is on")
        parameters["p3"] = 100

    # warunek na obecność siRNA, zmienia wartość paramteru p2
    if is_siRNA:
        print("Parameters: siRNA")
        parameters["p2"] = 0.02
    else:
        print("Parameters: no siRNA")
        parameters["p2"] = 440
  
    # warunek sprawdza czy w mamy uszkodzenia DNA i odpowiednio zmienia wartość parametru d2
    if DNA_damage:
        print("Parameters: DNA damage")
        parameters["d2"] = 1.375e-4
    else:
        print("Parameters: no DNA damage")
        parameters["d2"] = 0.1

    counter = 0
    p53_array = np.array([p53])
    p53_time_array = np.array([counter])
    while counter < time:

        new_hop =  calculate_hop_length_p53(max_percent_change=5, hop=10, p53=p53_array[-1], NDMm=NDMm, parameters=parameters)
        counter += new_hop
        p53_time_array = np.append(p53_time_array, counter)
        p53_array = np.append(p53_array, hop_result_p53(p53_array[-1], NDMm, new_hop, parameters))


    print(p53_array)
    print(p53_time_array)
    # rysowanie funkcji ilości białek w czasie
    plt.plot(p53_array, p53_time_array, label="p53", color="#0033cc")
   
    plt.ylabel("P53 value")
    plt.xlabel("Time")
    plt.legend(loc="upper left")
    plt.show()
  


# wywołanie funkcji
RK_method(hop=10, time=17280 * 3)
RK_method_V2(hop=10, time=17280 * 3)
