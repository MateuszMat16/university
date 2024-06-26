import numpy as np
import matplotlib.pyplot as plt
import random


# obliczanie skołonności poszczególnych reakcji

# tworzenie p53
def creation_p53(parameters):
    p1 = parameters["p1"]
    return p1


# degradacja p53
def destruction_p53(p53, NDMm, parameters):
    d1 = parameters["d1"]

    result = d1 * p53 * pow(NDMm, 2)
    return result


# tworzenie NDMst
def creation_NDMst(p53, parameters):
    p2 = parameters["p2"]
    k2 = parameters["k2"]

    result = p2 * ((p53**4) / (pow(p53, 4) + pow(k2, 4)))
    return result

# przemiana NDMst w NDMm
def change_NDMst_to_NDMm(NDMst, PTEN, parameters):
    k1 = parameters["k1"]
    k3 = parameters["k3"]

    result = k1 * (pow(k3, 2) /(pow(k3, 2) + pow(PTEN, 2))) * NDMst
    return result


# degradacja NDMst
def destruction_NDMst(NDMst, parameters):
    d2 = parameters["d2"]

    result = d2 * NDMst
    return result


# degradacja NDMm
def destruction_NDMm(NDMm, parameters):
    d2 = parameters["d2"]

    result = d2 * NDMm
    return result


# tworzenie PTEN
def creation_PTEN(p53, parameters):
    p3 = parameters["p3"]
    k2 = parameters["k2"]

    result = p3 * (pow(p53, 4) / (pow(p53, 4) + pow(k2, 4)))
    return result


# degradacja PTEN
def destruction_PTEN(PTEN, parameters):
    d3 = parameters["d3"]

    result = d3 * PTEN
    return result


# wyliczanie czasu tau
def calculate_tau(a_value):

    # wyjątek dzielenie przez zero
    if a_value == 0:
        return 10000
    
    # generates random number between 0.0 and 1
    random_number = random.random()
    result = (-1) * (1/a_value) * np.log(random_number)
    return round(result, 4) # zaokrąglenie do 4 miejsca po przecinku



# wyszukiwanie binarne najmniejszego tau
# funkcja zwraca dodatkowo pozycję w tablicy zwycięzcy
def find_min_tau(hop_arr):
    min = hop_arr[0]
    position = 0
    for i in range(1, len(hop_arr)):
        if min > hop_arr[i]:
            position = i
            min = hop_arr[i]
    
    return min, position


# Metoda Next Hop
# funkcja wykonuje jeden skok metodą pierwszej reakcji
def calculate_next_hop(p53, NDMm, NDMst, PTEN, parameters):
    
    a1 = creation_p53(parameters) 
    a2 = destruction_p53(p53, NDMm, parameters)
    a3 = creation_NDMst(p53, parameters)
    a4 = change_NDMst_to_NDMm(NDMst, PTEN, parameters)
    a5 = destruction_NDMst(NDMst, parameters)
    a6 = destruction_NDMm(NDMm, parameters)
    a7 = creation_PTEN(p53, parameters)
    a8 = destruction_PTEN(PTEN, parameters)

    a_values = np.array([a1, a2, a3, a4, a5, a6, a7, a8])
    tau_values = np.array([])
    
    for i in range(len(a_values)):
        tau_values = np.append(tau_values, calculate_tau(a_values[i]))

    tau, counter = find_min_tau(tau_values)


    # aktualnie counter jest pozycją w tablicy zaczynającej się od 0
    # uaktualniamy wartość w celu otrzymania numeru reakcji
    counter += 1

    if counter == 1:
        p53 += 1
    elif counter == 2:
        p53 -= 1
    elif counter == 3:
        NDMst += 1
    elif counter == 4:
        NDMst -= 1
        NDMm += 1
    elif counter == 5:
        NDMst -= 1
    elif counter == 6:
        NDMm -= 1
    elif counter == 7:
        PTEN += 1
    elif counter == 8:
        PTEN -= 1
    else:
        print("Something went wrong")

    return p53, NDMm, NDMst, PTEN, tau


# tworzenie symulacji dla metody następnej reakcji
def calculate_simulation_next_hop(p53=100, NDMm=100, NDMst=100, PTEN=100, time=17280, PTEN_off=False, is_siRNA=False, DNA_damage=False):

    # naprawia błąd przy którym liczba nie mieści się int32, więc zaczyna przyjmować ujemne wartości
    p53 = np.uint64(p53)
    
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

    # inicjowanie tablic z wartościami białek i czasu
    counter = 0
    p53_array = np.array([p53])
    NDMm_array = np.array([NDMm])
    NDMst_array = np.array([NDMst])
    PTEN_array = np.array([PTEN])
    time_array = np.array([counter])

    # warunek końca symulacji
    while time > counter:
        print(counter)
        new_p53, new_NDMm, new_NDMst, new_PTEN, tau = calculate_next_hop(p53_array[-1], NDMm_array[-1], NDMst_array[-1], PTEN_array[-1], parameters)
        print("tau", tau)
        counter += tau
        time_array = np.append(time_array, counter)
        p53_array = np.append(p53_array, new_p53)
        NDMm_array = np.append(NDMm_array, new_NDMm)
        NDMst_array = np.append(NDMst_array, new_NDMst)
        PTEN_array = np.append(PTEN_array, new_PTEN)
    
    

    # rysowanie funkcji ilości białek w czasie
    plt.plot(time_array, p53_array, label="p53", color="#0033cc")
    plt.ylabel("Protein value")
    plt.xlabel("Time [s]")
    plt.legend(loc="upper left")
    plt.savefig("p53_Gil_next_hop.png")
    plt.close()

    plt.plot(time_array, NDMm_array, label="NDMm", color="#ffff00")
    plt.ylabel("Protein value")
    plt.xlabel("Time [s]")
    plt.legend(loc="upper left")
    plt.savefig("NDMm_Gil_next_hop.png")
    plt.close()

    plt.plot(time_array,NDMst_array,  label="NDMst", color="#003300")
    plt.ylabel("Protein value")
    plt.xlabel("Time [s]")
    plt.legend(loc="upper left")
    plt.savefig("NDMst_Gil_next_hop.png")
    plt.close()

    plt.plot(time_array, PTEN_array, label="PTEN", color="#cc6699")
    plt.ylabel("Protein value")
    plt.xlabel("Time [s]")
    plt.legend(loc="upper left")
    plt.savefig("PTEN_Gil_next_hop.png")
    plt.close()
    # plt.savefig("RKIV_changing_step.png")
    # plt.close()

     # Zapis wyników do pliku
    with open("./result_next_hop.txt", "w") as file:
        line = "p53 \t NDMm \t NDMst \t PTEN \t time \n"
        file.write(line)
        for i in range(len(PTEN_array)):
            line = (str(p53_array[i]) + "\t" + str(NDMm_array[i]) + "\t" + str(NDMst_array[i]) + "\t" + str(PTEN_array[i]) + "\t" + str(time_array[i]) + "\n")
            file.write(line)


calculate_simulation_next_hop(time=2880)