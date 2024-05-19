
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


