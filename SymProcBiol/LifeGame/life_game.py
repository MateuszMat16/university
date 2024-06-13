import numpy as np
import pygame as pg
net = np.zeros([50, 50])

def get_neighbourhood_coordinates(x, y):
    coors = np.array([999, 999])
    x_coor = np.array([x - 1, x, x + 1])
    y_coor = np.array([y - 1, y, y + 1])

    for x_pos in x_coor:
        for y_pos in y_coor:

            if x == x_pos and y_pos == y:
                continue
            coors = np.vstack([coors, [x_pos, y_pos]])
    coors = np.delete(coors, 0, 0)
    return coors


def cnc_cell_status(x, y, net):
    
    neighbourhood = get_neighbourhood_coordinates(x, y)
    result = net[x, y]
    counter = 0
    for neigh in neighbourhood:
    
        if counter > 3:
            break
        counter = counter + net[neigh[0], neigh[1]]
    

    if net[x, y] == 1:
        
        if counter < 2:
            result = 0
        elif counter in (2, 3):
            pass
        elif counter > 3:
            result = 0
    
    elif net[x, y] == 0 and counter == 3:
        result = 1
    
    return result


def new_period(net, max_x, max_y):

    new_net = net
    for x in range (1 ,max_x - 1):
        for y in range(1, max_y - 1):
            new_net[x, y] = cnc_cell_status(x, y, net)
    return new_net

def calculate_simulation(iter_number, coordinates):
    net = np.zeros([52, 52])

    for coor in coordinates:
        if coor[0] < 51 and coor[1] < 51:
            net[coor[0], coor[1]] = 1
    
    net_arr = np.array(net)
    for i in range(iter_number):
        tmp = new_period(net_arr[-1], 51, 51)
        net_arr = np.append(net_arr, tmp)

    return net_arr    



def makegame(net_arr):
    print(net_arr)
    pg.init()

    # ustawienie kolorów
    white,black = (255,255,255),(0,0,0)

    #set display
    gameDisplay = pg.display.set_mode((1050,1050))

    #caption
    pg.display.set_caption("Life Game")

    #Size of squares
    size = 20


    gameDisplay.fill(white)


    for net in net_arr:
        print(net)
        for i in range(1,51):
            for z in range(1,51):
                if net[z, i] == 1:
                    pg.draw.rect(gameDisplay, black,[size*z,size*i,size,size])
                pg.draw.rect(gameDisplay, black,[size*z,size*i,size,size], 1)

        pg.display.update()


    #beginning of logic
    game_exit = False


    while not game_exit:
        for event in pg.event.get():
            if event.type == pg.QUIT:
                game_exit = True


    #quit from pygame & python
    pg.quit()


makegame(calculate_simulation(5, [[1, 2], [2, 2]]))