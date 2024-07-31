import sys
import pygame as pg
import time
import numpy as np
pg.init()

# utworzenie podstawowych zmiennych i inicjalizacji wartości związanych z grą
width, height = 1000, 920
fps = 10
fpsClock = pg.time.Clock()
screen = pg.display.set_mode((width, height))
font = pg.font.SysFont('Arial', 20)
objects = []

# utworzenie matrycy służącej do symulacji poszczególnych epok
game_playground = np.zeros([52, 52])

# funkcja zwracająca koordynaty sąsiednich komórek
def get_neighbourhood_coordinates(x, y, x_max, y_max):
    coors = np.array([999, 999])
    x_coor = np.array([x - 1, x, x + 1])
    y_coor = np.array([y - 1, y, y + 1])

    for x_pos in x_coor:
        for y_pos in y_coor:

            if x == x_pos and y_pos == y:
                continue
            if x_pos < 0 or y_pos < 0:
                continue
            if x_pos > x_max or y_pos > y_max:
                continue

            coors = np.vstack([coors, [x_pos, y_pos]])
    coors = np.delete(coors, 0, 0)
    return coors

# funkcja określająca stan komórki po przejściu w następną epokę
def cnc_cell_status(x, y, net):

    # inicjalizacja zmiennych
    x_max = np.shape(net)[0]
    y_max = np.shape(net)[1]
    neighbourhood = get_neighbourhood_coordinates(x, y, x_max -1, y_max- 1)
    result = net[x, y]
    counter = 0

    # sprawdzanie ilości żywych komórek w sąsiedzctwie
    for neigh in neighbourhood:
    
        if counter > 3:
            break
        # counter = counter + net[neigh[0], neigh[1]]
        if net[neigh[0], neigh[1]] == 1:

            counter += 1

    # warunki dla już żywej komórki
    if net[x, y] == 1:
        
        if counter < 2:
            result = 0
        elif counter in (2, 3):
            result = 1
        elif counter > 3:
            result = 0
    
    # warunek dla nieżywej komórki
    elif net[x, y] == 0 and counter == 3:
        result = 1
    
    return result

# funkcja określająca stan układu po przejściu do nowej epoki
def new_period(net, max_x, max_y):

    new_net = np.copy(net)
    for x in range (0 ,max_x):

        for y in range(0, max_y):
            new_net[x, y] = cnc_cell_status(x, y, net)
                
            if cnc_cell_status(x, y, net) == 1:
                print("coordinates: ", x, y)
    
    new_net = checkBorder(new_net, max_x, max_y)
    return new_net


# warunek brzegowy periodyczny
def checkBorder(net, max_x, max_y):
    new_net = np.copy(net)
    for i in range(0, max_x):
        for j in range(0, max_y):

            if new_net[i, j] == 1:
                
                if i == 0 and j == 0:
                    new_net[i, j] = 0
                    new_net[max_x - 2, max_y - 2] = 1
                    print("Found 1")
               
                elif i == max_x - 1 and j == max_y - 1:
                    new_net[i, j] = 0
                    new_net[1, 1] = 1
                    print("Found 2")

                elif i == 0:
                    new_net[i, j] = 0
                    new_net[max_x - 2, j] = 1
                    print("Found 3")

                elif i == max_x - 1:
                    new_net[i, j] = 0
                    new_net[1, j] = 1
                    print("Found 4")

                elif j == 0:
                    new_net[i, j] = 0
                    new_net[i, max_y - 2] = 1
                    print("Found 5")

                elif j == max_y - 1:
                    new_net[i, j] = 0
                    new_net[i, 1] = 1
                    print("Found 6")

    return new_net




# druga część

# klasa przycisku z podstawowymi wartościami
# na jej podstawie będą tworzone klasy przycisków funkcyjnych
class Button():

    # inicjalizowanie wartości w konstruktorze
    def __init__(self, x, y, width, height, 
                 buttonText = "Button", buttonColor="#ffffff"):
        print("Creating new buttom")
        self.x = x
        self.y = y
        self.witdth = width
        self.height = height
        self.buttonText = buttonText
        self.buttonColor = buttonColor
        self.buttonSurface = pg.Surface((self.witdth, self.height))
        self.buttonRect = pg.Rect(self.x, self.y, self.witdth, self.height)
        self.buttonSurf = font.render(self.buttonText, True, (20, 20, 20))
        self.buttonSurface.fill(self.buttonColor)
        # dodanie obiektu do listy obiektów do iteracji
        objects.append(self)

    
    # funkcja wykonywana w każdej iteracji podczas trwania programu
    def process(self):
        #self.buttonSurface.fill("#ffffff")
        self.buttonSurface.blit(self.buttonSurf, [
            self.buttonRect.width/2 - self.buttonSurf.get_rect().width/2,
            self.buttonRect.height/2 - self.buttonSurf.get_rect().height/2
        ])

        screen.blit(self.buttonSurface, self.buttonRect)
       

# klasa przycisku - komórki
# używana do symulacji komórki na planszy
class GameButton(Button):
    
    def __init__(self, x, y, width, height, buttonText="Button", 
                 buttonColor="#ffffff", coor1 = 0, coor2 = 0, value=False):
        super().__init__(x, y, width, height, buttonText, buttonColor)
        # dodatkowo dodane koordynaty związane z siatką
        self.coor1 = coor1
        self.coor2 = coor2
        self.value = value

        self.gameButtonColors = {
            "default": self.buttonColor,
            "pressed": "#0000ff"
        }


    def setButtonColor(self, newColor):
        self.buttonColor = newColor
    

    def getButtonColor(self):
        return self.buttonColor


    def getCoor1(self):
        return self.coor1
    

    def getCoor2(self):
        return self.coor2
    
    def getValue(self):
        return self.value
    
  

    # metoda zmiany wartości
    def setValue(self, newValue):
        self.value = newValue
        if newValue == False:
            self.setButtonColor(self.gameButtonColors["default"])
            self.buttonSurface.fill(self.buttonColor)
         
        else:
            self.setButtonColor(self.gameButtonColors["pressed"])
            self.buttonSurface.fill(self.buttonColor)

    # metoda iteracyjna
    def process(self):
      super().process()
      mousePos = pg.mouse.get_pos()
      if pg.mouse.get_pressed(num_buttons=3)[0] and self.buttonRect.collidepoint(mousePos):
        
        if (self.buttonColor == self.gameButtonColors["default"]):  
            self.setButtonColor(self.gameButtonColors["pressed"])
            self.buttonSurface.fill(self.buttonColor)
            self.value = not self.value


        else:
            self.setButtonColor(self.gameButtonColors["default"])
            self.buttonSurface.fill(self.buttonColor)
            self.value = not self.value
        
        time.sleep(0.2)

# klasa przycisku zaczynającego symulację  
class StartButton(Button):

    def __init__(self, x, y, width, height, buttonText="Button",
                  buttonColor="#ffffff", gameStarted=False):
        super().__init__(x, y, width, height, buttonText, buttonColor)
        # zmienna określająca, czy gram ma się rozpocząć
        self.gameStarted = gameStarted

    def setButtonColor(self, newColor):
        self.buttonColor = newColor

    def process(self):
        super().process()
        # iteracyjne sprawdzanie, czy gra się nie rozpoczęła
        if self.gameStarted:
            # rozpoczynanie gry
            self.processGame()
         

        mousePos = pg.mouse.get_pos()
        if pg.mouse.get_pressed(num_buttons=3)[0] and self.buttonRect.collidepoint(mousePos):
            self.startGame()
            self.setButtonColor("#0000ff")
            self.buttonSurface.fill(self.buttonColor)
            print("Game started!")
            time.sleep(0.2)

    def startGame(self):
        self.gameStarted = True

    # metoda zmieniająca pola komórek w grze
    def processGame(self):
        self.readGameStatus()
        self.changeGameStatus()
        time.sleep(0.5)

           
    # metoda zczytująca wartości z aktualnego stanu 
    def readGameStatus(self):
        for row in gameButtonArray:
            for button in row:
                tmp = button.getValue()
                if tmp:
                    game_playground[button.getCoor1() + 1, button.getCoor2() + 1] = 1

                else:
                    game_playground[button.getCoor1() + 1, button.getCoor2() + 1] = 0


    # metoda określająca i zmieniające epokę w grze
    def changeGameStatus(self):
        new_net = new_period(game_playground, 52, 52)
        myCounter.update()
        for i in range(1, 51):
                for j in range(1, 51):
                    if new_net[i, j] == 1:
                        gameButtonArray[i - 1][j - 1].setValue(True)
                        
                    else:
                        gameButtonArray[i - 1][j - 1].setValue(False)


# klasa licznika - przycisku
# główne zdanie - liczenie epok
class CounterButton(Button):
    def __init__(self, x, y, width, height, buttonText="Button", buttonColor="#ffffff", counter=0):
        super().__init__(x, y, width, height, buttonText, buttonColor)
        self.counter = counter
    

    def process(self):
        super().process()
    
    # metoda wywoływana przy zmianie epoki
    # zwiększa wartość licznika o 1
    def update(self):
        self.counter += 1
        self.buttonSurface = pg.Surface((self.witdth, self.height))
        self.buttonSurface.fill("#ffffff")
        self.buttonSurf = font.render(str(self.counter), True, (20, 20, 20))


# inicjowanie przycisków nie należących do komórek
StartButton(20, 20, 120, 60, "Start")
myCounter = CounterButton(200, 20, 80, 60, "")

# tworzenie macierzy dla przycisków komórek
gameButtonArray = []

# tworzenie i dodawanie do macierzy przycisków komórek
# kolumny
for i in range(50):
    gameButtomRow = []
    # rzędy
    for j in range(50):
        tmp = GameButton(100 + 16 * i, 100 + 16 * j, 15, 15, "", "#ffffff", i, j)
        gameButtomRow.append(tmp)
    
    gameButtonArray.append(gameButtomRow)



# Pętla gry
while True:
    # wypełnianie tła kolorem
    screen.fill((20, 20, 20))
    # zakończenie programu
    for event in pg.event.get():
        if event.type == pg.QUIT:
            pg.quit()
            sys.exit()

    # iteracyjne procesowanie obiektów
    for object in objects:
        object.process()

    pg.display.flip()
    fpsClock.tick(fps)