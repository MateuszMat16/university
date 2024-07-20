import sys
import pygame as pg

pg.init()

width, height = 850, 700
fps = 10
fpsClock = pg.time.Clock()
screen = pg.display.set_mode((width, height))
font = pg.font.SysFont('Arial', 20)
objects = []



class Button():

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
        objects.append(self)
       
    def process(self):
        #self.buttonSurface.fill("#ffffff")
        self.buttonSurface.blit(self.buttonSurf, [
            self.buttonRect.width/2 - self.buttonSurf.get_rect().width/2,
            self.buttonRect.height/2 - self.buttonSurf.get_rect().height/2
        ])

        screen.blit(self.buttonSurface, self.buttonRect)
       
        
class GameButton(Button):
    
    def __init__(self, x, y, width, height, buttonText="Button", 
                 buttonColor="#ffffff", coor1 = 0, coor2 = 0, value=False):
        super().__init__(x, y, width, height, buttonText, buttonColor)
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
    
  


    def setValue(self, newValue):
        self.value = newValue
        if newValue == False:
            self.setButtonColor(self.gameButtonColors["default"])
            self.buttonSurface.fill(self.buttonColor)
         
        else:
            print("eesfsf")
            self.setButtonColor(self.gameButtonColors["pressed"])
            self.buttonSurface.fill(self.buttonColor)

    def process(self):
      super().process()
      mousePos = pg.mouse.get_pos()
      if pg.mouse.get_pressed(num_buttons=3)[0] and self.buttonRect.collidepoint(mousePos):
        
        if (self.buttonColor == self.gameButtonColors["default"]):  
            self.setButtonColor(self.gameButtonColors["pressed"])
            self.buttonSurface.fill(self.buttonColor)
            self.value = not self.value
            print(self.value)

        else:
            self.setButtonColor(self.gameButtonColors["default"])
            self.buttonSurface.fill(self.buttonColor)
            self.value = not self.value
            print(self.value)
        




Button(20, 20, 120, 60, "My Button")
x = GameButton(200, 200, 140, 50, "Some")
x.setValue(True)

# Game loop.
while True:
    screen.fill((20, 20, 20))
    for event in pg.event.get():
        if event.type == pg.QUIT:
            pg.quit()
            sys.exit()

    for object in objects:
        object.process()

    pg.display.flip()
    fpsClock.tick(fps)