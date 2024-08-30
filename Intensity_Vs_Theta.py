# this will generate a plot showing the intensity of light vs the theta angle which the observer looks into the sky

import math
from PIL import Image, ImageDraw
from plane_parallel import intensity_at_ground
import numpy as np
import matplotlib.pyplot as plt
theta = range(181)
#x = [theta]
#y = [range(181)]
#im = Image.new('L', (256, 256), 255)
#draw = ImageDraw.Draw(im)
#extend this to phi 0, 45, 90
def main():
    tau_atm = 0.5
    theta_0 = 145
    
    #Theta_obs1 Intensity1 generation
    phi = 0
    theta_obs1 = list()
    intensity1 = list()
    for theta in range (90, 180):
        theta_obs1.append(-(180-theta))
        intensity1.append(intensity_at_ground(theta_0, float(theta), phi, tau_atm))
    phi = 180
    for theta in range (180, 90, -1):
        theta_obs1.append(180-theta)
        intensity1.append(intensity_at_ground(theta_0, float(theta), phi, tau_atm))
    
    #Theta_obs2 Intensity2 generation
    phi = 45
    theta_obs2 = list()
    intensity2 = list()
    for theta in range (90,180):
        theta_obs2.append(-(180-theta))
        intensity2.append(intensity_at_ground(theta_0, float(theta), phi, tau_atm))
    phi = 225
    for theta in range (180, 90, -1):
        theta_obs2.append(180-theta)
        intensity2.append(intensity_at_ground(theta_0, float(theta), phi, tau_atm))

    #Theta_obs3 Intensity3 generation
    phi = 90
    theta_obs3 = list()
    intensity3 = list()
    for theta in range (90, 180):
        theta_obs3.append(-(180-theta))
        intensity3.append(intensity_at_ground(theta_0, float(theta), phi, tau_atm))
    phi = 270
    for theta in range (180, 90, -1):
        theta_obs3.append(180-theta)
        intensity3.append(intensity_at_ground(theta_0, float(theta), phi, tau_atm))    

    plt.plot(theta_obs1,intensity1)
    plt.plot(theta_obs2,intensity2)
    plt.plot(theta_obs3,intensity3)
    plt.legend(['0-180','45-225','90-270'])
    plt.title('Intensity Vs. Theta')
    plt.xlabel('Theta')
    plt.ylabel('Intensity')
    plt.draw()
    plt.show()

if __name__ == '__main__':
    main()

#for i in range(190):
#    plt.title('Test Plot')
#    plt.xlabel('0')
#    plt.ylabel('I')
#    plt.plot(x,y)
#    plt.show()
#im.save('test.png')