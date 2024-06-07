"""CS_161_project.py - Simuluation of the time evolution of the solar system. 
 
Desciption:
   The simulation of the solar system requires the use of Newtonian Mechanics of planetary motion. This can be expressed through the Force Equation:
\vec{F_{i}} = \sum_{j=1}^{10} \frac{G m_{i} m_{j}}{|\vec{r_{i}}-\vec{r_{j}}|^{3}} (\vec{r_{i}}-\vec{r_{j}}).
This is Newton's law of universal gravitation. This expression is the total force that mass m_i feels from all the other masses in the system.
The other masses in the system are m_j and since there are ten other objects within the solar system, the total force is summed over ten masses. G
represents the Gravitational Constant, which equals 6.67e-11 [Nm^{2}/kg^2]. The vectors of r_i and r_j represent the positional vectors corresponding
to masses m_i and m_j respectively. In order to evolve the system, each object's position must be dependent on time. To relate a force being
exerted on an object and its position, numerical integration is required. The relationship stems from the velocity of an object being the time derivative
of its position. The acceleration of an object is the time derivative of its velocity. The acceleration of an object is equal to the total force exerted
on the object divided by its mass. In order to go from force to position, numerical integration is required.
   This calculation is performed in a method of each celestial object. In order to simulate the solar system's motion, each celestial object is an instance
of a class called Cele(). The celestial object's attributes include its name, mass, initial position, and initial velocity. From these quantities, its net force
can be calculated, and thus its new position per iteration of time. The iteration of time is referred to as dt. The naming of this quantity comes from what
It is represented in calculus, called a differential. All this quantity represents is a small increment in time. The smaller the increment is, the more accurate
the simulation becomes, but the more computationally intensive the calculation becomes. Since the solar system simulation is to scale, dt can be relatively
large. Common quantities of dt are usually values that approach zero. The numerical integration required for the update of each celestial object's new position
is performed in the method get_pos(). The last method __str__() returns the name of the celestial object.
   The attributes required for each celestial object is read data and can be found from this NASA resource: https://nssdc.gsfc.nasa.gov/planetary/factsheet/.
Note that for the initial position of each celestial object, its Perihelion was used and for its initial velocity, its maximal orbital velocity was used.
This comes from Kepler's laws which states that the closer an object is to its star, the faster it orbits it. By choosing the closest distance from the sun
(perihelion) and its maximum velocity, the simulated orbit is more accurate, and represents the eccentricity of the orbit. Along with the celestial object's
eccentric orbit, its inclination angle was implemented to increase its precision. This implementation was performed using trigonometry where the inclination angle
was based on the orbit's angle from the xy-plane. These calculations were performed beforehand, and the results of these calculations were stored in a .csv file.
In conclusion, all the attributes of the celestial object are stored in a .csv file.
   These parameters are read into the Python script using a library called Pandas. This occurs within the function data_load(). This function returns the list of all
celestial objects. The next process takes this list of celestial objects, and computes their new positions per time iteration. Now, In order to store all the
positions for each star, and not have their previous positions deleted, a dictionary is utilized. Here, the key of the dictionary is the name of the celestial object,
and the value is a list of numpy arrays, which is the list of all positions. Lastly, the visualization of each celestial object can be chosen. The chosen option are
to visualize each individual celestial object's orbit, the Earth and Moon's orbit, or the orbit of the entire solar system.
 
Kyle Gourlie
 
Design (pseudocode):

class Cele(object):
    def __init__(self,mass,pos,vel):
        self.mass,pos,vel = mass,pos,vel

    def get_pos(self):
        dt = t_j - t_{j-1} 
        t_j \in \{ t_0,t_1,t_2,...,t_n \} 
        \vec{F_{i}}(t_j) = \sum_{k \neq i}^{N} -\frac{G m_{i} m_{k}}{|\vec{r_{i}}(t_{j-1}) - \vec{r_{k}}(t_{j-1})|^3}(\vec{r_{i}}(t_{j-1}) - \vec{r_{k}}(t_{j-1})) 
        \vec{a_{i}}(t_j) = \frac{\vec{F_{i}}(t_j)}{m_i} 
        \vec{v_{i}}(t_j) = \vec{v_{i}}(t_{j-1}) + \vec{a_{i}}(t_j)dt 
        \vec{r_{i}}(t_j) = \vec{r_{i}}(t_{j-1}) + \vec{v_{i}}(t_j)dt 

    def __str__(self):
        return self.name

def data_load():
    celes = []
    with pandas.open(.csvfile) as file:
        for data in file:
        obj = Cele(data[0],data[1],data[2],data[3],data[4],data[5],data[6])
        celes.append(obj)
        del obj
        return celes

def plotter(dictionary of objects):
    decision = input('Which object would you like to see?')
    if decision = {S,Me,Ve,E,M,Ma,J,Sa,U,Ne,Pu}:
        plt.plot(chosen object)

def pos_calc(celes):
    dictionary of objects = {}
    for obj in celes:
        dictionary value = obj.mass
        dictioanry key = []
    
    yrs = input('enter number of years to evolve')
    n = 0
    while n<= yrs:
        for obj in celes:
            dictionary of objects[obj.name] = obj.get_pos()
        n = n + 1
    return dictionary of objects

if __name__ == "__main__":
    celes = data_load()
    pos_dict = pos_calc(celes)
    plotter(pos_dict)

Testing:
Doctest case is implemented in Cele class, testing the creation of the class, along with computing its position. The implementation of try/except statements 
were developed by entering invalid entry types into .csv file and viewing which errors occur. Invalid entries include non-float entries for Numpy arrays, and
invalid comma separated value formats. When exception occurs, the function and/or method returns None. Due to this choice of handling errors, additional try/except
statements are needed. There are two spots in the code where input from the user is required. To handle invalid inputs, these inputs are contained in loops, where
the loop breaks only when a valid entry is provided.
"""
#imports 
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import astropy.constants as c
import pandas as pd
import doctest

#Increment in time [s]. This value is small enough for good accuracy
dt = 1e5

class Cele(object):
    """Celestial object class for simulating the solar system.

    Parameters:
    name: Name of Celestial Object
    mass: Mass [kg] of Celestial Object
    pos: 3-Dimensional Numpy Array of the Position [m] of the Celestial Object
    vel: 3-Dimensional Numpy Array of the Velocity [m/s] of the Celestial Object

    Examples:
    >>> planet_x = Cele(name='planet_x',mass = 1e24, pos = np.array([0,5,0]), vel = np.array([0,0,2]))
    >>> print(planet_x)
    planet_x
    >>> print(planet_x.mass)
    1e+24
    >>> print(planet_x.pos)
    [0 5 0]
    >>> print(planet_x.vel)
    [0 0 2]
    >>> planet_y = Cele(name='planet_h',mass = 1e24, pos = 'hi', vel = np.array([0,0,2]))
    >>> print(planet_y)
    planet_h
    >>> print(planet_y.pos)
    hi
    >>> planet_y.get_pos()
    planet_h does not have attribute 'f_net' yet.
    >>> planet_y.f_net = 5
    >>> planet_y.get_pos()
    One of the vector quantities (f_net,vel,pos) is NOT a Numpy Array.
    """
    def __init__(self, name:str, mass:float, pos:np.array, vel:np.array):
        self.name = name
        self.mass = mass  
        self.pos = pos 
        self.vel = vel
        
    def get_pos(self):
        """Computes the new position of a celestial object based on numerical integration. It is important to note
        that this method should only be called up after the f_net attribute has been assigned to a celestial object.
        """
        try:
            #net acceleration
            #\vec{a_{i}}(t_j) = \frac{\vec{F_{i}}(t_j)}{m_i} 
            acc = self.f_net/self.mass
            #updated new velocity through euler method
            #\vec{v_{i}}(t_j) = \vec{v_{i}}(t_{j-1}) + \vec{a_{i}}(t_j)dt 
            self.vel = self.vel + acc*dt
            #updated new position through euler method
            #\vec{r_{i}}(t_j) = \vec{r_{i}}(t_{j-1}) + \vec{v_{i}}(t_j)dt 
            self.pos = self.pos + self.vel * dt

        except AttributeError:
             print(f'{self.name} does not have attribute \'f_net\' yet.')
             return None
        
        except np.core._exceptions._UFuncNoLoopError:
            print('One of the vector quantities (f_net,vel,pos) is NOT a Numpy Array.')
            return None
    
    def __str__(self):
        """Returns the name of the celestial object"""
        return f"{self.name}"
    


def data_load():
    """Loads parameters required for the Cele class including name, mass, initial position, and initial velocity

    Returns:
        cele_objs (list): list of instances of Cele
        None: if csv file was not in proper comma separated values format
    """
    cele_objs = []
    #try/except used to catch possible format errors in csv file or directory name
    try:
        file = pd.read_csv('CS_161_project_data.csv', header = 0)
        #list of names of the celestial objects
        names = [i for i in file.columns] 
        #deleting 'Labels' from name list
        del names[0]
        
        for name in names:
            #creates a list of parameters for each celestial object
            l_param = [param for param in file[name]]
            #creates celestial object based on all 7 parameters
            cele_obj = Cele(name = name,mass = l_param[0], pos = np.array([l_param[1],l_param[2],l_param[3]]),
                            vel = np.array([l_param[4],l_param[5],l_param[6]]))
            cele_objs.append(cele_obj)
            #object is no longer needed
            del cele_obj
        return cele_objs
    
    except pd.errors.ParserError:
        print('\ncsv file is not in comma separated values format\n')
        return None
    
    except FileNotFoundError:
        print('Directory of csv file is invalid.')
        return None



def force(obj_1: Cele, obj_list:list):
    """Computes the total force exerted upon obj_1 due to the list of celestial objects, obj_list

    Latex Expression:
    \vec{F_{i}}(t_j) = \sum_{k \neq i}^{N} -\frac{G m_{i} m_{k}}{|\vec{r_{i}}(t_{j-1}) - \vec{r_{k}}(t_{j-1})|^3}(\vec{r_{i}}(t_{j-1}) - \vec{r_{k}}(t_{j-1}))

    Args:
        obj_1 (Cele): Celestial object that the total force is calculated on
        obj_list (list): The total list of celestial objects.

    Returns:
        None: If the data used in the csv is incompatable with Numpy Arrays
    """
    f_net = []
    for obj in obj_list:
        #This ensures that it is not calculating the force exerted on itself
        if obj.name == obj_1.name:
            continue
        else:
            #try and except for issues caused by readings from csv file and numpy 
            try:
                #The benefit of numpy arrays is the ease of numpy arthmetic
                #displacement vector between the two objects
                r = obj_1.pos - obj.pos
                #the magnitude of the displacement vector
                r_mag = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
                #Newton's law of universal gravitational with a minus sign since gravity is attractive
                F = -(c.G.value * obj_1.mass * obj.mass)/r_mag**3 * r
                f_net.append(F)

            except np.core._exceptions._UFuncNoLoopError:
                print('\nThe numpy arrays used in computation is incompatable with data contained in csv.\n')
                return None
        #object is no longer needed
        del obj
    #sets the net force of the object
    obj_1.f_net = sum(f_net)



def pos_calc(celes):
     #dictionary of lists of numpy arrays
    pos_dict = {}
    n_vals = []
    for obj in celes:
        #sets each entry in the dictionary with a name and a empty list
        pos_dict[obj.name] = []

    n = 0
    while n <= it_er:
        for obj in celes:
            force(obj,celes)
            #try and except to catch issues with calculating forces. It is needed because it can't be calcualed, then NONE is returned
            try:
                obj.get_pos()
            except AttributeError:
                print()
                print('f_net does not exist.\nCheck to make sure data contained in csv are valid types.\n')
                exit()
            #stores the value in the list contained in dictionary
            pos_dict[obj.name].append(obj.pos)
        n_vals.append(n /3.154e7*dt) 
        n = n + 1
        
        #completion tab to show progress because sometimes it takes a long time
        print(f'{round(n/it_er*100,2)}% Complete',end='\r')
    return pos_dict, n_vals



def plotter(pos_dict,n_vals):
    """Plots the data of the solar system simulation based on choice of celestial object

    Args:
        pos_dict (dict): dictionary containing the name and all positions of each object
        n_vals (list): list of iterable time values
    """    """"""
    #variable used to kill the whileloop
    valve = True

    while valve:
        quest = input('Do you want to see plots [y/n]? ').lower()

        if quest == 'y':
            key_name = {'S':'Sun', 'Me':'Mercury', 'Ve':'Venus',
                        'E':'Earth', 'M': 'Moon', 'Ma': 'Mars', 
                        'J': 'Jupiter', 'Sa': 'Saturn', 'U': 'Uranus',
                        'Ne': 'Neptune', 'Pu': 'Pluto'}
            print("""\n
Keys:
        
Entire Solar System: D
The Sun: S
Mercury: Me
Venus: Ve
Earth: E
The Earth and Moon: EM
Mars: Ma
Jupiter: J
Saturn: Sa
Uranus: U
Neptune: Ne
Pluto: Pu\n""")
            valve_2 = True
            while valve_2:
                sel_name = input('Enter key for which celestial orbit you wish to view: ')
                
                if sel_name == 'D':
                    #creation of figure subplot
                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection='3d')
                    ax.set_xlabel('X [m]')
                    ax.set_ylabel('Y [m]')
                    ax.set_zlabel('Z [m]')
                    #iterates through the key and values of dictionary
                    for name, positions in pos_dict.items():
                        #unzips the Numpy Array as x,y,z coordinates
                        x, y, z = zip(*positions)
                        #plots the data
                        ax.plot(x, y, z, label=name)
                    #adds title of the plot
                    plt.title(f'{yrs} Orbit')
                    ax.legend()
                    plt.show()
                    valve_2 = False
                    
                elif sel_name == 'EM':
                    #creation of figure subplot
                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection='3d')
                    ax.set_xlabel('X [m]')
                    ax.set_ylabel('Y [m]')
                    ax.set_zlabel('Z [m]')
                    x_E, y_E, z_E = zip(*pos_dict['Earth'])
                    x_M, y_M, z_M = zip(*pos_dict['Moon'])
                    ax.plot(x_E, y_E, z_E, label='Earth')
                    ax.plot(x_M, y_M, z_M, label='Moon')
                    plt.title(f'{yrs} Orbit')
                    ax.legend()
                    plt.show()
                    valve_2 = False

                elif sel_name in key_name.keys():
                    #creation of figure subplot
                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection='3d')
                    ax.set_xlabel('X [m]')
                    ax.set_ylabel('Y [m]')
                    ax.set_zlabel('Z [m]')
                    x, y, z = zip(*pos_dict[key_name[sel_name]])
                    #scatter plot to use a colormap
                    mapp = ax.scatter(x, y, z, label=key_name[sel_name],c = n_vals, cmap=cm.hot)
                    plt.colorbar(mapp,label = 'yrs')
                    plt.title(f'{yrs} Orbit')
                    ax.legend()
                    plt.show()
                    valve_2 = False
                
                else:
                    print('\nInvalid Response\n')

        elif quest == 'n':
            valve = False
        else:
            print('\nInvalid Response\n')





if __name__ == "__main__":
    doctest.testmod()
    #loads all the celestrial objects as list of python objects
    while True:
        try:
            yrs = float(input('How many years do you want the solar system to time evolve?'))
            if yrs > 0:
                break
            else:
                print('\nTime only moves forward. Choose a positive number.\n')
        except ValueError:
            print('\nInvalid Type. Try Again.\n')
    
    #calculation of iteration number based on yrs entered and dt value
    it_er = int(yrs*3.154e7/dt)
    celes = data_load()
    #if None is returned due to errors
    if celes == None:
        print('\nNo celestial objects were created\n')
        exit()
    
    pos_dict, n_vals = pos_calc(celes)
    plotter(pos_dict,n_vals)


    
    

