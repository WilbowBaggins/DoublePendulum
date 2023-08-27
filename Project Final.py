#Double Pendulum Program

#~~~~~~~~~~~~~
#Key
#Comments with only one # are a comment on how the a section of the program works.
#Comments with two ## are additional pieces of code which were removed either for simplicity or as they were necessary for testing at some point during programming and may still be necessary in the future.
#~~~~~~~~~~~~~


import numpy as np
from scipy import integrate
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#The lines above import the libaries used in the rest of the program.



L1Valid = False
ValidInputs = np.array([False,False,False,False,False,False,False,False])

while ValidInputs.all() == False:
    try:
        if ValidInputs[0] == False:
            L1 = float(input("Enter length of rod 1: "))
            ValidInputs[0] = True
        if ValidInputs[1] == False:
            L2 = float(input("Enter length of rod 2: "))
            ValidInputs[1] = True
        if ValidInputs[2] == False:
            M1 = float(input("Enter mass of mass on the end of rod 1: "))
            ValidInputs[2] = True
        if ValidInputs[3] == False:
            M2 = float(input("Enter mass of mass on the end of rod 2: "))
            ValidInputs[3] = True
        if ValidInputs[4] == False:
            theta1 = np.radians(float(input("Enter theta of rod 1: ")))
            ValidInputs[4] = True
        if ValidInputs[5] == False:
            theta2 = np.radians(float(input("Enter theta of rod 2: ")))
            ValidInputs[5] = True
        if ValidInputs[6] == False:
            theta1_dot = float(input("Enter theta dot of rod 1: "))
            ValidInputs[6] = True
        if ValidInputs[7] == False:
            theta2_dot = float(input("Enter theta dot of rod 2: "))
            ValidInputs[7] = True
    except:
        print("Error")
#The section of code above takes all the inputs for the initial values of the double pendulum with a simple capture system.


g1 = 9.81
g2 = 9.81

InitialCon = [theta1,theta2,theta1_dot,theta2_dot]

time = [0,100,10001]

t_span =(time[0],time[1])
times = np.linspace(time[0],time[1],time[2])

#The section above set the constant of gravity used in the program and sets the time span and number of steps.

def Pendulum(t,InitialCon):

    theta1,theta2,theta1_dot,theta2_dot = InitialCon[0],InitialCon[1],InitialCon[2],InitialCon[3]

    
    theta1_ddot = -(g1/L1)*M2*np.sin(theta1 - 2.*theta2)/(2.*M1 + M2 - M2*np.cos(2.*theta1 - 2.*theta2))\
                                  -(g1/L1)*(2.*M1 + M2)*np.sin(theta1)/(2.*M1 + M2 - M2*np.cos(2.*theta1 - 2.*theta2))\
                                  -2.*M2*(L2/L1)*theta2_dot**2*np.sin(theta1 - theta2)/(2.*M1 + M2 - M2*np.cos(2.*theta1 - 2.*theta2))\
                                  -M2*theta1_dot**2*np.sin(2.*theta1 - 2.*theta2)/(2.*M1 + M2 - M2*np.cos(2.*theta1 - 2.*theta2))
    #This is the equation supplied for the theta 1 double dot value.
    
    theta2_ddot = (g2/L2)*(M1 + M2)*(np.sin(2.*theta1 - theta2) - np.sin(theta2))/(2.*M1 + M2 - M2*np.cos(2.*theta1 - 2.*theta2))\
                                  +M2*theta2_dot**2*np.sin(2.*theta1 - 2.*theta2)/(2.*M1 + M2 - M2*np.cos(2.*theta1 - 2.*theta2))\
                                  +2.*(M1 + M2)*(L1/L2)*theta1_dot**2*np.sin(theta1 - theta2)/(2.*M1 + M2 - M2*np.cos(2.*theta1 - 2.*theta2))
    #This is the equation supplied for the theta 2 double dot value.

    
    
    Output = theta1_dot,theta2_dot,theta1_ddot,theta2_ddot
    return Output
#This pedulum function contains the two formular for both theta 1 and theta 2 double dot. The function takes arguments for the time and all the initial conditions set then returns all of the new values for each of the inital conditions.

def Energy(theta1,theta2,theta1_dot,theta2_dot):
    return ((0.5*(M1+M2)*(L1**2)*(theta1_dot**2)) + (0.5*M2*(L2**2)*(theta2_dot**2)) + (M2*L1*L2*theta1_dot*theta2_dot*np.cos(theta1-theta2)) - ((M1+M2)*9.81*L1*np.cos(theta1)) - (M2*9.81*L2*np.cos(theta2)))
#This function calculates the mechanical energy of the system based on the given conditions.

print(Energy(theta1,theta2,theta1_dot,theta2_dot))
#This prints the mechanical energy of the system for the initial conditions of the program.

##Methods = ['RK45','RK23','DOP853','Radau','BDF','LSODA']
##for i in range (6):
##    sol= integrate.solve_ivp(Pendulum,t_span,InitialCon,method=Methods[i],t_eval=times,rtol = 1e-7,atol = 1e-8)
##    plt.plot(sol.t,Energy(sol.y[0],sol.y[1],sol.y[2],sol.y[3]))
#This loops through all the of the integration methods that can be used with solve_ivp and determines the mechanical energy after each step of the intergration. These values are then plotted onto a graph of energy against time.

plt.xlabel("Time(s)")
plt.ylabel("Energy(J)")
plt.legend(['RK45','RK23','DOP853','Radau','BDF','LSODA'])
plt.title("Plot of time against total mechanical energy")
plt.show()
#This labels and displays the graph above.


##sol= integrate.solve_ivp(Pendulum,t_span,InitialCon,method='Radau',t_eval=times,rtol = 1e-7,atol = 1e-8)
#This is the integrate that solves the pendulum function using the Radau method

##def Theta1Event(t,Theta):
##    return (Theta[0])

##Theta1Event.direction = 1

def Theta2Event(t,Theta):
    ThetaSine = np.sin(Theta[1])
    ThetaVals = np.arcsin(ThetaSine)
    return (ThetaVals)

##Theta2Event.direction = 1

#These two function simply return the value of theta 1 and theta 2. The angles have been corrected by passing them through sin and then back through arcsin to get valid values for theta. The lines Theta1Event.direction = 1 and Theta2Event.direction = 1 finds the function only when theta is moving the the positive direction, this can also be set to -1 for the negative direction.

##SolPoincareTheta1 = integrate.solve_ivp(Pendulum,t_span,InitialCon,method='Radau',t_eval=times,events=Theta1Event,rtol = 1e-7,atol = 1e-8)
##print(SolPoincareTheta1.y_events)

SolPoincareTheta2 = integrate.solve_ivp(Pendulum,t_span,InitialCon,method='Radau',t_eval=times,events=Theta2Event,rtol = 1e-7,atol = 1e-8)
print(SolPoincareTheta2.y_events)

#These two functions solve the integration and then return the y values for when the theta 1 or theta 2 are equal to 0

##if SolPoincareTheta2.y.all() == sol.y.all():
##    print("True")
##else:
##    print("False")

##plt.plot(SolPoincareTheta1.y_events[0][:,1],SolPoincareTheta1.y_events[0][:,3],".")
##plt.xlabel("Theta 2")
##plt.ylabel("Theta 2 dot")
##plt.axis("scaled")
##plt.show()

#This produces a Poincare plot for when the theta 1 value is equal to zero

SolSine = np.sin(SolPoincareTheta2.y[1])
SolYVals = np.arcsin(SolSine)
##print(SolSine)
##print(np.arcsin(SolSine))

plt.plot(SolPoincareTheta2.t,SolYVals)
plt.title("Graph of Theta 2 aginst Time")
plt.xlabel("Time (s)")
plt.ylabel("Theta 2")
plt.show()
#This converts the incorrect radian values to the correct ones by putting the value through sin and back through arcsin then plots a graph of them to see roughly how many points on the Poincare plot there should be.

SolPoincareTheta1Sine = np.sin(SolPoincareTheta2.y_events[0][:,0])
SolPoincareTheta1Vals = np.arcsin(SolPoincareTheta1Sine)

SolPoincareTheta1dotSine = np.sin(SolPoincareTheta2.y_events[0][:,2])
SolPoincareTheta1dotVals = np.arcsin(SolPoincareTheta1dotSine)

#These corrects the value of theta for the Poincare plot by passing them through sin and arcsin. 


plt.plot(SolPoincareTheta1Vals,SolPoincareTheta1dotVals,".")
plt.title("Poincare plot of Theta 1 against Theta 1 dot when Theta 2 is equal to 0")
plt.xlabel("Theta 1")
plt.ylabel("Theta 1 dot")
##plt.axis("scaled")
plt.show()
#This produces a Poincare plot for when the theta 2 value is equal to zero

##plt.scatter(sol.y[0],sol.y[2])
##plt.scatter(sol.y[1],sol.y[3])
##plt.show()


x1 = np.sin(SolPoincareTheta2.y[0])*L1
y1 = -np.cos(SolPoincareTheta2.y[0])*L1
x2 = (np.sin(SolPoincareTheta2.y[1])*L2)+x1
y2 = -(np.cos(SolPoincareTheta2.y[1])*L2)+y1
plt.plot(x2,y2)
plt.plot(x1,y1)
plt.title("Phase space diagram of double pendulum")
plt.legend(["Pendulum 2","Pendulum 1"])
plt.xlabel("X")
plt.ylabel("Y")
plt.axis("scaled")
plt.show()
#This produces a phase space diagram of the double pendulum to see the path taken by the double pendulum.

