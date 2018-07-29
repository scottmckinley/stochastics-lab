
'''
Class of OU SDEs
import os
os.chdir('/Users/mjensen1/Dropbox/SPT_codes/Python_codes')
os.getcwd

myou = OUs(dim = 1, parValues =[1,10,1], time_seq = np.arange(0,10, 0.01))
my2dou = OUs(dim = 1, parValues =[1,10,1, math.pi/4], time_seq = np.arange(0,10, 0.01))
'''

import numpy as np
import random
import math
import matplotlib.pyplot as plt

#%% OU Block


class OUs(object):
    '''
    Simulated OU with drift process, ie the solution to
        dX_t = -kg(X_t - nu*t)dt + sqrt(2D)*dW_t

    Inputs are:
        dim: 1 or 2d OU process
        parValues = (D, kg, nu)  if 1 d
                    (D, kg, nu, angle in radians)  if 2 d 
        time_seq = time observations are made

    Each process has the following attributes:
        parNames: Name of the parameters of the SDE
        parValues: values of the parameters used to simulated the SDE
        time: time sequence the SDE was generated for
        Xpos: X position of one instance of the SDE
        Ypos: Y position of one instance of the SDE
    Note that setting the nu = 0, a centered OU process is returned
    '''

    def __init__(self,dim, parValues, time_seq):
        self.time = time_seq
        self.parValues = parValues
        if dim ==1:
            self.parNames = ['D','kg','nu']
            self.Xpos = sim_OUposition(dim, time_seq, parValues)
        elif dim == 2:
            sde_traj = sim_OUposition(dim, time_seq, parValues)
            self.parNames = ['D','kg','nu', 'angle']
            self.Xpos = sde_traj['Xpos']
            self.Ypos = sde_traj['Ypos']

    def plot_traj(self):
        #min_dt = min(np.diff(self.time))
        #Xpos_lim = max(abs(max(self.Xpos)), abs(min(self.Xpos))) + max(abs(max(self.Xpos)), abs(min(self.Xpos)))*0.1
        if len(self.parNames) == 3:
            plot_path(self.time, self.Xpos)
        else:
            plot_path(self.time, self.Xpos, self.Ypos)


#%% BM Block
'''Brownian Motion'''

class BMs(object):
    '''    
    Simulated Brownian motion, ie the solution to
        dX_t = sqrt(2D)*dW_t

    Inputs are:
        dim: 1 or 2d BM process
        parValues = (D)  if 1 d or 2D
        time_seq = time observations are made

    Each process has the following attributes:
        parNames: Name of the parameters of the SDE
        parValues: values of the parameters used to simulated the SDE
        time: time sequence the SDE was generated for
        Xpos: X position of one instance of the SDE
        Ypos: Y position of one instance of the SDE
    Note that setting the nu = 0, a centered OU process is returned
    '''
    
    def __init__(self,dim, parValues, time_seq):
        self.time = time_seq
        self.parValues = parValues
        if dim == 1:
            self.parNames = ['D']
            self.Xpos = sim_BMposition(dim, time_seq, parValues)
            self.dim = 1
        elif dim == 2:
            sde_traj = sim_BMposition(dim, time_seq, parValues)
            self.parNames = ['D']
            self.Xpos = sde_traj['Xpos']
            self.Ypos = sde_traj['Ypos']
            self.dim = 2

    def plot_traj(self):
        #min_dt = min(np.diff(self.time))
        #Xpos_lim = max(abs(max(self.Xpos)), abs(min(self.Xpos))) + max(abs(max(self.Xpos)), abs(min(self.Xpos)))*0.1
        if self.dim == 1:
            plot_path(self.time, self.Xpos)
        else:
            plot_path(self.time, self.Xpos, self.Ypos)
    

#%% Switching block
            
class switchOUs(object):
       
    '''
    Simulated switching OU with drift process, ie the solution to
        dX_t = -kg(X_t - nu*t)dt + sqrt(2D)*dW_t for tau_i <=t < tau_i+1

    Inputs are:
        dim: 1 or 2d switching OU process
        parValues = vec{D} vector of diffusivities
                    vec{kg}  vector of k/g
                    vec{nu} vector of velocities
                    vec{tau} vector of switch points
                     
        time_seq = time observations are made

    Each process has the following attributes:
        parNames: Name of the parameters of the SDE
        parValues: values of the parameters used to simulated the SDE
        time: time sequence the SDE was generated for
        Xpos: X position of one instance of the SDE
        Ypos: Y position of one instance of the SDE
    Note that setting the nu = 0, a centered OU process is returned
    '''
    
    def __init__(self, dim,time_seq, D_vec, kg_vec, nu_vec, tau_vec = ['random', 2],angle = None):
        self.time = time_seq
        self.dim = dim
        
        if type(tau_vec) == int:
            tau_vec = [tau_vec]
            self.tau = tau_vec
        elif type(tau_vec[0]) == str:
            tau_vec = np.sort(np.random.randint(5,len(time_seq)-4, tau_vec[1] ))
            self.tau = tau_vec
        else:
            self.tau = tau_vec
        
        if type(D_vec) == int or type(D_vec) == float:
            self.D = np.repeat(float(D_vec), len(tau_vec)+1)
        else:
            self.D = D_vec
            
        if type(kg_vec) == int or type(kg_vec) == float:
            self.kg = np.repeat(float(kg_vec), len(tau_vec) + 1)
        else:
            self.kg = kg_vec
            
        if type(nu_vec) == int or type(nu_vec) == float:
            self.nu = np.repeat(float(nu_vec), len(tau_vec) + 1)
        else:
            self.nu = nu_vec
            
        if dim == 1:
            sde_traj= sim_OUswitch(dim, time_seq, self.D, self.kg, self.nu, self.tau)
            self.Xpos = sde_traj['Xpos']
        elif dim == 2:
            sde_traj = sim_OUswitch(dim, time_seq, self.D, self.kg, self.nu, self.tau)
            self.Xpos = sde_traj['Xpos']
            self.Ypos = sde_traj['Ypos']
            self.angle = angle  
    
    def plot_traj(self):
        #min_dt = min(np.diff(self.time))
        #Xpos_lim = max(abs(max(self.Xpos)), abs(min(self.Xpos))) + max(abs(max(self.Xpos)), abs(min(self.Xpos)))*0.1
        if self.dim == 1:
            plot_switch_path(self.tau, self.time, self.Xpos)
        else:
            plot_switch_path(self.tau, self.time, self.Xpos, self.Ypos)
    
            
    
            




                
 #%%              

'''Required Simulation Functions'''

def sim_BMposition(dim, time_seq, parValues):
    diff_constant = float(parValues[0])
    dt_seq = np.diff(time_seq)
    BM_var = np.multiply(np.repeat(2*diff_constant,len(dt_seq)),dt_seq)
    BM_sd = [np.sqrt(var)for var in BM_var]
    
    if dim == 1:
        x_inc = np.multiply(BM_sd,np.random.normal(loc=0,scale = 1, size = len(dt_seq)))
        x_inc = x_inc.tolist(); x_inc.append(0)
        Xn = np.cumsum(x_inc).tolist()
        return {'Xpos': Xn}
    if dim == 2:
        x_inc = np.multiply(BM_sd,np.random.normal(loc=0,scale = 1, size = len(dt_seq)))
        x_inc = x_inc.tolist();x_inc.append(0)
        Xn = np.cumsum(x_inc).tolist()
        y_inc = np.multiply(BM_sd,np.random.normal(loc=0,scale = 1, size = len(dt_seq)))
        y_inc = y_inc.tolist();y_inc.append(0)
        Yn = np.cumsum(y_inc).tolist()
        return {'Xpos': Xn, 'Ypos':Yn}


def sim_OUposition(dim,time_seq, parValues):
    diff_constant = float(parValues[0])
    spring_constant = float(parValues[1])
    velocity_constant = float(parValues[2])

    num_steps = time_seq.shape[0]
    dt_seq = np.diff(time_seq)
    rho_seq = math.e**(-1*dt_seq*spring_constant)
    OU_variance = (diff_constant/spring_constant)*(1 - math.e**(-2*dt_seq*spring_constant))

    if dim == 1:
        Xn = list(); Xn.append(0)
        for n in range(1, num_steps):
            An = (time_seq[n] - time_seq[n-1]*rho_seq[n-1]) - (1/spring_constant)*(1 - rho_seq[n-1])
            new_X = rho_seq[n-1]*Xn[n-1] + velocity_constant*An + math.sqrt(OU_variance[n-1])*random.gauss(0,1)
            Xn.append(new_X)
        return Xn
    if dim == 2:
        angle_constant = float(parValues[3])
        Xn = list(); Xn.append(0)
        Yn = list(); Yn.append(0)
        for n in range(1,num_steps):
            An = (time_seq[n] - time_seq[n-1]*rho_seq[n-1]) - (1/spring_constant)*(1 - rho_seq[n-1])
            mu_X = rho_seq[n-1]*Xn[n-1] + velocity_constant*An*math.cos(angle_constant)
            mu_Y = rho_seq[n-1]*Yn[n-1] + velocity_constant*An*math.sin(angle_constant)

            new_X = mu_X + math.sqrt(OU_variance[n-1])*random.gauss(0,1); Xn.append(new_X)
            new_Y = mu_Y + math.sqrt(OU_variance[n-1])*random.gauss(0,1); Yn.append(new_Y)
        return {'Xpos': Xn, 'Ypos':Yn}

                                                                    
def sim_OUswitch(dim, time_seq, D_vec, kg_vec, nu_vec, tau_vec, angle =  None):
    diff_constant = D_vec
    spring_constant = kg_vec
    velocity_constant = nu_vec
    switch_pts = tau_vec; switch_pts.insert(0,0); switch_pts.append(len(time_seq)-1)
    
    dt_seq = np.diff(time_seq)

    if dim == 1:
        Xn = list(); Xn.append(0)
        
        for j in range(0, len(spring_constant)):
            rho_seq = math.e**(-1*dt_seq*spring_constant[j])
            OU_variance = (diff_constant[j]/spring_constant[j])*(1 - math.e**(-2*dt_seq*spring_constant[j]))
            
            for n in range(switch_pts[j]+1, switch_pts[j+1]+1):
                An = (time_seq[n] - time_seq[n-1]*rho_seq[n-1]) - (1/spring_constant[j])*(1 - rho_seq[n-1])
                drift_mean = rho_seq[n-1]*Xn[n-1] + velocity_constant[j]*(
                        An - time_seq[switch_pts[j]]*(1- rho_seq[n-1])) + Xn[switch_pts[j]]*(1 - rho_seq[n-1])
                new_X = drift_mean+ math.sqrt(OU_variance[n-1])*random.gauss(0,1)
                Xn.append(new_X)
        return {'Xpos': Xn}
      
        
    if dim == 2:
        angle_constant = float(angle)
        Xn = list(); Xn.append(0)
        Yn = list(); Yn.append(0)
        for j in range(0, len(switch_pts)):
            rho_seq = math.e**(-1*dt_seq*spring_constant[j])
            OU_variance = (diff_constant[j]/spring_constant[j])*(1 - math.e**(-2*dt_seq*spring_constant[j]))
            
            for n in range(switch_pts[j-1]+1,switch_pts[j]+1):
                An = (time_seq[n] - time_seq[n-1]*rho_seq[n-1]) - (1/spring_constant[j])*(1 - rho_seq[n-1])
                mu_X = rho_seq[n-1]*Xn[n-1] + velocity_constant[j]*math.cos(angle_constant)*(
                        An - time_seq[switch_pts[j]]*(1 - rho_seq[n-1])) + Xn[switch_pts[j]]*(1 - rho_seq[n-1])
                mu_Y = rho_seq[n-1]*Yn[n-1] + velocity_constant[j]*math.sin(angle_constant)*(
                        An - time_seq[switch_pts[j]]*(1 - rho_seq[n-1])) + Yn[switch_pts[j]]*(1 - rho_seq[n-1])

                new_X = mu_X + math.sqrt(OU_variance[n-1])*random.gauss(0,1); Xn.append(new_X)
                new_Y = mu_Y + math.sqrt(OU_variance[n-1])*random.gauss(0,1); Yn.append(new_Y)
        return {'Xpos': Xn, 'Ypos':Yn}

        
       


#%% Required plotting functions 
def plot_path(time_seq, Xposition, Yposition = None):
    min_dt = min(np.diff(time_seq))
    Xpos_lim = max(abs(max(Xposition)), abs(min(Xposition))) + max(abs(max(Xposition)), abs(min(Xposition)))*0.1
    plt.figure()
    if Yposition == None:
        plt.plot(time_seq, Xposition)
        plt.axis([0-min_dt, max(time_seq) + min_dt, -1*Xpos_lim, Xpos_lim])
        plt.ylabel('Position (microns)')
        plt.xlabel('Time (seconds)')
        #plt.title(r'$D = $')
    else:
        Ypos_lim = max(abs(max(Yposition)), abs(min(Yposition))) + max(abs(max(Yposition)), abs(min(Yposition)))*0.1
        pos_lim = max(Xpos_lim, Ypos_lim)
        plt.subplot(131)
        plt.plot(Xposition, Yposition)
        plt.axis([-1*pos_lim, pos_lim, -1*pos_lim, pos_lim])
        plt.xlabel('X Position (microns)')
        plt.ylabel('Y Position (microns)')

        plt.subplot(132)
        plt.plot(time_seq, Xposition)
        plt.axis([0 - min_dt, max(time_seq) + min_dt ,-1*pos_lim, pos_lim])
        plt.ylabel('X Position (microns)')
        plt.xlabel('Time (seconds)')

        plt.subplot(133)
        plt.plot(time_seq, Yposition)
        plt.axis([0-min_dt, max(time_seq) + min_dt ,-1*pos_lim, pos_lim])
        plt.ylabel('Y Position (microns)')
        plt.xlabel('Time (seconds)')

        plt.show()

    
    
def plot_segments(x_coors, y_coors, breaks, subplot_index = [1,1,1], plot_labels = 
                  ["Time (second)", "Position (micron)"]):
    mycolors_all = ["#5E4FA2" ,"#466DB0", "#348BBB", "#50A9AF" ,"#6DC4A4", 
                "#91D3A4" ,"#B4E0A2","#D3ED9B" ,"#EBF7A0", "#F8FCB4", "#FEF6B1", 
                "#FEE695" ,"#FDD07D" , "#FDB567","#F99655" ,"#F47346", "#E65948", 
                "#D6404E", "#BA2148" ,"#9E0142"]
    
    breaks.append(len(x_coors)-1)
    breaks.insert(0, 0)

    x_seg = [x_coors[breaks[i-1]:breaks[i] +1]  for i in range(1, len(breaks))]
    y_seg = [y_coors[breaks[i-1]:breaks[i]+1]  for i in range(1, len(breaks))]
    
    color_index = int(np.floor(len(mycolors_all)/len(breaks)))

    ax = plt.subplot(subplot_index[0],subplot_index[1],subplot_index[2])
    
    for j in range(len(x_seg)):
        ax.plot(x_seg[j], y_seg[j], color = mycolors_all[j*color_index])
        ax.scatter(x_coors[breaks[j]],y_coors[breaks[j]], 
                   color = mycolors_all[j*color_index],marker ='o' )
    
        plt.xlabel(plot_labels[0]) #labels[0])
        plt.ylabel(plot_labels[1]) #labels[0])


    

def plot_switch_path(tau_vec, time_seq, Xposition, Yposition = None):

    plt.figure()
    if Yposition == None:
        plot_segments(time_seq, Xposition, breaks = tau_vec)
    else:
        plot_segments(Xposition, Yposition, breaks = tau_vec, subplot_index= [1,3,1],
                     plot_labels = ['X Position (microns)','Y Position (microns)'])

        plot_segments(time_seq, Xposition, breaks = tau_vec, subplot_index= [1,3,2])
       
        plot_segments(time_seq, Yposition, breaks = tau_vec, subplot_index= [1,3,3],
                     plot_labels = ['Time (seconds)','Y Position (microns)'])


#%%
   

#myous =OUs(dim = 2, parValues =[1,10,1, math.pi/4], time_seq = np.arange(0,10, 0.1))


#plot_segments(myous.time, myous.Ypos, [25,75])
#plot_switch_path([25,75],time_seq =  myous.time,Xposition = myous.Xpos, Yposition = myous.Ypos)


#my_switch = sim_OUswitch(dim= 1, time_seq = np.arange(0,10, 0.1), D_vec = [1.0,1.0,1.0], 
#                         kg_vec = [10.0,10.0, 10.0], nu_vec = [10.0,0.0, -10.0], tau_vec =[20, 30])


