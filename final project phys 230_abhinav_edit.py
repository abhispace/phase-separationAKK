#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 22:32:05 2021

@author: abhinav kevin and karina for PHYS 230 class
"""

import numpy as np
import matplotlib.pyplot as plt
#import scipy.stats
import matplotlib.cm as cm
import math
import time

# returns random number
def randfunc(xsize,ysize):
    return np.random.rand(xsize,ysize)

def randfunc1(xsize,ysize,zsize):
    return np.random.rand(xsize,ysize,zsize)

# program constants here
size = 32
x = np.arange(size)
y = np.arange(size)

# intializing the lattice grid
xlat1, ylat1 = np.meshgrid(x, y, sparse=False, indexing='ij')
xlat2, ylat2 = np.meshgrid(x, y, sparse=False, indexing='ij')

# simulation constants below
time_step = 10
snapshot_step = int(time_step/10)
rand_list = [-1,0,1]
kT = 1

# stores 0 (presence of solvent) or 1 (presence of tetramer): one is tertamer and the other is solvent
def setup():
    occ1 = np.empty([size, size])
    occ2 = np.empty([size, size, 2])
    init_rand1 = randfunc(np.shape(occ1)[0],np.shape(occ1)[1])
    init_rand2 = randfunc1(np.shape(occ2)[0],np.shape(occ2)[1], 2)
    
    for i in range(0,np.shape(occ1)[0]):
        for j in range(0,np.shape(occ1)[1]):
            if init_rand1[i,j] > prob_tet:
                occ1[i,j] = 0
            else:
                occ1[i,j] = 1

            for k in range(0,2):
                if init_rand2[i,j,k] > prob_dim:
                    occ2[i,j,k] = 0
                else:
                    occ2[i,j,k] = 1
                    
    return occ1,occ2


def plotfunc(x,y,xs,ys,xd1,yd1,xd2,yd2,xsd1,ysd1,xsd2,ysd2,i):
#    x = numpy.arange(0, 1, 0.05)
#    y = numpy.power(x, 2)
    
    fig = plt.figure(figsize=(10, 10), dpi=300)
    ax = fig.gca()
    ax.set_xticks(np.arange(0, size, 2))
    ax.set_yticks(np.arange(0, size, 2))
    plt.scatter(x, y, s=80, c = 'g', marker = 's', label='tetramer')
    plt.scatter(xs, ys, c = 'b', s=20, alpha = 0.7, label='solvent')
    plt.quiver(xd1, yd1, 1, 0, color= 'r', alpha = 0.8, label='dimer') # horizontal dimer
    plt.quiver(xd2, yd2, 0, 1, color= 'r', alpha = 0.8) # vertical dimer
    plt.quiver(xsd1, ysd1, 1, 0, color= 'b', alpha = 0.7, label='solvent') # horizontal solvent
    plt.quiver(xsd2, ysd2, 0, 1, color= 'b', alpha = 0.7) # vertical solvent
    plt.grid(linewidth = 0.1)
    plt.legend(bbox_to_anchor=(0.85,1.14), loc="upper left")
    ax.set_aspect('equal')
#    plt.show()
#    if anim_flag == 0:
    plt.savefig("snapshots/lattice"+"_"+str(prob_tet)+"_"+str(prob_dim)+"_step_"+str(i)+"_.png")
#    else:
#        plt.savefig("img000"+str(iterval)".png")
    plt.clf()
    plt.close()
    
def denplot(ratio,t,cutoff):    
    # Plot the density map using nearest-neighbor interpolation
    size_plt = size-2*cutoff
    x_plt = np.arange(size_plt)
    y_plt = np.arange(size_plt)
    
    # intializing the lattice grid
    xlat_plt, ylat_plt = np.meshgrid(x_plt, y_plt, sparse=False, indexing='ij')
    fig = plt.figure(figsize=(10,7), dpi=600)
    im = plt.pcolormesh(xlat_plt,ylat_plt,ratio)
    cmap = plt.get_cmap('PiYG')
    fig.colorbar(im)
    plt.savefig("snapshots/density"+"_"+str(prob_tet)+"_"+str(prob_dim)+"_step_"+str(t)+"_.png")
    plt.clf()
    plt.close()

def count(occ1):
    cutoff = 4
    num_count = np.zeros((size-2*cutoff,size-2*cutoff))

    for i in range(cutoff,np.shape(occ1)[0]-cutoff):
        for j in range(cutoff,np.shape(occ1)[1]-cutoff):
            for x in range(-cutoff,cutoff+1):
                for y in range(-cutoff,cutoff+1):
                    if occ1[i+x,j+y] == 1:
                        num_count[i-cutoff,j-cutoff] += 1
                        
                        
    box_size = 2*cutoff+1
    tot_nodes = box_size*box_size
    ratio = (num_count-(tot_nodes-num_count))/tot_nodes

    return num_count,ratio,cutoff,box_size,tot_nodes

# finding energy:
def enfunc(occ1,occ2):
    en_val = 0.0
    for i in range(0,np.shape(occ1)[0]):
        for j in range(0,np.shape(occ1)[1]):
            for k in range(0,np.shape(occ2[0,0])[0]):
                #case 1: tetramer-dimer
                if occ1[i,j] == 1 and occ2[i,j,k] == 1:
                    en_val += -2
                #case 2: tetramer-solvent
                if occ1[i,j] == 1 and occ2[i,j,k] == 0:
                    en_val += +1
                #case 3: solvent-dimer
                if occ1[i,j] == 0 and occ2[i,j,k] == 1:
                    en_val += +1
                #case 4: solvent-solvent
                if occ1[i,j] == 0 and occ2[i,j,k] == 0:
                    en_val += -1
                    
                    
    return en_val


# function to locate and plot tetramers and dimers
def locate(occ1,occ2):
    xt = np.where(occ1 == 1)[0]    # tetramer location
    yt = np.where(occ1 == 1)[1]    # tetramer location
    xs = np.where(occ1 == 0)[0]    # solvent location
    ys = np.where(occ1 == 0)[1]    # solvent location

    xd1 = np.where(occ2[:,:,0] == 1)[0]    # horizontal dimer location
    yd1 = np.where(occ2[:,:,0] == 1)[1]    # horizontal dimer location

    xd2 = np.where(occ2[:,:,1] == 1)[0]    # verticle dimer location
    yd2 = np.where(occ2[:,:,1] == 1)[1]    # verticle dimer location

    xsd1 = np.where(occ2[:,:,0] == 0)[0]    # horizontal solvent location
    ysd1 = np.where(occ2[:,:,0] == 0)[1]    # horizontal solvent location

    xsd2 = np.where(occ2[:,:,1] == 0)[0]    # verticle solvent location
    ysd2 = np.where(occ2[:,:,1] == 0)[1]    # verticle solvent location
    
    #xsd = np.where(occ2 == 0)[0]   # solvent at bond location
    #ysd = np.where(occ2 == 0)[1]   # solvent at bond location

    return xt,yt,xs,ys,xd1,yd1,xd2,yd2,xsd1,ysd1,xsd2,ysd2


def accept_reject(old_energy,new_energy):
    # condition to accept the new lattice
    accept = False
    
    if (new_energy <= old_energy):
        accept = True
    #   print(" 1 : ",accept)
    else:
        # Now apply the Monte Carlo test - compare
        # exp( -(E_new - E_old) / kT ) >= rand(0,1)
        x = math.exp( -(new_energy - old_energy) / kT )

        if (x >= np.random.uniform(0.0,1.0)):
            accept = True
#            print("accepting outlier case")
        else:
            accept = False
    #       print(" 2 : ",accept)

    return accept


def swap(occ_matrix, flag):
    # old occupancy lattice
    #old_occ1 = occ1
#    print(flag)

    boolval = True
    if flag == 0:
        while (boolval == True):
            #pick a random lattice node to swap
            rand_i = np.random.randint(1, size-1);
            rand_j = np.random.randint(1, size-1);
            #print(rand_i,rand_j)

            # pick a random neighbor to swap with
            rand_di = np.random.choice(rand_list)
            rand_dj = np.random.choice(rand_list)
            #print(rand_di,rand_dj)

            boolval = (occ_matrix[rand_i,rand_j] == occ_matrix[rand_i+rand_di,rand_j+rand_dj])
    
        #swapping the two nodes
        temp = occ_matrix[rand_i,rand_j]
        occ_matrix[rand_i,rand_j] = occ_matrix[rand_i+rand_di,rand_j+rand_dj]
        occ_matrix[rand_i+rand_di,rand_j+rand_dj] = temp
        return occ_matrix, rand_i, rand_j, rand_di, rand_dj, 2   # return 2 for the node swap case; the value is garbage it is needed so that     

    else: # swap bonds
        while (boolval == True):
            #pick a random lattice node that the bond is associated with
            rand_i = np.random.randint(1, size-1);
            rand_j = np.random.randint(1, size-1);
            
            rand_case = np.random.choice([0,1,2,3]) # random case for swapping
            rand_orient = np.random.choice([0,1])  # choose between horizontal or vertical position
            if rand_orient == 0: # horizontal case
                new_orient = 1
                if rand_case == 0:
                    rand_di = 0
                    rand_dj = -1
                if rand_case == 1:
                    rand_di = 0
                    rand_dj = 0
                if rand_case == 2:
                    rand_di = 1
                    rand_dj = 0
                if rand_case == 3:
                    rand_di = 1
                    rand_dj = -1
            else:
                new_orient = 0
                if rand_case == 0:
                    rand_di = 0
                    rand_dj = 0
                if rand_case == 1:
                    rand_di = -1
                    rand_dj = 0
                if rand_case == 2:
                    rand_di = -1
                    rand_dj = 1
                if rand_case == 3:
                    rand_di = 0
                    rand_dj = 1
                
            
            boolval = (occ_matrix[rand_i,rand_j,rand_orient] == occ_matrix[rand_i+rand_di,rand_j+rand_dj,new_orient])
    
        #swapping the two bonds
        temp = occ_matrix[rand_i,rand_j,rand_orient]
        occ_matrix[rand_i,rand_j,rand_orient] = occ_matrix[rand_i+rand_di,rand_j+rand_dj,new_orient]
        occ_matrix[rand_i+rand_di,rand_j+rand_dj,new_orient] = temp
#        print(occ_matrix[rand_i,rand_j,rand_orient],occ_matrix[rand_i+rand_di,rand_j+rand_dj,new_orient])
#        print("actual swap", occ_matrix[rand_i,rand_j, rand_orient], occ_matrix[rand_i+rand_di,rand_j+rand_dj, new_orient], rand_i,rand_j, rand_orient," -swap- ",rand_i+rand_di,rand_j+rand_dj, new_orient)


        
        return occ_matrix, rand_i, rand_j, rand_di, rand_dj, rand_orient
    
    
def reverse_swap(occ_matrix, rand_i, rand_j, rand_di, rand_dj, dk, flag):
    if flag == 0: # node case
        temp = occ_matrix[rand_i,rand_j]
        occ_matrix[rand_i,rand_j] = occ_matrix[rand_i+rand_di,rand_j+rand_dj]
        occ_matrix[rand_i+rand_di,rand_j+rand_dj] = temp 
    elif flag == 1: # bond case
        if dk == 0:
            new_dk = 1
        elif dk == 1:
            new_dk = 0
        temp = occ_matrix[rand_i,rand_j, dk]
        occ_matrix[rand_i,rand_j, dk] = occ_matrix[rand_i+rand_di,rand_j+rand_dj, new_dk]
        occ_matrix[rand_i+rand_di,rand_j+rand_dj, new_dk] = temp 
#        print("reverse swap", occ_matrix[rand_i,rand_j, dk], occ_matrix[rand_i+rand_di,rand_j+rand_dj, new_dk], rand_i,rand_j, dk," -reverse swap- ",rand_i+rand_di,rand_j+rand_dj, new_dk)

    return occ_matrix


def simulation(occ1, occ2):
    naccept = 0
    nreject = 0
    en_list = []  # list of energies after each time step

    for t in range(0, time_step):
    
    # swap random nodes ----------------------------------------------------------------
    # calculate the old energy
        old_energy = enfunc(occ1,occ2)
        occ1, rand_i, rand_j, rand_di, rand_dj, dk = swap(occ1,0)
    # new energy calculation
        new_energy = enfunc(occ1,occ2)
        accept = accept_reject(old_energy, new_energy)
        
        if accept:
            # accept the move
            naccept += 1
            total_energy = new_energy
        else:
            # reject the move - restore the old coordinates
            nreject += 1
            occ1 = reverse_swap(occ1, rand_i, rand_j, rand_di, rand_dj, dk, 0)

            total_energy = old_energy
#        print ("-----",old_energy, total_energy)
            
# moved from here 1
        # swap random bonds ----------------------------------------------------------------
        # calculate the old energy
#        old_energy = enfunc(occ1,occ2)
        old_energy = total_energy
        occ2, rand_i, rand_j, rand_di, rand_dj, dk = swap(occ2,1)
              
    # new energy calculation
        new_energy = enfunc(occ1,occ2)
        accept = accept_reject(old_energy, new_energy)
        
        if accept:
            # accept the move
            naccept += 1
            total_energy = new_energy
#            print (old_energy, total_energy)
        else:
            # reject the move - restore the old coordinates
            nreject += 1
            occ2 = reverse_swap(occ2, rand_i, rand_j, rand_di, rand_dj, dk, 1)
            total_energy = old_energy
            
        en_list.append(total_energy)

        if (t%snapshot_step == 0):
            print("\n prob dim: ", prob_dim, " prob tet: ", prob_tet, " percent: ", int(10*t/snapshot_step))
            xt,yt,xs,ys,xd1,yd1,xd2,yd2,xsd1,ysd1,xsd2,ysd2 = locate(occ1,occ2)
            plotfunc(xt,yt,xs,ys,xd1,yd1,xd2,yd2,xsd1,ysd1,xsd2,ysd2,t)
            num_count,ratio,cutoff,box_size,tot_nodes = count(occ1)
            denplot(ratio,t,cutoff)

    xt,yt,xs,ys,xd1,yd1,xd2,yd2,xsd1,ysd1,xsd2,ysd2 = locate(occ1,occ2)
    plotfunc(xt,yt,xs,ys,xd1,yd1,xd2,yd2,xsd1,ysd1,xsd2,ysd2,t+1)
    num_count,ratio,cutoff,box_size,tot_nodes = count(occ1)
    denplot(ratio,t+1,cutoff)

    AllEnergies.append(en_list)
    
    return num_count, ratio
    
    
#generating plots of all simulations on same plot using the AllEnergies list 
def plot_all():
    Time=np.arange(0,time_step)
    fig = plt.figure(figsize=(10,7), dpi=600) 
    for i in range(0,len(AllEnergies)):
        label = "Tet prob: " + str(prob_tet_list[i])
        plt.plot(Time,AllEnergies[i],label = label)
    plt.xlabel('Time step',size=15)
    plt.ylabel('Free energy',size=15)
    plt.title('Energy of System as function of time \n ',size=20)
    plt.legend(loc='upper right')
    plt.savefig('energy_plots/EnergyvsTime_prob_dim_const_'+str(prob_dim)+"_"+str(time_step)+'.png', format="png",bbox_inches='tight')
    plt.clf()
    plt.close()

#generating cluster size plot wrt tetramer-dimer ratio    
def plot_cluster():
    fig = plt.figure(figsize=(10,7), dpi=600) 
    plt.plot(r,max_count)
    plt.xlabel('tetramer-Dimer ratio',size=15)
    plt.ylabel('largest cluster size',size=15)
    plt.title('dimer prob: '+str(prob_dim),size=20)
    plt.legend(loc='upper right')
    plt.savefig('energy_plots/cluster_size_prob_dim_const_'+str(prob_dim)+"_"+str(time_step)+'.png', format="png",bbox_inches='tight')
    plt.clf()
    plt.close()    
    
#def cluster(occ):

    
#Runnung this at multiple tetramer probabilites and appending to the list generated in the previous cell
#The tetramer probabilities will be [0.2,0.4,0.6,0.8,1.0]


# filling the occupancy arrays
#prob_tet = 1.0   # probability of tetramers at nodes
#prob_dim = 0.5   # probability of dimers at bonds

prob_dim_list = [0.2,0.5,0.8]
#prob_dim_list = [0.8]
prob_tet_list = [0.2,0.4,0.5,0.6,0.8]
#prob_tet_list = [0.2]

#Running this cell only once to make an empty list 
#which we will use to add the list of energies at different tetramer probabilities

for prob_dim in prob_dim_list:
    AllEnergies=[]
    max_count = []
    max_ratio = []
    cl_count = []
    cl_ratio = []
    r=[]
    iterval = 0
    for prob_tet in prob_tet_list:
        occ1, occ2 = setup()
        countval,ratioval = simulation(occ1,occ2)
        cl_count.append(countval),cl_ratio.append(ratioval)
        max_count.append(np.amax(cl_count[iterval]))
        max_ratio.append(np.amax(cl_ratio[iterval]))
        r.append(prob_tet/prob_dim)
        iterval += 1
    plot_cluster()
    plot_all()

"""
max_count = []
max_ratio = []
for i in range(0,len(cl_count)):
    max_count.append(np.amax(cl_count[i]))
    max_ratio.append(np.amax(cl_ratio[i]))

r=[]
for prob_dim in prob_dim_list:
    for prob_tet in prob_tet_list:
        r.append(prob_tet/prob_dim)
plt.scatter(r,max_count)

#np.shape(AllEnergies)
"""

"""    
"""  