### FUNCTIONS FOR EFFICIENT RUNNING OF SIMULATIONS ###

def ol_rm(domain_start, domain_end, box_length, num_p, radius=1, mass=1, momentum_max=3):
    bl = box_length
    ds=domain_start
    de=domain_end
    N=num_p #number particles
    r=radius #radius particles
    m=mass #particle mass
    p_max=momentum_max #max momentum

    x_points = np.linspace(ds,de,int(N**0.5))
    y_points = np.linspace(0,bl,int(N**0.5))

    #particles
    init_random=[] #part
    for x in x_points:
        for y in y_points:
            pos_x = x
            pos_y = y
            pos=Vector(pos_x,pos_y)
            mom_x=(random()-0.5)*2*p_max                  #allow for positive and negative momenta
            mom_y=(random()-0.5)*2*p_max
            mom=Vector(mom_x,mom_y)
            init_random.append(Particle(pos,mom,radius=r,mass=m))
            
    return init_random


def test_pressure(number_particles):
    # parameters
    bl = 40
    timestep = 0.25
    steps = 1000
    p = ol_rm(domain_start=0, domain_end=bl, box_length=bl, num_p=number_particles, radius=0.1, mass=1, p_max=2)

    s_press=Simulation(particles=p,box_length=bl,dt=timestep)

    for step in range(steps):
        s_press.step()

    plt.plot(np.linspace(0,timestep*steps,steps),s_press.pressure)
    
    
    return np.mean(s_press.pressure), e_t(s_press)  


### OTHER FUNCTIONS FOR DATA ANALYSIS ###

def e_t(simulation):
    timestep = simulation.dt
    steps = len(simulation.trajectory)

    total_ke = []
    for t in simulation.trajectory:
        instant_ke = 0
        for p in t:
            instant_ke += p.momentum.norm()**2/(2*p.mass)
        total_ke.append(instant_ke)

    ke_mean = np.mean(total_ke)   
    ke_std = np.std(total_ke)   

    time = np.linspace(0,steps*timestep,steps)
    av_long = np.linspace(ke_mean,ke_mean,steps)

    plt.figure(figsize=[8,5])
    plt.title('Total Kinetic Energy of Simulation')
    plt.xlabel('Time')
    plt.ylabel('Kinetic Energy')
    plt.ylim(ke_mean-5*ke_std,ke_mean+5*ke_std)
    plt.plot(time,total_ke,linewidth=0.4, label='instantaneous w/ std',color='blue')
    plt.errorbar(time,total_ke,yerr=ke_std,marker='',linestyle='',elinewidth=0.02, color='blue')
    plt.plot(time,av_long,marker='',linestyle='--',linewidth=1,color='black', label='mean')
    plt.legend(loc='best')

    print('MEAN KE', ke_mean)
    print('STD KE', ke_std)
    
    return ke_mean


def speed_std(simulation):
    timestep = simulation.dt
    steps = len(simulation.trajectory)

    total_vstd = []
    for t in simulation.trajectory:
        instant_v = []
        for p in t:
            instant_v.append(p.velocity().norm())
        total_vstd.append(np.std(instant_v))

    time = np.linspace(0,steps*timestep,steps)

    plt.figure(figsize=[8,5])
    plt.title('Standard Deviation of Particle Speeds')
    plt.xlabel('Time')
    plt.ylabel('Standard Deviation')
    #plt.ylim(ke_mean-5*ke_std,ke_mean+5*ke_std)
    plt.plot(time,total_vstd,linewidth=0.4, label='',color='blue')
    #plt.legend(loc='best')
          
        
def pos_profile(simulation):
    bl = simulation.box_length # simulation box length
    bin_len = simulation.particles[0].radius/10

    bins = np.linspace(0,bl,int(bl/bin_len))
    bin_hist_tot = []

    for t in simulation.trajectory:
        bin_hist = []
        for i in bins:
            bin_hist.append(0)
        for p in t:
            for b in range(len(bins)):
                if p.position.x >= bins[b] and p.position.x <= bins[b]:
                    bin_hist[b] += 1 
        bin_hist_tot.append(bin_hist)


    bin_hist_av = []

    for i in range(len(bins)):
        total = 0
        for t in range(len(bin_hist_tot)):
            total += bin_hist_tot[t][i]
        bin_hist_av.append(total/len(bin_hist_tot))    

    plt.figure(figsize=[8,5])
    plt.title('Density Profile of Randomly Distributed Gas')
    plt.xlabel('x spacial axis')
    plt.ylabel('Number of Particles')
    plt.xlim(0,bl)
    plt.ylim(0,1)
    plt.bar(bins,bin_hist_av,width=1)
    #plt.plot(bins,bin_hist_av)
    #plt.legend(loc='best')
    
    
def expansion(simulation):
    timestep = simulation.dt
    steps = len(simulation.trajectory)
    bl=simulation.box_length #box_length

    total_lhs = []
    total_rhs = []
    for t in simulation.trajectory:
        instant_lhs = 0
        instant_rhs = 0
        for p in t:
            if  p.position.x <= 50:
                instant_lhs += 1
            elif p.position.x >= 50:
                instant_rhs += 1
        total_lhs.append(instant_lhs)
        total_rhs.append(instant_rhs)

    time = np.linspace(0,steps*timestep,steps+1)

    plt.figure(figsize=[8,5])
    plt.title('Diffusion into a Vacuum')
    plt.xlabel('Time')
    plt.ylabel('Number of Particles')
    plt.plot(time,total_lhs,label='Left Hand Side')
    plt.plot(time,total_rhs,label='Right Hand Side')
    plt.legend(loc='best')
    
    
def mixing(simulation):
    '''NOTE: currently set to the same number of particles of each type'''
    timestep = simulation.dt
    steps = len(simulation.trajectory)
    bl=simulation.box_length #box_length
    N1 = int(len(simulation.particles)/2)
    N2 = int(len(simulation.particles)/2)

    total1_lhs = []
    total1_rhs = []
    total2_lhs = []
    total2_rhs = []
    for t in simulation.trajectory:
        instant1_lhs = 0
        instant1_rhs = 0
        instant2_lhs = 0
        instant2_rhs = 0
        for p in t[:N1]:
            if  p.position.x <= 50:
                instant1_lhs += 1
            elif p.position.x >= 50:
                instant1_rhs += 1
        for p in t[N1:]:
            if  p.position.x <= 50:
                instant2_lhs += 1
            elif p.position.x >= 50:
                instant2_rhs += 1
        total1_lhs.append(instant1_lhs)
        total1_rhs.append(instant1_rhs)
        total2_lhs.append(instant2_lhs)
        total2_rhs.append(instant2_rhs)


    time = np.linspace(0,steps*timestep,steps+1)

    plt.figure(figsize=[8,5])
    plt.title('Diffusion into different Domain')
    plt.xlabel('Time [Timestep = 0.25]')
    plt.ylabel('Number of particles')
    plt.plot(time,total1_rhs,label='Particle 1 (Small) Right Hand Side')
    plt.plot(time,total2_lhs,label='Particle 2 (Large) Left Hand Side')
    plt.legend(loc='best')

    plt.figure(figsize=[8,5])
    plt.title('Diffusion out - into different Domain')
    plt.xlabel('Time [Timestep = 0.25]')
    plt.ylabel('Number of particles')
    plt.plot(time,total1_lhs,label='Particle 1 (Small) Left Hand Side')
    plt.plot(time,total2_rhs,label='Particle 2 (Large) Right Hand Side')
    plt.legend(loc='best')  

    
def test_pressure_vs_num(nump_list):
    press = [[],[],[]]
    for i in nump_list:
        print('for N =',i)
        cpi, et = test_pressure(i) 
        print('Calculated Pressure', cpi)
        press[0].append(cpi)
        ipi = et/40**2
        print('Pressure of Ideal Gas', ipi)
        press[1].append(ipi)
        vdwpi = et/(40**2-i*np.pi*0.1**2)
        print('Pressure of VDW Gas', vdwpi)
        press[2].append(vdwpi)
        print()
        
    plt.figure(figsize=[8,5])
    plt.title('Pressure as a Function of Particle Number')
    plt.xlabel('Number of Particles')
    plt.ylabel('Pressure')
    plt.plot(nump_list,press[0],label='Calculated Pressure')
    plt.plot(nump_list,press[1],label='Ideal Gas Pressure',linestyle='--')
    plt.legend(loc='best')
    
    plt.figure(figsize=[8,5])
    plt.title('Difference between Ideal and Simulated Gas')
    plt.xlabel('Number of Particles')
    plt.ylabel('Difference in Pressure')
    plt.plot(nump_list,np.array(press[0])-np.array(press[1]),label='Simulated Pressure minus Ideal')
    plt.legend(loc='best')
        
