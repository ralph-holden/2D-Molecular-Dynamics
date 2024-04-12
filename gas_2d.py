#Object-oriented implementation of a hard-disks molecular dynamics simulations

# imports
from itertools import combinations
import numpy as np
import matplotlib.pyplot as plt


class Vector():
    '''2D vectors'''

    def __init__(self, i1,i2):
        '''Initialise vectors with x and y coordinates'''
        self.x = i1
        self.y = i2

    def __add__(self,other):
        '''Use + sign to implement vector addition'''
        return (Vector(self.x+other.x,self.y+other.y))

    def __sub__(self,other):
        '''Use - sign to implement vector "subtraction"'''
        return (Vector(self.x-other.x,self.y-other.y))

    def __mul__(self,number):
        '''Use * sign to multiply a vector by a scaler on the left'''
        return Vector(self.x*number,self.y*number)

    def __rmul__(self,number):
        '''Use * sign to multiply a vector by a scaler on the right'''
        return Vector(self.x*number,self.y*number)

    def __truediv__(self,number):
        '''Use / to multiply a vector by the inverse of a number'''
        return Vector(self.x/number,self.y/number)

    def __repr__(self):
        '''Represent a vector by a string of 2 coordinates separated by a space'''
        return '{x} {y}'.format(x=self.x, y=self.y)

    def copy(self):
        '''Create a new object which is a copy of the current.'''
        return Vector(self.x,self.y)

    def dot(self,other):
        '''Calculate the dot product between two 2D vectors'''
        return self.x*other.x + self.y*other.y

    def norm(self):
        '''Calculate the norm of the 2D vector'''
        return (self.x**2+self.y**2)**0.5
    
    def orthogonal(self):
        '''Gives one of two orthogonal vectors to self '''
        return Vector(self.y,-1*self.x)
        

class Particle():
    def __init__(self, position, momentum, radius, mass):
        self.position = position
        self.momentum = momentum
        self.radius = radius
        self.mass = mass  

    def velocity(self):
        return self.momentum/self.mass # uses __truediv__ in Vector
        
    def copy(self):
        return Particle(self.position, self.momentum, self.radius, self.mass)
    
    def overlap(self, other):
        inter_vector = self.position-other.position
        min_approach = self.radius + other.radius
        return inter_vector.norm() <= min_approach
    
    def true_collision(self,other):
        ''' Checks if one of the three collision cases, thereby excluding the 
        fourth, 'non-collision' overlap case'''
        collision_vector = (self.position - other.position)/((self.position - other.position).norm())
        if self.momentum.dot(collision_vector) < 0 and other.momentum.dot(collision_vector) > 0:
            return True
        if self.momentum.dot(collision_vector)*other.momentum.dot(collision_vector) >= 0:
            if self.momentum.dot(collision_vector) >= 0 and other.momentum.dot(collision_vector) >= 0:
                return self.velocity().dot(collision_vector) < other.velocity().dot(collision_vector)
            if self.momentum.dot(collision_vector) <= 0 and other.momentum.dot(collision_vector) <= 0:
                return self.velocity().dot(collision_vector) < other.velocity().dot(collision_vector)   
            
    def apply_collision(self,p2):
        # collision parameters
        collision_vector = self.position - p2.position
        collision_vector_norm = collision_vector/(collision_vector.norm())
        # project momentum vectors onto collision_vector_norm & orthogonal 
        vec_proj_collision_p1 = collision_vector_norm * self.momentum.dot(collision_vector_norm)
        vec_proj_collision_p2 = collision_vector_norm * p2.momentum.dot(collision_vector_norm)
        vec_proj_orthogonal_p1 = self.momentum - vec_proj_collision_p1
        vec_proj_orthogonal_p2 = p2.momentum - vec_proj_collision_p2
        # rebuild momenta by addingf collision (<- swapped) & orthogonal vectors
        self.momentum = ((self.mass-p2.mass)*vec_proj_collision_p1 + 2*self.mass*vec_proj_collision_p2)/(self.mass+p2.mass) + vec_proj_orthogonal_p1
        p2.momentum = ((p2.mass-self.mass)*vec_proj_collision_p2 + 2*p2.mass*vec_proj_collision_p1)/(self.mass+p2.mass) + vec_proj_orthogonal_p2



class Simulation():
    ''' NOTE: Pressure & RDF calculations take a lot of time! 
    May want to turn off these features by commenting them inside the step method'''
    def __init__(self, particles, box_length, dt, thermostat=False, thermostat_val=0, periodic_boundary = False):
        self.particles = particles
        self.box_length = box_length
        self.dt = dt
        self.trajectory = []
        self.record_state()
        self.thermostat = thermostat
        if thermostat != False:
            self.thermostat_val = thermostat_val * 0.16 * (1/0.4*(2*np.pi)**0.5)*np.e**(-1/(2*thermostat_val**2) * np.linspace(-2,2,2000)**2)
        self.periodic_boundary = periodic_boundary
        
        total_particle_volume = 0
        for p in self.particles:
            total_particle_volume += np.pi * p.radius**2
        self.density = total_particle_volume/box_length**2
   
    def apply_box_collisions(self, particle):
        if particle.position.x+particle.radius >= self.box_length and particle.momentum.x >= 0:
            particle.momentum.x = -1*particle.momentum.x
        if particle.position.x-particle.radius <= 0 and particle.momentum.x <= 0:
            particle.momentum.x = -1*particle.momentum.x
        if particle.position.y+particle.radius >= self.box_length and particle.momentum.y >= 0:
            particle.momentum.y = -1*particle.momentum.y
        if particle.position.y-particle.radius <= 0 and particle.momentum.y <= 0:
            particle.momentum.y = -1*particle.momentum.y
        
    def apply_thermostat(self, particle):
        if particle.position.x+particle.radius >= self.box_length and particle.momentum.x >= 0:
            particle.momentum.x = -1*self.thermostat_val[np.random.randint(0,2000)]
        if particle.position.x-particle.radius <= 0 and particle.momentum.x <= 0:
            particle.momentum.x = self.thermostat_val[np.random.randint(0,2000)]
        if particle.position.y+particle.radius >= self.box_length and particle.momentum.y >= 0:
            particle.momentum.y = -1*self.thermostat_val[np.random.randint(0,2000)]
        if particle.position.y-particle.radius <= 0 and particle.momentum.y <= 0:
            particle.momentum.y = self.thermostat_val[np.random.randint(0,2000)]
            
    def apply_periodic_boundary(self,particle):
        if particle.position.x >= self.box_length and particle.momentum.x >= 0:
            particle.position.x = particle.position.x - self.box_length
        if particle.position.x <= 0 and particle.momentum.x <= 0:
            particle.position.x = particle.position.x + self.box_length
        if particle.position.y >= self.box_length and particle.momentum.y >= 0:
            particle.position.y = particle.position.y - self.box_length
        if particle.position.y <= 0 and particle.momentum.y <= 0:
            particle.position.y = particle.position.y + self.box_length
            
    def apply_particle_collisions(self):
        for p1, p2 in combinations(self.particles,2):
            if p1.overlap(p2) and p1.true_collision(p2):
                p1.apply_collision(p2)
                        
    def record_state(self):
        state = []
        for p in self.particles:
            state.append(p.copy())
        self.trajectory.append(state)

    def step(self):
        self.apply_particle_collisions()
        for p in self.particles:
            if self.thermostat == False and self.periodic_boundary == False:
                self.apply_box_collisions(p)
            if self.periodic_boundary == True and self.thermostat == False:
                self.apply_periodic_boundary(p)
            if self.periodic_boundary == False and self.thermostat != False:
                self.apply_thermostat(p)
            p.position = p.position + (self.dt*p.velocity())
        self.record_state()
      
    
    
### FUNCTIONS FOR DATA ANALYSIS ###

### SPECIALISED FUNCTIONS FOR PRESSURE & RDF ###

def calc_pressure(simulation):
    pressures = []
    for t in simulation.trajectory:
        total_dp = 0
        for p in t:
            if p.position.x+p.radius >= simulation.box_length or p.position.x-p.radius <= 0:
                total_dp += 2*(p.momentum.x**2)**0.5
            if p.position.y+p.radius >= simulation.box_length or p.position.y-p.radius <= 0:
                total_dp += 2*(p.momentum.y**2)**0.5 
        pressures.append(total_dp/(simulation.dt*simulation.box_length*4))
    return pressures
        
def calc_rdf(simulation):
    distances_total = []
    for t in simulation.trajectory:
        distances = []
        for p1, p2 in combinations(t,2):
            distances.append((p1.position-p2.position).norm())
        distances_total.append(distances)
    return distances_total
    
    
