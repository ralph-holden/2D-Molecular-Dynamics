            if self.periodic_boundary == True:
                shadow_p1 = Particle(Vector(-1,-1),Vector(0,0),0.1,1)
                shadow_p2 = Particle(Vector(-1,-1),Vector(0,0),0.1,1)
                if p1.position.x+p1.radius >= self.box_length:
                    shadow_p1 = p1.copy()
                    shadow_p1.position.x = shadow_p1.position.x - self.box_length
                if p2.position.x-p2.radius >= self.box_length:
                    shadow_p2 = p2.copy()
                    shadow_p2.position.x = shadow_p2.position.x - self.box_length
                if p1.position.x-p1.radius <= self.box_length:
                    shadow_p1 = p1.copy()
                    shadow_p1.position.x = shadow_p1.position.x + self.box_length         
                if p2.position.x-p2.radius <= self.box_length:
                    shadow_p2 = p2.copy()
                    shadow_p2.position.x = shadow_p2.position.x + self.box_length
                if p1.position.y+p1.radius >= self.box_length:
                    shadow_p1 = p1.copy()
                    shadow_p1.position.y = shadow_p1.position.y - self.box_length
                if p2.position.y-p2.radius >= self.box_length:
                    shadow_p2 = p2.copy()
                    shadow_p2.position.y = shadow_p2.position.y - self.box_length
                if p1.position.y-p1.radius <= self.box_length:
                    shadow_p1 = p1.copy()
                    shadow_p1.position.y = shadow_p1.position.y + self.box_length        
                if p2.position.y-p2.radius <= self.box_length:
                    shadow_p2 = p2.copy()
                    shadow_p2.position.y = shadow_p2.position.y + self.box_length
                #if shadow_p2.overlap(p1) and shadow_p2.true_collision(p1):
                #    shadow_p2.apply_collision(p1)
                #    p2.momentum = shadow_p2.momentum
                if shadow_p1.overlap(p2) and shadow_p1.true_collision(p2):
                    shadow_p1.apply_collision(p2)
                    p1.momentum = shadow_p1.momentum

### EXTRA CODE ###

# FOR COLLISIONS ACROSS PERIODIC BOUNDARY

    def create_shadows(self):
        pns_list = []
        for p in self.particles:
            pns_list.append(p.copy())
        for p2 in self.particles:
            count = 0
            if p2.position.x-p2.radius >= self.box_length:
                shadow_p2 = p2.copy()
                shadow_p2.position.x = shadow_p2.position.x - self.box_length
                count += 1
            if p2.position.x-p2.radius <= self.box_length:
                shadow_p2 = p2.copy()
                shadow_p2.position.x = shadow_p2.position.x + self.box_length
                count += 1
            if p2.position.y-p2.radius >= self.box_length:
                shadow_p2 = p2.copy()
                shadow_p2.position.y = shadow_p2.position.y - self.box_length
                count += 1
            if p2.position.y-p2.radius <= self.box_length:
                shadow_p2 = p2.copy()
                shadow_p2.position.y = shadow_p2.position.y + self.box_length
                count += 1
            if count >= 1:
                pns_list.append(p2)
        return pns_list
                
    def apply_boundary_collisions(self):
        for p1, p2 in combinations(self.create_shadows(),2):
            if p1.overlap(p2) and p1.true_collision(p2):
                p1.apply_collision(p2)

    def apply_collision(self,p1,p1_shadow,p2):
        if p1_shadow.overlap(p2) and p1_shadow.both_toward(p2):
            # collision parameters
            collision_vector = p1_shadow.position - p2.position
            collision_vector_norm = collision_vector/(collision_vector.norm())
            # project momentum vectors onto collision_vector_norm & orthogonal 
            vec_proj_collision_p1 = collision_vector_norm * p1.momentum.dot(collision_vector_norm)
            vec_proj_collision_p2 = collision_vector_norm * p2.momentum.dot(collision_vector_norm)
            vec_proj_orthogonal_p1 = p1.momentum - vec_proj_collision_p1
            vec_proj_orthogonal_p2 = p2.momentum - vec_proj_collision_p2
            # rebuild momenta by addingf collision (<- swapped) & orthogonal vectors
            p1.momentum = ((p1.mass-p2.mass)*vec_proj_collision_p1 + 2*p1.mass*vec_proj_collision_p2)/(p1.mass+p2.mass) + vec_proj_orthogonal_p1
            p2.momentum = ((p2.mass-p1.mass)*vec_proj_collision_p2 + 2*p2.mass*vec_proj_collision_p1)/(p1.mass+p2.mass) + vec_proj_orthogonal_p2
    
    def apply_periodic_boundary_collision(self,particle):
        for p1, p2 in combinations(self.particles,2):
            if p1.position.x + p1.radius >= self.box_length and p1.momentum.x >= 0:
                p1_shadow = p1.copy()
                p1_shadow.position.x = p1.position.x - self.box_length
                self.apply_collision(p1,p1_shadow,p2)
            if p1.position.x - p1.radius<= 0 and p1.momentum.x <= 0:
                p1_shadow = p1.copy()
                p1_shadow.position.x = p1.position.x + self.box_length
                self.apply_collision(p1,p1_shadow,p2)
            if p1.position.y + p1.radius >= self.box_length and p1.momentum.y >= 0:
                p1_shadow = p1.copy()
                p1_shadow.position.y = p1.position.y - self.box_length
                self.apply_collision(p1,p1_shadow,p2)
            if p1.position.y - p1.radius <= 0 and p1.momentum.y <= 0:
                p1_shadow = p1.copy()
                p1_shadow.position.x = p1.position.x + self.box_length 
                self.apply_collision(p1,p1_shadow,p2)
            
