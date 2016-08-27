"""
Rutherford Simulation
"""

import random
from math import sqrt, sin, cos, atan, radians, degrees

class Point(object):
    """Point class for points"""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __sub__(self, other):
        if(isinstance(other, Point)):
            return Vector(self.x - other.x,
                          self.y - other.y,
                          self.z - other.z)

    def __mul__(self, other):
        if(isinstance(other, Vector) or isinstance(other, Point)):
            return self.x * other.x + self.y * other.y + self.z * other.z

    def __add__(self, other):
        if(isinstance(other, Point)):
            return Point(self.x + other.x, self.y + other.y, self.z + other.z)
        elif(isinstance(other, Vector)):
            return Point(self.x + other.x, self.y + other.y, self.z + other.z)

    def distance(self, other):
        return sqrt((self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2)

    def __str__(self):
        return 'Coordinates are: (%f %f %f)' % (self.x, self.y, self.z)

class Vector(object):
    """Vector class for vector operations"""
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __add__(self, other):
        if(isinstance(other, Vector)):
            return Vector(self.x + other.x, self.y + other.y, self.z + other.z)

    def __mul__(self, other):
        if(isinstance(other, int) or isinstance(other, float)):
            return Vector(self.x * other, self.y * other, self.z * other)
            #Dot product
        if(isinstance(other, Vector)):
            return self.x * other.x + self.y * other.y + self.z * other.z

    def __sub__(self, other):
        return Vector(self.x - other.x,
                     self.y - other.y,
                     self.z - other.z)

    def calc_leng(self):
        return sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)

    def normalize(self):
        leng = sqrt(self.x ** 2 + self.y ** 2 + self.z ** 2)
        return Vector(self.x / leng, self.y / leng, self.z / leng)

    def cross_product(self, other):
        return Vector(self.y * other.z - self.z * other.y, self.z * other.x - self.x * other.z, self.x * other.y - self.y * other.x)

    def __str__(self):
        return 'The Vector is (%f %f %f)\n' % (self.x, self.y, self.z)

class Plane(object):
    def __init__(self, point, vector):
        self.point = point
        self.vector = vector

class Sphere(object):
    """Description of gold atoms"""
    radius = 1.4 * 10

    def __init__(self, position):
        self.position = position

    def __str__(self):
        return 'X:  %d  Y:  %d  Z:  %d' % (self.position.x, self.position.y, self.position.z)

class Particle(object):
    """Description of particles"""
    def __init__(self, point, vector):
        self.point = point
        self.vector = vector

    def __str__(self):
        return 'The particle is at the position:  %f  %f  %f\nAnd it has the direction vector:  %d  %d  %d\n' % (self.point.x, self.point.y, self.point.z, self.vector.x, self.vector.y, self.vector.z)

def compute_deviated_ray(b, A, B):
    #Compute the angle in radians. Don't forget to convert later.
    tg_half_angle = predefined_rutherford_const / b
    tg_angle = 2 * tg_half_angle / (1 - tg_half_angle * tg_half_angle)
    angle = atan(tg_angle)
    dev_angle = radians(90) - angle #1.57 = 90 degrees in radians

    dev_axis = A - B
    second_vector = Point(B.x, B.y, B.z + 10) - B

    dev_axis = dev_axis.normalize()
    second_vector = second_vector.normalize()

    normal_vector = dev_axis.cross_product(second_vector)
    normal_vector = normal_vector.normalize()

    #Compute the new ray of the  particle
    vrot = dev_axis * cos(dev_angle) + normal_vector.cross_product(dev_axis) * sin(dev_angle)
    vrot.normalize()

    return vrot

def intersect_with_plane(r, p1, plane, flag):
    global intersected_with_gold_atoms
    global went_trough

    if(flag == False):
         intersected_with_gold_atoms += 1
         weird_rays.append(r)
    else:
        went_trough +=1

    r = r.normalize()

    if(float(r * plane.vector) > 0.0000001):
        t = ((plane.point - p1) * plane.vector) / (r * plane.vector)
        plane_intersection_point = p1 + (r * t)
        plane_intersections.append(plane_intersection_point)
        if(flag == False):
            int_with_atoms.append(plane_intersection_point)
    else:
        global not_intersected_with_the_plane
        not_intersected_with_the_plane += 1



# Z * e * e / (2 * pi * epsilon0 * m * v * v)
predefined_rutherford_const = ((97.0 * 2.0 * (1.602 ** 2)) / (2.0 * 3.14 * 8.854 * 6.644 * (15 ** 2))) * (10 ** -11)
z_coord = 1000
n = 1000000
gold_atoms = []
plane_intersections = []
int_with_atoms = []
weird_rays = []
distance_between_atoms = 408
radius_of_atom = 136
radius_of_nucleus = 1.4 * (10 ** (-2))
max_size_of_atom_grid = 10 ** 4
plane = Plane(Point(0, 0, z_coord * 2), Vector(0, 0, 1))

intersected_with_gold_atoms = 0
went_trough = 0
not_intersected_with_the_plane = 0
number_of_golds = 0

for debug in xrange(0, n):
    particle = Particle(Point(random.randint(0,max_size_of_atom_grid), random.randint(0,max_size_of_atom_grid), 0), Vector(0, 0, 1))

    # Translated position of the particle
    t_x = particle.point.x - 68
    t_y = particle.point.y - 68
    i_x = t_x / 544
    i_y = t_y / 544

    #Create the neighbouring spheres
    gold_atoms = [Sphere(Point(i_x * 544 + 68, i_y * 544 + 68, z_coord)),
                  Sphere(Point(i_x * 544 + 68, (i_y + 1) * 544 + 68, z_coord)),
                  Sphere(Point((i_x + 1) * 544 + 68, i_y * 544 + 68, z_coord)),
                  Sphere(Point((i_x + 1) * 544 + 68, (i_y + 1) * 544 + 68, z_coord))]

    #Find the closest spheres
    minDist = 10 ** 100
    for atom in gold_atoms:
        dist = particle.point.distance(atom.position)
        if(dist < minDist):
            minDist = dist
            final_atom = atom

    #Calculate b for Rutherford formulae
    new_point = Point(particle.point.x, particle.point.y, z_coord - (predefined_rutherford_const * (10 ** 12)))
    b = new_point.distance(final_atom.position)

    #Check if the ray intersects the final_atom
    L = final_atom.position - particle.point
    tc = L * particle.vector

    #If ray does not intersect the sphere compute deviation
    if(tc < 0):
        dev_ray = compute_deviated_ray(b, new_point, final_atom.position)
        intersect_with_plane(dev_ray, new_point, plane, True)
    else:
        d = float(sqrt(L * L - tc * tc))
        #If ray does not intersect the sphere compute deviation
        if(d > final_atom.radius or d < 0.0):
            dev_ray = compute_deviated_ray(b, new_point, final_atom.position)
            intersect_with_plane(dev_ray, new_point, plane, True)
        else:
            #Compute the intersection point p1
            thc = sqrt(final_atom.radius * final_atom.radius - d * d)
            t1 = tc - thc
            t2 = tc + thc
            p1 = particle.point + particle.vector * t1

            #Calculate normal vector
            normal = p1 - final_atom.position
            normal = normal.normalize()

            #If the incident ray is paralel to the normal
            i = p1 - particle.point
            i = i.normalize()
            if(i.x == normal.x and i.y == normal.y and i.z == normal.z):
                r = normal
            else:
                #Calculate reflected vector
                r = i - normal * (i * normal) * 2
            """
            print'--> Start <--'
            print p1
            print '\n'
            print final_atom.position
            print '\n'
            print i
            print '\n'
            print normal
            print '\n'
            print r
            print '--> End <--' 
            """
            #Check if the reflected ray intersects the final plane
            intersect_with_plane(r, p1, plane, False)

foo = open('foo.txt', 'w')
for points in plane_intersections:
    print ('The particle intersected the plane at:  %f  %f  %f\n') % (points.x, points.y, points.z)
    foo.write('%f, %f\n' % (points.x, points.y))

for other_points in int_with_atoms:
    print ('The particle bounced and ended up at:  %f  %f  %f\n') % (other_points.x, other_points.y, other_points.z)

for ray in weird_rays:
    print ('Weird ray:  %f  %f  %f\n') % (ray.x, ray.y, ray.z)

print ('%d particles bounced of gold atoms') % (intersected_with_gold_atoms)
print ('%d particles went trough with little deviation') % (went_trough)
print ('%d particles not intersected with the plane') % (not_intersected_with_the_plane)
print (number_of_golds)
