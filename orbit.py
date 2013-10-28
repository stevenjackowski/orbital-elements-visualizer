__author__ = 'sajackow'

import numpy as np


class Orbit:
    """
    The Orbit class contains the orbital elements of that orbit and tools
    necessary to generate r-vectors of the orbit for plotting.
    """
    def __init__(self, ecc, inc, raan, aop, a, resolution=250):
        """
        Initializes Orbit with the given eccentricity, inclination,
        right ascension, argument of perigee, and semi-major axis.

        """
        # standard gravitational parameter for Earth
        self.mu = 398600.0

        # an array of angles that will be used to iterate over the entire
        # orbit, resolution controls how many data points will be created
        self.resolution = resolution
        self.thetas = np.linspace(0.0, 2*np.pi, self.resolution)

        # angles are converted to radians to avoid use of cosd/sind
        self.ecc = ecc
        self.inc = inc*np.pi/180
        self.raan = raan*np.pi/180
        self.aop = aop*np.pi/180
        self.a = a

        self.r_perifocal = self.update_perifocal(self.a, self.ecc)
        self.q_mat = self.calculate_transform(self.inc, self.raan, self.aop)
        self.r_geocentric = self.update_geocentric()

    def get_rvectors(self):
        """
        Returns the array of geocentric position vectors, which has a size of
        3 x self.resolution.
        """
        return self.r_geocentric

    def update(self, ecc, inc, raan, aop, a):
        """
        This method should be used to update the orbit once it has already
        been initialized. Does not have a return value.
        """

        # check for unneeded calculations of the perifocal frame vectors
        if not (self.a + 0.1 >= a >= self.a - 0.1) or not \
                (self.ecc + 0.001 >= ecc >= self.ecc - 0.001):
            self.r_perifocal = self.update_perifocal(a, ecc)

        inc = inc*np.pi/180
        raan = raan*np.pi/180
        aop = aop*np.pi/180

        # check for unneeded calculations of the transformation matrix
        if not (self.inc + 0.001 >= inc >= self.inc - 0.001) or \
            not (self.raan + 0.001 >= raan >= self.raan - 0.001) or \
                not (self.aop + 0.001 >= aop >= self.aop - 0.001):
            self.q_mat = self.calculate_transform(inc, raan, aop)

        self.a = a
        self.ecc = ecc
        self.inc = inc
        self.raan = raan
        self.aop = aop
        self.r_geocentric = self.update_geocentric()

    def update_perifocal(self, a, ecc):
        """
        Returns a new array containing r-vectors of the orbit in perifocal
        frame. This information is later transformed into the geocentric
        equatorial frame.
        """
        r_perifocal = np.zeros(self.resolution*3).reshape(3, self.resolution)
        for i, theta in enumerate(self.thetas):
            cos_theta = np.cos(theta) # avoid an extra cos calculation
            r_perifocal[:, i] = a*(1.0 - ecc**2)/(1 + ecc*cos_theta)*np.array(
                [cos_theta, np.sin(theta), 0])
        return r_perifocal

    def update_geocentric(self):
        """
        Returns a new array containing r-vectors of the orbit in a perifocal
        frame. This method should only be called after q_mat and r_perifocal
        have already been updated.
        """
        r_geocentric = np.zeros(self.resolution*3).reshape(3, self.resolution)
        for i in range(self.resolution):
            r_geocentric[:, i] = np.dot(self.q_mat, self.r_perifocal[:, i])
        return r_geocentric

    def calculate_transform(self, inc, raan, aop):
        """
        Calculates the transformation matrix, q_mat, which can transform
        vectors from the perifocal frame to the geocentric equatorial frame.
        """
        # calculate sin and cos values before hand to prevent unnecessary
        # calls to sin and cos
        sin_aop = np.sin(aop)
        cos_aop = np.cos(aop)
        sin_inc = np.sin(inc)
        cos_inc = np.cos(inc)
        sin_raan = np.sin(raan)
        cos_raan = np.cos(raan)

        r3_aop = np.array([[cos_aop, sin_aop, 0], [-sin_aop, cos_aop, 0],
                           [0, 0, 1]])
        r1_inc = np.array([[1, 0, 0], [0, cos_inc, sin_inc], [0, -sin_inc,
                                                              cos_inc]])
        r3_raan = np.array([[cos_raan, sin_raan, 0], [-sin_raan, cos_raan,
                                                      0], [0, 0, 1]])
        return np.dot(np.dot(r3_aop, r1_inc), r3_raan).T


