import numpy as np

class Vector:
    def ortho_cross(self, a, b):
        """
        Takes a cross of two vectors in orthonormal coordinates.

        Parameters
        ----------
        a : list of int
            The first vector
        b : list of int
            The second vector
        """
        [a1, a2, a3] = a
        [b1, b2, b3] = b
        return np.array([a2*b3 - a3*b2, a3*b1 - a1*b3, a1*b2 - a2*b1], dtype=np.float32)

    def ortho_dot(self, a, b):
        """
        Takes a dot of two vectors in orthonormal coordinates.

        Parameters
        ----------
        a : list of int
            The first vector
        b : list of int
            The second vector
        """
        [ar, at, ap] = a
        [br, bt, bp] = b
        return (ar*br + at*bt + ap*bp)

    def get_magnitude(self, a):
        """
        Takes the magnitude of a vector

        Parameters
        ----------
        a : list of int
            The vector
        """
        [ar, at, ap] = a
        return np.sqrt(ar*ar + at*at + ap*ap)

    def growth_prec(self, dvdt, v):
        """
        Calculates the growth and precession parts for an arbitrary vector given its vector time derivative.
        The precession part is calculated for the xy components only - Bryance
        """
        [ax, ay, az] = dvdt
        [vx, vy, vz] = v

        vmag = get_magnitude(v)

        # growth part
        growth = ortho_dot(v, dvdt) / vmag

        # precession part
        prec = (vx*ay - vy*ax) / (vx**2 + vy**2)

        return [growth, prec]


def wrap_phase(phi, positive=False):
    """
    Wraps phase of angular quantities, if positive is true they are confined to [0, 2pi] else they are confined to [-pi, pi]

    Parameters
    ----------
    positive : bool
        Confines angular quantities to [0, 2pi]
    """
    if positive:
        return (phi) % (2*np.pi)
    else:
        return (phi + np.pi) % (2*np.pi) - np.pi