__doc__ = """"""

from elastica.surface.surface_base import SurfaceBase
import numpy as np
from numpy.testing import assert_allclose
from elastica.utils import Tolerance


class Sleeve(SurfaceBase):
    def __init__(self, position_collection: np.ndarray, normal_collection: np.ndarray):
        """
        Sleeve initializer.

        Parameters
        ----------
        position_collection: np.ndarray
            position of the centerline of the sleeve.
            Expect (n_nodes,3)-shaped array.
        normal_collection: np.ndarray
            The normal vectors of the sleeve, must be normalized.
            Expect (n_nodes,3)-shaped array.
        """
        self.position_collection = position_collection
        self.normal_collection = normal_collection
        
        N = position_collection.shape[1]
        self.response_force_sleeve = np.zeros((3,N-1))
        self.displacement_sleeve= np.zeros((3,N-1))
        self.response_couple_sleeve= np.zeros((3,N-1))
        self.response_couple_local = np.zeros((3,N-1))
        self.rotation_sleeve= np.zeros((3,N-1))


