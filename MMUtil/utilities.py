import numpy as np
import math


def get_distance(at1, at2):

    """ Finds the distance between two atoms
    :param at1: (list) a list of xyz coordinates of atom1
    :param at2: (list) a list of xyz coordinates of atom2
    :return: (float) the distance between 2 atoms
    """

    return math.sqrt((at1[0]-at2[0])**2+(at1[1]-at2[1])**2+(at1[2]-at2[2])**2)

