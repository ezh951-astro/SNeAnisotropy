
import numpy as np


def load_scalefactors(snapshot=None):
    scalefactors = np.loadtxt(	'/n/home01/ptorrey/IllustrisAuxFiles/Illustris_scalefactors.txt' )

    if snapshot==None:
        val = scalefactors[:]
    else:
        val = scalefactors[snapshot]

    return val




