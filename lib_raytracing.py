import raytracing.thorlabs as thorlabs
from .imagingpath import *
import raytracing.eo as eo
import raytracing.olympus as olympus
path = ImagingPath()
path.label = "Demo #1: lens f = 5cm, infinite diameter"
path.append(Space(d=10))
path.append(Lens(f=5))
path.append(Space(d=10))
path.display()