#Cantilever model.
#-------------------------------------------------------------------------------------------
# Tenkte vi kunne bruke denne som et master script for å kjøre Cantiliver modellen vår.
# Da prøver vi å ikke endra sååå mye fra CorotBeam og TwoSimpleBeamModels men istede importere
# fra dem som pakker.
#-------------------------------------------------------------------------------------------
# Pakkene vi trenger
import math
import numpy as np
import matplotlib.pyplot as plt
import CorotBeam_with_TODO as CorotBeam
import matplotlib.animation as anm
from copy import deepcopy
import TwoSimpleBeamModels_with_TODO as sbeam

