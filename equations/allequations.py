from .pedro import Pedro
from .eqntest import EqKink
from .brusselator import Brusselator
from .fitzhughnagumo import FitzHughNagumo
from .duffing import Duffing
from .cgle import ComplexGinzburgLandau
#from .glvortex import GLVortex

ALL_EQUATIONS = {'Brusselator' : Brusselator, 'Pedro' : Pedro, 'Eqntest' : EqKink,
 'FitzHugh-Nagumo' : FitzHughNagumo, 'Duffing': Duffing, 'ComplexGinzburgLandau': ComplexGinzburgLandau}