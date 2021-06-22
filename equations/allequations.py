from equations.dim1.pedro import Pedro
from equations.dim1.eqntest import EqKink
from equations.dim1.brusselator import Brusselator
from equations.dim1.fitzhughnagumo import FitzHughNagumo
from equations.dim1.duffing import Duffing
from equations.dim1.cgle import ComplexGinzburgLandau
from equations.dim1.tsh import TuringSwiftHohenberg
from equations.dim0.duffing import Duffing
from equations.dim0.pendulum import Pendulum
from equations.dim0.lorenz import Lorenz
#from .glvortex import GLVortex

ALL_EQUATIONS = {0: {'Duffing': Duffing, 'Pendulum': Pendulum, 'Lorenz': Lorenz},
                1: {'Brusselator' : Brusselator, 
                    'Pedro' : Pedro, 'Eqntest' : EqKink,
                    'FitzHugh-Nagumo' : FitzHughNagumo, 
                    'Duffing': Duffing, 
                    'ComplexGinzburgLandau': ComplexGinzburgLandau,
                    'Turing-Swift-Hohenberg': TuringSwiftHohenberg}}