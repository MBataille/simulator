from equations.pedro import Pedro
from equations.eqntest import EqKink
from equations.brusselator import Brusselator
from equations.fitzhughnagumo import FitzHughNagumo
from equations.duffing import Duffing
from equations.cgle import ComplexGinzburgLandau
from equations.tsh import TuringSwiftHohenberg
from equations.dim0.duffing1d import Duffing1D
from equations.dim0.pendulum import Pendulum
from equations.dim0.lorenz import Lorenz
#from .glvortex import GLVortex

ALL_EQUATIONS = {'Brusselator' : Brusselator, 'Pedro' : Pedro, 'Eqntest' : EqKink,
 'FitzHugh-Nagumo' : FitzHughNagumo, 'Duffing': Duffing, 'ComplexGinzburgLandau': ComplexGinzburgLandau,
 'Turing-Swift-Hohenberg': TuringSwiftHohenberg, 'Duffing1D' : Duffing1D}

ALL_EQUATIONS = {0: {'Duffing': Duffing1D, 'Pendulum': Pendulum, 'Lorenz': Lorenz},
                1: {'Brusselator' : Brusselator, 
                    'Pedro' : Pedro, 'Eqntest' : EqKink,
                    'FitzHugh-Nagumo' : FitzHughNagumo, 
                    'Duffing': Duffing, 
                    'ComplexGinzburgLandau': ComplexGinzburgLandau,
                    'Turing-Swift-Hohenberg': TuringSwiftHohenberg}}