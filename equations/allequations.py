from equations.dim1.pedro import Pedro
from equations.dim1.eqntest import EqKink
from equations.dim1.brusselator import Brusselator
from equations.dim1.fitzhughnagumo import FitzHughNagumo
from equations.dim1.duffing import Duffing as Duffing1
from equations.dim1.cgle import ComplexGinzburgLandau
from equations.dim1.tsh import TuringSwiftHohenberg
from equations.dim0.duffing import Duffing as Duffing0
from equations.dim0.pendulum import Pendulum
from equations.dim0.lorenz import Lorenz
from equations.dim0.homo_bif import HomoclinicBifurcation
from equations.dim0.vdp import VanDerPol
from equations.dim1.rene import Rene
from equations.dim1.lle import LugiatoLefeverEquation
#from .glvortex import GLVortex

ALL_EQUATIONS = {0: {'Duffing': Duffing0, 'Pendulum': Pendulum, 'Lorenz': Lorenz,
                    'HomoclinicBifurcation': HomoclinicBifurcation,
                    'VanDerPol': VanDerPol},
                1: {'Brusselator' : Brusselator, 
                    'Pedro' : Pedro, 'Rene': Rene,'Eqntest' : EqKink,
                    'FitzHugh-Nagumo' : FitzHughNagumo, 
                    'Duffing': Duffing1, 
                    'ComplexGinzburgLandau': ComplexGinzburgLandau,
                    'Turing-Swift-Hohenberg': TuringSwiftHohenberg,
                    'LugiatoLefever': LugiatoLefeverEquation}}