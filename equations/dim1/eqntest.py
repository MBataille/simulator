from tkinter import *
from equations.equation import Equation, Parameter
import numpy as np

#######
# Eq:
#	\partial_t U = alpha + epsilon * U - U ^3 + Laplace(U)
#######
# In order to create an equation:
#
# 	1. Initialize params by specifying their names and
#  initial values in the 'initParams' dictionnary (see below).
#  Note: It must always contain parameters 'dt' and 'dx'.
#
#   2. Then call the Equation constructor 'Equation.__init__()'
#  and MUST include the following two arguments:
#
#		 a) Name wich will appear in the start screen when you
# 		select the eq.
#
#		 b) initParams dictionnary (see point 1)
#
#	   It MAY also include the following keyword arguments:
#	   if you're in a hurry you can skip directly to iv)
#
#		 i) dim (integer). Number of spatial dimensions, 
# 		the default is dim=1 and at the moment the code 
# 		doesn't support other values :( so you can leave it
# 		as it is.
#
#		 ii) N (integer). Number of points for discretization.
# 		The default is N=200, this value will
#		appear in the start screen when you select this equation,
# 		but it can be changed by the user in the start screen.
#
#		iii) isComplex (bool). Whether the data is complex or not. By
#		default it is set to False and at the moment this is not 
# 		correctly implemented :( so don't bother changing it.
#		
#		iv) n_fields (integer). Number of fields to simulate. By default
#		n_fields=1. If you change this number, you'll have to specify
#		the name of the fields (see v))
#
#		v) fieldNames (list of str). List of names for each field. By
#		default fieldNames=['u']
#	
#	3. Define the rhs (Right Hand Side) of the equation (see the right hand side
# 	of the equation at the top of this file for reference). It must be
#	called rhs and accept arguments (t, u). Where t is the current time
# 	and u is a vector containing the field(s) to integrate. This function
#	must return the derivative of u.
#
#	4. Add an image of the equation in the images folder. Note that its filename
#	must be 'name.jpg' where 'name' is the name you set for the equation (see 2.a)
#	
# 	5. (Optional) If you want to add specific initial conditions for an equation
#	You may define it as a method that begins with 'setInitialCondition' such as
#	the 'setInitialConditionKink' in this case. This way the code will understand
#	that this methods sets an initial condition and will display it as an option
#	on the start screen.
#
#	Note: If you want to add an equation with multiple fields, you can look at
#	brusselator.py to see an example.

class EqKink(Equation):
	''' In order to cr
	'''
	def __init__(self):
		initParams = { 'dt' : 0.1,
							'dx': 0.5, 
							'eps': 1, 
							'alpha': 0}

		Equation.__init__(self, 'Eqntest', initParams, dim=1)
		

	def rhs(self, t, u):
		v = self.getCurrentParams() # always get current values of parameters
		dx = v['dx']; eps = v['eps']; alpha = v['alpha'] # and extract them like this

		# calculate the derivative:
		dudt = alpha + eps*u - u**3  + self.Laplace1D(u)
		return dudt

	def setInitialConditionKink(self):
		x = self.getX() # get the current value of x
		eps = self.getParam('eps')

		# calculate the initial condition
		kink = np.sqrt(eps) * np.tanh(x * np.sqrt(eps/2))

		# set the initial condition.
		self.setInitialCondition(kink)