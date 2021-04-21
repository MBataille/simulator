import numpy as np

import os
import sys
import platform

import tkinter as tk
from tkinter import ttk

from windows import StartPage, MainPage

from equations.allequations import ALL_EQUATIONS


DATAFOLDER = 'data/'

class SimApp(tk.Tk):
	
	def __init__(self, *args, **kwargs):
		
		tk.Tk.__init__(self, *args, **kwargs)

		## icon
		# t.Tk.iconbitmap(self, default='clienticon.ico')
		tk.Tk.wm_title(self, 'Simulator v0.1')
		self.geometry('512x768')
		self.style = ttk.Style()
		system = platform.system()
		if system == 'Darwin': # mac
			theme = 'aqua'
		elif system == 'Linux':
			theme = 'clam'
		else:
			theme = 'default'
		self.style.theme_use(theme) 
		container = tk.Frame(self)

		container.pack(side='top', fill='both', expand=True)

		container.grid_rowconfigure(0, weight=1)
		container.grid_columnconfigure(0, weight=1)

		self.frames = {}

		pages = (StartPage, MainPage)

		self.initDataFolders()

		self.setEq('Brusselator')
		self.eq.setInitialConditionZero()

		for F in pages:
			frame = F(container, self)
			self.frames[F] = frame
			frame.grid(row=0, column=0, sticky='nsew')

		self.show_frame(StartPage)
		self.activePlot = None

		### load eq and init cond


	def initDataFolders(self):
		# check if data folder exists
		if not os.path.exists(DATAFOLDER):
			os.mkdir(DATAFOLDER)

		for eqname in ALL_EQUATIONS:
			if not os.path.exists(DATAFOLDER + eqname):
				os.mkdir(DATAFOLDER + eqname)

	def getEqInitConds(self):
		self.initConds =  self.eq.getInitialConditions(), self.eq.getSavedStatesNames()
		return self.initConds[0] + self.initConds[1]

	def initCond_isMethod(self, indx):
		n_methods = len(self.initConds[0])
		if indx < n_methods:
			return True
		return False

	def setEqNi(self, Ni): # eventually interpolate
		self.eq.setNi(Ni)

	def getEqInitCond_N(self): # assuming init cond is already loaded
		return self.eq.Ni

	def setEqInitCond(self, indx):
		if self.initCond_isMethod(indx):
			getattr(self.eq, 'setInitialCondition' + self.initConds[0][indx])()
		else:
			n_methods = len(self.initConds[0])
			self.eq.loadState(self.initConds[1][indx - n_methods])
	
		self.current_initcond_indx = indx

	def resetEqInitCond(self):
		self.setEqInitCond(self.current_initcond_indx)

	def setEq(self, name):
		### eq should be a string or sth
		self.eq = ALL_EQUATIONS[name]()

	def getBoundaryCondition(self):
		return self.eq.getBoundaryCondition()

	def setBoundaryCondition(self, indx):
		allbcs = self.getListBoundaryConditions()
		if indx < len(allbcs):
			self.eq.setBoundaryCondition(allbcs[indx])
			return True
		return False

	def getListBoundaryConditions(self):
		return self.eq.getAllBoundaryConditions()

	def show_frame(self, cont):
		# Raise frame to the top
		frame = self.frames[cont]
		frame.tkraise()
		frame.activate()


def main():
	app = SimApp()
	app.mainloop()

if __name__ == '__main__':
	main()