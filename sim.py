import numpy as np

import os
import sys

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
		tk.Tk.wm_title(self, 'Simulator')
		self.geometry('512x768')
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

	def setEqInitCond(self, indx):
		n_methods = len(self.initConds[0])
		if indx >= n_methods: # initial condition is a saved state
			self.eq.loadState(self.initConds[1][indx - n_methods])
#			initParams = self.eq.initParams
#			name = self.eq.name
#			self.eq = 
		else: # intial condition comes from a method
			getattr(self.eq, 'setInitialCondition' + self.initConds[0][indx])()


	def setEq(self, name):
		### eq should be a string or sth
		self.eq = ALL_EQUATIONS[name]()


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