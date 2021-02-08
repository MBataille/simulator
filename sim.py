import numpy as np

import os
import sys

import tkinter as tk

from windows import StartPage, MainPage

from eqntest import EqTest
from pedro import Pedro

class SimApp(tk.Tk):
	
	def __init__(self, *args, **kwargs):
		
		tk.Tk.__init__(self, *args, **kwargs)

		## icon
		# t.Tk.iconbitmap(self, default='clienticon.ico')
		tk.Tk.wm_title(self, 'Simulator')
		self.geometry('1024x768')
		container = tk.Frame(self)

		container.pack(side='top', fill='both', expand=True)

		container.grid_rowconfigure(0, weight=1)
		container.grid_columnconfigure(0, weight=1)

		self.frames = {}

		pages = (StartPage, MainPage)

		self.setEq(Pedro)
		self.eq.setInitialConditionKink()

		for F in pages:
			frame = F(container, self)
			self.frames[F] = frame
			frame.grid(row=0, column=0, sticky='nsew')

		self.show_frame(StartPage)

		### load eq and init cond


	def setEq(self, eq):
		### eq should be a string or sth
		self.eq = eq()


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