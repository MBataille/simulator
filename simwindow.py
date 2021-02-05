from tkinter import *
from tkinter import ttk
from paramwindow import Parameter,ParameterWindow, updateF
from eqntest import EqTest
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)

# Implement the default Matplotlib key bindings.
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import matplotlib.animation as animation


import numpy as np



class SimWindow:
	def __init__(self, root, eq):
		self.root = root
		self.eq = eq
		initCond = eq.initCond
		self.fig = Figure(figsize=(5, 4), dpi=100)


		### Draw inital condition
		# self.line, = self.fig.add_subplot(111).plot(self.x, 2 * np.sin(2 * np.pi * (self.x)))
		self.line, = self.fig.add_subplot(111).plot(self.eq.x[1:-1], self.eq.initCond[1:-1])

		self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
		self.canvas.draw()
		self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

		self.toolbar = NavigationToolbar2Tk(self.canvas, self.root)
		self.toolbar.update()
		self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

		self.canvas.mpl_connect('key_press_event', self.on_key_press)

		self.qbtn = ttk.Button(master=self.root, text='Quit', command=self._quit)
		self.qbtn.pack(side=BOTTOM)

		self.ani = animation.FuncAnimation(self.fig, self.animate, np.arange(1, 200), interval=25, blit=False)

	def animate(self, i):
		self.eq.updateX()
		k = (i-1) % 10
		if k == 0:
			#print('restart')
			self.eq.solve()
			self.eq.initCond = self.eq.sol[:, -1]
		self.line.set_ydata(self.eq.sol[:, k][1:-1])
		return self.line,

	def _quit(self):
		self.root.quit()
		self.root.destroy()

	def on_key_press(self, event):
		print('You pressed {}'.format(event.key))
		key_press_handler(event, self.canvas, self.toolbar)

def main():
	root = Tk()

	
	eq = EqTest()
	eq.setInitialConditionKink()
	params = eq.parameters
	sw = SimWindow(root, eq)
	newWindow = Toplevel(root)
	pw = ParameterWindow(newWindow, '400x500', params, updateF)

	root.mainloop()

if __name__ == '__main__':
	main()

