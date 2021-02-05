from tkinter import *
from tkinter import ttk
import tkinter.font as tkFont




class Parameter:
	def __init__(self, name, val):
		self.name = name
		self.val = val

def updateF(*args):
	print('### UPDATED ###')
	for p in params:
		print(p.val.get())

class ParameterWindow:
	def __init__(self, root, geometry, parameters, update):


		bigFontStyle = tkFont.Font(family='Helvetica', size='26')
		medFontStyle = tkFont.Font(family='Helvetica', size='20')

		self.root = root
		self.root.title('Parameters')
		self.root.geometry(geometry) # string: sizexsize
		self.parameters = parameters
		self.update = update

		self.titlelbl = ttk.Label(self.root, text='Parameters', width=15, font=bigFontStyle)
		self.paramslbls = [ttk.Label(self.root, text=self.parameters[p].name, width=6, font=medFontStyle) for p in self.parameters]
		self.paramsvals = [ttk.Entry(self.root, textvariable=self.parameters[p].val, width=4, font=medFontStyle) for p in self.parameters]

		self.titlelbl.grid(column=0, row=0, columnspan=3)
		for i in range(len(self.paramslbls)):
			self.paramslbls[i].grid(column=0, row=i+1, columnspan=2)
			self.paramsvals[i].grid(column=3, row=i+1)

		#self.root.bind('<Return>', self.update)


class ParameterPanel(ttk.Frame):
	def __init_(self, master, **kwargs):
		ttk.Frame.__init__(self, master, **kwargs, width=200, height=600)
		self.master = master
		self.params = params
		self.updateFunc = updateF
		self.initWindow()

	def setParameters(self, params):
		self.params = params

	def setUpdateFunc(self, updateFunc):
		self.updateFunc = updateFunc

	def initWindow(self):
		self.master.title('Parameters')
		#self.pack(fill='both', expand=1)

		self.titlelbl = ttk.Label(self, text='Parameters', width=30)
		self.paramslbls = [ttk.Label(self, text=p.name, width=20) for p in self.params]
		self.paramsvals= [ttk.Entry(self, textvariable=p.val, width=5) for p in self.params]

		self.titlelbl.grid(column=0, row=0, columnspan=3)
		for i in range(len(params)):
			self.paramslbls[i].grid(column=0, row=i+1, columnspan=2)
			self.paramsvals[i].grid(column=3, row=i+1)

		self.master.bind('<Return>', self.update)

#root.title("Parameters")


def main():
	root = Tk()

	params = [Parameter('alpha', DoubleVar()), Parameter('gamma', DoubleVar()), Parameter('sigma', DoubleVar())]
	pw = ParameterWindow(root, '400x720', params, updateF)
	root.mainloop()

if __name__ == '__main__':
	main()

# ppanel.setParameters(params)
# ppanel.setUpdateFunc(updateF)
# ppanel.initWindow()

# root.columnconfigure(0, weight=1)
# root.rowconfigure(0, weight=1)

# for i in range(len(params)+1):
# 	mainFrame.rowconfigure(i, weight=1)

# mainFrame.columnconfigure(0, weight=3)
# mainFrame.columnconfigure(1, weight=1)

