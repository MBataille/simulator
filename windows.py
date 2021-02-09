import tkinter as tk
from tkinter import ttk

from PIL import Image, ImageTk

import numpy as np

import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib import style

from scipy.interpolate import interp1d

LARGE_FONT = ('Times', 14)
MED_FONT = ('Times', 10)
style.use('ggplot')

ST_ROWS = 60
RAW_N = 30

class LineBuilder:
    def __init__(self, line):

        self.line = line
        self.xs = list(line.get_xdata())
        self.ys = list(line.get_ydata())
        self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

    def update(self):
        pass

    def __call__(self, event):

        print('click', event)
        if event.inaxes != self.line.axes: return
        self.xs.append(event.xdata)
        self.ys.append(event.ydata)
        self.line.set_data(self.xs, self.ys)
        self.line.figure.canvas.draw()

class ImageWindow(tk.Frame):
    def __init__(self, master, filename):
        tk.Frame.__init__(self, master)
        self.master = master
        
        load = Image.open(filename)
        w, h = load.size

        self.master.geometry(str(w)+ 'x' +str(h))

        render = ImageTk.PhotoImage(load)
        img = ttk.Label(self, image=render)
        img.image = render
        img.place(x=0, y=0)

        self.pack(fill=tk.BOTH, expand=1)



class StartPage(tk.Frame):

    def __init__(self, parent, controller):

        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text='Interactive simulation interface', font=LARGE_FONT)
        label.pack(pady=10, padx=10)

        btn_next = ttk.Button(self, text='Continue',
                        command=lambda: controller.show_frame(MainPage))
        btn_next.pack()

    def activate(self):
        pass

    def deactivate(self):
        pass


class MainPage(tk.Frame):

    def __init__(self, parent, controller):

        self.parent = parent
        self.controller = controller
        self.active = False
        self.anim_stopped = False

        tk.Frame.__init__(self, parent)

        btn1 = ttk.Button(self, text='Back home',
                        command=lambda: controller.show_frame(StartPage))
        btn1.grid(row=0, column=0)

        label = ttk.Label(self, text='Kit', font=LARGE_FONT)
        label.grid(row=0, column=1)

        btn2 = ttk.Button(self, text='Reset',
                        command=self.controller.eq.setInitialConditionKink)
        btn2.grid(row=0, column=2)

        ### Initialize plot
        #x = np.arange(10)
        #y = x**2

        # change line?
        #line, = a.plot(x, y)
        #lb = LineBuilder(line)
        #line, = fig.add_subplot(111).plot(x, y)


    def activate(self):

        self.active = True
        if self.anim_stopped:
            self.ani.event_source.start()
            return
        self.mustClear = False
        self.clicked = False
        self.released = False
        self.i_release = 0

        ## assuming eq and init cond are already defined
        eq = self.controller.eq


        ### create param window

        self.createParamWindow(eq)
        self.createImageWindow(eq)

        ### initialize plot and draw

        self.fig = Figure(figsize=(5,5), dpi=100)
        self.ax = self.fig.add_subplot(211)
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('u')
        self.line, = self.ax.plot(eq.x, eq.initCond)
        # make line interactive...
        #self.lb = LineBuilder(self.line)

        ## Spatiotemporal diagram
        spatiotemp = np.zeros((ST_ROWS, eq.N))
        self.ax2 = self.fig.add_subplot(212, sharex=self.ax)
        self.imvals = spatiotemp
        x0 = eq.x[0]
        xf = eq.x[-1]
        self.im = self.ax2.imshow(self.imvals, extent=[x0, xf, 0, 60], aspect='auto')
        #self.ax2.set_xlabel('x')
        self.ax2.set_ylabel('t')
        self.fig.colorbar(self.im, ax=self.ax2, orientation='horizontal')
        print(self.im)

        ## Draw canvas
        
        self.container_mpl = tk.Frame(self)

        self.canvas = FigureCanvasTkAgg(self.fig, self.container_mpl)
        self.canvas.draw()

        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('button_release_event', self.off_click)
        self.canvas.mpl_connect('motion_notify_event', self.move_click)


        self.toolbar = NavigationToolbar2Tk(self.canvas, self.container_mpl)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        #self.container_mpl.geometry('512x740')
        self.container_mpl.grid(row=1, column=0, columnspan=3, sticky='ns')


        #self.rowconfigure(0, weight=1)
        self.rowconfigure(1, weight=1)

        #self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        #self.columnconfigure(2, weight=1)

        #btn3 = ttk.Button(self, text='Clear', command=self.clear)

        #btn3.pack()

        ### Start animation

        self.ani = animation.FuncAnimation(self.fig, self.animate, frames=60, interval=100, blit=False)

    def on_click(self, event):
        if event.inaxes == self.ax:
            self.clicked = True
            self.updateY(event.xdata, event.ydata)
            self.animate(-1)
            #print(event.xdata, event.ydata)

    def off_click(self, event):
        self.clicked = False
        self.released = True
        self.animate(-1)

    def move_click(self, event):
        if event.inaxes == self.ax and self.clicked:
            self.updateY(event.xdata, event.ydata)
            self.animate(-1)

    def stop_animation(self):
        self.anim_stopped = True
        self.ani.event_source.stop()

    def deactivate(self):
        self.active = False

        # stop animation
        self.stop_animation

        self.killParamWindow()

    def createParamWindow(self, eq):

        # create new window
        self.parentParamWindow = tk.Toplevel(self)
        self.paramWindow = ParameterWindow(self.parentParamWindow, '400x500', eq)

    def createImageWindow(self, eq):
        self.parentImageWindow = tk.Toplevel(self)
        filename = eq.name + '.jpg'
        self.imageWindow = ImageWindow(self.parentImageWindow, filename)

    def killParamWindow(self):

        self.parentParamWindow.quit()
        self.parentParamWindow.destroy()

    def solve_cycle(self):
        eq = self.controller.eq
        eq.solve()
        eq.initCond = eq.sol[:, -1]

    def animate(self, i):
        eq = self.controller.eq
        eq.updateX()

        if i == -1:
            i = self.last_i
        else:
            self.last_i = i

        if self.released:
            eq.initCond = self.ys
            self.i_release = i
            self.solve_cycle()
            self.released = False
            return self.line,

        if self.clicked:
            self.line.set_data(eq.x, self.ys)
            #self.line.figure.canvas.draw()
            return self.line,

        k = (i-self.i_release) % 60
        print(i, k)
        if k == 0:
            self.solve_cycle()
        if self.mustClear:
            self.line, = self.ax.plot(eq.x, eq.sol[:,k])
            self.ax2.cla()
            self.im = self.ax2.imshow(np.zeros((ST_ROWS, eq.N)), extent=[x0, xf, 0, 60], aspect='auto')
            self.mustClear = False
        else:
            self.line.set_ydata(eq.sol[:, k])
            self.imvals[k, :] = eq.sol[:, k]
            self.im.set_data(self.imvals)

            xmin = np.min(eq.sol[:, k])
            xmax = np.max(eq.sol[:, k])

            vmin = np.min(self.imvals)
            vmax = np.max(self.imvals)

            self.im.set_clim(vmin, vmax)
            self.ax.set_ylim(vmin - abs(vmin)/10, vmax + abs(vmax)/10)


        self.ys = eq.sol[:, k]
        return self.line, self.im

    def raw_xtoi(self, x): # from raw_x to i
        return self.xtoi(x, isRaw=True)
        
    def xtoi(self, x, isRaw=False, truncate=False):
        x0 = self.controller.eq.x[0]
        xf = self.controller.eq.x[-1]
        _N = RAW_N if isRaw else self.controller.eq.N
        i = (x - x0)/(xf - x0) * _N
        if not truncate:
            i = round(i)
        return int(i)      

    def x0xf_to_xs(self, x0, xf):
        pass

    def xstois(self, xs):
        i_s = np.zeros(len(xs), dtype='int64')
        for i in range(len(xs)):
            i_s[i] = xtoi(xs[i])

    def updateY(self, xdata, ydata):
        xi = self.xtoi(xdata)
        if xi < len(self.ys) - 2 and xi >= 2:

            ### interpolate
            raw_dx = (self.controller.eq.x[-1] - self.controller.eq.x[0]) / (RAW_N - 1)
            x0 = xdata - raw_dx
            xf = xdata + raw_dx

            i_0 = self.xtoi(x0)
            i_f = self.xtoi(xf)
            i_s = np.array([i_0, xi, i_f])

            xs_interval = self.controller.eq.x[i_s]
            ys_interval = np.array([self.ys[i_0], ydata, self.ys[i_f]])
            ys_interp = interp1d(xs_interval, ys_interval)

            i_0 = self.xtoi(xdata - raw_dx/2)
            i_f = self.xtoi(xdata + raw_dx/2)

            # replace old values of ys
            self.ys[i_0:i_f+1] = ys_interp(self.controller.eq.x[i_0:i_f+1])
        elif xi < 2:
            self.ys[0:2] = ydata
        else:
            self.ys[-2:] = ydata


    def clear(self):
        self.ax.clear()
        self.mustClear = True
        #self.line, = self.ax.plot([], [])

class SpatioTemporalDiagram(tk.Frame):

    def __init__(self, root, eq):

        self.root = root
        self.root.title('Spatio-temporal Diagram')

        self.eq = eq
        self.parameters = self.eq.params

        # 60 rows
        self.ROWS = 60
        self.vals = np.zeros((self.ROWS, self.eq.N))

    def animate(self, i):
        pass


class ParameterWindow:

    def __init__(self, root, geometry, eq, update=None):

        self.root = root
        self.root.title('Parameters')
        #self.root.geometry(geometry) # string: sizexsize
        self.eq = eq
        self.parameters = self.eq.params
        self.update = update

        self.titlelbl = ttk.Label(self.root, text='Parameters', width=15, font=LARGE_FONT)
        self.paramslbls = [ttk.Label(self.root, text=self.parameters[p].name, width=6, font=MED_FONT) for p in self.parameters]
        self.paramsvals = [ttk.Entry(self.root, textvariable=self.parameters[p].val, width=4, font=MED_FONT) for p in self.parameters]

        self.titlelbl.grid(column=0, row=0, columnspan=3)
        for i in range(len(self.paramslbls)):
            self.paramslbls[i].grid(column=0, row=i+1, columnspan=2)
            self.paramsvals[i].grid(column=3, row=i+1)

        self.resetBtn = ttk.Button(self.root, text='Reset', command=self.reset)
        self.resetBtn.grid(column=0, row=len(self.paramslbls)+1)

    def reset(self):
        for p in self.eq.initParams:
            self.parameters[p].val.set(self.eq.initParams[p])
