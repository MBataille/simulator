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
from customtoolbar import CustomToolbar

from scipy.interpolate import interp1d

from solvers import INTEGRATION_METHODS

from equations.allequations import ALL_EQUATIONS


LARGE_FONT = ('Arial', 16)
MED_FONT = ('Arial', 12)
COLORS = 'none blue orange green red purple brown pink gray olive cyan'.split(' ')


PROFILE = 'profile'
SPATIOTEMPORAL = 'spatiotemporal'

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

        self.parent = parent
        self.controller = controller

        self.listbox = tk.Listbox(self)
        for eqname in ALL_EQUATIONS:
            self.listbox.insert(tk.END, eqname)

        self.listbox.bind('<<ListboxSelect>>', self.change_select)
        self.listbox.pack()

        btn_next = ttk.Button(self, text='Continue',
                        command=lambda: controller.show_frame(MainPage))
        btn_next.pack()

    def activate(self):
        pass

    def deactivate(self):
        pass

    def change_select(self, event):
        selection = event.widget.curselection()
        if selection:
            index = selection[0]
            data = event.widget.get(index)
            self.controller.setEq(ALL_EQUATIONS[data])
            self.controller.eq.setInitialConditionZero()

class InspectorProfile(tk.Frame):

    def __init__(self, parent):

        tk.Frame.__init__(self, parent)
        
        self.parent = parent

        self.titlelbl = ttk.Label(self, text='Inspector: Profile', font=LARGE_FONT)
        self.titlelbl.grid(column=0, row=0, columnspan=3, padx=20, pady=10)

        self.y_minlbl = ttk.Label(self, text='y_min')
        self.y_maxlbl = ttk.Label(self, text='y_max')

        self.y_minlbl.grid(column=0, row=2, columnspan=2)
        self.y_maxlbl.grid(column=0, row=1, columnspan=2)

        y_min, y_max = -1, 1
        
        self.y_minvar = tk.StringVar(value=str(y_min))

        self.y_maxvar = tk.StringVar(value=str(y_max))

        self.y_minval = ttk.Entry(self, textvariable=self.y_minvar)
        self.y_maxval = ttk.Entry(self, textvariable=self.y_maxvar)

        self.y_minval.grid(column=3, row=2)
        self.y_maxval.grid(column=3, row=1)

        self.auto_ylim = True
        self.autotxt = tk.StringVar()
        self.autotxt.set('Auto')
        self.autobtn = ttk.Button(self, textvariable=self.autotxt, command=self.set_auto)
        self.autobtn.grid(column=0, row=3, columnspan=3)


        self.intmethodlbl = ttk.Label(self, text='Integration method: ', font=MED_FONT)
        self.intmethodlbl.grid(column=0, row=4, columnspan=2, padx=10, pady=10)

        self.combobox = ttk.Combobox(self, state='readonly')
        self.combobox['values'] = [int_method.code for int_method in INTEGRATION_METHODS]
        self.combobox.bind('<<ComboboxSelected>>', self.change_solver)
        self.combobox.current(0)
        self.description_text = tk.StringVar(value=INTEGRATION_METHODS[0].description)

        self.combobox.grid(column=3, row=4, padx=10, pady=10)

        self.descriptionlbl = ttk.Label(self, textvariable=self.description_text)
        self.descriptionlbl.grid(row=5, column=0, columnspan=3, padx=5, pady=5, sticky='e')

        self.choosefieldbl = ttk.Label(self, text='Show field: ', font=MED_FONT)
        self.choosefieldbl.grid(column=0, row=6, padx=10, pady=10)

        self.choosefieldcbox = ttk.Combobox(self, state='readonly')
        self.choosefieldcbox['values'] = self.parent.controller.eq.fieldNames
        self.choosefieldcbox.bind('<<ComboboxSelected>>', self.change_field)
        self.choosefieldcbox.current(0)
        self.choosefieldcbox.grid(column=1, row=6)
        self.active_field_indx = 0

        self.choosecolorlbl = ttk.Label(self, text='Color: ', font=MED_FONT)
        self.choosecolorlbl.grid(column=2, row=6)

        self.current_colors = [_i + 1 for _i in range(self.parent.controller.eq.n_fields)]

        self.choosecolorbox = ttk.Combobox(self, state='readonly')
        self.choosecolorbox['values'] = COLORS
        self.choosecolorbox.bind('<<ComboboxSelected>>', self.change_color)
        self.choosecolorbox.current(self.current_colors[0])
        self.choosecolorbox.grid(column=3, row=6)

    def change_field(self, event):
        mainpg = self.parent.mainpg
        mainpg.activeFieldToCurrent()
        self.active_field_indx = self.choosefieldcbox.current()
        mainpg.resetActiveField()
        self.choosecolorbox.current(self.current_colors[self.active_field_indx])

    def change_color(self, event):
        new_color = self.choosecolorbox.current()
        mainpg = self.parent.mainpg
        if self.current_colors[self.active_field_indx] != new_color:
            self.current_colors[self.active_field_indx] = new_color
            self.parent.mainpg.redraw_fields()

    def change_solver(self, event):
        solver_index = self.combobox.current()
        int_method = INTEGRATION_METHODS[solver_index]
        print(f'Selected {int_method.name}')
        # change description
        self.description_text.set(int_method.description)

        ## change solver
        self.parent.controller.eq.setSolver(int_method)

    def set_auto(self):
        if self.auto_ylim:
            self.auto_ylim = False
            self.autotxt.set('Manual')
        else:
            self.auto_ylim = True
            self.autotxt.set('Auto')

    def set_ylim(self, ymin, ymax):
        self.y_minvar.set(str(ymin))
        self.y_maxvar.set(str(ymax))

    def get_ylim(self):
        try:
            y_min = float(self.y_minvar.get())
            y_max = float(self.y_maxvar.get())
            return y_min, y_max
        except ValueError as e:
            return None, None

    def activate(self):
        pass

    def deactivate(self):
        pass

    def setAx(self, ax):
        self.ax = ax
        ymin, ymax = self.ax.get_ylim()
        self.y_minvar.set(ymin)
        self.y_maxvar.set(ymax)

class InspectorSpatiotemporal(tk.Frame):

    def __init__(self, parent):

        tk.Frame.__init__(self, parent)
        #self.ax2 = ax2
        self.parent = parent

        self.titlelbl = ttk.Label(self, text='Inspector: Spatiotemporal diagram', font=LARGE_FONT)
        self.titlelbl.grid(column=0, row=0, columnspan=3, padx=20, pady=10)

    def activate(self):
        pass

    def deactivate(self):
        pass

class InspectorWindow(tk.Frame):

    def __init__(self, parent, controller, mainpg):
        tk.Frame.__init__(self, parent)

        self.parent = parent
        self.parent.title('Inspector')
        self.mainpg = mainpg
        self.controller = controller

        self.profile = InspectorProfile(self)
        self.spatiotemp = InspectorSpatiotemporal(self)

        self.profile.grid(row=0, column=0)
        self.spatiotemp.grid(row=0, column=0)

        self.showProfile()
        self.grid(row=0, column=0)

    def showProfile(self):
        print('showing profile')
        self.showFrame(self.profile)   

    def showSpatiotemp(self):
        print('showing spatiotemp')
        self.showFrame(self.spatiotemp)

    def showFrame(self, frame):
        frame.tkraise()
        # frame.activate()



class ParameterWindow:

    def __init__(self, root, geometry, eq, update=None):

        self.root = root
        self.root.title('Parameters')
        # self.root.geometry(geometry) # string: sizexsize
        self.eq = eq
        self.parameters = self.eq.parameters
        self.update = update

        self.titlelbl = ttk.Label(self.root, text='Parameters', width=15, font=LARGE_FONT)
        self.paramslbls = [ttk.Label(self.root, text=self.parameters[p].name, font=MED_FONT) for p in self.parameters]
        self.paramsvals = [ttk.Entry(self.root, textvariable=self.parameters[p].var, font=MED_FONT) for p in self.parameters]

        self.titlelbl.grid(column=0, row=0, columnspan=3)
        for i in range(len(self.paramslbls)):
            self.paramslbls[i].grid(column=0, row=i+1, columnspan=1)
            self.paramsvals[i].grid(column=1, row=i+1, columnspan=2)

        self.resetBtn = ttk.Button(self.root, text='Reset', command=self.reset)
        self.resetBtn.grid(column=0, row=len(self.paramslbls)+1)

    def reset(self):
        for p in self.eq.initParams:
            self.parameters[p].var.set(str(self.eq.initParams[p]))



class MainPage(tk.Frame):

    def __init__(self, parent, controller):

        self.parent = parent
        self.controller = controller
        self.active = False
        self.anim_stopped = False

        tk.Frame.__init__(self, parent)

        btn1 = ttk.Button(self, text='Back home',
                        command=self.back_home)
        btn1.grid(row=0, column=0)

        label = ttk.Label(self, text='Interactive simulation', font=LARGE_FONT)
        label.grid(row=0, column=1)

        btn2 = ttk.Button(self, text='Reset',
                        command=self.controller.eq.setInitialConditionZero)
        btn2.grid(row=0, column=2)

        ### Initialize plot
        #x = np.arange(10)
        #y = x**2

        # change line?
        #line, = a.plot(x, y)
        #lb = LineBuilder(line)
        #line, = fig.add_subplot(111).plot(x, y)

    def deactivate(self):

        self.active = False
        self.killParamWindow()
        self.killInspectorWindow()
        self.killImageWindow()
        self.stop_animation()

    def back_home(self):
        self.deactivate()
        self.controller.show_frame(StartPage)

    def draw_fields(self):
        self.lines = []
        print(len(self.Fields), len(self.inspector.profile.current_colors))
        for i in range(len(self.Fields)):
            curfield = self.Fields[i]
            color = self.inspector.profile.current_colors[i]
            if COLORS[color] == 'none':
                line = None
            else:
                line, = self.ax.plot(self.eqX, curfield, 'tab:' + COLORS[color])
            self.lines.append(line)

    def activate(self):

        self.active = True

        self.mustClear = False
        self.clicked = False
        self.released = False
        self.i_release = 0

        ## assuming eq and init cond are already defined
        eq = self.controller.eq
        eq.updateX()


        ### create param window

        self.createParamWindow(eq)
        self.createImageWindow(eq)
        self.createInspectorWindow(eq)

        ### initialize plot and draw

        self.fig = Figure(figsize=(5,5), dpi=100)
        self.ax = self.fig.add_subplot(211)
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('u')

        self.eqX = eq.x
        self.Fields = eq.getInitCondFields()

        #self.line, = self.ax.plot(eq.x, eq.initCond)
        self.draw_fields()

        ## Spatiotemporal diagram
        spatiotemp = np.zeros((ST_ROWS, eq.Ni))
        self.ax2 = self.fig.add_subplot(212, sharex=self.ax)
        self.imvals = spatiotemp
        x0 = eq.x[0]
        xf = eq.x[-1]
        self.im = self.ax2.imshow(self.imvals, extent=[x0, xf, 0, ST_ROWS], aspect='auto')
        #self.ax2.set_xlabel('x')
        self.ax2.set_ylabel('t')
        self.fig.colorbar(self.im, ax=self.ax2, orientation='horizontal')

        ## Draw canvas
        
        self.container_mpl = tk.Frame(self)

        self.canvas = FigureCanvasTkAgg(self.fig, self.container_mpl)
        self.canvas.draw()

        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('button_release_event', self.off_click)
        self.canvas.mpl_connect('motion_notify_event', self.move_click)


        self.toolbar = CustomToolbar(self.canvas, self.container_mpl, self)
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

        self.isEditing = False
        self.isPaused = False
        self.k_spatiotemp = 0
        self.k_sol = 0

        self.inspector.profile.setAx(self.ax)
        ### Start animation

        self.ani = animation.FuncAnimation(self.fig, self.animate, frames=60, interval=100, blit=False)

    def edit(self):
        print(f'Edit clicked!, isEditing = {not self.isEditing}')
        if self.isEditing:
            self.isEditing = False
        else:
            self.isEditing = True
            # self.animate(-1)

    def pause(self):
        if self.isPaused: return
        self.isPaused = True

    def play(self):
        if not self.isPaused: return
        self.isPaused = False
        #self.release = True
        # if self.isEditing:
        #     self.isEditing = False
        #     self.released = True
        self.animate(-1)

    def open_params(self):
        try:
            if self.parentParamWindow.state() == 'normal':
                self.parentParamWindow.focus_set()
        except Exception as e:
            self.createParamWindow(self.controller.eq)

    def open_equation(self):
        try:
            if self.parentImageWindow.state() == 'normal':
                self.parentImageWindow.focus_set()
        except Exception as e:
            self.createImageWindow(self.controller.eq)

    def open_inspector(self):
        try:
            if self.parentInspectorWindow.state() == 'normal':
                self.parentInspectorWindow.focus_set()
        except Exception as e:
            self.createInspectorWindow(self.controller.eq)

    def on_click(self, event):
        if event.inaxes == self.ax:
            self.controller.activePlot = (self.ax, PROFILE)
            self.inspector.showProfile()
            if self.isEditing:
                self.clicked = True
                self.updateY(event.xdata, event.ydata)
            #self.animate(-1)
            #print(event.xdata, event.ydata)
        elif event.inaxes == self.ax2:
            self.inspector.showSpatiotemp()
            self.controller.activePlot = (self.ax2, SPATIOTEMPORAL)


    def off_click(self, event):
        if self.clicked:
            self.clicked = False
            self.released = True
            self.isEditing = False
            self.animate(-1)
            # if not self.isPaused:
            #     self.play()

    def move_click(self, event):
        if event.inaxes == self.ax and self.clicked:
            self.updateY(event.xdata, event.ydata)
            #self.animate(-1)

    def stop_animation(self):
        self.anim_stopped = True
        self.ani.event_source.stop()

    def createParamWindow(self, eq):

        # create new window
        self.parentParamWindow = tk.Toplevel(self)
        self.paramWindow = ParameterWindow(self.parentParamWindow, '400x500', eq)

    def createImageWindow(self, eq):
        self.parentImageWindow = tk.Toplevel(self)
        filename = eq.name + '.jpg'
        self.imageWindow = ImageWindow(self.parentImageWindow, filename)

    def createInspectorWindow(self, eq):

        self.parentInspectorWindow = tk.Toplevel(self)
        self.inspector = InspectorWindow(self.parentInspectorWindow, self.controller, self)

    def killParamWindow(self):
        self.parentParamWindow.destroy()

    def killInspectorWindow(self):
        self.parentInspectorWindow.destroy()

    def killImageWindow(self):
        self.parentImageWindow.destroy()

    def solve_cycle(self):
        eq = self.controller.eq
        eq.solve()
        eq.initCond = eq.sol[:, -1]
        self.k_sol = 0

    def update_fields(self):
        for i in range(len(self.Fields)):
            if self.lines[i] is not None:
                self.lines[i].set_ydata(self.Fields[i])
                self.lines[i].set_xdata(self.controller.eq.x)

    def redraw_fields(self):
        self.ax.cla()
        self.draw_fields()

    def resetActiveField(self):
        self.active_Field = self.Fields[self.inspector.profile.active_field_indx]

    def activeFieldToCurrent(self):
        self.Fields[self.inspector.profile.active_field_indx] = self.active_Field

    def animate(self, i):
        if not self.active: return
        eq = self.controller.eq
        eq.updateX()

        if i == -1:
            i = self.last_i # no estuy usando
    
            if self.released:
                print('Released!')
                self.activeFieldToCurrent()
                eq.setInitCondFields(self.Fields)
                self.solve_cycle()
                self.released = False
                return self.lines

        elif not self.isPaused:
            self.last_i = i # tpco toy usando

        # if self.released:
        #     self.i_release = i
        #     self.released = False

        # if self.isEditing:
            # self.Fields[self.inspector.profile.active_field_indx] = self.active_Field
        #     self.update_fields()
        #     #self.line.figure.canvas.draw()
            # return self.lines

        if self.isPaused or self.isEditing:
            return self.lines + [self.im]

        k = (i-self.i_release) % ST_ROWS
        # k = i % ST_ROWS
        # print(f'i = {i}, i_r = {self.i_release}, i_s = {self.i_start_pause}, k = {k}, k_spatiotemp = {k_spatiotemp}.')
        # print(i, k)
        if self.k_sol == 0:
            self.solve_cycle()
        self.Fields = eq.getFields(self.k_sol)
        self.resetActiveField()

        if self.mustClear:
            self.redraw_fields()
            self.ax2.cla()
            self.im = self.ax2.imshow(np.zeros((ST_ROWS, eq.N)), extent=[x0, xf, 0, ST_ROWS], aspect='auto')
            self.mustClear = False
        else:
            self.update_fields() # draw fields
            self.imvals[self.k_spatiotemp, :] = self.Fields[0] # for the moment
            self.im.set_data(self.imvals)

            self.k_spatiotemp = (self.k_spatiotemp + 1)% 60

            xmin = np.min(self.Fields)
            xmax = np.max(self.Fields)

            vmin = np.min(self.imvals)
            vmax = np.max(self.imvals)

            self.im.set_clim(vmin, vmax)
            self.ax.set_xlim(eq.x[0], eq.x[-1])
            if self.inspector.profile.auto_ylim:
                ymin = vmin - abs(vmin)/10
                ymax = vmax + abs(vmax)/10

                self.ax.set_ylim(min(ymin, xmin - abs(xmin)/10), max(ymax, xmax + abs(xmax)/10))
                self.inspector.profile.set_ylim(min(ymin, xmin - abs(xmin)/10), max(ymax, xmax + abs(xmax)/10))
            else:
                ymin, ymax = self.inspector.profile.get_ylim()
                if ymin is not None:
                    self.ax.set_ylim(ymin, ymax)

        self.k_sol = (self.k_sol + 1) % ST_ROWS
        return self.lines + [self.im]

    def raw_xtoi(self, x): # from raw_x to i
        return self.xtoi(x, isRaw=True)
        
    def xtoi(self, x, isRaw=False, truncate=False):
        x0 = self.controller.eq.x[0]
        xf = self.controller.eq.x[-1]
        Ni = round(self.controller.eq.N / self.controller.eq.n_fields)
        _N = RAW_N if isRaw else Ni
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
        #self.active_Field = self.Fields[self.inspector.profile.active_field_indx]
        self.resetActiveField()
        if xi < len(self.active_Field) - 4 and xi >= 4:

            ### interpolate
            raw_dx = (self.controller.eq.x[-1] - self.controller.eq.x[0]) / (RAW_N - 1)
            x0 = xdata - raw_dx/2
            xf = xdata + raw_dx/2

            i_0 = self.xtoi(x0)
            i_f = self.xtoi(xf)
            i_s = np.array([i_0, xi, i_f])

            xs_interval = self.controller.eq.x[i_s]
            ys_interval = np.array([self.active_Field[i_0], ydata, self.active_Field[i_f]])
            ys_interp = interp1d(xs_interval, ys_interval)

            i_0 = self.xtoi(xdata - raw_dx/2)
            i_f = self.xtoi(xdata + raw_dx/2)

            # replace old values of ys
            self.active_Field[i_0:i_f+1] = ys_interp(self.controller.eq.x[i_0:i_f+1])
        elif xi < 4:
            self.active_Field[0:4] = ydata
        else:
            self.active_Field[-4:] = ydata
        self.activeFieldToCurrent()
        self.update_fields()

    def clear(self):
        self.ax.clear()
        self.mustClear = True
        #self.line, = self.ax.plot([], [])

