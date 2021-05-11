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
#from customtoolbar import CustomToolbar

from scipy.interpolate import interp1d

from solvers import INTEGRATION_METHODS

from equations.allequations import ALL_EQUATIONS

LARGE_FONT = ('Helvetica', 16)
MED_FONT = ('Helvetica', 11)
COLORS = 'black none blue orange green red purple brown pink gray olive cyan'.split(' ')

COLORMAPS = 'viridis plasma hot cool coolwarm hsv jet'.split(' ')

PROFILE = 'profile'
SPATIOTEMPORAL = 'spatiotemporal'
IMAGESFOLDER = 'images/'

ICONSFOLDER = 'icons/'

style.use('ggplot')

ST_ROWS = 60
RAW_N = 30

def i10p(val):
    """ Increase val in 10%  """
    return val + 0.1 * abs(val)

def d10p(val):
    """ Decrease val in 10% """
    return val - 0.1 * abs(val)
class ImageWindow(tk.Frame):
    def __init__(self, master, filename):
        tk.Frame.__init__(self, master)
        self.master = master
        self.master.title('Equation')

        load = Image.open(IMAGESFOLDER + filename)
        w, h = load.size

        self.master.geometry(str(w) + 'x' + str(h))

        render = ImageTk.PhotoImage(load)
        img = ttk.Label(self, image=render)
        img.image = render
        img.place(x=0, y=0)

        self.pack(fill=tk.BOTH, expand=1)


class StartPage(tk.Frame):

    def __init__(self, parent, controller):

        tk.Frame.__init__(self, parent)
        label = ttk.Label(self, text='Interactive simulation interface', font=LARGE_FONT)
        label.grid(row=0, column=0, columnspan=2,  padx=10, pady=5)

        self.parent = parent
        self.controller = controller

        self.eq_listbox = tk.Listbox(self)
        for eqname in ALL_EQUATIONS:
            self.eq_listbox.insert(tk.END, eqname)
        self.eq_listbox.selection_set(0)

        self.eq_lbl = ttk.Label(self, text='Choose equation', font=MED_FONT)
        self.ic_lbl = ttk.Label(self, text='Choose initial condition', font=MED_FONT)

        self.eq_lbl.grid(row=1, column=0, padx=10, pady=5)
        self.ic_lbl.grid(row=1, column=1, sticky='ew', padx=5, pady=5)

        self.eq_listbox.bind('<<ListboxSelect>>', self.change_select_eq)
        self.eq_listbox.grid(row=2, column=0, padx=5, pady=5, sticky='e')

        self.initcond_listbox = tk.Listbox(self)
        self.initcond_listbox.bind('<<ListboxSelect>>', self.change_select_ic)
        self.fill_initcond_listbox()       

        self.initcond_listbox.grid(row=2, column=1, padx=5, pady=5, sticky='w')

        self.n_lbl = ttk.Label(self, text='Number of pts =  ', font=MED_FONT)
        self.n_lbl.grid(row=3, column=0, padx=5, pady=5, sticky='e')

        self.n_var = tk.StringVar(value='200')
        self.n_ent = tk.Entry(self, textvariable=self.n_var)
        self.n_ent.grid(row=3, column=1, padx=5, pady=5, sticky='w')

        self.statustext = tk.StringVar()
        self.statuslbl = ttk.Label(self, textvariable=self.statustext)
        self.statuslbl.grid(row=4, column=0, columnspan=2)

        btn_next = ttk.Button(self, text='Continue',
                              command=self.cont)
        btn_next.grid(row=5, column=0, columnspan=2)

        #self.rowconfigure(0, weight=1)
        #self.rowconfigure(1, weight=1)
        #   self.rowconfigure(2, weight=1)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)

    def activate(self):
        pass

    def deactivate(self):
        pass

    def clear_initcond_listbox(self):
        last = self.initcond_listbox.size() - 1
        self.initcond_listbox.delete(0, last=last)

    def fill_initcond_listbox(self):
        self.clear_initcond_listbox()
        initconds = self.controller.getEqInitConds()
        for initcond in initconds:
            self.initcond_listbox.insert(tk.END, initcond)
        self.initcond_listbox.selection_set(0)

    def change_select_eq(self, event):
        selection = event.widget.curselection()
        if selection:
            index = selection[0]
            name = event.widget.get(index)
            self.controller.setEq(name)
            self.fill_initcond_listbox()

    def change_select_ic(self, event):
        selection = event.widget.curselection()
        if selection:
            index = selection[0]
            if not self.controller.initCond_isMethod(index):
                self.controller.setEqInitCond(index)
                Ni = self.controller.getEqInitCond_N()
                self.n_var.set(str(Ni))
                self.n_ent.configure(state='readonly')
            else:
                self.n_ent.configure(state=tk.NORMAL)

    def cont(self):
        try:
            index = self.initcond_listbox.curselection()[0]
        except IndexError:
            self.statustext.set('Please enter an initial condition to continue.')
        print(f'selected {self.initcond_listbox.get(index)} initcond')

        try:
            Ni = int(self.n_var.get())
        except ValueError:
            self.statustext.set('Please enter a valid value for N')
            return

        if self.controller.initCond_isMethod(index):
            self.controller.setEqNi(Ni)
            self.controller.setEqInitCond(index)

        else: # interpolate
            pass
        self.controller.show_frame(MainPage)


class InspectorProfile(tk.Frame):

    def __init__(self, parent):

        tk.Frame.__init__(self, parent)

        self.parent = parent

        self.titlelbl = ttk.Label(self, text='Inspector: Profile', font=LARGE_FONT)
        self.titlelbl.grid(column=0, row=0, columnspan=3, padx=20, pady=10)

        self.y_minlbl = ttk.Label(self, text='y_min', font=MED_FONT)
        self.y_maxlbl = ttk.Label(self, text='y_max', font=MED_FONT)

        self.y_minlbl.grid(column=0, row=2)
        self.y_maxlbl.grid(column=0, row=1)

        y_min, y_max = -1, 1

        self.y_minvar = tk.StringVar(value=str(y_min))

        self.y_maxvar = tk.StringVar(value=str(y_max))

        self.y_minval = ttk.Entry(self, textvariable=self.y_minvar, font=MED_FONT, width=15)
        self.y_maxval = ttk.Entry(self, textvariable=self.y_maxvar, font=MED_FONT, width=15)

        self.y_minval.grid(column=1, row=2, padx=10)
        self.y_maxval.grid(column=1, row=1, padx=10)

        self.x_minlbl = ttk.Label(self, text='x_min', font=MED_FONT)
        self.x_maxlbl = ttk.Label(self, text='x_max', font=MED_FONT)

        self.x_minlbl.grid(column=2, row=2)
        self.x_maxlbl.grid(column=2, row=1)

        x_min, x_max = -1, 1

        self.x_minvar = tk.StringVar(value=str(x_min))

        self.x_maxvar = tk.StringVar(value=str(x_max))

        self.x_minval = ttk.Entry(self, textvariable=self.x_minvar, font=MED_FONT, width=15)
        self.x_maxval = ttk.Entry(self, textvariable=self.x_maxvar, font=MED_FONT, width=15)

        self.x_minval.grid(column=3, row=2, padx=10)
        self.x_maxval.grid(column=3, row=1, padx=10)

        self.autotxt = tk.StringVar()
        self.autotxt.set('Auto')
        self.autobtn = ttk.Button(self, textvariable=self.autotxt, command=self.set_auto)
        self.autobtn.grid(column=3, row=3, padx=10, pady=10)

        self.intmethodlbl = ttk.Label(self, text='Integration method: ', font=MED_FONT)
        self.intmethodlbl.grid(column=0, row=4, padx=10, pady=10, sticky='ew', columnspan=2)

        self.combobox = ttk.Combobox(self, state='readonly', width=10)
        self.combobox['values'] = [int_method.code for int_method in INTEGRATION_METHODS]
        self.combobox.bind('<<ComboboxSelected>>', self.change_solver)
        self.combobox.current(0)
        self.description_text = tk.StringVar(value=INTEGRATION_METHODS[0].description)

        self.combobox.grid(column=3, row=4, padx=10, pady=10)

        self.descriptionlbl = ttk.Label(self, textvariable=self.description_text, font=MED_FONT)
        self.descriptionlbl.grid(row=5, column=1, columnspan=3, padx=5, pady=5, sticky='ew')

        self.bc_lbl = ttk.Label(self, text='Boundary Condition: ', font=MED_FONT)
        self.bc_lbl.grid(row=6, column=0, padx=10, pady=10, sticky='ew', columnspan=2)

        self.bc_combobox = ttk.Combobox(self, state='readonly', width=10)
        self.bc_combobox['values'] = self.parent.mainpg.controller.getListBoundaryConditions()
        self.bc_combobox.bind('<<ComboboxSelected>>', self.change_bc)
        self.bc_combobox.current(0) # neumann

        self.bc_combobox.grid(row=6, column=3)

        self.ellapsed_time_desc_lbl = ttk.Label(self, text='Ellapsed time: ', font=MED_FONT)
        self.ellapsed_time_desc_lbl.grid(row=3, column=0, padx=10, pady=5)

        self.time_container = tk.Frame(self)

        self.resetimg = tk.PhotoImage(file=ICONSFOLDER + 'reset20.png')

        self.ellapsed_time_var = tk.StringVar(value='0')
        self.ellapsed_time_val_lbl = ttk.Label(self.time_container, textvariable=self.ellapsed_time_var)
        self.ellapsed_time_val_lbl.grid(row=0, column=1, padx=5, pady=5)

        self.elapsed_time_reset_btn = ttk.Button(self.time_container, text='Reset', command=self.reset_time, image=self.resetimg)
        self.elapsed_time_reset_btn.grid(row=0, column=2, padx=5, pady=5)

        self.time_container.grid(row=3, column=1)

        self.choosefieldbl = ttk.Label(self, text='Show field: ', font=MED_FONT)
        self.choosefieldbl.grid(column=0, row=8, padx=10, pady=10)

        self.choosefieldcbox = ttk.Combobox(self, state='readonly', width=10)
        self.choosefieldcbox['values'] = self.parent.controller.eq.fieldNames + self.parent.controller.eq.auxFieldNames
        self.choosefieldcbox.bind('<<ComboboxSelected>>', self.change_field)
        self.choosefieldcbox.current(0)
        self.choosefieldcbox.grid(column=1, row=8)

        self.choosecolorlbl = ttk.Label(self, text='Color: ', font=MED_FONT)
        self.choosecolorlbl.grid(column=2, row=8)

        self.choosecolorbox = ttk.Combobox(self, state='readonly', width=10)
        self.choosecolorbox['values'] = COLORS
        self.choosecolorbox.bind('<<ComboboxSelected>>', self.change_color)
        self.choosecolorbox.current(2)
        self.choosecolorbox.grid(column=3, row=8)

    def initProfilePlot(self, pplot):
        self.pplot = pplot
        self.setAx(pplot.ax)

    def changeProfilePlot(self, pplot):
        self.pplot = pplot

        self.choosefieldcbox.current(self.pplot.active_field_indx)
        self.choosecolorbox.current(self.pplot.getCurrentColorActiveField())
        self.autotxt.set("Auto" if self.pplot.auto_lim else "Manual")

    def change_bc(self, event):
        indx = self.bc_combobox.current() # store controller directly?
        self.parent.mainpg.controller.setBoundaryCondition(indx)
        self.parent.mainpg.controller.solve_cycle()

    def reset_time(self): # interact w P plot
        self.parent.mainpg.t = 0

    def change_field(self, event):
        self.pplot.activeFieldToCurrent()
        self.pplot.active_field_indx = self.choosefieldcbox.current()
        self.pplot.resetActiveField()
        self.choosecolorbox.current(self.pplot.getCurrentColorActiveField())

    def change_color(self, event): # interact w P Plot
        new_color = self.choosecolorbox.current()
        curcolor = self.pplot.getCurrentColorActiveField()
        if curcolor != new_color:
            self.pplot.setCurrentColorActiveField(new_color)

    def change_solver(self, event):
        solver_index = self.combobox.current()
        int_method = INTEGRATION_METHODS[solver_index]
        print(f'Selected {int_method.name}')
        # change description
        self.description_text.set(int_method.description)

        # change solver
        self.parent.controller.eq.setSolver(int_method) # call controller


    def set_auto(self):
        if self.pplot.auto_lim:
            self.pplot.auto_lim = False
            self.autotxt.set('Manual')
        else:
            self.pplot.auto_lim = True
            self.autotxt.set('Auto')

    def set_ylim(self, ymin, ymax):
        self.y_minvar.set('{:.5f}'.format(ymin))
        self.y_maxvar.set('{:.5f}'.format(ymax))

    def set_xlim(self, xmin, xmax):
        self.x_minvar.set('{:.5f}'.format(xmin))
        self.x_maxvar.set('{:.5f}'.format(xmax))

    def get_ylims(self):
        try:
            y_min = float(self.y_minvar.get())
            y_max = float(self.y_maxvar.get())
            return y_min, y_max
        except ValueError:
            return None, None

    def get_xlims(self):
        try:
            x_min = float(self.x_minvar.get())
            x_max = float(self.x_maxvar.get())
            return x_min, x_max
        except ValueError:
            return None, None

    def activate(self):
        self.grid(row=0, column=0)

    def deactivate(self):
        self.grid_forget()

    def setAx(self, ax): # P Plot
        self.ax = ax
        xlims, ylims = self.pplot.get_ax_lims()
        self.set_ylim(*ylims)
        self.set_xlim(*xlims)

    def setTime(self, t):
        self.ellapsed_time_var.set('{:.2f}'.format(t))


class InspectorSpatiotemporal(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)
        self.parent = parent

        self.titlelbl = ttk.Label(self, text='Inspector: Spatiotemporal diagram', font=LARGE_FONT)
        self.titlelbl.grid(column=0, row=0, columnspan=3, padx=20, pady=10)

        self.v_minlbl = ttk.Label(self, text='v_min', font=MED_FONT)
        self.v_maxlbl = ttk.Label(self, text='v_max', font=MED_FONT)

        self.v_minlbl.grid(column=4, row=2)
        self.v_maxlbl.grid(column=4, row=1)

        v_min, v_max = -1, 1

        self.v_minvar = tk.StringVar(value=str(v_min))

        self.v_maxvar = tk.StringVar(value=str(v_max))

        self.v_minval = ttk.Entry(self, textvariable=self.v_minvar, width=15, font=MED_FONT)
        self.v_maxval = ttk.Entry(self, textvariable=self.v_maxvar, width=15, font=MED_FONT)

        self.v_minval.grid(column=5, row=2)
        self.v_maxval.grid(column=5, row=1)

        self.x_minlbl = ttk.Label(self, text='x_min', font=MED_FONT)
        self.x_maxlbl = ttk.Label(self, text='x_max', font=MED_FONT)

        self.x_minlbl.grid(column=0, row=2)
        self.x_maxlbl.grid(column=0, row=1)

        x_min, x_max = -1, 1

        self.x_minvar = tk.StringVar(value=str(x_min))

        self.x_maxvar = tk.StringVar(value=str(x_max))

        self.x_minval = ttk.Entry(self, textvariable=self.x_minvar, width=15, font=MED_FONT)
        self.x_maxval = ttk.Entry(self, textvariable=self.x_maxvar, width=15, font=MED_FONT)

        self.x_minval.grid(column=1, row=2)
        self.x_maxval.grid(column=1, row=1)

        ####  ####

        self.y_minlbl = ttk.Label(self, text='y_min', font=MED_FONT)
        self.y_maxlbl = ttk.Label(self, text='y_max', font=MED_FONT)

        self.y_minlbl.grid(column=2, row=2)
        self.y_maxlbl.grid(column=2, row=1)

        y_min, y_max = -1, 1

        self.y_minvar = tk.StringVar(value=str(y_min))

        self.y_maxvar = tk.StringVar(value=str(y_max))

        self.y_minval = ttk.Entry(self, textvariable=self.y_minvar, width=15, font=MED_FONT)
        self.y_maxval = ttk.Entry(self, textvariable=self.y_maxvar, width=15, font=MED_FONT)

        self.y_minval.grid(column=3, row=2)
        self.y_maxval.grid(column=3, row=1)

        self.autotxt = tk.StringVar()
        self.autotxt.set('Auto')
        self.autobtn = ttk.Button(self, textvariable=self.autotxt, command=self.set_auto)
        self.autobtn.grid(column=0, row=3, columnspan=3)


        self.choosefieldbl = ttk.Label(self, text='Show field: ', font=MED_FONT)
        self.choosefieldbl.grid(column=0, row=4, padx=10, pady=10)

        self.choosefieldcbox = ttk.Combobox(self, state='readonly', width=10)
        self.choosefieldcbox['values'] = self.parent.controller.eq.fieldNames + self.parent.controller.eq.auxFieldNames
        self.choosefieldcbox.bind('<<ComboboxSelected>>', self.change_field)
        self.choosefieldcbox.current(0)
        self.choosefieldcbox.grid(column=1, row=4, padx=10, pady=10)

        self.choosecolorlbl = ttk.Label(self, text='Colormap: ', font=MED_FONT)
        self.choosecolorlbl.grid(column=0, row=5, padx=10, pady=10)

        self.choosecolorbox = ttk.Combobox(self, state='readonly', width=10)
        self.choosecolorbox['values'] = COLORMAPS
        self.choosecolorbox.bind('<<ComboboxSelected>>', self.change_cmap)
        self.choosecolorbox.current(0)
        self.choosecolorbox.grid(column=1, row=5, padx=10, pady=10)

    def initSpatioTempPlot(self, stplot):
        self.stplot = stplot
        self.setIm(stplot.im)

    def changeSTPlot(self, stplot):
        self.stplot = stplot

        for limtype in ('x', 'y', 'v'):
            self.set_lim(*stplot.get_ax_lim(limtype), limtype)

        self.choosefieldcbox.current(stplot.active_field_indx)
        self.choosecolorbox.current(stplot.getCurrentColormap())
        self.autotxt.set("Auto" if self.stplot.auto_lim else "Manual")

    def setIm(self, im): # also interact with S-T Plot
        self.im = im
        vmin, vmax = self.im.get_clim()
        self.set_vlim(vmin, vmax)

    def get_minmax_type(self, type):
        if type == 'x':
            var_min, var_max = self.x_minvar, self.x_maxvar
        elif type == 'y':
            var_min, var_max = self.y_minvar, self.y_maxvar
        elif type == 'v':
            var_min, var_max = self.v_minvar, self.v_maxvar
        return var_min, var_max

    def get_lim(self, type):
        minvar, maxvar = self.get_minmax_type(type)
        try:
            _min = float(minvar.get())
            _max = float(maxvar.get())
            return _min, _max
        except ValueError:
            return None, None

    def set_lim(self, _min, _max, _type):
        var_min, var_max = self.get_minmax_type(_type)
        var_min.set('{:.5f}'.format(_min))
        var_max.set('{:.5f}'.format(_max))

    def set_xlim(self, xmin, xmax):
        self.set_lim(xmin, xmax, 'x')
    
    def set_ylim(self, ymin, ymax):
        self.set_lim(ymin, ymax, 'y')

    def set_vlim(self, vmin, vmax):
        self.set_lim(vmin, vmax, 'v')

    def get_xlim(self):
        return self.get_lim('x')

    def get_ylim(self):
        return self.get_lim('y')

    def get_vlim(self):
        return self.get_lim('v')

    def set_auto(self):
        if self.stplot.auto_lim:
            self.stplot.auto_lim = False
            self.autotxt.set('Manual')
        else:
            self.stplot.auto_lim = True
            self.autotxt.set('Auto')

    def change_field(self, event): # this should interact with S-T Plot // prep
        new_active_field_indx = self.choosefieldcbox.current()
        self.stplot.setActiveField(new_active_field_indx)
        self.choosecolorbox.current(self.stplot.getCurrentColormap())

    def change_cmap(self, event): # same sis // prep
        new_active_cmap = self.choosecolorbox.current()
        self.stplot.setCurrentColormap(new_active_cmap)

    def activate(self):
        self.grid(row=0, column=0)

    def deactivate(self):
        self.grid_forget()


class InspectorWindow(tk.Frame):

    def __init__(self, parent, controller, mainpg):
        tk.Frame.__init__(self, parent)

        self.parent = parent
        self.parent.title('Inspector')
        self.mainpg = mainpg
        self.controller = controller

        self.profile = InspectorProfile(self)
        self.spatiotemp = InspectorSpatiotemporal(self)

        #self.profile.grid(row=0, column=0)
        #self.spatiotemp.grid(row=0, column=0)

        self.activeFrame = None
        self.showProfile()
        self.grid(row=0, column=0)

    def showProfile(self):
        print('showing profile')
        self.showFrame(self.profile)
        if self.controller.activePlot is not None:
            pplot = self.controller.activePlot[-1]
            self.profile.changeProfilePlot(pplot)
        self.activeFrame = self.profile

    def showSpatiotemp(self):
        print('showing spatiotemp')
        self.showFrame(self.spatiotemp)
        if self.controller.activePlot is not None:
            self.spatiotemp.changeSTPlot(self.controller.activePlot[-1])
        self.activeFrame = self.spatiotemp

    def showFrame(self, frame):
        if self.activeFrame:
            self.activeFrame.deactivate()
        #frame.tkraise()
        frame.activate()

class EntryWindow:

    def __init__(self, root, mainpg, eq, k_sol, record=False):

        self.root = root
        self.eq = eq
        self.k_sol = k_sol
        self.mainpg = mainpg
        self.record = record

        self.txtlbl = ttk.Label(self.root, text='Insert name', font=MED_FONT)
        self.entry = ttk.Entry(self.root)
        self.entry.bind('<Return>', self.save)
        self.savebtn = ttk.Button(self.root, text='save', command=self.save)

        self.statustext = tk.StringVar()
        self.status = ttk.Label(self.root, textvariable=self.statustext, font=MED_FONT)

        self.warned = False

        self.txtlbl.grid(row=0, column=0, columnspan=2)
        self.entry.grid(row=1, column=0, padx=5, pady=5)
        self.savebtn.grid(row=1, column=1, padx=5, pady=5)
        self.status.grid(row=2, column=0, columnspan=2)

    def save(self, *args):
        name = self.entry.get()
        cond = self.eq.isFolder(name) if self.record else self.eq.isState(name)
        print(f'cond is {cond}')
        if (not cond) or self.warned:
            if self.record:
                self.mainpg.startRecord(name)
            else:
                self.eq.saveState(self.k_sol, name)
            self.die()
        else:
            self.statustext.set('Name already exists. Are you sure you want to overwrite it?')
            self.warned = True


    def die(self):
        self.mainpg.play()
        self.mainpg.kill_entry_window()


class ParameterWindow:

    def __init__(self, root, geometry, eq, update=None):

        self.root = root
        # self.root.geometry(geometry) # string: sizexsize
        self.eq = eq
        self.parameters = self.eq.parameters
        self.update = update

        self.titlelbl = ttk.Label(self.root, text='Parameters', width=15, font=LARGE_FONT)
        self.paramsvals = []
        self.paramslbls = []
        self.titlelbl.grid(column=0, row=0, columnspan=3)

        i = 0
        for p in self.parameters:
            lbl = ttk.Label(self.root, text=self.parameters[p].name, font=MED_FONT)
            lbl.grid(column=0, row=i + 1)
            self.paramslbls.append(lbl)

            vals = ttk.Entry(self.root, textvariable=self.parameters[p].var, font=MED_FONT, width=10)
            vals.grid(column=2, row=i + 1, padx=5)
            self.paramsvals.append(vals)

            if eq.initRange is not None:
                from_, to = eq.initRange[self.parameters[p].name]
            else:
                from_, to = 0, 1
            scale = ttk.Scale(self.root, orient=tk.HORIZONTAL, variable=self.parameters[p].doublevar, command=self.parameters[p].udv, from_=from_, to=to)
            scale.grid(column=1, row=i+1)
            i += 1

        self.resetBtn = ttk.Button(self.root, text='Reset', command=self.reset)
        self.resetBtn.grid(column=0, row=i+1)

    def reset(self):
        for p in self.eq.initParams:
            self.parameters[p].var.set(str(self.eq.initParams[p]))


class ProfilePlot:
    def __init__(self, fig, ax, profile_inspector, controller):
        self.fig = fig
        self.ax = ax
        self.profile_inspector = profile_inspector
        self.controller = controller
        self.eq = controller.eq

        self.profile_inspector.initProfilePlot(self)

        self.Fields = self.eq.getInitCondFields()
        self.auxFields = self.eq.getAuxFields(*self.Fields)

        self.active_field_indx = 0
        self.current_colors = [_i + 2 for _i in range(self.eq.n_fields + self.eq.n_aux_fields)]
        
        self.set_label('x', 'u')
        self.draw_fields()

        self.ymin_global = np.inf
        self.ymax_global = -np.inf

        self.auto_lim = True

    def set_label(self, xlabel, ylabel):
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)

    def get_shown_field_indx(self):
        """ Returns a list of indices of all shown fields """
        return [i for i, color in enumerate(self.current_colors) if color != 0]

    def getCurrentColorActiveField(self):
        return self.current_colors[self.active_field_indx]

    def setCurrentColorActiveField(self, color):
        self.current_colors[self.active_field_indx] = color
        self.redraw_fields()


    def auto_update_lims(self):
        xmin, xmax = self.eq.get_xlims()
        # get min and max between all currently shown fields
        ymin, ymax = self.eq.get_fields_lims(indices=self.get_shown_field_indx())
        ymin, ymax = d10p(ymin), i10p(ymax) # decrease and increase in 10%

        self.ymin_global = min(ymin, self.ymin_global)
        self.ymax_global = max(ymax, self.ymax_global)

        if self.isProfileInspectorActive(): # only if this plot is selected, update lims
            self.profile_inspector.set_xlim(xmin, xmax)
            self.profile_inspector.set_ylim(self.ymin_global, self.ymax_global)

        self.set_lims((xmin, xmax), (self.ymin_global, self.ymax_global))

    def set_lims(self, xlims, ylims):
        xmin, xmax = xlims
        ymin, ymax = ylims

        if xmin is not None:
            self.ax.set_xlim(xmin, xmax)
        
        if ymin is not None:
            self.ax.set_ylim(ymin, ymax)

    def get_ax_lims(self):
        return self.ax.get_xlim(), self.ax.get_ylim()

    def isProfileInspectorActive(self):
        return self.profile_inspector.pplot == self

    def manual_update_lims(self):
        if self.isProfileInspectorActive():
            self.set_lims(self.profile_inspector.get_xlims(), self.profile_inspector.get_ylims())

    def update_fields(self): # P Plot // feed fields?
        for i in range(len(self.Fields)):
            if self.lines[i] is not None:
                self.lines[i].set_ydata(self.Fields[i])
                self.lines[i].set_xdata(self.eq.x)
        for j in range(len(self.auxFields)):
            i = j + self.eq.n_fields 
            if self.lines[i] is not None:
                self.lines[i].set_ydata(self.auxFields[j])
                self.lines[i].set_xdata(self.eq.x)


    def draw_fields(self):
        self.lines = []
        for i in range(len(self.Fields)):
            curfield = self.Fields[i]
            color = self.current_colors[i]
            if COLORS[color] == 'none':
                line = None
            else:
                line, = self.ax.plot(self.eq.x, curfield, 'tab:' + COLORS[color])
            self.lines.append(line)
        for i in range(len(self.auxFields)):
            curfield = self.auxFields[i]
            color = self.current_colors[i+self.eq.n_fields]
            if COLORS[color] == 'none':
                line = None
            else:
                line, = self.ax.plot(self.eq.x, curfield, 'tab:' + COLORS[color])
            self.lines.append(line)

    def redraw_fields(self): # P Plot /// ok
        self.ax.cla()
        self.draw_fields()

    def resetActiveField(self):
        if not self.isActiveFieldAux():
            self.active_Field = self.Fields[self.active_field_indx]

    def activeFieldToCurrent(self):
        if not self.isActiveFieldAux():
            self.Fields[self.active_field_indx] = self.active_Field

    def saveEditandTick(self):
        self.activeFieldToCurrent()
        self.eq.tick(newInitCondFields=self.Fields)

    def update(self):
        self.Fields = self.eq.getCurrentFields()
        self.auxFields = self.eq.getCurrentAuxFields()
        self.resetActiveField()
        self.update_fields()

        if self.auto_lim:
            self.auto_update_lims()
        else:
            self.manual_update_lims()
        return self.lines

    def raw_xtoi(self, x):  # from raw_x to i
        return self.xtoi(x, isRaw=True)

    def xtoi(self, x, isRaw=False, truncate=False):
        x0 = self.controller.eq.x[0]
        xf = self.controller.eq.x[-1]
        Ni = round(self.controller.eq.N / self.controller.eq.n_fields)
        _N = RAW_N if isRaw else Ni
        i = (x - x0) / (xf - x0) * _N
        if not truncate:
            i = round(i)
        return int(i)

    def isActiveFieldAux(self):
        return self.active_field_indx >= self.eq.n_fields

    def updateY(self, xdata, ydata): # improve this pls
        if self.isActiveFieldAux(): return 
        xi = self.xtoi(xdata)
        self.resetActiveField()
        if xi < len(self.active_Field) - 4 and xi >= 4:

            ### interpolate
            raw_dx = (self.eq.x[-1] - self.eq.x[0]) / (RAW_N - 1)
            x0 = xdata - raw_dx / 2
            xf = xdata + raw_dx / 2

            i_0 = self.xtoi(x0)
            i_f = self.xtoi(xf)
            i_s = np.array([i_0, xi, i_f])

            xs_interval = self.eq.x[i_s]
            ys_interval = np.array([self.active_Field[i_0], ydata, self.active_Field[i_f]])
            ys_interp = interp1d(xs_interval, ys_interval)

            i_0 = self.xtoi(xdata - raw_dx / 2)
            i_f = self.xtoi(xdata + raw_dx / 2)

            # replace old values of ys
            self.active_Field[i_0:i_f + 1] = ys_interp(self.eq.x[i_0:i_f + 1])
        elif xi < 4:
            self.active_Field[0:4] = ydata
        else:
            self.active_Field[-4:] = ydata
        self.activeFieldToCurrent()
        self.update_fields()

class SpatioTemporalPlot():

    def __init__(self, fig, ax, spatiotemp_inspector, controller):
        self.fig = fig
        self.ax = ax
        self.spatiotemp_inspector = spatiotemp_inspector
        self.controller = controller
        self.eq = controller.eq
        self.k_st = 0 # k spatiotemp
        self.st_rows = ST_ROWS

        self.Fields = self.eq.getInitCondFields()
        self.auxFields = self.eq.getAuxFields(*self.Fields)

        self.auto_lim = True
        self.active_field_indx = 0
        self.current_cmap = [0] * (self.eq.n_fields + self.eq.n_aux_fields)

        self.imvals = np.zeros((self.st_rows, self.eq.getN()))
        self.plot()

        self.spatiotemp_inspector.initSpatioTempPlot(self)

    def plot(self):
        active_cmap = COLORMAPS[self.getCurrentColormap()]
        self.im = self.ax.imshow(self.imvals, cmap=active_cmap, extent=[self.eq.x[0], self.eq.x[-1], 0, self.st_rows], aspect='auto', origin='lower')
        self.colorbar = self.fig.colorbar(self.im, ax=self.ax, orientation='horizontal')
        self.set_label('x', 't (time steps)')

    def set_label(self, xlabel, ylabel):
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)

    def getCurrentColormap(self):
        return self.current_cmap[self.active_field_indx]

    def setCurrentColormap(self, color_indx):
        if color_indx != self.current_cmap[self.active_field_indx]:
            self.current_cmap[self.active_field_indx] = color_indx
            self.replot()

    def setActiveField(self, indx):
        if indx != self.active_field_indx:
            self.active_field_indx = indx
            self.reset()
            self.replot()

    def getActiveField(self):
        return self.eq.getCurrentField(self.active_field_indx)

    def replot(self):
        self.ax.clear()
        self.colorbar.remove()
        self.plot()

    def isSpatioTempInspectorActive(self):
        return self.spatiotemp_inspector.stplot == self

    def reset(self):
        self.imvals = np.zeros((self.st_rows, self.eq.getN()))
        self.k_st = 0
        self.imvals[self.k_st, :] = self.getActiveField()
        self.replot()
    
    def auto_update(self):

        #  calc lims        
        vmin = np.min(self.imvals)
        vmax = np.max(self.imvals)

        # update axes lims
        self.im.set_extent([self.eq.x[0], self.eq.x[-1], 0, self.st_rows])
        self.im.set_clim(vmin, vmax)

        # tell inspector
        if self.isSpatioTempInspectorActive():
            self.spatiotemp_inspector.set_xlim(self.eq.x[0], self.eq.x[-1])
            self.spatiotemp_inspector.set_ylim(0, self.st_rows)
            self.spatiotemp_inspector.set_vlim(vmin, vmax)

    def set_lim(self, limtype, _min, _max):
        if _min is None: return
        if limtype == 'x':
            self.ax.set_xlim(_min, _max)
        elif limtype == 'v':
            self.im.set_clim(_min, _max)
        elif limtype == 'y':
            if _max > self.st_rows:
                new_rows = int(_max) - ST_ROWS
                self.imvals = np.append(self.imvals, np.zeros((new_rows, self.eq.getN())), axis=0)
                self.st_rows += new_rows
                self.replot()
            elif int(_max) != self.st_rows and int(_max) !=  0:
                self.st_rows = int(_max)
                self.imvals = self.imvals[:self.st_rows, :]
                self.k_st = self.k_st % self.st_rows
                self.replot()
            self.ax.set_ylim(_min, _max) 

    def get_ax_lim(self, limtype):
        if limtype == 'x':
            return self.ax.get_xlim()
        elif limtype == 'y':
            return self.ax.get_ylim()
        elif limtype == 'v':
            return self.im.get_clim()

    def manual_update(self):
        if not self.isSpatioTempInspectorActive(): return
        for limtype in ('x', 'y', 'v'):
            _min, _max = self.spatiotemp_inspector.get_lim(limtype)
            self.set_lim(limtype, _min, _max)


    def update(self):
        self.imvals[self.k_st, :] = self.getActiveField()
        self.im.set_data(self.imvals)

        if self.auto_lim:
            self.auto_update()
        else:
            self.manual_update()

        self.k_st = (self.k_st + 1) % self.st_rows
        return self.im


class PlotWindow(tk.Frame):

    def __init__(self, parent, controller, mainpg):
        self.parent = parent
        self.controller = controller
        self.mainpg = mainpg
        self.inspector = mainpg.inspector

        self.clicked = False
        self.released = False
        self.anim_stopped = False
        self.isEditing = False
        self.isPaused = False

        tk.Frame.__init__(self, parent)

        self.fig = Figure(figsize=(10, 10), dpi=100)
        self.ax = self.fig.add_subplot(111)

        ## draw canvas

        self.container_mpl = tk.Frame(self)

        self.canvas = FigureCanvasTkAgg(self.fig, self.container_mpl)
        self.canvas.draw()

        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('button_release_event', self.off_click)
        self.canvas.mpl_connect('motion_notify_event', self.move_click)

        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # self.container_mpl.geometry('512x740')
        self.container_mpl.grid(row=0, column=0, sticky='ns')
        self.columnconfigure(0, weight=1)
        self.rowconfigure(0, weight=1)

    def start_anim(self):
        event_source = self.mainpg.ani.event_source
        self.ani = animation.FuncAnimation(self.fig, self.animate, frames=60, interval=100, blit=False, event_source=event_source)

    def on_click(self, event):
        pass

    def off_click(self, event):
        pass

    def move_click(self, event):
        pass

    def animate(self, i):
        pass

class SpatiotemporalWindow(PlotWindow):

    def __init__(self, parent, controller, mainpg, spatiotemporal_inspector):

        self.spatiotemporal_inspector = spatiotemporal_inspector
        PlotWindow.__init__(self, parent, controller, mainpg)

        self.stplot = SpatioTemporalPlot(self.fig, self.ax, spatiotemporal_inspector, controller)
        self.im  = self.stplot.im

    def animate(self, i):
        if self.isPaused or self.isEditing:
            return self.im,
        return self.stplot.update(),

    def on_click(self, event):
        if event.inaxes == self.ax:
            self.controller.activePlot = (self.ax, SPATIOTEMPORAL, self.stplot)
            self.inspector.showSpatiotemp()

    def off_click(self, event):
        pass

    def move_click(self, event):
        pass


class ProfileWindow(PlotWindow):

    def __init__(self, parent, controller, mainpg, profile_inspector):
        self.profile_inspector = profile_inspector
        PlotWindow.__init__(self, parent, controller, mainpg)
        self.pplot = ProfilePlot(self.fig, self.ax, self.profile_inspector, self.controller)
        self.lines = self.pplot.lines

    def animate(self, i):
        if self.isPaused or self.isEditing:
            return self.lines, 
        self.lines = self.pplot.update()
        return self.lines,

    def on_click(self, event):
        pplot = self.pplot
        if event.inaxes == pplot.ax: # call p plot
            self.controller.activePlot = (pplot.ax, PROFILE, pplot)
            self.inspector.showProfile()
            if self.isEditing:
                self.clicked = True
                pplot.updateY(event.xdata, event.ydata)

    def off_click(self, event):
        if self.clicked:
            self.clicked = False
            self.released = True
            self.isEditing = False
            # change this later
            pplot = self.controller.activePlot[-1]
            pplot.saveEditandTick()

    def move_click(self, event):
        if self.clicked:
            if event.inaxes == self.pplot.ax:
                self.pplot.updateY(event.xdata, event.ydata)

class MainPage(tk.Frame):

    def __init__(self, parent, controller):

        self.parent = parent
        self.controller = controller
        self.active = False
        self.anim_stopped = False

        tk.Frame.__init__(self, parent)

        self.btn1img = tk.PhotoImage(file=ICONSFOLDER + 'home.png')
        btn1 = ttk.Button(self, text='Back home',
                          command=self.back_home, image=self.btn1img)
        btn1.grid(row=0, column=0)

        initcond_container = tk.Frame(self)

        label = ttk.Label(self, text='Interactive PDE simulation', font=LARGE_FONT)
        label.grid(row=0, column=1)

        iclabel = ttk.Label(initcond_container, text='Current initial condition: ', font=MED_FONT)
        iclabel.grid(row=0, column=0, padx=5, pady=5)

        self.initcond_name = tk.StringVar(value='')
        ciclabel = ttk.Label(initcond_container, textvariable=self.initcond_name, font=MED_FONT)
        ciclabel.grid(row=0, column=1, padx=5, pady=5)

        btn_reset = ttk.Button(initcond_container, text='Reset',
                          command=self.controller.resetEqInitCond)

        btn_zero = ttk.Button(initcond_container, text='Zero', 
                            command=self.controller.setEqInitCondZero)
        
        btn_reset.grid(row=0, column=2, pady=5)
        btn_zero.grid(row=0, column=3, pady=5)
        
        initcond_container.grid(row=1, column=0, columnspan=3)

    def deactivate(self):

        self.active = False
        self.controller.activePlot = None
        self.btn_container.grid_forget()
        self.killParamWindow()
        self.killInspectorWindow()
        self.killImageWindow()
        self.killAuxWindows()
        self.stop_animation()

    def killAuxWindows(self):
        for parent_aw in self.parentNewWindows:
            parent_aw.destroy()

    def back_home(self):
        self.deactivate()
        self.controller.show_frame(StartPage)

    def activate(self):

        self.active = True

        self.mustClear = False # who are you
        self.clicked = False
        self.released = False
        self.i_release = 0 

        self.btn_container = tk.Frame(self)
        self.btn_container.grid(row=2, columnspan=3)

        self.playimg = tk.PhotoImage(file=ICONSFOLDER + 'play.png')
        self.editimg = tk.PhotoImage(file=ICONSFOLDER + 'pen2.png')
        self.pauseimg = tk.PhotoImage(file=ICONSFOLDER + 'pause.png')
        self.inspimg = tk.PhotoImage(file=ICONSFOLDER + 'info.png')
        self.eqimg = tk.PhotoImage(file=ICONSFOLDER + 'equation.png')
        self.paramimg = tk.PhotoImage(file=ICONSFOLDER + 'params.png')
        self.saveimg = tk.PhotoImage(file=ICONSFOLDER + 'save.png')
        self.recordimg = tk.PhotoImage(file=ICONSFOLDER + 'record1.png')
        self.recordingimg = tk.PhotoImage(file=ICONSFOLDER + 'record2.png')
        self.recordoffimg = tk.PhotoImage(file=ICONSFOLDER + 'record3.png')
        self.addimg = tk.PhotoImage(file = ICONSFOLDER + 'add.png')

        self.btn_e = ttk.Button(self.btn_container,
                                command=self.edit, image=self.editimg)
        self.btn_e.grid(row=0, column=0)

        self.btn_play = ttk.Button(self.btn_container,
                                   command=self.play, image=self.playimg)
        self.btn_play.grid(row=0, column=1)

        self.btn_pause = ttk.Button(self.btn_container,
                                    command=self.pause, image=self.pauseimg)
        self.btn_pause.grid(row=0, column=2)

        self.btn_pause = ttk.Button(self.btn_container,
                                    command=self.save, image=self.saveimg)
        self.btn_pause.grid(row=0, column=3)

        self.btn_record = ttk.Button(self.btn_container,
                                   command=self.record, image=self.recordimg)
        self.btn_record.grid(row=0, column=4)

        self.btn_params = ttk.Button(self.btn_container, 
                                     command=self.open_params, image=self.paramimg)
        self.btn_params.grid(row=0, column=5)

        self.btn_eq = ttk.Button(self.btn_container,
                                 command=self.open_equation, image=self.eqimg)
        self.btn_eq.grid(row=0, column=6)

        self.btn_insp = ttk.Button(self.btn_container,
                                   command=self.open_inspector, image=self.inspimg)
        self.btn_insp.grid(row=0, column=7)

        self.btn_add = ttk.Button(self.btn_container, command=self.add_window, image=self.addimg)
        self.btn_add.grid(row=0, column=8)

        self.recording = False

        ## assuming eq and init cond are already defined
        eq = self.controller.eq
        eq.updateX()

        ### create param window

        self.createParamWindow(eq)
        self.createImageWindow(eq)
        self.createInspectorWindow(eq)

        ### initialize plot and draw

        self.fig = Figure(figsize=(10, 10), dpi=100)
        self.ax = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212) # STW

        self.pplots = [ProfilePlot(self.fig, self.ax, self.inspector.profile, self.controller) ]
        self.stplot = SpatioTemporalPlot(self.fig, self.ax2, self.inspector.spatiotemp, self.controller)

        self.eqX = eq.x
        self.Fields = eq.getInitCondFields()

        ## Draw canvas

        self.container_mpl = tk.Frame(self)

        self.canvas = FigureCanvasTkAgg(self.fig, self.container_mpl)
        self.canvas.draw()

        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('button_release_event', self.off_click)
        self.canvas.mpl_connect('motion_notify_event', self.move_click)

        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        # self.container_mpl.geometry('512x740')
        self.container_mpl.grid(row=3, column=0, columnspan=3, sticky='ns')


        # self.rowconfigure(0, weight=1)
        self.rowconfigure(3, weight=1)

        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.columnconfigure(2, weight=1)

        self.isEditing = False
        self.isPaused = False
        self.t = 0 # same sis

        self.parentNewWindows = []
        self.newWindows = []
        ### Start animation

        # save event_source
        self.ani = animation.FuncAnimation(self.fig, self.animate, frames=60, interval=100, blit=False)

    def set_init_cond_zero(self):
        self.controller.setEqInitCondZero()

    def add_window(self):
        parent_new_window = tk.Toplevel(self)
        parent_new_window.geometry('512x512')
        self.parentNewWindows.append(parent_new_window)
        if self.controller.activePlot[1] == PROFILE:
            pw = ProfileWindow(parent_new_window, self.controller, self, self.inspector.profile)
        elif self.controller.activePlot[1] == SPATIOTEMPORAL:
            pw = SpatiotemporalWindow(parent_new_window, self.controller, self, self.inspector.spatiotemp)
        pw.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.newWindows.append(pw)
        pw.start_anim()

    def record(self):
        if self.recording:
            self.btn_record.configure(image=self.recordimg)
            self.controller.eq.stopRecording()
            self.recording = False
            print('Stop Recording...')
        else:
            self.pause()
            self.create_entry_window(record=True)

    def animate_record(self):
        self.i_record += 1
        if self.i_record >= 4:
            img = self.recordoffimg if self.img_recording else self.recordingimg
            self.btn_record.configure(image=img)
            self.img_recording = not self.img_recording
            self.i_record = 0

    def startRecord(self, recordname):
        self.controller.eq.startRecording(recordname)
        self.recording = True
        self.i_record = 0
        self.img_recording = True
        self.btn_record.configure(image=self.recordingimg)

    def edit(self): # send signal to P Plot that is editing
        print(f'Edit clicked!, isEditing = {not self.isEditing}')
        if self.isEditing:
            self.isEditing = False
        else:
            self.isEditing = True

    def pause(self):
        if self.isPaused: return
        self.isPaused = True

    def play(self):
        if not self.isPaused: return
        self.isPaused = False

    def save(self):
        self.pause()
        self.create_entry_window()
        #self.controller.eq.saveState(self.k_sol, f'k_{self.k_sol}')

    def open_params(self):
        try:
            if self.parentParamWindow.state() == 'normal':
                self.parentParamWindow.focus_set()
        except Exception:
            self.createParamWindow(self.controller.eq)

    def open_equation(self):
        try:
            if self.parentImageWindow.state() == 'normal':
                self.parentImageWindow.focus_set()
        except Exception:
            self.createImageWindow(self.controller.eq)

    def open_inspector(self):
        try:
            if self.parentInspectorWindow.state() == 'normal':
                self.parentInspectorWindow.focus_set()
        except Exception:
            self.createInspectorWindow(self.controller.eq)

    def create_entry_window(self, record=False):
        self.parentEntryWindow = tk.Toplevel(self)
        self.entrywindow = EntryWindow(self.parentEntryWindow, self, self.controller.eq, self.controller.eq.k_sol, record=record)

    def kill_entry_window(self):
        self.parentEntryWindow.destroy()

    def on_click(self, event):
        for pplot in self.pplots:
            if event.inaxes == pplot.ax: # call p plot
                self.controller.activePlot = (pplot.ax, PROFILE, pplot)
                self.inspector.showProfile()
                if self.isEditing:
                    self.clicked = True
                    pplot.updateY(event.xdata, event.ydata)
        if event.inaxes == self.ax2:
            self.controller.activePlot = (self.ax2, SPATIOTEMPORAL, self.stplot)
            self.inspector.showSpatiotemp()

    def off_click(self, event):
        if self.clicked:
            self.clicked = False
            self.released = True
            self.isEditing = False
            # change this later
            pplot = self.controller.activePlot[-1]
            pplot.saveEditandTick()

    def move_click(self, event):
        if self.clicked:
            for pplot in self.pplots:
                if event.inaxes == pplot.ax:
                    pplot.updateY(event.xdata, event.ydata)

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

    def getAnimIterables(self): # this will change
        lines = []
        for pplot in self.pplots:
            lines += pplot.lines
        return lines + [self.stplot.im]

    def animate(self, i):
        if not self.active: return
        eq = self.controller.eq
        eq.updateX()

        # Paused part should stay here and isEditing go to PPlot
        if self.isPaused or self.isEditing:
            return self.getAnimIterables()

        eq.tick()
        self.Fields = eq.getCurrentFields() # equation tick /// ok
        self.inspector.profile.setTime(self.t) # PPlot + equation i guesss

        for pplot in self.pplots:
            pplot.update()

        self.stplot.update()
        if self.recording: # in equation tick() /// ok
            self.animate_record() # except for this 

        self.t += eq.getParam('dt') # PPlot
        return self.getAnimIterables()


