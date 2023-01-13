from numpy.lib.arraysetops import isin
from equations.equation import SOLVE_EVERY_TI
import tkinter as tk
from tkinter import StringVar, ttk
from tkinter.constants import ACTIVE

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
PHASEPLOT = 'phaseplot'
TEMPORALPLOT = 'temporalplot'
IMAGESFOLDER = 'images/'

ICONSFOLDER = 'icons/'

style.use('ggplot')

ST_ROWS = 200
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
        label = ttk.Label(self, text='Interactive Differential Equations Analysis Simulator', font=LARGE_FONT)
        label.grid(row=0, column=0, columnspan=3,  padx=10, pady=5)

        self.parent = parent
        self.controller = controller

        self.dim_listbox = tk.Listbox(self, width=15)
        dims = list(ALL_EQUATIONS)
        dims.sort()
        for dim in dims:
            self.dim_listbox.insert(tk.END, f'Dim {dim}')

        self.dim_listbox.bind('<<ListboxSelect>>', self.change_select_dim)
        self.dim_listbox.selection_set(0); self.dim=0
        self.dim_listbox.grid(row=2, column=0)

        self.n_var = tk.StringVar(value='200')
        self.n_ent = tk.Entry(self, textvariable=self.n_var)
        self.n_ent.grid(row=3, column=1, padx=5, pady=5, sticky='w')

        self.search_var = tk.StringVar(value='search')
        self.search_var.trace_add('write', self.search_callback)
        self.search_ent = tk.Entry(self, textvariable=self.search_var)
        self.search_ent.grid(row=3, column=2, padx=5, pady=5, sticky='w')

        self.eq_listbox = tk.Listbox(self)
        self.fill_eq_listbox(0, init=True)

        self.dim_lbl = ttk.Label(self, text='Choose dimension', font=MED_FONT)
        self.eq_lbl = ttk.Label(self, text='Choose equation', font=MED_FONT)
        self.ic_lbl = ttk.Label(self, text='Choose initial condition', font=MED_FONT)

        self.dim_lbl.grid(row=1, column=0, padx=5, pady=5)
        self.eq_lbl.grid(row=1, column=1, padx=5, pady=5)
        self.ic_lbl.grid(row=1, column=2, sticky='ew', padx=5, pady=5)

        self.eq_listbox.bind('<<ListboxSelect>>', self.change_select_eq)
        self.eq_listbox.grid(row=2, column=1, padx=5, pady=5)

        self.initcond_listbox = tk.Listbox(self, width=30)
        self.initcond_listbox.bind('<<ListboxSelect>>', self.change_select_ic)
        self.fill_initcond_listbox()       

        self.initcond_listbox.grid(row=2, column=2, padx=5, pady=5, sticky='w')

        self.n_lbl = ttk.Label(self, text='Number of pts =  ', font=MED_FONT)
        self.n_lbl.grid(row=3, column=0, padx=5, pady=5, sticky='e')

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

    def init_eq(self):
        name = self.eq_listbox.get(0)
        self.controller.setEq(name, dim=0)
        print(name)
        self.fill_initcond_listbox()

    def clear_listbox(self, lbox):
        last = lbox.size() - 1
        if last >= 0:
            lbox.delete(0, last=last)

    def search_callback(self, *args):
        search_string = self.search_var.get()
        if search_string == 'search' or search_string == '': return
        self.fill_initcond_listbox(filtr=search_string)

    def fill_initcond_listbox(self, filtr=None):
        self.clear_listbox(self.initcond_listbox)
        initconds = self.controller.getEqInitConds()
        if filtr is not None:
            initconds = [[initcond, i] for i, initcond in enumerate(initconds) if filtr in initcond]
        else:
            initconds = [[initcond, i] for i, initcond in enumerate(initconds)]
        self.list_initconds = initconds
        for initcond, i in initconds:
            self.initcond_listbox.insert(tk.END, initcond)
        self.initcond_listbox.selection_set(0)

    def fill_eq_listbox(self, dim, init=False):
        self.clear_listbox(self.eq_listbox)
        eqs = list(ALL_EQUATIONS[dim])
        eqs.sort()
        for eqname in eqs:
            self.eq_listbox.insert(tk.END, eqname)
        self.eq_listbox.selection_set(0)
        if init:
            self.controller.setEq(eqs[0], dim=dim)
            self.n_var.set(self.controller.eq.getN())

    def change_select_dim(self, event):
        selection = event.widget.curselection()
        if selection:
            index = selection[0]
            self.dim = int(event.widget.get(index)[-1])
            self.fill_eq_listbox(self.dim, init=True)
            self.fill_initcond_listbox()


    def change_select_eq(self, event):
        selection = event.widget.curselection()
        if selection:
            index = selection[0]
            name = event.widget.get(index)
            self.controller.setEq(name, dim=self.dim)
            self.n_var.set(self.controller.eq.getN())
            self.fill_initcond_listbox()

    def change_select_ic(self, event):
        selection = event.widget.curselection()
        if selection:
            i = selection[0]
            index = self.list_initconds[i][1]
            if not self.controller.initCond_isMethod(index):
                self.controller.setEqInitCond(index)
                Ni = self.controller.getEqInitCond_N()
                self.n_var.set(str(Ni))
                self.n_ent.configure(state='readonly')
                print(f'Selected {self.controller.current_initcond_name}, ')
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

        self.chooseupdtlbl = ttk.Label(self, text='Update: ', font=MED_FONT)
        self.chooseupdtlbl.grid(column=2, row=4, padx=10, pady=10)

        self.chooseupdtcbox = ttk.Combobox(self, state='readonly', width=10)
        self.chooseupdtcbox['values'] = ['always'] + self.parent.controller.eq.st_update_optns
        self.chooseupdtcbox.bind('<<ComboboxSelected>>', self.change_update_optn)
        self.chooseupdtcbox.current(0)
        self.chooseupdtcbox.grid(column=3, row=4, padx=10, pady=10)

        self.Tlbl = ttk.Label(self, text="T", font=MED_FONT)
        self.Tlbl.grid(column=4, row=4, padx=10, pady=10)

        self.Tvar = tk.StringVar()
        self.Tentry = ttk.Entry(self, textvariable=self.Tvar, width=10, font=MED_FONT, state=tk.DISABLED)
        self.Tentry.grid(column=5, row=4, padx=10, pady=10)

        
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

    def get_T_update(self):
        T = self.Tvar.get()
        try:
            T = float(T)
            return T
        except ValueError:
            return None

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

    def change_update_optn(self, event):
        new_optn = self.chooseupdtcbox.current()
        if new_optn == 0:
            self.stplot.init_default()
        elif new_optn == 1: # 'strobe'
            self.stplot.init_strobe()
        elif new_optn == 2:
            self.stplot.init_poincare()

    def activate(self):
        self.grid(row=0, column=0)

    def deactivate(self):
        self.grid_forget()

class InspectorTemporal(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)

        self.parent = parent

        self.titlelbl = ttk.Label(self, text='Inspector: Temporal Plot (Oscilloscope)', font=LARGE_FONT)
        self.titlelbl.grid(column=0, row=0, columnspan=3, padx=20, pady=10)

        self.y_minlbl = ttk.Label(self, text='u_min', font=MED_FONT)
        self.y_maxlbl = ttk.Label(self, text='u_max', font=MED_FONT)

        self.y_minlbl.grid(column=2, row=2)
        self.y_maxlbl.grid(column=2, row=1)

        y_min, y_max = -1, 1

        self.y_minvar = tk.StringVar(value=str(y_min))

        self.y_maxvar = tk.StringVar(value=str(y_max))

        self.y_minval = ttk.Entry(self, textvariable=self.y_minvar, font=MED_FONT, width=15)
        self.y_maxval = ttk.Entry(self, textvariable=self.y_maxvar, font=MED_FONT, width=15)

        self.y_minval.grid(column=3, row=2, padx=10)
        self.y_maxval.grid(column=3, row=1, padx=10)

        self.x_minlbl = ttk.Label(self, text='t_min', font=MED_FONT)
        self.x_maxlbl = ttk.Label(self, text='t_max', font=MED_FONT)

        self.x_minlbl.grid(column=0, row=2)
        self.x_maxlbl.grid(column=0, row=1)

        x_min, x_max = -1, 1

        self.x_minvar = tk.StringVar(value=str(x_min))

        self.x_maxvar = tk.StringVar(value=str(x_max))

        self.x_minval = ttk.Entry(self, textvariable=self.x_minvar, font=MED_FONT, width=15)
        self.x_maxval = ttk.Entry(self, textvariable=self.x_maxvar, font=MED_FONT, width=15)

        self.x_minval.grid(column=1, row=2, padx=10)
        self.x_maxval.grid(column=1, row=1, padx=10)

        self.autotxt = tk.StringVar()
        self.autotxt.set('Auto')
        self.autobtn = ttk.Button(self, textvariable=self.autotxt, command=self.set_auto)
        self.autobtn.grid(column=3, row=3, padx=10, pady=10)

        self.Nvar = tk.StringVar(value=f'Show unit: (N = {str(self.parent.controller.eq.getN())})')
        self.Nval = ttk.Label(self, textvariable=self.Nvar)

        self.Nval.grid(row=3, column=0, padx=10, pady=10)

#        self.unitlbl = ttk.Label(self, text='Show unit: ', font=MED_FONT)
#        self.unitlbl.grid(row=4, column=2, padx=10, pady=10)

        self.unitvar = StringVar(value='0')
        self.unitentry = ttk.Entry(self, textvariable=self.unitvar, font=MED_FONT, width=10)
        self.unitentry.grid(row=3, column=1, padx=10, pady=10)

        self.choosefieldbl = ttk.Label(self, text='Show field: ', font=MED_FONT)
        self.choosefieldbl.grid(column=0, row=4, padx=10, pady=10)

        self.active_field_indx = 0

        self.choosefieldcbox = ttk.Combobox(self, state='readonly', width=10)
        self.choosefieldcbox['values'] = self.parent.controller.eq.fieldNames
        self.choosefieldcbox.bind('<<ComboboxSelected>>', self.change_field)
        self.choosefieldcbox.current(0)
        self.choosefieldcbox.grid(column=1, row=4)

        self.choosecolorlbl = ttk.Label(self, text='Color: ', font=MED_FONT)
        self.choosecolorlbl.grid(column=2, row=4)

        self.choosecolorbox = ttk.Combobox(self, state='readonly', width=10)
        self.choosecolorbox['values'] = COLORS
        self.choosecolorbox.bind('<<ComboboxSelected>>', self.change_color)
        self.choosecolorbox.current(2)
        self.choosecolorbox.grid(column=3, row=4)

    def initTemporalPlot(self, tplot):
        self.tplot = tplot
        xlims, ylims = self.tplot.get_ax_lims()
        self.set_xlim(*xlims)
        self.set_ylim(*ylims)
        self.last_unit = self.tplot.unit
        self.unitvar.set(self.last_unit)


    def updateN(self):
        self.Nvar.set(f'Show unit: (N = {str(self.parent.controller.eq.getN())})')

    def getUnit(self):
        try:
            unit = int(self.unitvar.get())
        except ValueError:
            return self.last_unit
        if unit >= self.parent.controller.eq.getN():
            return self.last_unit
        self.last_unit = unit
        return unit

    def get_kmax(self):
        try:
            tmax = float(self.x_maxvar.get())
            return self.tmax_to_kmax(tmax)
        except ValueError:
            return None

    def kmax_to_tmax(self, kmax):
        dt = self.parent.controller.eq.getParam('dt')
        return round(dt * kmax, 4)

    def tmax_to_kmax(self, tmax):
        dt = self.parent.controller.eq.getParam('dt')
        return round(tmax/dt)

    def changePlot(self, tplot):
        self.tplot = tplot

        i = self.active_field_indx
        self.choosefieldcbox.current(i)
        self.choosecolorbox.current(self.tplot.getColorField(i))
        self.autotxt.set('Auto' if self.tplot.auto_lim else 'Manual')
        self.initTemporalPlot(tplot)

    def change_field(self, event):
        self.active_field_indx = self.choosefieldcbox.current()
        activeColor = self.tplot.getColorField(self.active_field_indx)
        self.choosecolorbox.current(activeColor)

    def change_color(self, event):
        new_color = self.choosecolorbox.current()
        cur_color = self.tplot.getColorField(self.active_field_indx)
        if cur_color != new_color:
            self.tplot.setColorField(new_color, self.active_field_indx)
        
    def set_auto(self):
        if self.tplot.auto_lim:
            self.tplot.auto_lim = False
            self.autotxt.set('Manual')
        else:
            self.tplot.auto_lim = True
            self.autotxt.set('Auto')
    
    def get_xlims(self):
        try:
            y_min = float(self.x_minvar.get())
            y_max = float(self.x_maxvar.get())
            return y_min, y_max
        except ValueError:
            return None, None

    def get_ylims(self):
        try:
            x_min = float(self.y_minvar.get())
            x_max = float(self.y_maxvar.get())
            return x_min, x_max
        except ValueError:
            return None, None

    def set_ylim(self, ymin, ymax):
        self.y_minvar.set('{:.5f}'.format(ymin))
        self.y_maxvar.set('{:.5f}'.format(ymax))

    def set_xlim(self, xmin, xmax):
        self.x_minvar.set('{:.5f}'.format(xmin))
        self.x_maxvar.set('{:.5f}'.format(xmax))

    def activate(self):
        self.grid(row=0, column=0)

    def deactivate(self):
        self.grid_forget()       

class InspectorPhase(tk.Frame):

    def __init__(self, parent):
        tk.Frame.__init__(self, parent)

        self.parent = parent

        self.titlelbl = ttk.Label(self, text='Inspector: Phase Plot', font=LARGE_FONT)
        self.titlelbl.grid(column=0, row=0, columnspan=3, padx=10, pady=10)

        self.is3D = False
        if self.parent.controller.eq.n_fields >= 3:
            self.val3D = tk.StringVar(value='3D')
            self.btn3D = ttk.Button(self, textvariable=self.val3D, command=self.plot3D)
            self.btn3D.grid(column=3, row=0, padx=5, sticky='e')

        self.y_minlbl = ttk.Label(self, text='y_min', font=MED_FONT)
        self.y_maxlbl = ttk.Label(self, text='y_max', font=MED_FONT)

        self.y_minlbl.grid(column=2, row=2)
        self.y_maxlbl.grid(column=2, row=1)

        y_min, y_max = -1, 1

        self.y_minvar = tk.StringVar(value=str(y_min))

        self.y_maxvar = tk.StringVar(value=str(y_max))

        self.y_minval = ttk.Entry(self, textvariable=self.y_minvar, font=MED_FONT, width=15)
        self.y_maxval = ttk.Entry(self, textvariable=self.y_maxvar, font=MED_FONT, width=15)

        self.y_minval.grid(column=3, row=2, padx=10)
        self.y_maxval.grid(column=3, row=1, padx=10)

        self.x_minlbl = ttk.Label(self, text='x_min', font=MED_FONT)
        self.x_maxlbl = ttk.Label(self, text='x_max', font=MED_FONT)

        self.x_minlbl.grid(column=0, row=2)
        self.x_maxlbl.grid(column=0, row=1)

        x_min, x_max = -1, 1

        self.x_minvar = tk.StringVar(value=str(x_min))

        self.x_maxvar = tk.StringVar(value=str(x_max))

        self.x_minval = ttk.Entry(self, textvariable=self.x_minvar, font=MED_FONT, width=15)
        self.x_maxval = ttk.Entry(self, textvariable=self.x_maxvar, font=MED_FONT, width=15)

        self.x_minval.grid(column=1, row=2, padx=10)
        self.x_maxval.grid(column=1, row=1, padx=10)

        self.autotxt = tk.StringVar()
        self.autotxt.set('Auto')
        self.autobtn = ttk.Button(self, textvariable=self.autotxt, command=self.set_auto)
        self.autobtn.grid(column=3, row=3, padx=10, pady=10)

        self.clearbtn = ttk.Button(self, text='Clear', command=self.clear)
        self.clearbtn.grid(column=0, row=3, padx=10, pady=10)

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



        self.ellapsed_time_desc_lbl = ttk.Label(self, text='Ellapsed time: ', font=MED_FONT)
        self.ellapsed_time_desc_lbl.grid(row=3, column=1, padx=10, pady=5)

        self.time_container = tk.Frame(self)

        self.resetimg = tk.PhotoImage(file=ICONSFOLDER + 'reset20.png')

        self.ellapsed_time_var = tk.StringVar(value='0')
        self.ellapsed_time_val_lbl = ttk.Label(self.time_container, textvariable=self.ellapsed_time_var)
        self.ellapsed_time_val_lbl.grid(row=0, column=1, padx=5, pady=5)

        self.elapsed_time_reset_btn = ttk.Button(self.time_container, text='Reset', command=self.reset_time, image=self.resetimg)
        self.elapsed_time_reset_btn.grid(row=0, column=2, padx=5, pady=5)

        self.time_container.grid(row=3, column=2)

        self.xfieldlbl = ttk.Label(self, text='x axis:', font=MED_FONT)
        self.yfieldlbl = ttk.Label(self, text='y axis:', font=MED_FONT)
        self.xfieldlbl.grid(row=6, column=0, padx=5, pady=5)
        self.yfieldlbl.grid(row=6, column=2, padx=5, pady=5)
        
        self.xfieldcbox = ttk.Combobox(self, state='readonly', width=10)
        fieldNames = self.parent.controller.eq.fieldNames
        self.xfieldcbox['values'] = fieldNames
        self.xfieldcbox.bind('<<ComboboxSelected>>', self.change_x_field)
        self.xfieldcbox.current(0)
        self.xfieldcbox.grid(row=6, column=1, padx=5, pady=5)

        self.yfieldcbox = ttk.Combobox(self, state='readonly', width=10)
        fieldNames = self.parent.controller.eq.fieldNames
        self.yfieldcbox['values'] = fieldNames
        self.yfieldcbox.bind('<<ComboboxSelected>>', self.change_y_field)
        self.yfieldcbox.current(1)
        self.yfieldcbox.grid(row=6, column=3, padx=5, pady=5)

    def inspector_3D_to_2D(self):
        self.val3D.set('3D')
        self.is3D = False

        # Remove zmin zmax
        for widget in (self.z_maxlbl, self.z_minlbl, self.z_minval, self.z_maxval, self.zfieldcbox, self.zfieldlbl):
            widget.grid_forget()

    def inspector_2D_to_3D(self):
        self.val3D.set('2D')
        self.is3D = True

        # Add zmin zmax to inspector
        self.z_minlbl = ttk.Label(self, text='z_min', font=MED_FONT)
        self.z_maxlbl = ttk.Label(self, text='z_max', font=MED_FONT)

        self.z_minlbl.grid(column=4, row=2)
        self.z_maxlbl.grid(column=4, row=1)

        z_min, z_max = -1, 1

        self.z_minvar = tk.StringVar(value=str(z_min))

        self.z_maxvar = tk.StringVar(value=str(z_max))

        self.z_minval = ttk.Entry(self, textvariable=self.z_minvar, font=MED_FONT, width=15)
        self.z_maxval = ttk.Entry(self, textvariable=self.z_maxvar, font=MED_FONT, width=15)

        self.z_minval.grid(column=5, row=2, padx=10)
        self.z_maxval.grid(column=5, row=1, padx=10)

        self.zfieldlbl = ttk.Label(self, text='z axis:', font=MED_FONT)
        self.zfieldlbl.grid(row=6, column=4, padx=5, pady=5)
        
        self.zfieldcbox = ttk.Combobox(self, state='readonly', width=10)
        fieldNames = self.parent.controller.eq.fieldNames
        self.zfieldcbox['values'] = fieldNames
        self.zfieldcbox.bind('<<ComboboxSelected>>', self.change_z_field)
        self.zfieldcbox.current(2)
        self.zfieldcbox.grid(row=6, column=5, padx=5, pady=5)

    def plot3D(self):
        if self.is3D:
        
            self.inspector_3D_to_2D()

            # Remove plot
            fig = self.pplot.fig
            ax = self.pplot.ax
            fig.delaxes(ax)

            if isinstance(self.pplot.parent, MainPage):
                new_ax = fig.add_subplot(211)
            else:
                new_ax = fig.add_subplot(111)
            
            new_pplot = PhasePlot(fig, new_ax, self.pplot.mainpg, parent=self.pplot.parent)

            if isinstance(self.pplot.parent, MainPage):
                self.pplot.mainpg.plots[0] = new_pplot
            else:
                self.pplot.parent.plot = new_pplot
            
            self.changePlot(new_pplot)

        else:

            # Improve 2D plot to 3D plot
            fig = self.pplot.fig
            ax = self.pplot.ax
            fig.delaxes(ax)

            if isinstance(self.pplot.parent, MainPage):
                ax3D = fig.add_subplot(211, projection='3d')
            else:
                ax3D = fig.add_subplot(111, projection='3d')
            
            new_pplot = Phase3DPlot(fig, ax3D, self.pplot.mainpg, parent=self.pplot.parent)

            self.inspector_2D_to_3D()

            if isinstance(self.pplot.parent, MainPage):
                self.pplot.mainpg.plots[0] = new_pplot
            else:
                self.pplot.parent.plot = new_pplot

            self.changePlot(new_pplot)


    def change_x_field(self, event):
        new_x_field = self.xfieldcbox.current()
        if new_x_field != self.pplot.x_field:
            self.pplot.changeShownFields(new_x_field=new_x_field)
    
    def change_y_field(self, event):
        new_y_field = self.yfieldcbox.current()
        if new_y_field != self.pplot.y_field:
            self.pplot.changeShownFields(new_y_field=new_y_field)
    
    def change_z_field(self, event):
        new_z_field = self.zfieldcbox.current()
        if new_z_field != self.pplot.z_field:
            self.pplot.changeShownFields(new_z_field=new_z_field)

    def initPhasePlot(self, pplot):
        self.pplot = pplot
        self.autotxt.set('Auto' if self.pplot.auto_lim else 'Manual')
        lims = self.pplot.get_ax_lims()
        self.set_xlim(*lims[0])
        self.set_ylim(*lims[1])
        if self.is3D:
            self.set_zlim(*lims[2])

    def update_xy_fields(self):
        x_field = self.pplot.x_field
        y_field = self.pplot.y_field
        self.xfieldcbox.current(x_field)
        self.yfieldcbox.current(y_field)

        if self.is3D:
            z_field = self.pplot.z_field
            self.zfieldcbox.current(z_field)

    def changePlot(self, pplot):
        if isinstance(pplot, PhasePlot) and self.is3D:
            self.inspector_3D_to_2D()

        if isinstance(pplot, Phase3DPlot) and not self.is3D:
            self.inspector_2D_to_3D()

        self.initPhasePlot(pplot)
        self.update_xy_fields()

    def change_solver(self, event):
        solver_index = self.combobox.current()
        int_method = INTEGRATION_METHODS[solver_index]
        print(f'Selected {int_method.name}')
        # change description
        self.description_text.set(int_method.description)

        # change solver
        self.parent.controller.eq.setSolver(int_method) # call controller

    def reset_time(self): # interact w P plot
        self.parent.mainpg.t = 0

    def setTime(self, t):
        self.ellapsed_time_var.set('{:.2f}'.format(t))

    def set_auto(self):
        if self.pplot.auto_lim:
            self.pplot.auto_lim = False
            self.autotxt.set('Manual')
        else:
            self.pplot.auto_lim = True
            self.autotxt.set('Auto')
    
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

    def get_zlims(self):
        if not self.is3D: return
        try:
            z_min = float(self.z_minvar.get())
            z_max = float(self.z_maxvar.get())
            return z_min, z_max
        except ValueError:
            return None, None

    def set_ylim(self, ymin, ymax):
        self.y_minvar.set('{:.5f}'.format(ymin))
        self.y_maxvar.set('{:.5f}'.format(ymax))

    def set_xlim(self, xmin, xmax):
        self.x_minvar.set('{:.5f}'.format(xmin))
        self.x_maxvar.set('{:.5f}'.format(xmax))

    def set_zlim(self, zmin, zmax):
        if not self.is3D: return
        self.z_minvar.set('{:.5f}'.format(zmin))
        self.z_maxvar.set('{:.5f}'.format(zmax))

    def activate(self):
        self.grid(row=0, column=0)

    def deactivate(self):
        self.grid_forget()

    def clear(self):
        self.pplot.clear()

class InspectorWindow(tk.Frame):

    def __init__(self, parent, controller, mainpg):
        tk.Frame.__init__(self, parent)

        self.parent = parent
        self.parent.title('Inspector')
        self.mainpg = mainpg
        self.controller = controller

        #self.profile.grid(row=0, column=0)
        #self.spatiotemp.grid(row=0, column=0)

        self.activeFrame = None
        if self.controller.eq.dim == 1:
            self.profile = InspectorProfile(self)
            self.spatiotemp = InspectorSpatiotemporal(self)
            self.showProfile()
        elif self.controller.eq.dim == 0:
            self.temporal = InspectorTemporal(self)
            self.phase = InspectorPhase(self)
            self.showPhase()
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

    def showTemporal(self):
        print('showing temp')
        self.showFrame(self.temporal)
        if self.controller.activePlot is not None:
            self.temporal.changePlot(self.controller.activePlot[-1])
        self.activeFrame = self.temporal

    def showPhase(self):
        print('showing phase')
        self.showFrame(self.phase)
        if self.controller.activePlot is not None:
            self.phase.changePlot(self.controller.activePlot[-1])
        self.activeFrame = self.phase

    def show(self, kind):
        if kind == PROFILE:
            self.showProfile()
        elif kind == SPATIOTEMPORAL:
            self.showSpatiotemp()
        elif kind == TEMPORALPLOT:
            self.showTemporal()
        elif kind == PHASEPLOT:
            self.showPhase()
        # elif Phase

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

        self.txtlbl = ttk.Label(self.root, text='Enter name', font=MED_FONT)
        self.entry = ttk.Entry(self.root, width=15)
        self.entry.bind('<Return>', self.save)
        self.savebtn = ttk.Button(self.root, text='save', command=self.save)
        
        add_row = 0
        if record:
            add_row += 1
            self.interval_lbl = ttk.Label(self.root, text='Save every X time steps')

            self.interval_val = tk.StringVar(value='1')
            self.interval_entry = ttk.Entry(self.root, textvariable=self.interval_val, width=15)
            self.interval_lbl.grid(row=1, column=0)
            self.interval_entry.grid(row=1, column=1)

        self.statustext = tk.StringVar()
        self.status = ttk.Label(self.root, textvariable=self.statustext, font=MED_FONT)

        self.savefig_val = tk.BooleanVar()
        self.savefig_checkbox = ttk.Checkbutton(self.root, text='Save figure?',
                                    variable=self.savefig_val)
        self.savefig_checkbox.grid(row=1+add_row, column=1, columnspan=1)

        self.savestate_val = tk.BooleanVar(value=True)
        self.savestate_checkbox = ttk.Checkbutton(self.root, text='Save State?',
                                    variable=self.savestate_val)
        self.savestate_checkbox.grid(row=1+add_row, column=0, columnspan=1)

        self.warned = False

        self.txtlbl.grid(row=0, column=0, padx=5, pady=5)
        self.entry.grid(row=0, column=1, padx=5, pady=5)
        
        self.savebtn.grid(row=3+add_row, column=1, padx=5, pady=5)
        self.status.grid(row=4+add_row, column=0, columnspan=2)

    def save(self, *args):
        name = self.entry.get()
        alreadyExists = self.eq.isFolder(name) if self.record else self.eq.isState(name)
        try:
            if self.record:
                interval = int(self.interval_val.get())
            isIntervalOk = True
        except ValueError:
            isIntervalOk = False

        print(f'alreadyExists is {alreadyExists}')
        if ((not alreadyExists) or self.warned) and isIntervalOk:
            savefig = self.savefig_val.get()
            savestate = self.savestate_val.get()
            if self.record:
                self.mainpg.startRecord(name, savefig=savefig, savestate=savestate, interval=interval)
            else:
                if savefig:
                    self.mainpg.savefig(name)
                if savestate:
                    self.eq.saveState(self.k_sol, name)
            self.die()
        else:
            if alreadyExists:
                self.statustext.set('Name already exists. Are you sure you want to overwrite it?')
                self.warned = True
            if not isIntervalOk:
                self.statustext.set('Please enter a valid number (integer) for the interval.')

    def record_figure(self):
        pass

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
    def __init__(self, fig, ax, mainpg):
        self.fig = fig
        self.ax = ax
        self.mainpg = mainpg
        self.profile_inspector = self.mainpg.inspector.profile
        self.controller = self.mainpg.controller
        self.eq = self.controller.eq
        self.kind = PROFILE
        self.clicked = False

        self.profile_inspector.initProfilePlot(self)

        self.Fields = self.eq.getInitCondFields()
        self.auxFields = self.eq.getAuxFields(*self.Fields)
        self.markers = self.eq.getMarkers(*self.Fields, *self.auxFields)

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
        for k in range(len(self.markers)):
            i = k + self.eq.n_fields + self.eq.n_aux_fields
            _, (ym, yM) = self.get_ax_lims()
            self.lines[i].set_xdata([self.markers[k]]*2)
            self.lines[i].set_ydata([ym, yM])


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
        for m in self.markers:
            _, ylims = self.get_ax_lims()
            line,  = self.ax.plot([m,m], [ylims[0], ylims[1]], linestyle='dashed', color='tab:red')
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

    def getAnimIterables(self):
        return self.lines

    def update(self):
        self.Fields = self.eq.getCurrentFields()
        self.auxFields = self.eq.getCurrentAuxFields()
        self.markers = self.eq.getCurrentMarkers()
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

    def on_click(self, event):
        if event.inaxes == self.ax: # call p plot
            self.mainpg.setActivePlot(self)
            if self.mainpg.isEditing:
                self.clicked = True
                self.updateY(event.xdata, event.ydata)

    def off_click(self, event):
        if self.clicked:
            self.clicked = False
            self.released = True
            self.mainpg.isEditing = False
            self.saveEditandTick()

    def move_click(self, event):
        if self.clicked:
            if event.inaxes == self.ax:
                self.updateY(event.xdata, event.ydata)

class SpatioTemporalPlot:

    def __init__(self, fig, ax, mainpg):
        self.fig = fig
        self.ax = ax
        self.mainpg = mainpg
        self.spatiotemp_inspector = mainpg.inspector.spatiotemp
        self.controller = mainpg.controller
        self.eq = self.controller.eq
        self.k_st = 0 # k spatiotemp
        self.st_rows = ST_ROWS
        self.kind = SPATIOTEMPORAL

        self.Fields = self.eq.getInitCondFields()
        self.auxFields = self.eq.getAuxFields(*self.Fields)

        self.auto_lim = True
        self.active_field_indx = 0
        self.current_cmap = [0] * (self.eq.n_fields + self.eq.n_aux_fields)

        self.update_st_func = self.update_st_default

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
                new_rows = int(_max) - self.st_rows
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

    def on_click(self, event):
        if event.inaxes == self.ax:
            self.mainpg.setActivePlot(self)

    def off_click(self, event):
        pass

    def move_click(self, event):
        pass

    def getAnimIterables(self):
        return [self.im]

    def init_default(self):
        self.reset()
        v = self.eq.getCurrentParams()
        self.spatiotemp_inspector.Tentry.configure(state='disabled')
        self.spatiotemp_inspector.Tvar.set(v['dt'])
        self.update_st_func = self.update_st_default

    def init_strobe(self):
        self.reset()
        v = self.eq.getCurrentParams()
        self.k_strobe = 0
        if 'omega' in v:
            T_strobe = 2 * np.pi / v['omega']
        else:
            T_strobe = 20
        self.spatiotemp_inspector.Tentry.configure(state='enabled')
        self.spatiotemp_inspector.Tvar.set(T_strobe)
        self.N_strobe = int(round(T_strobe / v['dt']))
        self.update_st_func = self.update_strobe

    def update_strobe(self):
        v = self.eq.getCurrentParams()
        T_strobe = self.spatiotemp_inspector.get_T_update()
        if T_strobe is None: 
            N_strobe = self.N_strobe
        else:
            N_strobe = int(round(T_strobe / v['dt']))
        if N_strobe != self.N_strobe:
            self.k_strobe = 0
            self.N_strobe = N_strobe
            self.reset()
        if self.k_strobe == 0:
            self.k_strobe += 1
            return True
        self.k_strobe = (self.k_strobe + 1) % self.N_strobe
        return False

    def init_poincare(self):
        self.reset()
        last_field = self.getActiveField()
        cent, left, right = self.eq.getCurrentMarkers(indices=True) # centroid, left, right
        print(left, right)
        self.last_avg = np.mean(np.append(last_field[:left], last_field[right:]))
        self.update_st_func = self.update_poincare

    def update_poincare(self):
        new_field = self.getActiveField()
        cent, left, right = self.eq.getCurrentMarkers(indices=True)
        new_avg = np.mean(np.append(new_field[:left], new_field[right:]))

        if new_avg * self.last_avg < 0:
            self.last_avg = new_avg
            return True

        self.last_avg = new_avg
        return False

    def update_st_default(self):
        return True

    def update(self):
        if self.update_st_func():
            self.imvals[self.k_st, :] = self.getActiveField()
            self.im.set_data(self.imvals)
            self.k_st = (self.k_st + 1) % self.st_rows

        if self.auto_lim:
            self.auto_update()
        else:
            self.manual_update()

        return self.im


class PlotWindow(tk.Frame):

    def __init__(self, parent, mainpg, kind):
        self.parent = parent
        self.controller = mainpg.controller
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

        self.createPlot(kind)
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

    def createPlot(self, kind):
        if kind == TEMPORALPLOT:
            PlotType = TemporalPlot
        elif kind == SPATIOTEMPORAL:
            PlotType = SpatioTemporalPlot
        elif kind == PHASEPLOT:
            PlotType = PhasePlot
        elif kind == PROFILE:
            PlotType = ProfilePlot
        self.plot = PlotType(self.fig, self.ax, self.mainpg)
        if kind == PHASEPLOT:
            self.plot.parent = self

    def start_anim(self):
        event_source = self.mainpg.ani.event_source
        self.ani = animation.FuncAnimation(self.fig, self.animate, frames=60, interval=100, blit=False, event_source=event_source)

    def on_click(self, event):
        self.plot.on_click(event)

    def off_click(self, event):
        self.plot.off_click(event)

    def move_click(self, event):
        self.plot.move_click(event)

    def animate(self, i):
        if self.mainpg.isPaused or self.mainpg.isEditing:
            return self.plot.getAnimIterables()
        self.plot.update()
        return self.plot.getAnimIterables()

class TemporalPlot:
    def __init__(self, fig, ax, mainpg):
        self.fig = fig
        self.ax = ax
        self.temp_inspector = mainpg.inspector.temporal
        self.controller = mainpg.controller
        self.mainpg = mainpg
        self.eq = self.controller.eq

        self.kind = TEMPORALPLOT

        self.k_max = 4*SOLVE_EVERY_TI
        self.lines = []
        self.last_dt = self.eq.getParam('dt')
        self.ts = np.arange(self.k_max) * self.last_dt
        self.ys = np.zeros((self.eq.n_fields, len(self.ts)))
        self.current_colors = [_i + 2 for _i in range(self.eq.n_fields)]
        self.auto_lim = True

        self.unit = 0

        self.k_tp = 0

    def update_kmax(self, new_k_max):
        if new_k_max == self.k_max or new_k_max == 0: return
        if new_k_max > self.k_max:
            diff = new_k_max - self.k_max
            self.ys = np.append(self.ys, np.zeros((self.eq.n_fields, diff)), axis=1)
        else:
            self.ys = self.ys[:, :new_k_max]
            self.k_tp = self.k_tp % new_k_max
        self.k_max = new_k_max
        self.ts = np.arange(self.k_max) * self.eq.getParam('dt')

    def updateT(self):
        dt = self.eq.getParam('dt')
        if dt == self.last_dt: return
        self.ts = np.arange(self.k_max) * dt
        self.last_dt = dt

    def plot_fields(self):
        self.ax.clear()
        self.lines = []
        k_sol = self.k_tp
        for i in range(self.ys.shape[0]):
            color = COLORS[self.current_colors[i]]
            if color == 'none': 
                line = None
            else:
                line, = self.ax.plot(self.ts[:k_sol], self.ys[i, :k_sol], color='tab:' + color)
            self.lines.append(line)

    def getColorField(self, i):
        return self.current_colors[i]
    
    def setColorField(self, color, i):
        self.current_colors[i] = color
        self.plot_fields()

    def update_fields(self):
        k_sol = self.k_tp
        for i in range(self.ys.shape[0]):
            if self.lines[i] is not None:
                self.lines[i].set_xdata(self.ts[:k_sol])
                self.lines[i].set_ydata(self.ys[i, :k_sol])

    def auto_update_lim(self):
        xmin, xmax = self.ts[0], self.ts[-1]
        ymin, ymax = d10p(self.ys.min()), i10p(self.ys.max())
        self.set_lims((xmin, xmax), (ymin, ymax))

        if self.isInspectorActive():
            self.temp_inspector.set_xlim(xmin, xmax)
            self.temp_inspector.set_ylim(ymin, ymax)
        # call inspector

    def get_ax_lims(self):
        return self.ax.get_xlim(), self.ax.get_ylim()

    def isInspectorActive(self):
        return self.temp_inspector.tplot == self

    def set_lims(self, xlims, ylims):
        xmin, xmax = xlims
        ymin, ymax = ylims

        if xmin is not None:
            self.ax.set_xlim(xmin, xmax)
        
        if ymin is not None:
            self.ax.set_ylim(ymin, ymax)

    def manual_update_lim(self):
        if self.isInspectorActive():
            self.update_kmax(self.temp_inspector.get_kmax())
            self.set_lims(self.temp_inspector.get_xlims(), self.temp_inspector.get_ylims())

    def changeUnit(self, new_unit):
        if new_unit == self.unit: return
        self.unit = new_unit
        self.k_tp = 0

    def update(self):
        self.updateT()
        self.changeUnit(self.temp_inspector.getUnit())        
        Fields = self.eq.getCurrentFields()
        if self.k_tp == 0:
            for i in range(len(Fields)):
                self.ys[i, 0] = Fields[i][self.unit] # change this later
            self.k_tp = (self.k_tp + 1) % len(self.ts)
            self.plot_fields()
            return self.lines
        
        for i in range(len(Fields)):
            self.ys[i, self.k_tp] = Fields[i][self.unit]
        
        self.k_tp = (self.k_tp + 1) % len(self.ts)
        if self.auto_lim:
            self.auto_update_lim()
        else:
            self.manual_update_lim()
        self.update_fields()
        return self.lines
    
    def on_click(self, event):
        if event.inaxes == self.ax:
            self.mainpg.setActivePlot(self)

    def move_click(self, event):
        pass

    def off_click(self, event):
        pass

    def getAnimIterables(self):
        return self.lines

class Phase3DPlot:
    def __init__(self, fig, ax, mainpg, parent=None):
        self.fig = fig
        self.ax = ax
        self.mainpg = mainpg
        self.phase_inspector = mainpg.inspector.phase
        self.controller = mainpg.controller
        self.eq = self.controller.eq
        self.clicked = False
        if parent is not None:
            self.parent = parent

        self.x_field = 0
        self.y_field = 1
        self.z_field = 2

        self.kind = PHASEPLOT

        self.k_reset = 8*SOLVE_EVERY_TI

        self.resetXsandYs()
        self.lines = []
        self.pts = []

        self.auto_lim = True

        self.k_pp = 0
        self.set_field_label()

    def set_field_label(self):
        xlabel = self.eq.fieldNames[self.x_field]
        ylabel = self.eq.fieldNames[self.y_field]
        zlabel = self.eq.fieldNames[self.z_field]
        self.set_label(xlabel, ylabel, zlabel)

    def changeShownFields(self, new_x_field=None, new_y_field=None, new_z_field=None):
        if (new_x_field is None) and (new_y_field is None) and (new_z_field is None): return
        if new_x_field is not None:
            self.x_field = new_x_field
        if new_y_field is not None:
            self.y_field = new_y_field
        if new_z_field is not None:
            self.z_field = new_z_field
        self.set_field_label()
        self.resetXsandYs()
        self.clear()

    def getSelectedUnit(self):
        return self.phase_inspector.parent.temporal.getUnit()

    def set_label(self, xlabel, ylabel, zlabel):
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)
        self.ax.set_zlabel(zlabel)

    def move_click(self, event):
        pass

    def off_click(self, event):
        pass

    def on_click(self, event):
        if event.inaxes != self.ax: return
        self.mainpg.setActivePlot(self)

    def plot_fields(self):
        if len(self.xs) == 0: return
        self.ax.clear()
        self.lines = []
        self.pts = []
        for i in range(self.xs.shape[0]):
            line, = self.ax.plot(self.xs[i, :self.k_pp], self.ys[i, :self.k_pp], self.zs[i, :self.k_pp])
            self.lines.append(line)
        unit = self.getSelectedUnit()
        line, = self.ax.plot([self.xs[unit, self.k_pp-1]], [self.ys[unit, self.k_pp-1]], [self.zs[unit, self.k_pp-1]], 'ko')
        self.pts.append(line)

    def update_fields(self):
        if len(self.lines) == 0: return
        for i in range(self.xs.shape[0]):
            self.lines[i].set_data_3d(self.xs[i, :self.k_pp],
                                    self.ys[i, :self.k_pp],
                                    self.zs[i, :self.k_pp])

        unit = self.getSelectedUnit()

        self.pts[0].set_data_3d([self.xs[unit, self.k_pp-1]], 
                                [self.ys[unit, self.k_pp-1]],
                                [self.zs[unit, self.k_pp-1]])


    def resetXsandYs(self):
        self.xs = np.zeros((self.eq.getN(), self.k_reset))
        self.ys = np.zeros((self.eq.getN(), self.k_reset))
        self.zs = np.zeros((self.eq.getN(), self.k_reset))

    def get_ax_lims(self):
        return self.ax.get_xlim(), self.ax.get_ylim(), self.ax.get_zlim()

    def isInspectorActive(self):
        return self.phase_inspector.pplot == self

    def set_lims(self, xlims, ylims, zlims):
        xmin, xmax = xlims
        ymin, ymax = ylims
        zmin, zmax = zlims

        if xmin is not None:
            self.ax.set_xlim(xmin, xmax)
        
        if ymin is not None:
            self.ax.set_ylim(ymin, ymax)

        if zmin is not None:
            self.ax.set_zlim(zmin, zmax)

    def manual_update_lim(self):
        if self.isInspectorActive():
            self.set_lims(self.phase_inspector.get_xlims(), self.phase_inspector.get_ylims(), self.phase_inspector.get_zlims())
   
    def auto_update_lim(self):
        if len(self.xs) == 0: return
        xmin, xmax = self.xs.min(), self.xs.max()
        ymin, ymax = self.ys.min(), self.ys.max()
        zmin, zmax = self.zs.min(), self.zs.max()

        self.ax.set_xlim(d10p(xmin), i10p(xmax))
        self.ax.set_ylim(d10p(ymin), i10p(ymax))
        self.ax.set_zlim(d10p(zmin), i10p(zmax))

    def update(self):
        Fields = self.eq.getCurrentFields()
        # print(len(Fields), len(Fields[0]), self.eq.getN())
        if len(Fields[self.x_field]) != self.xs.shape[0]:
            self.clear() # resets Xs and Ys
        self.xs[:, self.k_pp] = Fields[self.x_field]
        self.ys[:, self.k_pp] = Fields[self.y_field]
        self.zs[:, self.k_pp] = Fields[self.z_field]
        self.k_pp = (self.k_pp + 1) % self.k_reset
        if self.k_pp == 1:
            self.plot_fields()
        else:
            self.update_fields()
        if self.auto_lim:
            self.auto_update_lim()
        else:
            self.manual_update_lim()
        # update lines
        return self.lines + self.pts

    def clear(self):
        self.resetXsandYs()
        self.k_pp = 0

    def getAnimIterables(self):
        return self.lines + self.pts


class PhasePlot:
    def __init__(self, fig, ax, mainpg, parent=None):
        self.fig = fig
        self.ax = ax
        self.mainpg = mainpg
        self.phase_inspector = mainpg.inspector.phase
        self.controller = mainpg.controller
        self.eq = self.controller.eq
        self.clicked = False

        if parent is not None:
            self.parent=parent

        self.x_field = 0
        self.y_field = 1

        self.kind = PHASEPLOT

        self.k_reset = 8*SOLVE_EVERY_TI

        self.resetXsandYs()
        self.lines = []
        self.pts = []

        self.auto_lim = True

        self.k_pp = 0
        self.set_field_label()

    def set_field_label(self):
        xlabel = self.eq.fieldNames[self.x_field]
        ylabel = self.eq.fieldNames[self.y_field]
        self.set_label(xlabel, ylabel)

    def changeShownFields(self, new_x_field=None, new_y_field=None):
        if (new_x_field is None) and (new_y_field is None): return
        if new_x_field is not None:
            self.x_field = new_x_field
        if new_y_field is not None:
            self.y_field = new_y_field
        self.set_field_label()
        self.resetXsandYs()
        self.clear()

    def getSelectedUnit(self):
        return self.phase_inspector.parent.temporal.getUnit()

    def set_label(self, xlabel, ylabel):
        self.ax.set_xlabel(xlabel)
        self.ax.set_ylabel(ylabel)

    def draw_user_line(self):
        # interpolate
        if len(self.xdata) != 0:
            self.user_line, = self.ax.plot(self.xdata, self.ydata, 'k-')

    def update_user_line(self):
        self.user_line.set_xdata(self.xdata)
        self.user_line.set_ydata(self.ydata)

    def on_click(self, event):
        if event.inaxes != self.ax: return
        self.mainpg.setActivePlot(self)
        if self.mainpg.isEditing:
            self.xdata = [event.xdata]
            self.ydata = [event.ydata]
            self.draw_user_line()
            self.clicked = True

    def move_click(self, event):
        if self.clicked and event.inaxes == self.ax:
            self.xdata.append(event.xdata)
            self.ydata.append(event.ydata)
            self.update_user_line()

    def off_click(self, event):
        if self.clicked:
            self.clicked = False
            self.mainpg.isEditing = False
            self.saveEditandTick()
    
    def data_to_initCond(self, xdata, ydata):
        initCond = []
        for i in range(self.eq.n_fields):
            if i == self.x_field:
                field = xdata
            elif i == self.y_field:
                field = ydata
            else:
                field = np.zeros(len(xdata))
            initCond.append(field)
        return initCond

    def saveEditandTick(self):
        newInitCond = self.data_to_initCond(self.xdata, self.ydata)
        Ni = len(self.xdata)
        self.eq.setNi(Ni)
        self.eq.tick(newInitCondFields=newInitCond)
        self.clear() # resets Xs and Ys
        self.phase_inspector.parent.temporal.updateN()
        self.k_pp = 0

    def plot_fields(self):
        if len(self.xs) == 0: return
        self.ax.clear()
        self.lines = []
        self.pts = []
        for i in range(self.xs.shape[0]):
            line, = self.ax.plot(self.xs[i, :self.k_pp], self.ys[i, :self.k_pp])
            self.lines.append(line)
        unit = self.getSelectedUnit()
        line, = self.ax.plot([self.xs[unit, self.k_pp-1]], [self.ys[unit, self.k_pp-1]], 'ko')
        self.pts.append(line)

    def update_fields(self):
        if len(self.lines) == 0: return
        for i in range(self.xs.shape[0]):
            self.lines[i].set_xdata(self.xs[i, :self.k_pp])
            self.lines[i].set_ydata(self.ys[i, :self.k_pp])
        unit = self.getSelectedUnit()
        # print(unit, self.k_pp, self.xs.shape)
        self.pts[0].set_xdata([self.xs[unit, self.k_pp-1]])
        self.pts[0].set_ydata([self.ys[unit, self.k_pp-1]])

    def resetXsandYs(self):
        self.xs = np.zeros((self.eq.getN(), self.k_reset))
        self.ys = np.zeros((self.eq.getN(), self.k_reset))

    def get_ax_lims(self):
        return self.ax.get_xlim(), self.ax.get_ylim()

    def isInspectorActive(self):
        return self.phase_inspector.pplot == self

    def set_lims(self, xlims, ylims):
        xmin, xmax = xlims
        ymin, ymax = ylims

        if xmin is not None:
            self.ax.set_xlim(xmin, xmax)
        
        if ymin is not None:
            self.ax.set_ylim(ymin, ymax)

    def manual_update_lim(self):
        if self.isInspectorActive():
            self.set_lims(self.phase_inspector.get_xlims(), self.phase_inspector.get_ylims())
   
    def auto_update_lim(self):
        if len(self.xs) == 0: return
        xmin, xmax = self.xs.min(), self.xs.max()
        ymin, ymax = self.ys.min(), self.ys.max()
        self.ax.set_xlim(d10p(xmin), i10p(xmax))
        self.ax.set_ylim(d10p(ymin), i10p(ymax))

    def update(self):
        Fields = self.eq.getCurrentFields()
        # print(len(Fields), len(Fields[0]), self.eq.getN())
        if len(Fields[self.x_field]) != self.xs.shape[0]:
            self.clear() # resets Xs and Ys
        self.xs[:, self.k_pp] = Fields[self.x_field]
        self.ys[:, self.k_pp] = Fields[self.y_field]
        self.k_pp = (self.k_pp + 1) % self.k_reset
        if self.k_pp == 1:
            self.plot_fields()
        else:
            self.update_fields()
        if self.auto_lim:
            self.auto_update_lim()
        else:
            self.manual_update_lim()
        # update lines
        return self.lines + self.pts

    def clear(self):
        self.resetXsandYs()
        self.k_pp = 0

    def getAnimIterables(self):
        return self.lines + self.pts


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

        label = ttk.Label(self, text='Interactive Differential Equations Analysis Simulator', font=LARGE_FONT)
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

    def resetParameterWindow(self):
        self.killParamWindow()
        self.createParamWindow(self.controller.eq)

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

    def setDim(self, dim):
        self.dim = dim

    def activate(self):

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

        if self.controller.eq.dim == 1:
            self.plots = [ProfilePlot(self.fig, self.ax, self), \
                SpatioTemporalPlot(self.fig, self.ax2, self)]
            self.eqX = eq.x

        elif self.controller.eq.dim == 0:
            self.plots = [PhasePlot(self.fig, self.ax, self, parent=self), \
                TemporalPlot(self.fig, self.ax2, self)]
            self.inspector.phase.initPhasePlot(self.plots[0])
            self.inspector.temporal.initTemporalPlot(self.plots[1])
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
        self.active = True
        self.anim_stopped = False

        self.t = eq.t0

        self.parentNewWindows = []
        self.newWindows = []
        ### Start animation

        # save event_source
        self.ani = animation.FuncAnimation(self.fig, self.animate, frames=60, interval=10, blit=False)

    def set_init_cond_zero(self):
        self.controller.setEqInitCondZero()

    def add_window(self):
        parent_new_window = tk.Toplevel(self)
        parent_new_window.geometry('512x512')

        self.parentNewWindows.append(parent_new_window)

        kind = self.controller.activePlot[1]
        pw = PlotWindow(parent_new_window, self, kind)

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

    def startRecord(self, recordname, interval=1, savefig=False, savestate=True):
        self.recording = True
        self.recording_fig = savefig
        self.recording_state = savestate
        self.recording_interval = interval

        callback = self.recordfig if savefig else None
        self.controller.eq.startRecording(recordname, interval=interval, callback=callback, savestate=savestate)

        self.i_record = 0
        self.img_recording = True
        self.btn_record.configure(image=self.recordingimg)

    def recordfig(self, name):
        self.savefig(name)

    def savefig(self, name):
        # figure folder = fig/dim/eqname/
        path = self.controller.eq.getFigureFolder() + name
        print(f'Saved figure in {path}')
        self.fig.savefig(path)

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
        for plot in self.plots:
            if event.inaxes == plot.ax:
                plot.on_click(event)

    def setActivePlot(self, aplot):
        self.controller.activePlot = (aplot.ax, aplot.kind, aplot)
        self.inspector.show(aplot.kind)

    def off_click(self, event):
        for plot in self.plots:
            plot.off_click(event)

    def move_click(self, event):
        for plot in self.plots:
            if event.inaxes == plot.ax:
                plot.move_click(event)

    def stop_animation(self):
        self.anim_stopped = True
        self.ani.event_source.stop()

    def createParamWindow(self, eq):
        # create new window
        self.parentParamWindow = tk.Toplevel(self)
        self.paramWindow = ParameterWindow(self.parentParamWindow, '400x500', eq)

    def createImageWindow(self, eq):
        self.parentImageWindow = tk.Toplevel(self)
        filename = eq.getDimFolder() + eq.name + '.jpg'
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
        iterables = []
        for plot in self.plots:
            iterables += plot.getAnimIterables()
        return iterables

    def animate(self, i):
        if not self.active or self.anim_stopped: return
        eq = self.controller.eq
        eq.updateX()

        if self.isPaused or self.isEditing:
            return self.getAnimIterables()

        eq.tick()
        self.Fields = eq.getCurrentFields()
        if self.controller.eq.dim == 1: 
            self.inspector.profile.setTime(self.t) # PPlot + equation i guesss
        elif self.controller.eq.dim == 0:
            self.inspector.phase.setTime(self.t)

        for plot in self.plots:
            plot.update()

        if self.recording:
            self.animate_record()

        self.t += eq.getParam('dt') # PPlot
        return self.getAnimIterables()


