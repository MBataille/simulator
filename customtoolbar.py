import numpy as np
import tkinter as tk
import matplotlib as mpl
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)

class CustomToolbar(NavigationToolbar2Tk):

    def __init__(self, canvas_, parent_, main_window):
        self.parent_frame = parent_
        self.main_window = main_window
        self.toolitems = (
            ('Home', 'Home view', 'home', 'home'),
            ('Back', 'Previous view', 'back', 'back'),
            ('Forward', 'Next view', 'forward', 'forward'),
            (None, None, None, None),
            ('Pan', 'Activate pan', 'move', 'pan'),
            ('Zoom', 'Activate zoom', 'zoom_to_rect', 'zoom'),
            ('Save', 'Save figure', 'filesave', 'save_figure'),
            (None, None, None, None),
            ('Edit', 'Edit state', 'sim-edit', 'edit'),
            ('Play', 'Continue simulation', 'sim-play', 'play'),
            ('Pause', 'Pause simulation', 'sim-pause', 'pause'),
            (None, None, None, None),
            ('Params', 'Open parameter window', 'subplots', 'open_params'),
            ('Equation', 'Open equation window', 'sim-equal', 'open_equation'),
            ('Inspector', 'Open inspector window', 'sim-info', 'open_inspector'),
            )
        NavigationToolbar2Tk.__init__(self,canvas_,parent_)

    def edit(self):
        self.main_window.edit()

    def play(self):
        self.main_window.play()

    def pause(self):
        self.main_window.pause()

    def open_params(self):
        self.main_window.open_params()

    def open_equation(self):
        self.main_window.open_equation()

    def open_inspector(self):
        self.main_window.open_inspector()
