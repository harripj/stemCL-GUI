# This Python file uses the following encoding: utf-8
from PySide2 import QtWidgets, QtGui
from matplotlib.backends.backend_qt5cairo import FigureCanvas
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class MPLWidget(FigureCanvas):
    '''Plotting widget that can be embedded in a PySide GUI.'''
    def __init__(self, figure=None, projection=None):
        if figure is None:
            figure = plt.Figure(tight_layout=True)
        super().__init__(figure)

        self.axes = self.figure.add_subplot(111, projection=projection)

        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
    # end class
