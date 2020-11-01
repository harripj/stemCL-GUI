# This Python file uses the following encoding: utf-8

import sys
import os

from PySide2.QtWidgets import (
    QApplication,
    QMainWindow,
    QLineEdit,
    QPushButton,
    QTableWidget,
    QWidget,
    QTableWidgetItem,
    QFileDialog,
    QBoxLayout,
    QComboBox,
    QSizePolicy,
    QTreeView,
    QCheckBox,
    QTextEdit,
    QProgressBar,
)
from PySide2.QtCore import QFile, QCoreApplication, Qt, QObject, SIGNAL
from PySide2.QtUiTools import QUiLoader

# from PySide2.QtCharts import QtCharts
from PySide2.QtGui import QVector3D, QStandardItemModel, QStandardItem
from PySide2.QtDataVisualization import QtDataVisualization

import math
import os
import pandas as pd
import numpy as np
from matplotlib.patches import Rectangle, Circle
from matplotlib.transforms import Bbox
from ase import Atoms
from ase.utils import rotate
from ase.data import atomic_names, covalent_radii
from ase.data.colors import jmol_colors

from .functions import stemcl_format_xyz

from matplotlib.backends.backend_qt5cairo import FigureCanvas
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from datetime import datetime

import pyopencl
import subprocess


class MPLWidget(FigureCanvas):
    """Plotting widget that can be embedded in a PySide GUI."""

    def __init__(self, figure=None, projection=None):
        if figure is None:
            figure = plt.Figure(tight_layout=True)
        super().__init__(figure)

        self.axes = self.figure.add_subplot(111, projection=projection)

        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

    # end class


class stemCL(QMainWindow):
    def __init__(self):
        super(stemCL, self).__init__()
        self.load_ui()
        self.setup_openCL()

        # setup class attr with no loaded data
        self.file_name_xyz = ""
        self.atoms = None
        self.data = None

    def setup_openCL(self):
        self.platforms = pyopencl.get_platforms()
        self.devices = [platform.get_devices() for platform in self.platforms]

        # populate treeview
        self.treeView_model = QStandardItemModel()
        self.treeView_model.setHorizontalHeaderLabels(["Platform"])  # , 'Device'])
        self.treeView_openCL.setModel(self.treeView_model)
        self.treeView_openCL.setUniformRowHeights(True)

        # connect command
        # self.treeView_model.itemChanged.connect(self.create_run_command)

        for i, platform in enumerate(self.platforms):
            parent = QStandardItem(platform.name)
            for d in self.devices[i]:
                parent.appendRow(QStandardItem(d.name))
            self.treeView_model.appendRow(parent)
            # span container columns
            self.treeView_openCL.setFirstColumnSpanned(
                i, self.treeView_openCL.rootIndex(), True
            )

        self.treeView_openCL.expandAll()

    def create_run_command(self):
        self.command = []

        if self.checkBox_caffeinate.checkState():
            self.command.append("caffeinate")

        self.command.append("stemcl")

        if not self.treeView_openCL.selectedIndexes():
            return
        else:
            indexes = self.treeView_openCL.selectedIndexes()[
                0
            ]  # only care about first selection (if somehow more than one)

        self.command.append(self.platforms[indexes.column()].name)
        self.command.append(str(indexes.row()))

        #        self.command.append('.')

        if not self.file_name_xyz:
            return
        else:
            # escape spaces in path
            self.command.append(
                # "{}".format(os.path.dirname(self.file_name_xyz).replace(" ", "\\ "))
                os.path.dirname(self.file_name_xyz).replace(" ", "\\ ")
            )

        self.lineEdit_run_command.setText(" ".join([str(i) for i in self.command]))

    def run_command(self):
        self.textEdit_console.clear()
        self.textEdit_console.insertHtml(str(datetime.now()))

        if self.checkBox_remove_progress_pgm.checkState():
            # remove progress.dat in directory so simulation runs again from scratch
            fname = os.path.join(os.path.dirname(self.file_name_xyz), "progress.pgm")
            if os.path.exists(fname):
                os.remove(fname)

        # print(" ".join(self.command))
        out = subprocess.Popen(
            " ".join(self.command),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            shell=True,
            universal_newlines=True,
        )

        lines = []
        while True:
            line = out.stdout.readline()
            if not line:
                break
            lines.append(line)

            # print to terminal
            print(line, end="")
            try:
                index = int(line.strip().index("%"))
                self.progressBar.setValue(int(line.strip()[index - 3 : index]))
            except ValueError:
                # simulation preparing, not running so % will not be found in string
                pass

            self.textEdit_console.append(line.strip())

        # self.textEdit_console.toMarkdown()

    def load_ui(self):
        loader = QUiLoader()
        path = os.path.join(os.path.dirname(__file__), "form.ui")
        ui_file = QFile(path)
        ui_file.open(QFile.ReadOnly)
        self.ui = loader.load(ui_file, self)
        ui_file.close()

        #        self.chartView = QtCharts.QChartView()
        #        self.chartView.setMinimumWidth(450)
        #        self.graphicsView = self.ui.findChild(QGraphicsView, "graphicsView")

        # graphics
        self.tab_summary = self.ui.findChild(QWidget, "tab_summary")

        #        QtDataVisualization.QAbstract3DGraph

        self.scatter3d = QtDataVisualization.Q3DScatter()
        self.container = QWidget.createWindowContainer(self.scatter3d)
        self.container.setMinimumWidth(450)

        self.scatterProxy = QtDataVisualization.QScatterDataProxy()
        self.scatterSeries = QtDataVisualization.QScatter3DSeries(self.scatterProxy)

        self.scatter3d.axisX().setTitle("x")
        self.scatter3d.axisY().setTitle("y")
        self.scatter3d.axisZ().setTitle("z")

        self.scatter3d.addSeries(self.scatterSeries)

        self.find_widgets()

    def find_widgets(self):
        # buttons
        self.pushButton_create_dat = self.ui.findChild(
            QPushButton, "pushButton_create_dat"
        )
        self.pushButton_create_dat.clicked.connect(self.onclick_create_dat)

        self.pushButton_load_xyz = self.ui.findChild(QPushButton, "pushButton_load_xyz")
        self.pushButton_load_xyz.clicked.connect(self.onclick_load_xyz)

        self.pushButton_view_3d = self.ui.findChild(QPushButton, "pushButton_view3d")
        self.pushButton_view_3d.clicked.connect(self.onclick_view_3d)

        self.pushButton_save_fig = self.ui.findChild(QPushButton, "pushButton_save_fig")
        self.pushButton_save_fig.clicked.connect(self.onclick_save_fig)

        self.pushButton_scherzer_defocus = self.ui.findChild(
            QPushButton, "pushButton_scherzer_defocus"
        )

        self.pushButton_update_scan_slices = self.ui.findChild(
            QPushButton, "pushButton_update_scan_slices"
        )
        self.pushButton_update_scan_slices.clicked.connect(self.plot_slices)

        # lineedits
        # tab1
        self.lineEdit_beam_voltage = self.ui.findChild(
            QLineEdit, "lineEdit_beam_voltage"
        )
        self.lineEdit_energy_spread = self.ui.findChild(
            QLineEdit, "lineEdit_energy_spread"
        )
        self.lineEdit_cs3 = self.ui.findChild(QLineEdit, "lineEdit_cs3")
        self.lineEdit_cs5 = self.ui.findChild(QLineEdit, "lineEdit_cs5")
        self.lineEdit_cc = self.ui.findChild(QLineEdit, "lineEdit_cc")
        self.lineEdit_defocus = self.ui.findChild(QLineEdit, "lineEdit_defocus")
        self.lineEdit_aperture_inner = self.ui.findChild(
            QLineEdit, "lineEdit_aperture_inner"
        )
        self.lineEdit_aperture_outer = self.ui.findChild(
            QLineEdit, "lineEdit_aperture_outer"
        )
        self.lineEdit_astig2_mag = self.ui.findChild(QLineEdit, "lineEdit_astig2_mag")
        self.lineEdit_astig2_angle = self.ui.findChild(
            QLineEdit, "lineEdit_astig2_angle"
        )
        self.lineEdit_astig3_mag = self.ui.findChild(QLineEdit, "lineEdit_astig3_mag")
        self.lineEdit_astig3_angle = self.ui.findChild(
            QLineEdit, "lineEdit_astig3_angle"
        )

        # tab2
        self.lineEdit_scan_x_min = self.ui.findChild(QLineEdit, "lineEdit_scan_x_min")
        self.lineEdit_scan_x_max = self.ui.findChild(QLineEdit, "lineEdit_scan_x_max")
        self.lineEdit_scan_y_min = self.ui.findChild(QLineEdit, "lineEdit_scan_y_min")
        self.lineEdit_scan_y_max = self.ui.findChild(QLineEdit, "lineEdit_scan_y_max")
        self.lineEdit_slice_thickness = self.ui.findChild(
            QLineEdit, "lineEdit_slice_thickness"
        )
        self.lineEdit_tilt_x = self.ui.findChild(QLineEdit, "lineEdit_tilt_x")
        self.lineEdit_tilt_y = self.ui.findChild(QLineEdit, "lineEdit_tilt_y")
        self.lineEdit_probe_wf_x = self.ui.findChild(QLineEdit, "lineEdit_probe_wf_x")
        self.lineEdit_probe_wf_y = self.ui.findChild(QLineEdit, "lineEdit_probe_wf_y")
        self.lineEdit_transmission_func_x = self.ui.findChild(
            QLineEdit, "lineEdit_transmission_func_x"
        )
        self.lineEdit_transmission_func_y = self.ui.findChild(
            QLineEdit, "lineEdit_transmission_func_y"
        )
        self.lineEdit_output_size_x = self.ui.findChild(
            QLineEdit, "lineEdit_output_size_x"
        )
        self.lineEdit_output_size_y = self.ui.findChild(
            QLineEdit, "lineEdit_output_size_y"
        )
        self.lineEdit_diffraction_n = self.ui.findChild(
            QLineEdit, "lineEdit_diffraction_n"
        )

        # tab run
        self.treeView_openCL = self.ui.findChild(QTreeView, "treeView_openCL")

        self.checkBox_caffeinate = self.ui.findChild(QCheckBox, "checkBox_caffeinate")
        self.checkBox_remove_progress_pgm = self.ui.findChild(
            QCheckBox, "checkBox_remove_progress_pgm"
        )

        self.lineEdit_run_command = self.ui.findChild(QLineEdit, "lineEdit_run_command")

        self.pushButton_run_command = self.ui.findChild(
            QPushButton, "pushButton_run_command"
        )
        self.pushButton_run_command.clicked.connect(self.run_command)

        self.pushButton_generate_command = self.ui.findChild(
            QPushButton, "pushButton_generate_command"
        )
        self.pushButton_generate_command.clicked.connect(self.create_run_command)

        self.progressBar = self.ui.findChild(QProgressBar, "progressBar")

        self.textEdit_console = self.ui.findChild(QTextEdit, "textEdit_console")
        #        This functionality is slow and has been replaced by a button
        #        for l in (
        #            self.lineEdit_scan_x_min,
        #            self.lineEdit_scan_x_max,
        #            self.lineEdit_scan_y_min,
        #            self.lineEdit_scan_y_max,
        #            self.lineEdit_slice_thickness,
        #        ):
        #            l.textEdited.connect(self.plot_slices)

        # tab3
        self.lineEdit_GPU_PROBE = self.ui.findChild(QLineEdit, "lineEdit_GPU_PROBE")
        self.lineEdit_SMOOTH_PROBE = self.ui.findChild(
            QLineEdit, "lineEdit_SMOOTH_PROBE"
        )
        self.lineEdit_num_parallel = self.ui.findChild(
            QLineEdit, "lineEdit_num_parallel"
        )
        self.lineEdit_use_hdd = self.ui.findChild(QLineEdit, "lineEdit_use_hdd")

        # tables
        self.tableWidget_detectors = self.ui.findChild(
            QTableWidget, "tableWidget_detectors"
        )
        self.tableWidget_xyz_atoms = self.ui.findChild(
            QTableWidget, "tableWidget_xyz_atoms"
        )

        # comboBox
        self.comboBox_plt = self.ui.findChild(QComboBox, "comboBox_plt")
        self.comboBox_plt.addItems(("yz", "xz", "xy"))  # 3 axes
        self.comboBox_plt.currentIndexChanged.connect(self.plot_atoms)

        self.tab_summary = self.ui.findChild(QWidget, "tab_summary")

        self.mplWidget = MPLWidget(projection=None)
        self.tab_summary.layout().addWidget(self.mplWidget, 0, 1, 2, 2)

    def set_scherzer_defocus(self):
        # source: JEOL https://www.jeol.co.jp/en/words/emterms/search_result.html?keyword=wavelength%20of%20electron
        voltage = float(self.lineEdit_beam_voltage.text()) * 1e3  # kV
        wavelength = 1.23e3 / np.sqrt(voltage * (1 + 9.78e-7 * voltage))  # pm

        defocus = (self.lineEdit_cs3.text() * 1e-3 * wavelength * 1e-12) ** (1.0 / 2)
        self.lineEdit_defocus.setText(f"{defocus:.3f}")

    def get_data_from_table(self, col_start=0, dataframe=False):
        # get data from table
        data = []
        for row in range(self.tableWidget_xyz_atoms.rowCount()):
            # ignore symbols column
            data_row = []
            for col in range(col_start, self.tableWidget_xyz_atoms.columnCount()):
                data_row.append(float(self.tableWidget_xyz_atoms.item(row, col).text()))
            data.append(data_row)

        if dataframe:
            header = []
            for col in range(col_start, self.tableWidget_xyz_atoms.columnCount()):
                header.append(
                    self.tableWidget_xyz_atoms.horizontalHeaderItem(col).text()
                )

            data = pd.DataFrame(data, columns=header)

        return data

    def add_data_to_table(self, data):
        self.tableWidget_xyz_atoms.setRowCount(len(data))

        for i in range(len(data)):
            for j in range(self.tableWidget_xyz_atoms.columnCount()):
                tableItem = self.tableWidget_xyz_atoms.item(i, j)
                # no previous data
                if tableItem is None:
                    tableItem = QTableWidgetItem(str(data[i, j]))
                    self.tableWidget_xyz_atoms.setItem(i, j, tableItem)
                else:
                    self.tableWidget_xyz_atoms.item(i, j).setText(str(data[i, j]))

    def update_scatter(self):
        data = self.get_data_from_table(col_start=1)
        # reset series
        self.scatterSeries.dataProxy().removeItems(
            0, self.scatterSeries.dataProxy().itemCount()
        )
        # add to series
        self.scatterSeries.dataProxy().addItems(
            [QtDataVisualization.QScatterDataItem(QVector3D(*i)) for i in data]
        )

    def dataframe_as_atoms(self, df):
        return Atoms(
            symbols=self.data["Symbol"], positions=self.data[["x", "y", "z"]].values
        )
        # if symbol is str

    #        return Atoms(symbols=[atomic_names[i] for i in self.data['Symbol']],
    #                     positions=self.data[['x', 'y', 'z']].values)

    def onclick_load_xyz(self, button):
        # self.tr for translate -> avoids filedialog crash
        self.file_name_xyz, filter = QFileDialog.getOpenFileName(
            filter=self.tr("xyz (*.xyz)"),
            dir=self.tr("~"),
            caption=self.tr("Open .xyz file"),
        )

        # file was selected, returns empty string if cancelled
        if self.file_name_xyz:
            self.data = pd.read_csv(
                self.file_name_xyz,
                skiprows=4,
                delimiter=" ",
                names=("Symbol", "x", "y", "z", "Occupancy", "Debye-Waller"),
            )
            self.atoms = self.dataframe_as_atoms(self.data)

            self.add_data_to_table(self.data.values)

            # update scan area
            self.lineEdit_scan_x_min.setText("{:.1f}".format(self.data["x"].min()))
            self.lineEdit_scan_x_max.setText("{:.1f}".format(self.data["x"].max()))

            self.lineEdit_scan_y_min.setText("{:.1f}".format(self.data["y"].min()))
            self.lineEdit_scan_y_max.setText("{:.1f}".format(self.data["y"].max()))

            # initial plot
            self.plot_atoms()
            self.update_scatter()

    def calculate_plotting_offset(self, rotation, scale, maxwidth):
        # from ase.visualize.plot and ase.io.utils
        R = np.dot(self.atoms.get_positions(), rotate(rotation))
        radii = covalent_radii[self.atoms.numbers]

        X1 = (R - radii[:, None]).min(0)
        X2 = (R + radii[:, None]).max(0)
        if show_unit_cell == 2:
            X1 = np.minimum(X1, cell_vertices.min(0))
            X2 = np.maximum(X2, cell_vertices.max(0))
        M = (X1 + X2) / 2
        S = 1.05 * (X2 - X1)
        w = scale * S[0]
        if w > maxwidth:
            w = maxwidth
            scale = w / S[0]
        h = scale * S[1]
        offset = np.array([scale * M[0] - w / 2, scale * M[1] - h / 2, 0])

        return offset

    def plot_slices(self):
        # remove all slice patches
        for i in range(len(self.mplWidget.axes.patches))[::-1]:
            patch = self.mplWidget.axes.patches[i]
            if patch.get_label() == "Slice":
                patch.remove()

        axes = self.comboBox_plt.currentText()
        # for overlays
        alpha = 0.3
        colors = ("r", "b", "g")

        # add scan area to xy plot
        patches = []
        labels = []

        if axes == "xy":
            x_min = float(self.lineEdit_scan_x_min.text())
            x_max = float(self.lineEdit_scan_x_max.text())
            y_min = float(self.lineEdit_scan_y_min.text())
            y_max = float(self.lineEdit_scan_y_max.text())

            rectangle = Rectangle(
                (x_min, y_min),
                x_max - x_min,
                y_max - y_min,
                color="r",
                alpha=alpha,
                ec="k",
                label="Slice",
            )
            self.mplWidget.axes.add_patch(rectangle)

            patches.append(rectangle)
            labels.append("Scan Area")

            self.mplWidget.axes.legend(patches, labels)

        # plot multislice layers
        elif axes == "xz" or axes == "yz":
            positions = self.atoms.get_positions()

            if axes == "xz":
                x_min = float(self.lineEdit_scan_x_min.text())
                x_max = float(self.lineEdit_scan_x_max.text())
                z_min = positions[:, 2].min()
                z_max = positions[:, 2].max()

            if axes == "yz":
                x_min = float(self.lineEdit_scan_y_min.text())
                x_max = float(self.lineEdit_scan_y_max.text())
                z_min = positions[:, 2].min()
                z_max = positions[:, 2].max()

            # in Angstrom
            slice_thickness = float(self.lineEdit_slice_thickness.text())

            # from stemCL code
            nslices = math.ceil((z_max - z_min) / slice_thickness)

            patches = []
            labels = []
            for i in range(nslices):
                rectangle = Rectangle(
                    (x_min, z_min + slice_thickness * i),
                    x_max - x_min,
                    slice_thickness,
                    color=colors[i % len(colors)],
                    alpha=alpha,
                    ec="k",
                    label="Slice",
                )
                self.mplWidget.axes.add_patch(rectangle)

                self.mplWidget.axes.text(
                    self.mplWidget.axes.get_xlim()[0],
                    z_min + slice_thickness * i,
                    str(i),
                    size=12,
                    color=rectangle.get_facecolor(),
                    bbox=dict(boxstyle="round", color="w", alpha=0.5),
                    ha="left",
                    va="bottom",
                )

                patches.append(rectangle)
                labels.append(f"Slice {i}")

        self.mplWidget.draw()

    def plot_atoms(self):
        # no data set, don't run function
        if self.atoms is None:
            return

        self.mplWidget.axes.cla()

        axes = self.comboBox_plt.currentText()

        #        if axes == "xy":
        #            rotation = ('0x,0y,0z')
        #        elif axes == "xz":
        #            rotation = ('90x,0y,0z')
        #        elif axes == "yz":
        #            rotation = ('0x,90y,90z')

        #        # default from ase.io.utils
        #        scale = 1
        #        maxwidth = 500
        #        offset = self.calculate_plotting_offset(rotation, scale, maxwidth)

        radii = covalent_radii[self.atoms.numbers]
        colors = jmol_colors[self.atoms.numbers.clip(max=len(jmol_colors) - 1)]
        positions = self.atoms.get_positions()

        # apply rotation
        #        positions = np.dot(self.atoms.get_positions(), rotate(rotation))

        pad = 2  # angstrom

        if axes == "xy":
            # add atoms from lowest z up
            sort_index = 2
            x_index = 0
            y_index = 1
        elif axes == "xz":
            # add atoms from lowest y up
            sort_index = 1
            x_index = 0
            y_index = 2
        elif axes == "yz":
            # add atoms from lowest x up
            sort_index = 0
            x_index = 1
            y_index = 2
        else:
            print("Axes wrong? This should not happen.")

        # add atoms to plot
        for i in np.argsort(positions[:, sort_index]):
            self.mplWidget.axes.add_patch(
                Circle(
                    positions[i, [x_index, y_index]],
                    radii[i],
                    ec="k",
                    color=colors[i],
                    label="Atom",
                )
            )

        self.mplWidget.axes.set_xlim(
            (positions[:, x_index] - radii).min() - pad,
            (positions[:, x_index] + radii).max() + pad,
        )
        self.mplWidget.axes.set_ylim(
            (positions[:, y_index] - radii).min() - pad,
            (positions[:, y_index] + radii).max() + pad,
        )

        self.mplWidget.axes.set_aspect("equal")

        self.mplWidget.axes.set_xlabel(axes[0])
        self.mplWidget.axes.set_ylabel(axes[1])

        self.mplWidget.figure.tight_layout()
        self.plot_slices()

    def onclick_view_3d(self, button):
        self.popup_window = QWidget()
        self.popup_window.setMinimumWidth(450)
        self.popup_window.setMinimumHeight(450)
        self.popup_window.setLayout(QBoxLayout(QBoxLayout.TopToBottom))

        self.popup_window.layout().addWidget(self.container)
        self.popup_window.show()

        self.update_scatter()

    def onclick_create_dat(self, button):
        _dir = (
            os.path.join(os.path.dirname(self.file_name_xyz), "parameter.dat")
            if self.file_name_xyz
            else "parameter.dat"
        )
        out, _filter = QFileDialog.getSaveFileName(
            self,
            dir=self.tr(_dir),
            caption=self.tr("Save .dat file"),
            filter=self.tr("dat (*.dat)"),
        )

        if out:
            with open(out, "w") as f:

                f.write(
                    "# STEM probe parameters, V0(kv), Cs3(mm), Cs5(mm), df(Angstroms), apert1,2(mrad), Cc (mm), energy spread (eV)\n"
                )
                line1 = " ".join(
                    [
                        f"{float(self.lineEdit_beam_voltage.text()):.3f}",
                        f"{float(self.lineEdit_cs3.text()):.3f}",
                        f"{float(self.lineEdit_cs5.text()):.3f}",
                        f"{float(self.lineEdit_defocus.text()):.3f}",
                        f"{float(self.lineEdit_aperture_inner.text()):.3f}",
                        f"{float(self.lineEdit_aperture_outer.text()):.3f}",
                        f"{float(self.lineEdit_cc.text()):.3f}",
                        f"{float(self.lineEdit_energy_spread.text()):.3f}",
                    ]
                )
                f.write(f"{line1}\n")

                f.write(
                    "# Magnitude and angle of 2-fold astigmatism (in Ang. and degrees)\n"
                )
                f.write(
                    f"{float(self.lineEdit_astig2_mag.text()):.3f} {float(self.lineEdit_astig2_angle.text()):.3f}\n"
                )

                f.write(
                    "# Magnitude and angle of 3-fold astigmatism (in Ang. and degrees)\n"
                )
                f.write(
                    f"{float(self.lineEdit_astig3_mag.text()):.3f} {float(self.lineEdit_astig3_angle.text()):.3f}\n"
                )

                f.write("# Crystal tilt x,y in mrad\n")
                f.write(
                    f"{float(self.lineEdit_tilt_x.text()):.3f} {float(self.lineEdit_tilt_y.text()):.3f}\n"
                )

                f.write("# Number of detectors\n")
                f.write(f"{self.tableWidget_detectors.rowCount()}\n")

                f.write(
                    "# Detectors min,max angles(mrad) of collector (must match number of detectors)\n"
                )
                detectors = []
                for i in range(self.tableWidget_detectors.rowCount()):
                    detector = []
                    for j in range(1, self.tableWidget_detectors.columnCount()):
                        detector.append(
                            float(self.tableWidget_detectors.item(i, j).text())
                        )
                    detectors.append(" ".join([f"{d:.3f}" for d in detector]) + "\n")

                f.write("".join(detectors))

                f.write("# scan area xi, xf, yi, yf in angstrom\n")
                line15 = " ".join(
                    [
                        f"{float(self.lineEdit_scan_x_min.text()):.3f}",
                        f"{float(self.lineEdit_scan_x_max.text()):.3f}",
                        f"{float(self.lineEdit_scan_y_min.text()):.3f}",
                        f"{float(self.lineEdit_scan_y_max.text()):.3f}",
                    ]
                )
                f.write(f"{line15}\n")

                f.write("# slice thickness (A)\n")
                f.write(f"{float(self.lineEdit_slice_thickness.text()):.3f}\n")

                f.write("# Size of probe wavefunction Nx,Ny in pixels\n")
                f.write(
                    f"{int(self.lineEdit_probe_wf_x.text())} {int(self.lineEdit_probe_wf_y.text())}\n"
                )

                f.write("# size of transmission function Nx, Ny in pixels\n")
                f.write(
                    f"{int(self.lineEdit_transmission_func_x.text())} {int(self.lineEdit_transmission_func_y.text())}\n"
                )

                f.write("# size of ouput image nx, ny in pixel\n")
                f.write(
                    f"{int(self.lineEdit_output_size_x.text())} {int(self.lineEdit_output_size_y.text())}\n"
                )

                f.write("# input filename\n")
                f.write(f"{os.path.basename(self.file_name_xyz)}\n")

                f.write(
                    "# optional calculation: save STEM diffraction pattern every N steps\n"
                )
                f.write(f"{int(self.lineEdit_diffraction_n.text())}\n")

                f.write(
                    "# do not change settings below this line, unless you know what you are doing (GPU_PROBE, SMOOTH_PROBE, num_parallel, use_hdd)\n"
                )
                line29 = " ".join(
                    [
                        f"{int(self.lineEdit_GPU_PROBE.text())}",
                        f"{int(self.lineEdit_SMOOTH_PROBE.text())}",
                        f"{int(self.lineEdit_num_parallel.text())}",
                        f"{int(self.lineEdit_use_hdd.text())}",
                    ]
                )
                f.write(f"{line29}\n")

    def onclick_save_fig(self):
        basename, fname = os.path.split(self.file_name_xyz)
        fname, _ext = os.path.splitext(fname)

        out = QFileDialog.getExistingDirectory(
            self, dir=self.tr(basename), caption=self.tr("Open directory"),
        )
        if out:
            for i in range(self.comboBox_plt.count()):
                self.comboBox_plt.setCurrentIndex(i)
                axes = self.comboBox_plt.currentText()

                # temporarily set figure size for save
                width = self.mplWidget.axes.figure.get_figwidth()
                height = self.mplWidget.axes.figure.get_figheight()

                bbox = Bbox([[0, 0], [width * 2, height * 2]])
                self.mplWidget.axes.figure.savefig(
                    os.path.join(basename, f"{fname}_{axes}.png"),
                    dpi=120,
                    bbox_inches=bbox,
                )


#    def plot_slices(self):
#        self.surface = QtDataVisualization.Q3DSurface()
#        self.surfaceProxy = QtDataVisualization.QSurfaceDataProxy()
#        self.surfaceSeries = QtDataVisualization.QSurface3DSeries(self.surfaceProxy)
#        self.surface.addSeries(self.surfaceSeries)

#        self.scatter3d.addCustomItem(QtDataVisualization.QCustom3DItem(QtDataVisualization.QScatterDataItem(QVector3D(*[100,100,100]))))

if __name__ == "__main__":
    QCoreApplication.setAttribute(Qt.AA_ShareOpenGLContexts)
    app = QApplication([])
    app.setQuitOnLastWindowClosed(True)
    #    app.aboutToQuit.connect(exiting)
    widget = stemCL()
    widget.ui.show()
    sys.exit(app.exec_())
