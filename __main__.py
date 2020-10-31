import sys
from PySide2.QtWidgets import QApplication
from PySide2.QtCore import QCoreApplication, Qt
from .interface import stemCL

QCoreApplication.setAttribute(Qt.AA_ShareOpenGLContexts)
app = QApplication([])
app.setQuitOnLastWindowClosed(False)
#    app.aboutToQuit.connect(exiting)
widget = stemCL()
widget.ui.show()
sys.exit(app.exec_())
