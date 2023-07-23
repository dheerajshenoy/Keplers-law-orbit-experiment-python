import sys
from PyQt6.QtWidgets import (QPushButton, QLineEdit, QColorDialog, QLabel,
                             QSplitter, QMainWindow, QApplication, QWidget,
                             QVBoxLayout, QHBoxLayout)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

class MainApplication(QMainWindow):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.mainWidget = QWidget()
        self.mainLayout = QVBoxLayout()
        self.setCentralWidget(self.mainWidget)
        self.mainWidget.setLayout(self.mainLayout)
        self.show()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    win = MainApplication()

    sys.exit(app.exec())
