import sys
import os
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QApplication, QPushButton, QLabel, QLineEdit,QFileDialog

from PyQt5 import QtCore, QtGui, QtWidgets
import pyqtgraph as pg
from pyqtgraph import PlotWidget
import pyvisa as visa
import numpy as np
import pandas as pd
from scipy import optimize
from Ui_pyODmeter import Ui_MainWindow
from osc import Oscilloscope


rigol = Oscilloscope('TCPIP::192.168.3.213::INSTR')
rigol.open()


def threeLevelAbsorption(x, offset, amp, OD, Delta_p, gamma_12, Omega_c, Delta_c):
    return offset+amp*np.exp(-OD*(3*((x-Delta_p-Delta_c)*(3*(x-Delta_p-Delta_c)+gamma_12*(x-Delta_p))+gamma_12*((Omega_c**2)/4-(x-Delta_p)*(x-Delta_p-Delta_c)+3*gamma_12))/(((Omega_c**2)/4-(x-Delta_p)*(x-Delta_p-Delta_c)+3*gamma_12)**2+(3*(x-Delta_p-Delta_c)+gamma_12*(x-Delta_p))**2)+0.2857/(1+((x+362)/3)**2)))

def twoLevelAbsorption(x, offset, amp, OD, Delta_p):
    return offset+amp*np.exp(-OD*(1/(1+((x-Delta_p)/3)**2)+0.2857/(1+((x+362)/3)**2)))

class PyODMeter(QtWidgets.QMainWindow,Ui_MainWindow):

    def __init__(self,oscilloscope = None):

        super(PyODMeter,self).__init__()
        self.setupUi(self)
        self.graph = self.set_graph_ui()
        self.bg.clicked.connect(self.takeBackground)
        self.data.clicked.connect(self.takeData)
        self.single.clicked.connect(self.singleFrame)
        self.Save.clicked.connect(self.saveData)
        self.STOP.clicked.connect(self.stopAcquire)
        self.clearBg.clicked.connect(self.clearBackground)

        self.osc = oscilloscope
        self.xbg=np.ndarray([1,1])
    
        self.ybg=0
        self.xdata=0
        self.ydata=0
        self.guess_value = [0.00, 1.00, 40, 0.00, 1.00, 1.00, 1.00]

        self.Offset.setText('%.2f'%self.guess_value[0])
        self.Amplitude.setText('%.2f'%self.guess_value[1])
        self.OD_3.setText('%.2f'%self.guess_value[2])
        self.Delta_p.setText('%.2f'%self.guess_value[3])
        self.gamma_12.setText('%.2f'%self.guess_value[4])
        self.Omega_c.setText('%.2f'%self.guess_value[5])
        self.Delta_c.setText('%.2f'%self.guess_value[6])
        self.fileIndex.setText('%o'%1)
     
        self.Offset.editingFinished.connect(lambda:self.setGuessValue(self.Offset.text(), 0))
        self.Amplitude.editingFinished.connect(lambda:self.setGuessValue(self.Amplitude.text(), 1))
        self.OD_3.editingFinished.connect(lambda:self.setGuessValue(self.OD_3.text(), 2))
        self.Delta_p.editingFinished.connect(lambda:self.setGuessValue(self.Delta_p.text(), 3))
        self.gamma_12.editingFinished.connect(lambda:self.setGuessValue(self.gamma_12.text(), 4))
        self.Omega_c.editingFinished.connect(lambda:self.setGuessValue(self.Omega_c.text(), 5))
        self.Delta_c.editingFinished.connect(lambda:self.setGuessValue(self.Delta_c.text(), 6))


        self.update_timer = QtCore.QTimer()
        self.update_timer.timeout.connect(self.singleFrame)
     #   self.update_timer.timeout.connect(self.updateText)

    def set_graph_ui(self):
        pg.setConfigOptions(antialias=True)
        win = pg.PlotWidget(background = 'w') 
        self.graph_layout.addWidget(win) 
        return win
        #graph = win.addPlot()

    def getOscData(self):
        xorigin = self.osc.getParam()['xorigin']
        xincrement = self.osc.getParam()['xincrement']
        xreference = self.osc.getParam()['xreference']

        tbase = self.osc.query(':TIMebase:OFFSet?')
        trig_delay = -11e-6
        tbase_start = int(np.floor((float(tbase)-trig_delay-xorigin)/xincrement))
        tbase_pt = int(100e-6/xincrement)
        tbase_end = tbase_start+tbase_pt

        xdata = np.linspace(-40, 40, tbase_pt)
        ydata = self.osc.getWave()
        return xdata,ydata[tbase_start:tbase_end]

    def takeBackground(self):
        if self.xbg.shape == (1,1):
            [self.xbg,self.ybg] = self.getOscData()
            self.ybg_ave=self.ybg
        else:
            [self.xbg,ybg] = self.getOscData()
            if self.ybg.ndim == 1:
                self.ybg = self.ybg[:,np.newaxis]
            ybg = ybg[:,np.newaxis]
            self.ybg = np.concatenate((ybg,self.ybg),axis=1)            
            self.ybg_ave=ybg.mean(axis = 1)
        self.graph.clear()
        self.graph.plot(self.xbg,self.ybg_ave,pen="b")
        self.graph.enableAutoRange(axis='y')


    def plotShow(self,x,y,ybg):
        self.graph.clear()
        self.graph.setYRange(-0.2,1.2,padding=0)
        if ybg.ndim == 1:
            transmission = np.divide(y,ybg)
        else:
            transmission = np.divide(y,np.squeeze(ybg.mean(axis = 1)))
        transmission[np.isnan(transmission)]=1;
        transmission[np.isinf(transmission)]=1;
        self.transmission = transmission;
        self.graph.plot(x,transmission,pen="b")
    

    def clearBackground(self):
        self.xbg = np.ndarray([1,1]);
        self.ybg = 0;

    def takeData(self):
        self.update_timer.start()

    def singleFrame(self):
        [self.xdata,self.ydata]=self.getOscData()
        self.graph.clear()
        if self.xbg.shape != (1,1):
            
            self.plotShow(self.xdata,self.ydata,self.ybg)
            if self.enableFitting.isChecked():
                self.FittingCurve()
        else:
            self.graph.plot(self.xdata,self.ydata,pen="b")
            self.graph.enableAutoRange(axis='y')

    def saveData(self):
        if self.fileDirectory.text() == "":
            self.fileDirectory.setText(QFileDialog.getExistingDirectory(self))
        folder = self.fileDirectory.text()  
        if not os.path.isdir(folder):
               os.makedirs(folder)   
        
        if self.fileName.text() == "":
            self.fileName.setText('Untitled')

        df = pd.DataFrame([self.xdata,self.ydata,self.xbg,self.ybg])

        if self.autoIndex.isChecked():
            fileName = folder+'/'+self.fileName.text()+'_'+self.fileIndex.text()+'.txt'
            self.fileIndex.setText('%o'%(int(self.fileIndex.text())+1));
        else:
            fileName = folder+'/'+self.fileName.text()+'.txt'
        if self.enableFitting.isChecked():
            comments = 'Offset,'+self.Offset_2.text()+'; Amplitude, '+self.Amplitude_2.text()+'; OD, '+self.OD.text()+'; Delta_p, '+self.Delta_p_2.text()+'; Omega_c, '+self.Omega_c_2.text()+'; Delta_c, '+self.Delta_c_2.text()
        else:
            comments = ''
        data= np.array((self.xbg,self.ybg_ave,self.ydata))
        np.savetxt(fileName,data.T,fmt='%.4f',header='x background data',newline='\n',footer=comments)
        self.saveStatus.setText("Successfully save data to "+fileName)
        
    def startAcquire(self): 
        self.update_timer.start()

    def stopAcquire(self):
        self.update_timer.stop()

    def FittingCurve(self):
        if self.comboBox.currentText() == "EIT":
            popt, pcov = optimize.curve_fit(threeLevelAbsorption, self.xdata, self.transmission, p0=self.guess_value, bounds=([-10, 0, 0, 0, 0, 0, -1], [0, 5, 100, 2, 1, 20, 1]))
            self.graph.plot(self.xdata, threeLevelAbsorption(self.xdata, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6]), pen=pg.mkPen(color=(255, 0, 0)), name='Fitting curve')
            self.coef = popt
            self.gamma12.setText('%.2f'%self.coef[4])
            self.gamma12_2.setText('%.2f'%self.coef[4])
            self.Omega_c_2.setText('%.2f'%self.coef[5])
            self.Delta_c_2.setText('%.2f'%self.coef[6])
        else:
            popt, pcov = optimize.curve_fit(twoLevelAbsorption, self.xdata, self.transmission, p0=self.guess_value[0:4], bounds=([-10, 0, 0, 0], [0, 5, 100, 2]))
            self.graph.plot(self.xdata, twoLevelAbsorption(self.xdata, popt[0], popt[1], popt[2], popt[3]), pen=pg.mkPen(color=(255, 0, 0)), name='Fitting curve')
            self.coef = popt
            self.gamma12.setText('N.A.')
            self.gamma12_2.setText('N.A.')
            self.Omega_c_2.setText('N.A.')
            self.Delta_c_2.setText('N.A.')

        self.err = np.sqrt(np.diag(pcov))

        self.OD.setText('%.2f'%self.coef[2])
        

        self.Offset_2.setText('%.2f'%self.coef[0])
        self.Amplitude_2.setText('%.2f'%self.coef[1])
        self.OD_2.setText('%.2f'%self.coef[2])
        self.Delta_p_2.setText('%.2f'%self.coef[3])

    def setGuessValue(self, value, index):
        self.guess_value[index] = float(value)

if __name__ == '__main__':
    app = QApplication([])
    window = PyODMeter(rigol)
#   w=QMainWindow()
#   window.setupUi(w)
    window.show()
    app.exit(app.exec_())
