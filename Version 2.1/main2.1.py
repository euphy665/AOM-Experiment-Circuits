import sys
import os
#from typing_extensions import ParamSpec
from PyQt5.QtWidgets import QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QApplication, QPushButton, QLabel, QLineEdit,QFileDialog
from PyQt5 import QtCore, QtGui, QtWidgets, QtChart
import pyqtgraph as pg
from pyqtgraph import PlotWidget
from pyqtgraph.graphicsItems.GridItem import GridItem
import pyvisa as visa
import numpy as np
import pandas as pd
from scipy import optimize
from Ui_pyODmeter import Ui_MainWindow
from osc import Oscilloscope
from lmfit import Model

rigol = Oscilloscope('TCPIP::192.168.3.243::INSTR')
rigol.open()


def threeLevelAbsorption(x, offset, amp, OD, Delta_p, gamma_12, Omega_c, Delta_c):
    temp = x-Delta_p-Delta_c
    return offset+amp*np.exp(-OD*(3*(temp*(3*temp+gamma_12*(x-Delta_p))+gamma_12*((Omega_c**2)/4-(x-Delta_p)*temp+3*gamma_12))/(((Omega_c**2)/4-(x-Delta_p)*temp+3*gamma_12)**2+(3*temp+gamma_12*(x-Delta_p))**2)+0.2857/(1+((x+362)/3)**2)))

def twoLevelAbsorption(x, offset, amp, OD, Delta_p):
    return offset+amp*np.exp(-OD*(1/(1+((x-Delta_p)/3)**2)+0.2857/(1+((x+362)/3)**2)))

class PyODMeter(QtWidgets.QMainWindow,Ui_MainWindow):

    def __init__(self,oscilloscope = None):

        super(PyODMeter,self).__init__()
        self.setupUi(self)
        self.graph,self.graph_2,self.graph_3 = self.set_graph_ui()
        self.bg.clicked.connect(self.takeBackground)
        self.data.clicked.connect(self.takeData)
        self.single.clicked.connect(self.singleFrame)
        self.Save.clicked.connect(self.saveData)
        self.STOP.clicked.connect(self.stopAcquire)
        self.clearBg.clicked.connect(self.clearBackground)

        self.osc = oscilloscope
        self.xbg=np.ndarray([1,1])
        self.statusPara=[True,True,True,True,True,True,True];

        self.ybg=0
        self.xdata=0
        self.ydata=0
        self.guess_value = [0.00, 1.00, 40, 0.00, 1.00, 1.00, 1.00]
        self.bounds = ([-10, 0, 0, 0, 0, 0, -1], [0, 5, 700, 2, 1, 20, 1]);
        self.gamma12.setText('N.A.')
        self.OD.setText('0')

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

        self.enableOffset.stateChanged.connect(lambda: self.validParams(0,self.enableOffset.isChecked()))
        self.enableAmp.stateChanged.connect(lambda: self.validParams(1,self.enableAmp.isChecked()))
        self.enableOD.stateChanged.connect(lambda: self.validParams(2,self.enableOD.isChecked()))
        self.enableDelta_p.stateChanged.connect(lambda: self.validParams(3,self.enableDelta_p.isChecked()))
        self.enableGamma12.stateChanged.connect(lambda: self.validParams(4,self.enableGamma12.isChecked()))
        self.enableOmega_c.stateChanged.connect(lambda: self.validParams(5,self.enableOmega_c.isChecked()))
        self.enableDelta_c.stateChanged.connect(lambda: self.validParams(6,self.enableDelta_c.isChecked()))

        self.update_timer = QtCore.QTimer()
        self.update_timer.timeout.connect(self.singleFrame)
     #   self.update_timer.timeout.connect(self.updateText)

    def set_graph_ui(self):
        pg.setConfigOptions(antialias=True)
        win = PlotViewOsc()
        win_2 = PlotViewPara()#pg.PlotWidget(background = 'w')
        win_3 = PlotViewPara()
        self.graph_layout.addWidget(win) 
        self.graph_layout_2.addWidget(win_2)
        self.graph_layout_3.addWidget(win_3)
        return [win, win_2, win_3]

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
        self.graph.plotOsc(self.xbg,self.ybg_ave,type='bg')

    def plotShow(self,x,y,ybg):
        
        if ybg.ndim == 1:
            transmission = np.divide(y,ybg)
        else:
            transmission = np.divide(y,np.squeeze(ybg.mean(axis = 1)))
        transmission[np.isnan(transmission)]=1;
        transmission[np.isinf(transmission)]=1;
        self.transmission = transmission;
        self.graph.plotOsc(x,transmission)  

    def clearBackground(self):
        self.xbg = np.ndarray([1,1]);
        self.ybg = 0;

    def takeData(self):
        self.update_timer.start()

    def singleFrame(self):
        [xdata,ydata]=self.getOscData()
        if self.smoothFitting.isChecked():
            n = self.smoothPara.value()
            for i in range(1,n):
                temp=self.getOscData()
                xdata = xdata+temp[0]
                ydata = ydata+temp[1]
            xdata = xdata/n
            ydata = ydata/n           
        self.xdata = xdata
        self.ydata = ydata

        if self.xbg.shape != (1,1):
            self.plotShow(self.xdata,self.ydata,self.ybg)
            if self.enableFitting.isChecked():
                self.FittingCurve(self.xdata,self.transmission)
        else:
            self.graph.plotOsc(self.xdata,self.ydata,type='rawdata')

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

    def FittingCurve(self, fit_x, fit_y):
        len_x = len(fit_x)
        if self.comboBox.currentText() == "EIT": 
            if self.enhancedFitting.isChecked():
            # enhanced fitting, make the EIT transmission more important in fitting.
                bw=20;
                index = np.where(fit_y[int(len_x/2-bw):int(len_x/2+bw)]>0.1);
                EITpeak_x=fit_x[int(len_x/2-bw):int(len_x/2+bw)][index]
                EITpeak_y=fit_y[int(len_x/2-bw):int(len_x/2+bw)][index]           
                fit_x=np.concatenate([fit_x, EITpeak_x.repeat(self.enhancedPara.value()-1)])
                fit_y=np.concatenate([fit_y, EITpeak_y.repeat(self.enhancedPara.value()-1)])

            popt, pcov = optimize.curve_fit(threeLevelAbsorption, fit_x, fit_y, p0=self.guess_value, bounds=self.bounds)     
            print(popt)
            self.graph.plotOsc(fit_x[0:len_x], threeLevelAbsorption(fit_x[0:len_x], popt[0], popt[1], popt[2], popt[3], popt[4], popt[5], popt[6]), type='fitting')
            self.coef = popt
            self.gamma12.setText('%.2f'%self.coef[4])
            self.gamma12_2.setText('%.2f'%self.coef[4])
            self.Omega_c_2.setText('%.2f'%self.coef[5])
            self.Delta_c_2.setText('%.2f'%self.coef[6])
            self.graph_3.addPoint(self.coef[4])
        else:
            popt, pcov = optimize.curve_fit(twoLevelAbsorption, fit_x, fit_y, p0=self.guess_value[0:4], bounds=(self.bounds[0][0:4],self.bounds[1][0:4]))
            self.graph.plotOsc(fit_x, twoLevelAbsorption(fit_x, popt[0], popt[1], popt[2], popt[3]), type='fitting')
            self.coef = popt
            self.gamma12.setText('N.A.')
            self.gamma12_2.setText('N.A.')
            self.Omega_c_2.setText('N.A.')
            self.Delta_c_2.setText('N.A.')
            self.graph_3.addPoint(0)
            

        self.err = np.sqrt(np.diag(pcov))
        self.OD.setText('%.2f'%self.coef[2])
        self.Offset_2.setText('%.2f'%self.coef[0])
        self.Amplitude_2.setText('%.2f'%self.coef[1])
        self.OD_2.setText('%.2f'%self.coef[2])
        self.Delta_p_2.setText('%.2f'%self.coef[3])
        self.graph_2.addPoint(self.coef[2])

    def setGuessValue(self, value, index):
        self.guess_value[index] = float(value)

    def validParams(self,index,status):
        bounds_ref = ([-10, 0, 0, 0, 0, 0, -1], [0, 5, 700, 2, 1, 20, 1]);
        epsilon = 0.00001
        if status :
            self.bounds[1][index] = bounds_ref[1][index]
            self.bounds[0][index] = bounds_ref[0][index]
        else :
            self.bounds[0][index] = self.guess_value[index]-epsilon
            self.bounds[1][index] = self.guess_value[index]+epsilon


class PlotViewOsc(PlotWidget):
    def __init__(self):
        super().__init__()
        self.setBackground('w')
        self.showGrid(x = True, y = True, alpha = 1.0)
        self.getPlotItem().setYRange(-0.2,1)
    #    self.getPlotItem().hideAxis('left')
    
    def plotOsc(self,x,y,type='transmission'):
        if type == 'rawdata':
            self.clear()
            self.plot(x,y,pen='b')
            self.enableAutoRange(axis='y') 
        elif type == 'bg':
            self.clear()
            self.plot(x,y,pen='b')
            self.enableAutoRange(axis='y') 
        elif type == 'fitting':
            self.plot(x,y,pen='r')
            
        elif type == 'transmission':
            self.clear()
            self.plot(x,y,pen='b')
            self.getPlotItem().setYRange(-0.2,1)






class PlotViewPara(PlotWidget):
   def __init__(self):
        super().__init__()        
        self.x = np.arange(20)
        self.y = np.zeros(20)
        self.setBackground('w')
        self.showGrid(x = True, y = True, alpha = 1.0)
        self.getPlotItem().hideAxis('left')
        self.getPlotItem().hideAxis('bottom')
        self.plot(self.x,self.y)
        self.show()

   def addPoint(self, y):
        self.clear()
        self.y = np.append(self.y[1:],y)
        self.plot(self.x,self.y,pen='k')

    

if __name__ == '__main__':
  #  QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling)
    app = QApplication([])
    window = PyODMeter(rigol)
    window.show()
    app.exit(app.exec_())