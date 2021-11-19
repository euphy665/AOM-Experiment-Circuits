# OD Meter
This python file is for OD and EIT real time fitting. 

How to use the program?

1. Make sure that the python environment is already set. 
3. Install PyQt5 and pyvisa through cmd pip installation. Some package version are listed as follow:
    PyQt5                              5.15.6 
    PyQt5-Qt5                          5.15.2 //
    PyQt5-sip                          12.9.0 //
    pyqtgraph                          0.12.3 //
4. Install NI Visa from https://www.ni.com/zh-cn/support/downloads/drivers/download.ni-visa.html
5. Connect the oscillascope through Ethernet to the same network as the computer to run the program
6. Get the IP address of the oscillascope and change the 15 line of main.py address.
7. Run the main.py address

Click on the main.py file and run.

Oscillascope type: Rigol. DS2202A, if other oscillascope is used, the osc file need to be changed. 

When measuring:

1. Click 'take background' for several times to average
2. Click 'take data' to get the data flow constantly
3. Click the 'enable fitting' if it is not checked.



--Not finished yet.
