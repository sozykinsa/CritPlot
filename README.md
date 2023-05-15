<a href="https://codecov.io/gh/sozykinsa/CritPlot" >
<img src="https://codecov.io/gh/sozykinsa/CritPlot/branch/main/graph/badge.svg?token=MS7qX21EMW"/>
</a>

# EDTop-Visual

EDTop-Visual (Electron Density Topology Visualizer) - graphical user interface for Critic2 and TOPOND programms. It is a cross-platform program. 

## Install
EDTop-Visual program is written in Python 3 (version >= 3.8). It has some dependences. To install the necessary modules, run in the terminal (command line):

pip3 install -r ./requirements.txt

To run the program, type

python3 edtop-visual.py

## Problems with Install?

Some operating systems may require additional packages to be installed:

Xubuntu: sudo apt-get install qtbase5-dev

You have to set the variable QT_API:

export QT_API=pyside2 (in linux, or https://www.architectryan.com/2018/08/31/how-to-change-environment-variables-on-windows-10/ in Windows)