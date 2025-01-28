#!/bin/bash
#FreeCAD path
export FREECADPATH=/Users/mayar/Desktop/ICL/4Y/Thesis/Codes/Contents/Resources/bin #eddychanged bin pathname to mac version

main_LOD3_file=$1 #eddychanged $

PYTHONPATH=${LODPATH} PATH=$PATH:${FREECADPATH} freecadcmd $main_LOD3_file $@

#From terminal, being inside src folder, activate environment and run in terminal as:  
#./LOD3.sh ../examples/p2_LOD3_00_School_main_LOD3.py
