#!/usr/bin/env python

import os
import sys
import argparse
from pathlib import Path
from sys import argv
from glob import glob
import operator
import logging

from frameIterations import *

from datetime import datetime



__author__ = 'Jas Kalayan'
__copyright__ = 'Copyright 2020, Poseidon Project'
__credits__ = ['Jas Kalayan', 'Jon Higham', 'Richard Henchman']
__license__ = 'GPL'
__version__ = '2.0.0'
__maintainer__ = 'Jas Kalayan'
__email__ = 'jkalayan@gmail.com'
__status__ = 'Development'





def poseidon(start='start', end='end', 
        step='step', lammps='lammps', amber='amber', gromacs='gromacs', 
        pdb='pdb', pureAtomNum='pureAtomNum', cutShell='cutShell', 
        excludedResnames='excludedResnames',perFrameObj='perFrameObj',
        water='water', verbose='verbose'):


    startTime = datetime.now()

    verbosePrint = print if verbose else lambda *a, **k: None

    print(startTime)

    initiateAnalyses(start, end, step, lammps, amber, gromacs, pdb, 
            pureAtomNum, cutShell, excludedResnames, perFrameObj, 
            water, verbosePrint, startTime)

    print(datetime.now() - startTime)


def main():

    try:
        usage = 'runPoseidon_beta.py [-h]'
        parser = argparse.ArgumentParser(description='Program for reading '\
                'in Molecular Dynamics Simulation files for: '\
                'Prediction Of a System\'s Entropy Including a '\
                'Determination Of its Nature - POSEIDON beta v2', usage=usage, 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        group = parser.add_argument_group('Options')
        group = parser.add_argument('-s', '--start', action='store', 
                default='1', help='starting frame number')
        group = parser.add_argument('-e', '--end', action='store', 
                default='1', help='end frame number')
        group = parser.add_argument('-dt', '--step', action='store', 
                default='1', help='steps between frames')
        group = parser.add_argument('-l', '--lammps', nargs='+', 
                metavar='file', help='list of lammps input files,'
                ' up to 5 files: amber_top, amber_crd, lammps_traj, '
                'lammps_force and lammps_energy files.')
        group = parser.add_argument('-a', '--amber', nargs='+', 
                metavar='file', help='list of amber input files,'
                ' up to 3 files: amber_top, amber_crd, amber_forces')
        group = parser.add_argument('-g', '--gromacs', nargs='+', 
                metavar='file', help='list of gromacs input files,'
                ' up to 2 files: gro_tpr, gro_trr')
        group = parser.add_argument('-p', '--pdb', nargs='+', 
                metavar='file', help='pdb file with multiple frames'
                ' up to 1 file: frames separated by ENDMDL')
        group = parser.add_argument('-pn', '--pureAtomNum', action='store', 
                default='1', help='reference molecule resid for pure liquid')
        group = parser.add_argument('-cs', '--cutShell', action='store', 
                default=None, help='include cutoff shell analysis, '\
                'add cutoff distance in angstrom')
        group = parser.add_argument('-ex', '--excludedResnames', 
                action='store', nargs='+',
                default=None, help='exclude a list of molecule names '\
                'from nearest non-like analysis')
        group = parser.add_argument('-obj', '--perFrameObj', 
                action='store_true', 
                help='output object every frame, use this for large systems')
        group = parser.add_argument('-wat', '--water', action='store', 
                default='WAT', help='resname for water molecules')
        group = parser.add_argument('-v', '--verbose', 
                action='store_true', 
                help='print out progress of each analysis step')

        op = parser.parse_args()
    except argparse.ArgumentError:
        logging.error('Command line arguments are ill-defined, '
        'please check the arguments.')
        raise
        sys.exit(1)



    poseidon(start=op.start, end=op.end, 
            step=op.step, lammps=op.lammps, amber=op.amber, 
            gromacs=op.gromacs, pdb=op.pdb,
            pureAtomNum=op.pureAtomNum, cutShell=op.cutShell,
            excludedResnames=op.excludedResnames, perFrameObj=op.perFrameObj, 
            water=op.water, verbose=op.verbose)

if __name__ == '__main__':
    main()


