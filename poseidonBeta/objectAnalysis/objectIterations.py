#!/usr/bin/env python

import sys
import logging
import argparse
from glob import glob
from pathlib2 import Path
import gc
import gzip

'''
try:
    import cPickle as pickle
except ModuleNotFoundError:
    import pickle
'''
import dill as pickle

from datetime import datetime
startTime = datetime.now()

from populateClasses import * 
from EECalculation import processEE


def memoryInfo(verbosePrint):
    '''
    '''
    process = psutil.Process(os.getpid())
    bytes_info = float(process.memory_info().rss)
    gb = bytes_info * 1e-9
    print('memory use:', round(gb, 3), 'GB')  # in gbytes




def weightingPopulation(weighting):
    '''
    Save a list of weightings per frame from biased simulations in the
    EEclass
    '''
    dataFile = weighting[0]
    try:
        frame_chunks = weighting[1]
    except IndexError:
        frame_chunks = 1
    weighting_list = []
    with open(dataFile) as data:
        for line in data:
            weighting_list.append(line)
    data.close()

    weighting_info = [weighting_list, frame_chunks]

    return weighting_info


def readAllObjects(inp='inp', paths='paths', 
        pathClasses='pathClasses', 
        temperature='temperature', 
        entropyEnergy='entropyEnergy', 
        level_list='level_list', 
        solvent='solvent', water='water',
        verbose='verbose', weighting='weighting', 
        forceUnits='forceUnits'):
    '''
    Read in list of objects and run analyses. 
    '''

    verbosePrint = print if verbose else lambda *a, **k: None
    print(startTime)



    waterTuple = ['WAT', 'wat', 'SOL', 'H2O', 'h2o', 'WAT_O', 'TIP3']
    if water != 'WAT':
        waterTuple = [water]
    if solvent == None: ##when solvent is NOT water
        solvent = waterTuple




    print('\nsolvent: %s' % (solvent))
    print('\nwater: %s' % (waterTuple))
    print('\n1. Populate Dictionaries\n')

    count = 1
    totAtoms = 0
    list_len = None
    atom_nums = []
    atom_count = 0
    objectIteration = 0

    EEclass = None 
    EEclass_residLevel_resname = None
    EEclass_soluteContacts = None
    EEclass_atomLevel = None


    class_str_list = ['EEclass', 'EEclass_residLevel_resname', 
            'EEclass_soluteContacts', 'EEclass_atomLevel', 
            'totAtoms', 'atom_nums', 'atom_count', 'objectIteration']


    ###if class objects already exist, then populate
    for path in pathClasses:

        def readObjs(Aclass, Aclass_str):
            '''
            read in each class separately
            '''
            pathDataFile = glob(path+'/'+Aclass_str+'.obj')
            for pd in pathDataFile:
                p = Path(pd)
                if p.is_file():
                    #globals()[Aclass] = None
                    print(p)
                    with gzip.open(str(pd), 'rb') as data:
                        gc.disable()
                        try:
                            Aclass = pickle.load(data)
                            #print(Aclass_str, Aclass)
                        except EOFError:
                            pass
                    data.close()
            return Aclass

        
        EEclass = readObjs(EEclass, 'EEclass')
        EEclass_residLevel_resname = readObjs(EEclass_residLevel_resname, 
                'EEclass_residLevel_resname')
        EEclass_soluteContacts = readObjs(EEclass_soluteContacts, 
                'EEclass_soluteContacts')
        EEclass_atomLevel = readObjs(EEclass_atomLevel, 
                'EEclass_atomLevel')
        totAtoms = readObjs(totAtoms, 'totAtoms')
        atom_nums = readObjs(atom_nums, 'atom_nums')
        atom_count = readObjs(atom_count, 'atom_count')
        objectIteration = readObjs(atom_count, 'objectIteration')



    memoryInfo(verbosePrint)
    print(datetime.now() - startTime)
    sys.stdout.flush()
    if totAtoms != 0:
        count = 2



    #sortFileNumerically
    def key_func(path):
        '''hard-coded to sort paths to match with weighting.csv'''
        path = path.split('/')

        sim_digits = ''
        array_digits = ''
        n_digits = 0
        ml_digits = ''
        for p in path:
            if 'sim' in p:
                for s in p:
                    if s.isdigit():
                        sim_digits += str(s)
                    else:
                        continue
            elif 'moleculeList' in p:
                for s in p:
                    if s.isdigit():
                        ml_digits += str(s)
                    else:
                        continue
            elif 'job_array' in p:
                for s in p:
                    if s.isdigit():
                        array_digits += str(s)
                    else:
                        continue
            else:
                continue

        try:
            n_digits = int(path[-2])
        except ValueError:
            pass

        if len(sim_digits) == 0:
            sim_digits = 0

        if len(array_digits) == 0:
            array_digits = 0

        if len(ml_digits) == 0:
            ml_digits = 0


        return int(sim_digits), int(array_digits), int(n_digits), int(ml_digits)


    weighting_info = None
    ##Add new information into exisiting classes
    for path in paths:
        pathDataFile = sorted(glob(path+'/'+inp), key=key_func)
        for pd in pathDataFile:
            p = Path(pd)
            if p.is_file():
                print(p)
                with gzip.open(str(pd), 'rb') as data:
                    gc.disable()
                    try:
                        atomList = pickle.load(data)
                        list_len = len(atomList[0])
                    except EOFError:
                        pass
                    totAtoms += len(atomList)
                    if count == 1:
                        for i in atomList:
                            atom_num = i[2]
                            if atom_num not in atom_nums:
                                atom_nums.append(atom_num)
                            else:
                                break
                        if weighting != None:
                            weighting_info = weightingPopulation(weighting)

                EEclass, EEclass_residLevel_resname, \
                        EEclass_soluteContacts, EEclass_atomLevel, \
                        atom_count = \
                            classPopulation(
                        atomList, entropyEnergy, 
                        level_list, count, 
                        atom_count, len(atom_nums), 
                        waterTuple, temperature, 
                        solvent, EEclass, 
                        EEclass_residLevel_resname, EEclass_soluteContacts, 
                        EEclass_atomLevel, 
                        weighting_info, p, verbosePrint)
                count += 1 #needed for class initiation
                gc.enable()
                data.close()
                memoryInfo(verbosePrint)
                print(datetime.now() - startTime)
                sys.stdout.flush()


    print(datetime.now() - startTime)
    memoryInfo(verbosePrint)
    sys.stdout.flush()

    
    try:
        totFrames = float(totAtoms) / float(len(atom_nums))
        print('\nTotal number of frames: %s' % (totFrames))
        print('Number of atoms in each frame: %s' % (len(atom_nums)))
        print('Number of variables in each list: %s' % (list_len))
    except ZeroDivisionError:
        logging.error('No frames to analyse, please chose correct '\
                'path to .obj files')
        sys.exit()

    ###once all the classes have been populated, calculate properties
    #and output to files
    num_frames = None
    print('\n2. Process Dictionaries')


    for level in level_list:

        print('---level:', level)
        EEclass2, DSclass2 = None, None

        if level in [None, 'moleculeLevel']:
            EEclass2 = EEclass
        if level in ['residLevel_resname']:
            EEclass2 = EEclass_residLevel_resname
        if level in ['soluteContacts']:
            EEclass2 = EEclass_soluteContacts
        if level in ['atomLevel']:
            EEclass2 = EEclass_atomLevel

        
        if entropyEnergy:
            name = 'EE'
            processEE(num_frames, totFrames, EEclass2, 
                    solvent, waterTuple, 
                    temperature, level, name, forceUnits, verbosePrint)


    #'''
    ##Save each class as an object so that we can continue populating
    #in stages.
    if len(paths) != 0:
        objectIteration += 1
        if pathClasses:
            print('\n3. Save Dictionaries')
            for Aclass in class_str_list:
                with gzip.GzipFile('%s.obj' % (Aclass), 'wb') as pickleFile:
                    pickle.dump((locals()[Aclass]), pickleFile, protocol=2)
                    pickleFile.close()
                    #locals()[Aclass] = None
            print('Number of objectIteration cycles saved: %s' % 
                    (objectIteration))
            print('Total atoms processed: %s' % (totAtoms))

    #'''


    sys.stdout.flush()
    print('\n')
    print(datetime.now() - startTime)

def main():

    try:
        usage = 'objectIterations.py [-h]'
        parser = argparse.ArgumentParser(description=\
                'Analysis of saved object from entropyAnalysis.py', 
                usage=usage, 
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        group = parser.add_argument_group('Options')
        group = parser.add_argument('-i', '--inp', action='store', 
                default='moleculeList*.obj', help='input object file name')
        group = parser.add_argument('-p', '--paths', action='store', 
                nargs='+', default=[], 
                help='list of paths to directories containing '
                'allMoleculeList.obj file.')
        group = parser.add_argument('-pc', '--pathClasses', action='store', 
                nargs='+', default=[], 
                help='list of paths to directories containing '
                '\'class\'.obj files to load and carry on populatng.')
        group = parser.add_argument('-T', '--temperature', 
                action='store', default=298, 
                help='include temperature (K) of the system, '\
                'default temperature used is 298 K')
        group = parser.add_argument('-SE', '--entropyEnergy', 
                action='store_true', 
                help='run entropy and energy analysis '\
                '(Smix, Sor, Svib, PE, KE)')
        group = parser.add_argument('-l', '--level_list', 
                action='store', 
                nargs='+', default=['moleculeLevel'], 
                help='choose refinement level of analysis: '\
                'moleculeLevel, residLevel_resname, '\
                'residLevel_atom, atomLevel, soluteContacts, res_atomLevel')
        group = parser.add_argument('-sol', '--solvent', 
                action='store', nargs='+',
                default=None, help='include resname of solvent molecules'\
                ' (case-sensitive)')
        group = parser.add_argument('-wat', '--water', action='store', 
                default='WAT', help='resname for water molecules')
        group = parser.add_argument('-v', '--verbose', 
                action='store_true', 
                help='print out progress of each analysis step')
        group = parser.add_argument('-w', '--weighting', 
                action='store', default=None, nargs='+',
                help='get weighting for each frame if simulations are '\
                'biased, currently the number on .obj files is used to '\
                'find the index in a list of all weightings. '\
                'The file containing weights should contain a list '\
                'of all weights in one column. The first argument '\
                'should be the file path of the weight file, '\
                'the second argument is the max number on each .obj file.')
        group = parser.add_argument('-fu', '--forceUnits', 
                action='store', default='kcal',
                help='provide units that forces are in - kJ or kcal.')

        op = parser.parse_args()
    except argparse.ArgumentError:
        logging.error('command line arguments are ill-defined, '
        'please check the arguments')
        raise
        sys.exit(1)

    readAllObjects(inp=op.inp, paths=op.paths, 
            pathClasses=op.pathClasses, 
            temperature=op.temperature, 
            entropyEnergy=op.entropyEnergy, 
            level_list=op.level_list,
            solvent=op.solvent, water=op.water, verbose=op.verbose,
            weighting=op.weighting, forceUnits=op.forceUnits)

if __name__ == '__main__':
    main()


