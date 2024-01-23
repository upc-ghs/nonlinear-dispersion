'''
General configuration for running models 
'''

import os 

BASEDIR = os.getcwd()
FOLDERS = {
        'sims'        : '../sims/',
        'lib'         : '../lib/',
    }
LIBFILE = 'libmf6.so'  # linux
#LIBFILE = 'libmf6.dll' # wos
