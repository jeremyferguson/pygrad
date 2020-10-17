from setuptools import setup, Extension
from Cython.Distutils import build_ext
import os
from Cython.Build import cythonize

import Cython
import glob
import numpy
from io import open


def returnPyxFiles(path):
    l = []
    for i in os.listdir(path):
        if i.endswith(".pyx") or i.endswith(".c"):
            l.append(path+i)
    return l


def returnPxdFiles(path):
    l = []
    for i in os.listdir(path):
        if i.endswith(".pxd"):
            l.append(path+i)
    return l
def makeExtensions(path):
    extensions = []
    for root,dirs,files in os.walk(path):
        for file in files:
            if file.endswith(".pyx"):
                moduleFiles = []
                moduleFiles.append(root+'/'+file)
                pathWithFile = root+'/'+file.split('.')[0]
                moduleName = pathWithFile.replace('/','.')
                extensions.append(Extension(moduleName,moduleFiles,include_dirs=[numpy.get_include(),'./Pygrad/','./Pygrad/C/']))
    return extensions
extensions = makeExtensions('Pygrad')
print(extensions)
setup(ext_modules=cythonize(extensions))
