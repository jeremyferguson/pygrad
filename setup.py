from setuptools import setup, Extension
from Cython.Distutils import build_ext
import os
import Cython
import numpy
from io import open
import glob

def returnPyxFiles(path):
    l = []
    for i in os.listdir(path):
        if i.endswith(".pyx") or i.endswith(".c") or i.endswith(".pxd"):
            l.append(path+i)
    return l


def returnPxdFiles(path):
    l = []
    for i in os.listdir(path):
        if i.endswith(".pxd") or i.endswith(".c"):
            l.append(path+i)
    return l

def makeExtensions(path):
    extensions = []
    for root,dirs,files in os.walk(path):
        for file in files:
            if file.endswith(".pyx"):
                moduleFiles = [] 
                moduleFiles.append(root+'/'+file)
                if root+'/'+(file.split('.')[0]+'.pxd') in glob.glob(root+"/*.pxd"):
                    moduleFiles.append(root+'/'+(file.split('.')[0]+'.pxd'))
                pathWithFile = root+'/'+file.split('.')[0]
                moduleName = pathWithFile.replace('/','.')
                extensions.append(Extension(moduleName,moduleFiles,include_dirs=[numpy.get_include(),'./Pygrad/','./Pygrad/C/']))
    return extensions
extensions = makeExtensions('Pygrad')
setup(
    setup_requires=[
        'cython>=0.2',
    ],
    zip_safe=False,
    name='Pygrad',  # Required
    packages=['Pygrad'],

    version='1.1.0',  # Required
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*, <4',
    package_dir={'Pygrad': 'Pygrad'},
    package_data={  # Optional
        'Pygrad': returnPxdFiles("./Pygrad/"),
    },

    install_requires=['numpy','cython','PyGasMix @ git+https://github.com/UTA-REST/PyGasMix.git#egg=PyGasMix-1.2.0'],  # Optional
    include_package_data = True,
    ext_modules = extensions,
    cmdclass={'build_ext': build_ext},
)

