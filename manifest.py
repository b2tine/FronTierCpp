#!/usr/bin/python

import os
import sys

excluded= tuple('config.h')

docfiles = ('README','LICENSE')

infile_prefixes = ['in-']
script_prefixes = ['build','mrun','crun','frun','qrun-inputdir']
file_prefixes = tuple(infile_prefixes + script_prefixes)
directory_prefixes = tuple(['refsol-'])

script_suffixes = ['.py','sh']
autotools_suffixes = ['.ac','.am','.m4']
src_suffixes = ['.c','.cpp','.h','.hpp','.cu','.cuh','.F90']
file_suffixes = tuple(script_suffixes + autotools_suffixes + src_suffixes)


def createFileList(projdir):
    filelist = list()
    for (root, dirs, files) in os.walk(projdir):
        for filename in files:
            if (filename.endswith(file_suffixes)
                    or filename.startswith(file_prefixes)
                    or filename in docfiles):
                if (filename not in excluded):
                    filelist.append(os.path.join(root,filename))
        for dirname in dirs:
            if (dirname.startswith(directory_prefixes)):
                dirpath = os.path.join(root,dirname)
                for f in os.listdir(dirpath):
                    filelist.append(os.path.join(dirpath,f))
    return filelist


def writeManifest(projdir):
    manifestFile = open('Manifest.txt','w')
    for file in createFileList(projdir):
        manifestFile.write(str(file) + "\n")
    manifestFile.close()
    return



#begin script
topdir = str(sys.argv[1])
writeManifest(topdir)



