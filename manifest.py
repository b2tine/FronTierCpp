#!/usr/bin/python

import os
import sys


docfiles = ('README','LICENSE')

infile_prefixes = ['in-']
script_prefixes = ['build','mrun','qrun-inputdir']
prefixes = tuple(infile_prefixes + script_prefixes)

autotools_suffixes = ['.ac','.am','.m4']
src_suffixes = ['.c','.cpp','.h','.hpp','.cu','.cuh']
suffixes = tuple(autotools_suffixes + src_suffixes)


def createFileList(projdir):
    filelist = list()
    for (root, dirs, files) in os.walk(projdir):
        for filename in files:
            if( filename.endswith(suffixes)
                    or filename.startswith(prefixes)
                    or filename in docfiles ):
                filelist.append(os.path.join(root,filename))
    return filelist


def writeManifest(projdir):
    manifestFile = open('manifest.txt','w')
    for file in createFileList(projdir):
        manifestFile.write(str(file) + "\n")
    manifestFile.close()
    return



topdir = str(sys.argv[1])
writeManifest(topdir)



