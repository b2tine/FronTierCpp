#!/usr/bin/python

import os
import sys

excludedfiles = tuple('config.h')

docfiles = ('README','LICENSE')

infile_prefixes = ['in-']
script_prefixes = ['build','mrun','qrun-inputdir']
prefixes = tuple(infile_prefixes + script_prefixes)

script_suffixes = ['.py']
autotools_suffixes = ['.ac','.am','.m4']
src_suffixes = ['.c','.cpp','.h','.hpp','.cu','.cuh']
suffixes = tuple(script_suffixes + autotools_suffixes + src_suffixes)


def createFileList(projdir):
    filelist = list()
    for (root, dirs, files) in os.walk(projdir):
        for filename in files:
            if( filename.endswith(suffixes)
                    or filename.startswith(prefixes)
                    or filename in docfiles ):
                if( !filename in excludedfiles ):
                    filelist.append(os.path.join(root,filename))
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



