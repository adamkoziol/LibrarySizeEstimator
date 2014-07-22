#!/usr/bin/env python
__author__ = 'akoziol'

from multiprocessing import Pool
# Shutil is useful for file moving/copying
import shutil
# Subprocess->call is used for making system calls
import subprocess
import sys
# Errno is used in the file creation command  - I think it's similar to the $! variable in Perl
import errno
# OS is used for file/folder manipulations
import os
import time
# Glob finds all the path names matching a specified pattern according to the rules used by the Unix shell
import glob

path = os.getcwd()

os.chdir("%s/Best_Assemblies" % path)

referenceFile = glob.glob("*.fa*")
references = ["%s/Best_Assemblies/%s" % (path, fastaFile) for fastaFile in referenceFile]

targets = [reference.split('.')[0] for reference in referenceFile]

# Create a dictionary of sorted tuples using zip
inputData = dict(zip(sorted(references), sorted(targets)))


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


dotcount = 0


def dotter():
    """This function is borrowed from Mike Knowles. It allows for the addition of pretty
    dots after every pool is finished its task. Additionally, it formats the dots such that
    there are only 80 dots per line, and the date is added at the start of each line"""
    global dotcount
    if dotcount <= 80:
        sys.stdout.write('.')
        # I added this flush command, as the dots were not being printed until the script
        # finished processing
        sys.stdout.flush()
        dotcount += 1
    else:
        sys.stdout.write('\n[%s].' % (time.strftime("%H:%M:%S")))
        dotcount = 0


def indexTargetsProcesses():
    sys.stdout.write('Indexing targets\n')
    indexTargetArgs = []
    if __name__ == '__main__':
        indexTargetsPool = Pool()
        # Initialise the pool of processes - it defaults to the number of processors
        for reference, target in inputData.iteritems():
            indexTargetArgs.append((reference, target))
        indexTargetsPool.map(indexTargets, indexTargetArgs)



def indexTargets((reference, target)):
    """Performs smalt index on the targets using the range of k-mers stored in the variable kmer"""
    filename = target.split('.')[0]
    # Create a new path to be created (if necessary) for the generation of the range of k-mers
    indexPath = "%s/targets/%s" % (path, filename)
    # Call the make_path function to make folders as necessary
    make_path(indexPath)
    shutil.copy(reference, indexPath)
    os.chdir(indexPath)
    indexFileSMI = "%s.smi" % filename
    if not os.path.isfile("%s/%s" % (indexPath, indexFileSMI)):
        indexCommand = "smalt index -k 20 -s 10 %s %s" % (target, reference)
        subprocess.call(indexCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        dotter()
    else:
        dotter()


def mappingProcesses():
    """Mapping threads!"""
    os.chdir(path)
    print '\nPerforming reference mapping'
    mappingProcessesArgs = []
    if __name__ == '__main__':
        mappingProcessesPool = Pool()
        # uses target
        for reference, target in inputData.iteritems():
            mappingProcessesArgs.append((reference, target))
        mappingProcessesPool.map(mapping, mappingProcessesArgs)


def mapping((reference, target)):
    """Performs the mapping of the simulated reads to the targets"""
    filename = target.split('.')[0]
    os.chdir("%s/%s" % (path, target))
    fastq1 = str("%s/%s/%s" % (path, target, glob.glob("*R1_001.fastq")[0]))
    fastq2 = str("%s/%s/%s" % (path, target, glob.glob("*R2_001.fastq")[0]))
    filePath = "%s/tmp/%s" % (path, target)
    make_path(filePath)
    targetPath = "%s/targets/%s/%s" % (path, filename, filename)
    if not os.path.isfile("%s/%s.bam" % (filePath, target)):
        smaltMap = "smalt map -o %s/%s.bam -f bam -x %s %s %s" \
                   % (filePath, target, targetPath, fastq1, fastq2)
        subprocess.call(smaltMap, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        dotter()
    else:
        dotter()


def extractingProcesses():
    """Mapping threads!"""
    os.chdir(path)
    print '\nExtracting insert sizes'
    extractingProcessesArgs = []
    if __name__ == '__main__':
        extractingProcessesPool = Pool()
        # uses target
        for reference, target in inputData.iteritems():
            extractingProcessesArgs.append((reference, target))
        extractingProcessesPool.map(extractInsertSize, extractingProcessesArgs)


def extractInsertSize((reference, target)):
    """Uses samtools view and Linux cut to extract the column of interest (column 9), which contains the distance between
    mapped paired reads"""
    # samtools view HG00418_A.bam | cut -f9 > HG00418_A.insertsizes.txt
    filePath = "%s/tmp/%s" % (path, target)
    if not os.path.isfile("%s/%s_insertsizes.csv" % (filePath, target)):
        extractCommand = "samtools view %s/%s.bam | cut -f9 > %s/%s_insertsizes.csv" % (filePath, target, filePath, target)
        subprocess.call(extractCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        dotter()
    else:
        dotter()


def graphingProcesses():
    """Mapping threads!"""
    print '\nGraphing results'
    graphingProcessesArgs = []
    if __name__ == '__main__':
        graphingProcessesPool = Pool()
        # uses target
        for reference, target in inputData.iteritems():
            graphingProcessesArgs.append(target)
        graphingProcessesPool.map(graphing, graphingProcessesArgs)


def graphing(target):
    """Uses samtools view and Linux cut to extract the column of interest (column 9), which contains the distance between
    mapped paired reads"""
    # samtools view HG00418_A.bam | cut -f9 > HG00418_A.insertsizes.txt
    filePath = "%s/tmp/%s" % (path, target)
    newPath = "%s/insertSizes" % (path)
    make_path(newPath)
    os.chdir(newPath)
    if not os.path.isfile("%s/%s_insert_sizes.pdf" % (newPath, target)):
        graphingCommand = "Rscript /home/blais/PycharmProjects/LibrarySizeEstimator/insertsizes.R %s %s" % (filePath, target)
        subprocess.call(graphingCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        dotter()
    else:
        dotter()


def formatOutput():
    os.chdir("%s/insertSizes" % path)
    print("\nFormatting Outputs")
    # Determine the folder name by taking the last folder name from the path
    folderName = path.split('/')[-1]
    if not os.path.isfile("%s_insertSizes.csv" % folderName):
        textfiles = glob.glob("*.txt")
        with open("%s_insertSizes.csv" % folderName, "a") as outputFile:
            outputFile.write("Strain\tMedian Insert Size\tStandard Deviation\n")
            for files in textfiles:
                infile = open(files, "r")
                inData = infile.read()
                infile.close()
                outputFile.write("%s\n" % inData)
                os.remove(files)
                dotter()
        print "\nFormatting complete"
    else:
        print "Formatting complete"


def pipeline():
    """Calls all the functions in a way that they can be multi-processed"""
    indexTargetsProcesses()
    #Start the mapping operations
    mappingProcesses()
    extractingProcesses()
    graphingProcesses()
    formatOutput()

start = time.time()
pipeline()
print "\nElapsed Time: %s seconds" % (time.time() - start)