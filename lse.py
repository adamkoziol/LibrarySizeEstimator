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

# The path is still hardcoded as, most of the time, this script is run from within Pycharm.
os.chdir("/media/nas/akoziol/Pipeline_development/LibrarySizeEstimation")
path = os.getcwd()

targets = [name for name in os.listdir(".") if os.path.isdir(name) and name != "Best_Assemblies"]

os.chdir("%s/Best_Assemblies" % path)

referenceFile = glob.glob("*.fa*")
references = ["%s/Best_Assemblies/" % path + fastaFile for fastaFile in referenceFile]

inputData = dict(zip(sorted(references), sorted(targets)))

count = 0

def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def indexTargetsProcesses():
    global count
    sys.stdout.write('Indexing targets\n') if count == 0 else sys.stdout.write('')
    count = 1
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
    indexFileSMI = "%s.smi" % filename
    indexFileSMA = "%s.sma" % filename
    if not os.path.isfile("%s/%s" % (indexPath, indexFileSMI)):
        indexCommand = "smalt index -k 20 -s 10 %s %s" % (target, reference)
        subprocess.call(indexCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        shutil.move(indexFileSMI, indexPath)
        shutil.move(indexFileSMA, indexPath)
        shutil.copy(reference, indexPath)
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')


def mappingProcesses():
    """Mapping threads!"""
    os.chdir(path)
    print '\nPerforming reference mapping'
    mappingProcessesArgs = []
    if __name__ == '__main__':
        mappingProcessesPool = Pool()
        # uses target
        for reference, target in inputData.iteritems():
            mappingProcessesArgs.append(target)
        mappingProcessesPool.map(mapping, mappingProcessesArgs)


def mapping(target):
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
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')


def extractingProcesses():
    """Mapping threads!"""
    os.chdir(path)
    print '\nExtracting insert sizes'
    extractingProcessesArgs = []
    if __name__ == '__main__':
        extractingProcessesPool = Pool()
        # uses target
        for reference, target in inputData.iteritems():
            extractingProcessesArgs.append(target)
        extractingProcessesPool.map(extractInsertSize, extractingProcessesArgs)


def extractInsertSize(target):
    """Uses samtools view and Linux cut to extract the column of interest (column 9), which contains the distance between
    mapped paired reads"""
    # samtools view HG00418_A.bam | cut -f9 > HG00418_A.insertsizes.txt
    filePath = "%s/tmp/%s" % (path, target)
    if not os.path.isfile("%s/%s_insertsizes.csv" % (filePath, target)):
        extractCommand = "samtools view %s/%s.bam | cut -f9 > %s/%s_insertsizes.csv" % (filePath, target, filePath, target)
        #print extractCommand
        subprocess.call(extractCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')

def pipeline():
    """Calls all the functions in a way that they can be multi-processed"""
    indexTargetsProcesses()
    #Start the mapping operations
    mappingProcesses()
    extractingProcesses()


start = time.time()
pipeline()
print "\nElapsed Time: %s seconds" % (time.time() - start)