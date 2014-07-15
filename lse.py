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

# The path is still hardcoded as, most of the time, this script is run from within Pycharm.
os.chdir("/media/nas/akoziol/Pipeline_development/SipprModelling/rMLST")
path = os.getcwd()

os.chdir("%s/reference" % path)

referenceFile = glob.glob("*.fa*")
references = ["%s/reference/" % path + fastaFile for fastaFile in referenceFile]
#reference = "Escherichia_coli_O157_H7_str_Sakai.fas"

os.chdir("%s/targets" % path)
targets = glob.glob("*.fa")

outPath = "%s/outputs" % path


def make_path(inPath):
    """from: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary \
    does what is indicated by the URL"""
    try:
        os.makedirs(inPath)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def indexTargetsProcesses():
    print '\nIndexing targets'
    indexTargetArgs = []
    if __name__ == '__main__':
        indexTargetsPool = Pool()
        indexTargetArgs.append((target, size))
        # Initialise the pool of processes - it defaults to the number of processors
        indexTargetsPool.map(indexTargets, indexTargetArgs)


def indexTargets():
    """Performs smalt index on the targets using the range of k-mers stored in the variable kmer"""
    print '\nIndexing targets'
    for target in targets:
        for size in kmer:
            filename = target.split('.')[0]
            # Create a new path to be created (if necessary) for the generation of the range of k-mers
            indexPath = "%s/targets/%s/%s_%s" % (path, filename, filename, size)
            # Call the make_path function to make folders as necessary
            make_path(indexPath)
            indexFileSMI = "%s.smi" % filename
            indexFileSMA = "%s.sma" % filename
            if not os.path.isfile("%s/%s" % (indexPath, indexFileSMI)):
                indexCommand = "smalt index -k %s -s 1 %s %s/targets/%s" % (size, filename, path, target)
                subprocess.call(indexCommand, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
                shutil.move(indexFileSMI, indexPath)
                shutil.move(indexFileSMA, indexPath)
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
        # uses kmer, targets, readLength, foldCoverage
        for rLength in readLength:
            for fCov in foldCoverage:
                for target in targets:
                    for size in kmer:
                        mappingProcessesArgs.append((rLength, fCov, target, size))
        mappingProcessesPool.map(mapping, mappingProcessesArgs)



def mapping((rLength, fCov, target, size)):
    """Performs the mapping of the simulated reads to the targets"""
    filename = target.split('.')[0]
    megaName = "rL%s_fC%s_%s_kmer%s" % (rLength, fCov, filename, size)
    filePath = "%s/tmp/rL%s/rL%s_fC%s" % (path, rLength, rLength, fCov)
    newPath = "%s/%s" % (filePath, megaName)
    make_path(newPath)
    targetPath = "%s/targets/%s/%s_%s" % (path, filename, filename, size)
    if not os.path.isfile("%s/%s.bam" % (newPath, megaName)):
        smaltMap = "smalt map -o %s/%s.bam -f bam -x %s/%s %s/%s_%s.fq" \
                   % (newPath, megaName, targetPath, filename, filePath, rLength, fCov)
        subprocess.call(smaltMap, shell=True, stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
        sys.stdout.write('.')
    else:
        sys.stdout.write('.')