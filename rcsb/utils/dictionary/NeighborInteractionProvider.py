##
#  File:           NeighborInteractionProvider.py
#  Date:           17-Feb-2021 jdw
#
#  Updated:
#
##
"""
Generators and accessors for target and ligand neighbor interactions

"""

import logging
import os.path
import time

from rcsb.utils.dictionary import __version__
from rcsb.utils.dictionary.DictMethodCommonUtils import DictMethodCommonUtils, LigandTargetInstance
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashableBase import StashableBase
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil
from rcsb.utils.repository.RepositoryProvider import RepositoryProvider

logger = logging.getLogger(__name__)


class TargetInteractionWorker(object):
    """A skeleton class that implements the interface expected by the multiprocessing
    for calculating non-polymer instance ligand and target neighbor interactions --
    """

    def __init__(self, repositoryProviderObj, **kwargs):
        self.__rpP = repositoryProviderObj
        _ = kwargs
        self.__commonU = DictMethodCommonUtils()

    def build(self, dataList, procName, optionsD, workingDir):
        """Enumerate non-polymer instance target interactions."""
        _ = workingDir
        successList = []
        failList = []
        retList = []
        diagList = []
        #
        try:
            distLimit = optionsD.get("distLimit", 6.0)
            for locatorObj in dataList:
                dataContainerList = self.__rpP.getContainerList([locatorObj])
                for dataContainer in dataContainerList:
                    entryId = dataContainer.getName()
                    rD = self.__getNeighborInfo(procName, dataContainer, distLimit)
                    retList.append((entryId, rD))
            #
            successList = list(set(dataList) - set(failList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)

            logger.debug("%s built target interactions for %d/%d entries failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))
        #
        return successList, retList, diagList

    def __getNeighborInfo(self, procName, dataContainer, distLimit):
        """Internal method return a dictionary target interactions."""
        rD = {}
        try:
            rD = self.__commonU.getNeighborInfo(dataContainer, distLimit=distLimit)

        except Exception as e:
            logger.exception("%s failing with %s", procName, str(e))
        return rD


class NeighborInteractionProvider(StashableBase):
    """Generators and accessors for non-polymer instance target interactions."""

    def __init__(self, cachePath, useCache=False, cfgOb=None, configName=None, **kwargs):
        #
        self.__version = __version__
        self.__cfgOb = cfgOb
        self.__configName = configName
        self.__cachePath = cachePath
        self.__fileLimit = kwargs.get("fileLimit", None)
        dirName = "neighbor-interactions"
        self.__dirPath = os.path.join(cachePath, dirName)
        super(NeighborInteractionProvider, self).__init__(self.__cachePath, [dirName])
        self.__numProc = kwargs.get("numProc", 2)
        self.__chunkSize = kwargs.get("chunkSize", 10)
        # useCache = kwargs.get("useCache", True)
        #
        #  - Configuration for stash services -
        #    Local target directory name to be stashed.  (subdir of dirPath)
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__rpP = RepositoryProvider(cfgOb=self.__cfgOb, numProc=self.__numProc, fileLimit=self.__fileLimit, cachePath=self.__cachePath)
        self.__neighborD = self.__reload(fmt="pickle", useCache=useCache)
        #

    def testCache(self, minCount=1):
        try:
            if minCount == 0:
                return True
            if self.__neighborD and minCount and "entries" in self.__neighborD and len(self.__neighborD["entries"]) >= minCount:
                logger.info("Target neighbor data for (%d) entries created %r version %r", len(self.__neighborD["entries"]), self.__neighborD["created"], self.__neighborD["version"])
                return True
        except Exception:
            pass
        return False

    def getLigandNeighborIndex(self, entryId):
        """Return the target neighbors for the non-polymer instances for the input entry.

        Args:
            entryId (str): entry identifier

        Returns:
            (dict): {ligandAsymId: {(targetAsymId, targetAuthSeqId): nnIndex1, (): nnIndex2}
        """
        try:
            return self.__neighborD["entries"][entryId.upper()]["ligandNeighborIndexD"]
        except Exception:
            pass
        return {}

    def getTargetNeighborIndex(self, entryId):
        """Return the ligand neighbors for the polymer or branched entity instances in the input entry.

        Args:
            entryId (str): entry identifier

        Returns:
            (dict): {(targetAsymId, targetAuthSeqId): {(ligandAsymId): nnIndex1, (): nnIndex2}

        """
        try:
            return self.__neighborD["entries"][entryId.upper()]["targetNeighborIndexD"]
        except Exception:
            pass
        return {}

    def getNearestNeighborList(self, entryId):
        """Return the list of nearest neighbors for the entry.

        Args:
            entryId (str): entry identifier

        Returns:
            list: [LigandTargetInstance(), ...]

        """
        try:
            return self.__neighborD["entries"][entryId.upper()]["nearestNeighbors"]
        except Exception:
            pass
        return []

    def getLigandNeighborBoundState(self, entryId):
        """Return the dictionary of ligand instances with isBound boolean status.

        Args:
            entryId (str): entry identifier

        Returns:
            (dict): {ligandAsymId: True if isBound,  ...  }
        """
        try:
            return self.__neighborD["entries"][entryId.upper()]["ligandIsBoundD"]
        except Exception:
            pass
        return {}

    def getAtomCounts(self, entryId):
        """Return the non-polymer instance atom counts for the input entry (all reported atoms).

        Args:
            entryId (str): entry identifier

        Returns:
            (dict): {asymId: {'FL': count, 'altA': count, 'altB': count, ... }}
        """
        try:
            return self.__neighborD["entries"][entryId.upper()]["ligandAtomCountD"]
        except Exception:
            pass
        return {}

    def getHydrogenAtomCounts(self, entryId):
        """Return the non-polymer instance hydrogen atom counts for the input entry.

        Args:
            entryId (str): entry identifier

        Returns:
            (dict): {asymId: {'FL': count, 'altA': count, 'altB': count, ... }}
        """
        try:
            return self.__neighborD["entries"][entryId.upper()]["ligandHydrogenAtomCountD"]
        except Exception:
            pass
        return {}

    def __getInstanceOccupancySumD(self, entryId):
        """Return instance heavy atom occupancy for the input entry (all reported heavy atoms).

        Args:
            entryId (str): entry identifier

        Returns:
            (dict): {asymId: {'FL': count, 'altA': count, 'altB': count, ... }}
        """
        try:
            return self.__neighborD["entries"][entryId.upper()]["occupancySumD"]
        except Exception:
            pass
        return {}

    def hasEntry(self, entryId):
        """Return if the input entry is stored in the cache of non-polymer instance target interactions.

        Args:
            entryId (str): entry identifier

        Returns:
            (bool): True if entry is in the cache or False otherwise
        """
        try:
            return entryId in self.__neighborD["entries"]
        except Exception:
            pass
        return False

    def getEntries(self):
        """Return a list of entry identifier for which non-polymer instance target interactions are stored.

        Returns:
            (list): [entryId, entryId, ... ]
        """
        try:
            return list(self.__neighborD["entries"].keys())
        except Exception:
            pass
        return []

    def generate(self, distLimit=5.0, updateOnly=False, fmt="pickle", indent=0):
        """Generate and export non-polymer target interactions for all of the structures in the repository.

        Args:
            distLimit (float, optional): interaction distance. Defaults to 5.0.
            updateOnly (bool):  only calculate interactions for new entries.  Defaults to False.
            fmt (str, optional): export file format. Defaults to "pickle".
            indent (int, optional): json format indent. Defaults to 0.

        Returns:
            bool: True for success or False otherwise
        """
        ok = False
        try:
            tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
            tD = self.__calculateNeighbors(distLimit=distLimit, numProc=self.__numProc, chunkSize=self.__chunkSize, updateOnly=updateOnly)
            self.__neighborD = {"version": self.__version, "created": tS, "entries": tD}
            kwargs = {"indent": indent} if fmt == "json" else {"pickleProtocol": 4}
            targetFilePath = self.__getTargetFilePath(fmt=fmt)
            ok = self.__mU.doExport(targetFilePath, self.__neighborD, fmt=fmt, **kwargs)
            logger.info("Wrote %r status %r", targetFilePath, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def reload(self, fmt="pickle"):
        self.__neighborD = self.__reload(fmt=fmt, useCache=True)
        return self.__neighborD is not None

    def __reload(self, fmt="pickle", useCache=True):
        """Reload from the current cache file."""
        try:
            targetFilePath = self.__getTargetFilePath(fmt=fmt)
            tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
            neighborD = {"version": self.__version, "created": tS, "entries": {}}
            logger.debug("useCache %r targetFilePath %r", useCache, targetFilePath)
            #
            if useCache and self.__mU.exists(targetFilePath):
                neighborD = self.__mU.doImport(targetFilePath, fmt=fmt)
                if fmt != "pickle":
                    for _, nD in neighborD["entries"].items():
                        nD["nearestNeighbors"] = [LigandTargetInstance(*neighbor) for neighbor in nD["nearestNeighbors"]]
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return neighborD

    def __getTargetFilePath(self, fmt="pickle"):
        ext = "pic" if fmt == "pickle" else "json"
        pth = os.path.join(self.__dirPath, "neighbor-data." + ext)
        return pth

    def __calculateNeighbors(self, distLimit=5.0, numProc=2, chunkSize=10, updateOnly=False):
        """Calculate non-polymer target interactions for all repository structure files.

        Args:
            distLimit (float, optional): interaction distance limit. Defaults to 5.0.
            numProc (int, optional): number of processes to use. Defaults to 2.
            chunkSize (int, optional): incremental chunk size used for distribute work processes. Defaults to 10.

        Returns:
            (dict): {entryId: {asymId: [TargetLigandInteraction()], ...}, ...}
        """
        contentType = "pdbx"
        mergeContent = None
        rD = {}
        exD = {}
        #
        # updateOnly - will reuse any existing data loaded when this is instantiated
        #              otherwise the cache context is cleared before the calculation.
        if updateOnly:
            exD = {k: True for k in self.getEntries()}
            logger.info("Reusing (%d) entries", len(exD))
            rD = self.__neighborD["entries"] if "entries" in self.__neighborD else {}
        #
        locatorObjList = self.__rpP.getLocatorObjList(contentType=contentType, mergeContentTypes=mergeContent, excludeIds=exD)
        logger.info("Starting with %d entries numProc %d updateOnly (%r)", len(locatorObjList), self.__numProc, updateOnly)
        #
        rWorker = TargetInteractionWorker(self.__rpP)
        mpu = MultiProcUtil(verbose=True)
        optD = {"distLimit": distLimit}
        mpu.setOptions(optD)
        mpu.set(workerObj=rWorker, workerMethod="build")
        ok, failList, resultList, _ = mpu.runMulti(dataList=locatorObjList, numProc=numProc, numResults=1, chunkSize=chunkSize)
        if failList:
            logger.info("Target interaction build failures (%d): %r", len(failList), failList)
        #
        for (entryId, nD) in resultList[0]:
            rD[entryId] = nD
        #
        logger.info("Completed with multi-proc status %r failures %r total entries with data (%d)", ok, len(failList), len(rD))
        return rD

    def convert(self, fmt1="json", fmt2="pickle"):
        #
        targetFilePath = self.__getTargetFilePath(fmt=fmt1)
        self.__neighborD = self.__mU.doImport(targetFilePath, fmt=fmt1)
        #
        targetFilePath = self.__getTargetFilePath(fmt=fmt2)
        ok = self.__mU.doExport(targetFilePath, self.__neighborD, fmt=fmt2, pickleProtocol=4)
        return ok
