##
#  File:           TargetInteractionProvider.py
#  Date:           17-Feb-2021 jdw
#
#  Updated:
#
##
"""
Generators and accessors for non-polymer instance target interactions

"""

import logging
import os.path
import time

from rcsb.utils.dictionary import __version__
from rcsb.utils.dictionary.DictMethodCommonUtils import DictMethodCommonUtils, LigandTargetInstance
from rcsb.utils.io.MarshalUtil import MarshalUtil
from rcsb.utils.io.StashUtil import StashUtil
from rcsb.utils.multiproc.MultiProcUtil import MultiProcUtil
from rcsb.utils.repository.RepositoryProvider import RepositoryProvider

logger = logging.getLogger(__name__)


class TargetInteractionWorker(object):
    """A skeleton class that implements the interface expected by the multiprocessing
    for calculating non-polymer instance target interactions --
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
                    instD, countD = self.__buildTargetInteractions(procName, dataContainer, distLimit)
                    retList.append((entryId, instD, countD))
            #
            successList = sorted(set(dataList) - set(failList))
            if failList:
                logger.info("%s returns %d definitions with failures: %r", procName, len(failList), failList)

            logger.debug("%s built target interactions for %d/%d entries failures %d", procName, len(retList), len(dataList), len(failList))
        except Exception as e:
            logger.exception("Failing %s for %d data items %s", procName, len(dataList), str(e))
        #
        return successList, retList, diagList

    def __buildTargetInteractions(self, procName, dataContainer, distLimit):
        """Internal method return a dictionary target interactions."""
        rD = {}
        try:
            rD, countD = self.__commonU.getNonpolymerInstanceNeighborInfo(dataContainer, distLimit=distLimit)
        except Exception as e:
            logger.exception("%s failing with %s", procName, str(e))
        return rD, countD


class TargetInteractionProvider(object):
    """Generators and accessors for non-polymer instance target interactions."""

    def __init__(self, cfgOb, configName, cachePath, **kwargs):
        #
        self.__version = __version__
        self.__cfgOb = cfgOb
        self.__configName = configName
        self.__cachePath = cachePath
        self.__fileLimit = kwargs.get("fileLimit", None)
        self.__dirPath = os.path.join(cachePath, "target-interactions")
        self.__numProc = kwargs.get("numProc", 2)
        self.__chunkSize = kwargs.get("chunkSize", 10)
        useCache = kwargs.get("useCache", True)
        #
        #  - Configuration for stash services -
        #    Local target directory name to be stashed.  (subdir of dirPath)
        #
        self.__stashDir = "nonpolymer-target-interactions"
        #
        self.__mU = MarshalUtil(workPath=self.__dirPath)
        self.__rpP = RepositoryProvider(cfgOb=self.__cfgOb, numProc=self.__numProc, fileLimit=self.__fileLimit, cachePath=self.__cachePath)
        self.__targetD = self.__reload(fmt="json", useCache=useCache)
        #

    def testCache(self, minCount=0):
        try:
            if minCount == 0:
                return True
            if self.__targetD and minCount and len(self.__targetD["interactions"]) >= minCount:
                logger.info("Target interactions (%d) created %r version %r", len(self.__targetD["interactions"]), self.__targetD["created"], self.__targetD["version"])
                return True
        except Exception:
            pass
        return False

    def getInteractions(self, entryId):
        """Return the target interactions for the non-polymer instances for the input entry.

        Args:
            entryId (str): entry identifier

        Returns:
            (dict): {asymId: [LigandTargetInstance(), LigandTargetInstance(), ...]}
        """
        try:
            return self.__targetD["interactions"][entryId.upper()]
        except Exception:
            pass
        return {}

    def getAtomCounts(self, entryId):
        """Return the non-polymer instance occupancy weighted atome counts for the input entry.

        Args:
            entryId (str): entry identifier

        Returns:
            (dict): {asymId: {'FL': count, 'altA': occ*count, 'altB': occ*count, ... }}
        """
        try:
            return self.__targetD["atomCounts"][entryId.upper()]
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
            return entryId in self.__targetD["interactions"]
        except Exception:
            pass
        return False

    def getEntries(self):
        """Return a list of entry identifier for which non-polymer instance target interactions are stored.

        Returns:
            (list): [entryId, entryId, ... ]
        """
        try:
            return list(self.__targetD["interactions"].keys())
        except Exception:
            pass
        return []

    def generate(self, distLimit=5.0, updateOnly=False, fmt="json", indent=0):
        """Generate and export non-polymer target interactions for all of the structures in the repository.

        Args:
            distLimit (float, optional): interaction distance. Defaults to 5.0.
            updateOnly (bool):  only calculate interactions for new entries.  Defaults to False.
            fmt (str, optional): export file format. Defaults to "json".
            indent (int, optional): json format indent. Defaults to 0.

        Returns:
            bool: True for success or False otherwise
        """
        ok = False
        try:
            tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
            tD, cD = self.__calculateNeighbors(distLimit=distLimit, numProc=self.__numProc, chunkSize=self.__chunkSize, updateOnly=updateOnly)
            self.__targetD = {"version": self.__version, "created": tS, "interactions": tD, "atomCounts": cD}
            kwargs = {"indent": indent} if fmt == "json" else {}
            targetFilePath = self.__getTargetFilePath(fmt=fmt)
            ok = self.__mU.doExport(targetFilePath, self.__targetD, fmt=fmt, **kwargs)
            logger.info("Wrote %r status %r", targetFilePath, ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def reload(self):
        self.__targetD = self.__reload(fmt="json", useCache=True)
        return self.__targetD is not None

    def __reload(self, fmt="json", useCache=True):
        """Reload from the current cache file."""
        try:
            targetFilePath = self.__getTargetFilePath(fmt=fmt)
            tS = time.strftime("%Y %m %d %H:%M:%S", time.localtime())
            targetD = {"version": self.__version, "created": tS, "interactions": {}, "atomCounts": {}}
            logger.debug("useCache %r targetFilePath %r", useCache, targetFilePath)
            #
            if useCache and self.__mU.exists(targetFilePath):
                targetD = self.__mU.doImport(targetFilePath, fmt=fmt)
                for _, neighborD in targetD["interactions"].items():
                    for asymId in neighborD:
                        neighborD[asymId] = [LigandTargetInstance(*neighbor) for neighbor in neighborD[asymId]]
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return targetD

    def __getTargetFilePath(self, fmt="json"):
        pth = os.path.join(self.__dirPath, "nonpolymer-target-interactions", "inst-target-interactions." + fmt)
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

        rD = {}
        cD = {}
        contentType = "pdbx"
        mergeContent = None
        exD = None
        if updateOnly:
            if not self.testCache(minCount=10):
                ok = self.reload()
            exD = {k: True for k in self.getEntries()}
            rD = self.__targetD["interactions"] if "interactions" in self.__targetD else {}
            cD = self.__targetD["atomCounts"] if "atomCounts" in self.__targetD else {}
        #
        locatorObjList = self.__rpP.getLocatorObjList(contentType=contentType, mergeContentTypes=mergeContent, excludeIds=exD)
        logger.info("Starting with %d numProc %d updateOnly (%r)", len(locatorObjList), self.__numProc, updateOnly)
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
        for (entryId, tD, qD) in resultList[0]:
            rD[entryId] = tD
            cD[entryId] = qD
        #
        logger.info("Completed with multi-proc status %r failures %r total entry interactions (%d) atom counts (%d)", ok, len(failList), len(rD), len(cD))
        return rD, cD

    def toStash(self):
        ok = False
        try:
            userName = self.__cfgOb.get("_STASH_AUTH_USERNAME", sectionName=self.__configName)
            password = self.__cfgOb.get("_STASH_AUTH_PASSWORD", sectionName=self.__configName)
            basePath = self.__cfgOb.get("_STASH_SERVER_BASE_PATH", sectionName=self.__configName)
            url = self.__cfgOb.get("STASH_SERVER_URL", sectionName=self.__configName)
            urlFallBack = self.__cfgOb.get("STASH_SERVER_FALLBACK_URL", sectionName=self.__configName)
            ok = self.__toStash(url, basePath, userName=userName, password=password)
            ok = self.__toStash(urlFallBack, basePath, userName=userName, password=password)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return ok

    def __toStash(self, url, stashRemoteDirPath, userName=None, password=None, remoteStashPrefix=None):
        """Copy tar and gzipped bundled cache data to remote server/location.

        Args:
            url (str): server URL (e.g. sftp://hostname.domain) None for local host
            stashRemoteDirPath (str): path to target directory on remote server
            userName (str, optional): server username. Defaults to None.
            password (str, optional): server password. Defaults to None.
            remoteStashPrefix (str, optional): channel prefix. Defaults to None.

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            stU = StashUtil(os.path.join(self.__dirPath, "stash"), "nonpolymer-target-interactions")
            ok = stU.makeBundle(self.__dirPath, [self.__stashDir])
            if ok:
                ok = stU.storeBundle(url, stashRemoteDirPath, remoteStashPrefix=remoteStashPrefix, userName=userName, password=password)
        except Exception as e:
            logger.error("Failing with url %r stashDirPath %r: %s", url, stashRemoteDirPath, str(e))
        return ok

    def fromStash(self):
        try:
            minCount = 10
            userName = self.__cfgOb.get("_STASH_AUTH_USERNAME", sectionName=self.__configName)
            password = self.__cfgOb.get("_STASH_AUTH_PASSWORD", sectionName=self.__configName)
            basePath = self.__cfgOb.get("_STASH_SERVER_BASE_PATH", sectionName=self.__configName)
            url = self.__cfgOb.get("STASH_SERVER_URL", sectionName=self.__configName)
            #
            ok = self.__fromStash(url, basePath, userName=userName, password=password)
            ok = self.reload()
            ok = self.testCache(minCount=minCount)
            if not ok:
                urlFallBack = self.__cfgOb.get("STASH_SERVER_FALLBACK_URL", sectionName=self.__configName)
                ok = self.__fromStash(urlFallBack, basePath, userName=userName, password=password)
                ok = self.testCache(minCount=minCount)
                ok = self.reload()
        except Exception as e:
            logger.exception("Failing with %s", str(e))

        return ok

    def __fromStash(self, url, stashRemoteDirPath, userName=None, password=None, remoteStashPrefix=None):
        """Restore local cache from a tar and gzipped bundle to fetched from a remote server/location.

        Args:
            url (str): server URL (e.g. sftp://hostname.domain) None for local host
            stashRemoteDirPath (str): path to target directory on remote server
            userName (str, optional): server username. Defaults to None.
            password (str, optional): server password. Defaults to None.
            remoteStashPrefix (str, optional): channel prefix. Defaults to None.

        Returns:
            (bool): True for success or False otherwise
        """
        ok = False
        try:
            stU = StashUtil(os.path.join(self.__dirPath, "stash"), "nonpolymer-target-interactions")
            ok = stU.fetchBundle(self.__dirPath, url, stashRemoteDirPath, remoteStashPrefix=remoteStashPrefix, userName=userName, password=password)
        except Exception as e:
            logger.error("Failing with url %r stashDirPath %r: %s", url, stashRemoteDirPath, str(e))
        return ok
