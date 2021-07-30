# File:    DictmethodResourceCacheWorkflow.py
# Author:  J. Westbrook
# Date:    28-Jul-2021
# Version: 0.001
#
# Update:
#
##
"""
Workflow to rebuild and stash "buildable" cache resources.
"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time

from rcsb.utils.dictionary.DictMethodResourceProvider import DictMethodResourceProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.io.FileUtil import FileUtil

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)


class DictMethodResourceCacheWorkflow(object):
    def __init__(self, workPath):
        self.__cachePath = os.path.join(workPath, "workflow-output", "CACHE")
        self.__configPath = os.path.join(workPath, "exdb-config-example.yml")
        configName = "site_info_remote_configuration"
        self.__configName = configName
        self.__cfgOb = ConfigUtil(configPath=self.__configPath, defaultSectionName=configName)
        self.__startTime = time.time()
        logger.debug("Starting at %s", time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def reportUsage(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 1.0e6, unitS)
        endTime = time.time()
        logger.info("Completed at %s (%.4f seconds)", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def buildResourceCache(self):
        rP = DictMethodResourceProvider(
            self.__cfgOb,
            configName=self.__configName,
            cachePath=self.__cachePath,
            restoreUseStash=True,
            restoreUseGit=True,
            providerTypeExclude=None,
        )
        ok = rP.cacheResources(useCache=False, doBackup=True, useStash=True, useGit=True)
        logger.info("Cache rebuild status (%r)", ok)

    def testRecoverCacheFromStash(self):
        # remove any cache directory

        fU = FileUtil()
        cPth = os.path.abspath(self.__cachePath)
        if fU.exists(cPth):
            fU.remove(cPth)
        else:
            logger.info("Bad cache path %r", cPth)
        #
        rP = DictMethodResourceProvider(
            self.__cfgOb,
            configName=self.__configName,
            cachePath=self.__cachePath,
            restoreUseStash=True,
            restoreUseGit=False,
            providerTypeExclude=None,
        )
        ok = rP.cacheResources(useCache=True, doRestore=True, doBackup=False)
        logger.info(">>> Stash recovery test status (%r)", ok)

    def testRecoverCacheFromGit(self):
        # remove any cache directory
        fU = FileUtil()
        cPth = os.path.abspath(self.__cachePath)
        if fU.exists(cPth):
            fU.remove(cPth)
        else:
            logger.info("Bad cache path %r", cPth)
        #
        rP = DictMethodResourceProvider(
            self.__cfgOb,
            configName=self.__configName,
            cachePath=self.__cachePath,
            restoreUseStash=False,
            restoreUseGit=True,
            providerTypeExclude=None,
        )
        ok = rP.cacheResources(useCache=True, doRestore=True, doBackup=False)
        logger.info(">>> Git recovery test status (%r)", ok)

    def syncResourceCache(self):
        for providerName in [
            "NeighborInteractionProvider instance",
            "GlycanProvider instance",
            "DrugBankTargetCofactorProvider instance",
            "ChEMBLTargetCofactorProvider instance",
            "PharosTargetCofactorProvider instance",
            "CARDTargetFeatureProvider instance",
            "IMGTTargetFeatureProvider instance",
            "SAbDabTargetFeatureProvider instance",
        ]:
            rP = DictMethodResourceProvider(self.__cfgOb, configName=self.__configName, cachePath=self.__cachePath)
            ok = rP.syncCache(providerName, self.__cfgOb, self.__configName, self.__cachePath, remotePrefix=None, sourceCache="stash")
            logger.info("Sync %r status (%r)", providerName, ok)


if __name__ == "__main__":
    dmrWf = DictMethodResourceCacheWorkflow(workPath="./")
    dmrWf.buildResourceCache()
    dmrWf.syncResourceCache()
