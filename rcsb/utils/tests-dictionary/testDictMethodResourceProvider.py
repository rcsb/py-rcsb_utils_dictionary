# File:    DictmethodResourceProviderTests.py
# Author:  J. Westbrook
# Date:    17-Jul-2021
# Version: 0.001
#
# Update:
#
##
"""
Tests for caching resources used by dictionary methods.

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
import unittest

from rcsb.utils.dictionary.DictMethodResourceProvider import DictMethodResourceProvider
from rcsb.utils.dictionary.NeighborInteractionProvider import NeighborInteractionProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


class DictmethodResourceProviderTests(unittest.TestCase):
    skipFull = platform.system() != "Darwin"
    buildTestingCache = False

    def setUp(self):
        isMac = platform.system() == "Darwin"
        self.__excludeType = None if isMac else "optional"
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__configPath = os.path.join(self.__mockTopPath, "config", "dbload-setup-example.yml")
        configName = "site_info_configuration"
        self.__configName = configName
        self.__cfgOb = ConfigUtil(configPath=self.__configPath, defaultSectionName=configName, mockTopPath=self.__mockTopPath)
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        #
        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 1.0e6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testGetResource(self):
        """Test restore a single packages from  git stash storage"""
        resourceName = "SiftsSummaryProvider instance"
        rP = DictMethodResourceProvider(
            self.__cfgOb,
            configName=self.__configName,
            cachePath=self.__cachePath,
            restoreUseStash=False,
            restoreUseGit=True,
            providerTypeExclude=self.__excludeType,
        )
        obj = rP.getResource(resourceName, useCache=False, default=None, doBackup=False, useStash=False, useGit=True)
        self.assertTrue(obj is not None)

    def testResourceCache(self):
        """Test restore cache from git stash storage"""
        rP = DictMethodResourceProvider(
            self.__cfgOb,
            configName=self.__configName,
            cachePath=self.__cachePath,
            restoreUseStash=False,
            restoreUseGit=True,
            providerTypeExclude=self.__excludeType,
        )
        ok = rP.cacheResources(useCache=True, doRestore=True)
        self.assertTrue(ok)

    # ---- Maintenance tests -----
    @unittest.skipUnless(buildTestingCache, "Maintenance task to construct testing cache")
    def testBuildResourceCache(self):
        """Build the testing cache resources from scratch and update git storage"""
        rP = DictMethodResourceProvider(
            self.__cfgOb,
            configName=self.__configName,
            cachePath=self.__cachePath,
            restoreUseStash=False,
            restoreUseGit=False,
            providerTypeExclude=None,
        )
        ok = rP.cacheResources(useCache=False, doBackup=True, useStash=False, useGit=True, clearCache=True, providerSelect="buildable")
        self.assertTrue(ok)
        # ---
        niP = NeighborInteractionProvider(self.__cachePath, useCache=False, cfgOb=self.__cfgOb, configName=self.__configName, numProc=2, fileLimit=50)
        ok = niP.generate(distLimit=5.0, updateOnly=False, fmt="pickle")
        ok = niP.reload()
        self.assertTrue(ok)
        ok = niP.testCache(minCount=25)
        self.assertTrue(ok)
        if ok:
            okB = niP.backup(self.__cfgOb, self.__configName, remotePrefix=None, useGit=True, useStash=False)
            self.assertTrue(okB)
        # --- Sync the other non-buildable resource data
        configName = "site_info_remote_configuration"
        cfgOb = ConfigUtil(configPath=self.__configPath, defaultSectionName=configName, mockTopPath=self.__mockTopPath)
        rP = DictMethodResourceProvider(cfgOb, configName=configName, cachePath=self.__cachePath, restoreUseStash=True, restoreUseGit=True)
        for providerName in [
            "CARDTargetFeatureProvider instance",
            "ChEMBLTargetCofactorProvider instance",
            "DrugBankTargetCofactorProvider instance",
            "GlycanProvider instance",
            "IMGTTargetFeatureProvider instance",
            # "NeighborInteractionProvider instance",
            "PharosProvider instance",
            "PharosTargetCofactorProvider instance",
            "PubChemProvider instance",
            "SAbDabTargetFeatureProvider instance",
        ]:
            # ok = rP.syncCache(providerName, cfgOb, configName, self.__cachePath, remotePrefix=None, sourceCache="stash")
            obj = rP.getResource(providerName, useCache=False, default=None, doRestore=True, doBackup=False, useStash=False, useGit=False)
            self.assertTrue(obj is not None)

    @unittest.skip("Maintenance test")
    def testCacheResourceToGit(self):
        """Update the cache for an individual buildable resource"""
        resourceName = "Scop2Provider instance"
        rP = DictMethodResourceProvider(
            self.__cfgOb,
            configName=self.__configName,
            cachePath=self.__cachePath,
            restoreUseStash=True,
            restoreUseGit=False,
            providerTypeExclude=self.__excludeType,
        )
        obj = rP.getResource(resourceName, useCache=False, default=None, doBackup=True, useStash=False, useGit=True)
        self.assertTrue(obj is not None)


def dictResourceCacheSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DictmethodResourceProviderTests("testResourceCache"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = dictResourceCacheSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
