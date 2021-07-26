# File:    DictmethodResourceProviderTests.py
# Author:  J. Westbrook
# Date:    17-Jul-2021
# Version: 0.001
#
# Update:
#
##
"""
Tests for applying dictionary methods defined as references to helper plugin methods .

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
from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


class DictmethodResourceProviderTests(unittest.TestCase):
    skipFull = True

    def setUp(self):
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

    def testResourceCache(self):
        rP = DictMethodResourceProvider(self.__cfgOb, configName=self.__configName, cachePath=self.__cachePath)
        ok = rP.cacheResources(useCache=True, doRestore=True, useStash=False, useGit=True)
        self.assertTrue(ok)

    @unittest.skipIf(skipFull, "Large test")
    def testBuildResourceCache(self):
        rP = DictMethodResourceProvider(self.__cfgOb, configName=self.__configName, cachePath=self.__cachePath)
        ok = rP.cacheResources(useCache=False, doBackup=False, useStash=False, useGit=True)
        self.assertTrue(ok)

    def testGetResource(self):
        resourceName = "SiftsSummaryProvider instance"
        rP = DictMethodResourceProvider(self.__cfgOb, configName=self.__configName, cachePath=self.__cachePath)
        obj = rP.getResource(resourceName, useCache=False, default=None, doBackup=False, useStash=False, useGit=True)
        self.assertTrue(obj is not None)

    @unittest.skip("Troubleshooting test")
    def testSyncResourceCache(self):
        configName = "site_info_remote_configuration"
        cfgOb = ConfigUtil(configPath=self.__configPath, defaultSectionName=configName, mockTopPath=self.__mockTopPath)
        for providerName in [
            "GlycanProvider instance",
            "DrugBankTargetCofactorProvider instance",
            "ChEMBLTargetCofactorProvider instance",
            "PharosTargetCofactorProvider instance",
            "CARDTargetFeatureProvider instance",
            "IMGTTargetFeatureProvider instance",
            "SAbDabTargetFeatureProvider instance",
        ]:
            rP = DictMethodResourceProvider(self.__cfgOb, configName=self.__configName, cachePath=self.__cachePath)
            ok = rP.syncCache(providerName, cfgOb, configName, self.__cachePath, remotePrefix=None, sourceCache="stash")
            self.assertTrue(ok)


def dictResourceCacheSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DictmethodResourceProviderTests("testResourceCache"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = dictResourceCacheSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
