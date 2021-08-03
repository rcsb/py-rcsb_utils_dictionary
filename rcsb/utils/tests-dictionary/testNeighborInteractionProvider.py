##
# File:    NeighborInteractionProviderTests.py
# Author:  J. Westbrook
# Date:    18-Feb-2021
#
# Update:
#
#
##
"""
Tests for generators and accessors for non-polymer instance target interactions

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


from rcsb.utils.dictionary.NeighborInteractionProvider import NeighborInteractionProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class NeighborInteractionProviderTests(unittest.TestCase):
    skipFlag = platform.system() != "Darwin"
    buildTestingCache = True

    def setUp(self):
        mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        configPath = os.path.join(mockTopPath, "config", "dbload-setup-example.yml")
        self.__configName = "site_info_configuration"
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        self.__fileLimit = None if platform.system() == "Darwin" else 5
        #
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=self.__configName, mockTopPath=mockTopPath)

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 1.0e6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testNeighborInteractionProviderBootstrap(self):
        """Test case: generate and load neighbor and occupancy data"""
        try:
            niP = NeighborInteractionProvider(self.__cachePath, useCache=False, cfgOb=self.__cfgOb, configName=self.__configName, numProc=2, fileLimit=self.__fileLimit)
            ok = niP.generate(distLimit=5.0, updateOnly=False, fmt="pickle")
            self.assertTrue(ok)
            ok = niP.generate(distLimit=5.0, updateOnly=True, fmt="pickle")
            self.assertTrue(ok)
            ok = niP.reload()
            self.assertTrue(ok)
            ok = niP.testCache(minCount=self.__fileLimit if self.__fileLimit else 30)
            self.assertTrue(ok)
            #
            niP = NeighborInteractionProvider(self.__cachePath, useCache=True, cfgOb=self.__cfgOb, configName=self.__configName, numProc=2)
            ok = niP.testCache(minCount=self.__fileLimit if self.__fileLimit else 30)
            self.assertTrue(ok)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skipIf(skipFlag, "Long test")
    def testStashRemote(self):
        try:
            niP = NeighborInteractionProvider(self.__cachePath, useCache=True, cfgOb=self.__cfgOb, configName=self.__configName, numProc=2, fileLimit=self.__fileLimit)
            ok = niP.testCache()
            self.assertTrue(ok)
            ok = niP.restore(self.__cfgOb, self.__configName, remotePrefix=None, useGit=True, useStash=True)
            self.assertTrue(ok)
            ok = niP.reload()
            self.assertTrue(ok)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def targetInteractionSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(NeighborInteractionProviderTests("testNeighborInteractionProviderBootstrap"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = targetInteractionSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
