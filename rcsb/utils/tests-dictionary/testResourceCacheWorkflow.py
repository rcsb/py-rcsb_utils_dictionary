##
# File:    testResourceCacheWorkflow.py
# Author:  A. Evans
# Date:    3-Mar-2023
#
# Update:
#
#
##
"""
Tests for resource cache workflow for stash
"""

__docformat__ = "google en"
__author__ = "Alicia Evans"
__email__ = "alicia.evans@rcsb.org"
__license__ = "Apache 2.0"

import logging
import os
import platform
import resource
import time
import unittest

from rcsb.utils.dictionary.DictMethodResourceCacheWorkflow import DictMethodResourceCacheWorkflow

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class DictMethodResourceCacheWorkflowTests(unittest.TestCase):
    skipFlag = True

    def setUp(self):
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__configPath = os.path.join(self.__mockTopPath, "config", "dbload-setup-example.yml")
        self.__configName = "site_info_configuration"
        self.__useCache = False
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 1.0e6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    @unittest.skipIf(skipFlag, "Test requires files that are too large for Docker")
    def testBuildResourceCacheStash(self):
        try:
            tiWf = DictMethodResourceCacheWorkflow(
                configPath=self.__configPath,
                configName=self.__configName,
                cachePath=self.__cachePath,
                mockTopPath=self.__mockTopPath,
            )
            ok = tiWf.buildResourceCache()
            self.assertTrue(ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def stashResourcesSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DictMethodResourceCacheWorkflowTests("testBuildResourceCacheStash"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = stashResourcesSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
