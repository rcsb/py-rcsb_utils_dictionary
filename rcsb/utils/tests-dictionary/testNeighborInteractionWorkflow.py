##
# File:    NeighborInteractionWorkflowTests.py
# Author:  J. Westbrook
# Date:    21-Feb-2021
#
# Update:
#
#
##
"""
Tests for ligand and target interaction generation/update workflow
"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import time
import unittest

from rcsb.utils.dictionary.NeighborInteractionWorkflow import NeighborInteractionWorkflow

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class NeighborInteractionWorkflowTests(unittest.TestCase):
    skipFlag = True

    def setUp(self):
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__configPath = os.path.join(HERE, "test-data", "stash-config-example.yml")
        self.__configName = "site_info_configuration"
        self.__useCache = False
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testTargetInteractionUpdate(self):
        try:
            tiWf = NeighborInteractionWorkflow(
                configPath=self.__configPath,
                configName=self.__configName,
                cachePath=self.__cachePath,
                useCache=self.__useCache,
                mockTopPath=self.__mockTopPath,
                numProc=2,
            )
            ok = tiWf.update(incremental=False)
            self.assertTrue(ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skipIf(skipFlag, "Internal test")
    def testTargetInteractionRestore(self):
        try:
            tiWf = NeighborInteractionWorkflow(
                configPath=self.__configPath,
                configName=self.__configName,
                cachePath=self.__cachePath,
                useCache=False,
                mockTopPath=self.__mockTopPath,
                numProc=2,
            )
            ok = tiWf.restore()
            self.assertTrue(ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def targetInteractionSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(NeighborInteractionWorkflowTests("testTargetInteractionUpdate"))
    suiteSelect.addTest(NeighborInteractionWorkflowTests("testTargetInteractionRestore"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = targetInteractionSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
