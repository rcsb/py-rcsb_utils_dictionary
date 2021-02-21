##
# File:    TargetInteractionProviderTests.py
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

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import time
import unittest


from rcsb.utils.dictionary.TargetInteractionProvider import TargetInteractionProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class TargetInteractionProviderTests(unittest.TestCase):
    skipFlag = True

    def setUp(self):
        mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__cachePath = os.path.join(TOPDIR, "CACHE")
        configPath = os.path.join(HERE, "test-data", "stash-config-example.yml")
        self.__configName = "site_info_configuration"
        self.__startTime = time.time()
        logger.info("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        #
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=self.__configName, mockTopPath=mockTopPath)

    def tearDown(self):
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testTargetInteractionProviderBootstrap(self):
        try:
            tiP = TargetInteractionProvider(self.__cfgOb, self.__configName, self.__cachePath, useCache=False)
            ok = tiP.generate(distLimit=5.0, updateOnly=False, fmt="json", indent=3)
            self.assertTrue(ok)
            ok = tiP.generate(distLimit=5.0, updateOnly=True, fmt="json", indent=3)
            self.assertTrue(ok)
            ok = tiP.reload()
            self.assertTrue(ok)
            ok = tiP.testCache(minCount=85)
            self.assertTrue(ok)
            #
            tiP = TargetInteractionProvider(self.__cfgOb, self.__configName, self.__cachePath, useCache=True)
            ok = tiP.testCache(minCount=85)
            self.assertTrue(ok)
            for entryId in tiP.getEntries():
                intD = tiP.getInteractions(entryId)
                for asymId, neighborL in intD.items():
                    tD = {}
                    for neighbor in neighborL:
                        logger.debug("%s %s %r", entryId, asymId, neighbor)
                        tD.setdefault(asymId, set()).add((neighbor.partnerEntityId, neighbor.partnerAsymId, neighbor.connectType))
                    if len(tD) > 1:
                        logger.info("%s %s (%d) %r", entryId, asymId, len(tD), tD)
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skipIf(skipFlag, "Requires internal configuration")
    def testStashRemote(self):
        try:
            #
            tiP = TargetInteractionProvider(self.__cfgOb, self.__configName, self.__cachePath, useCache=True)
            ok = tiP.testCache()
            self.assertTrue(ok)
            ok = tiP.toStash()
            self.assertTrue(ok)
            #
            ok = tiP.fromStash()
            self.assertTrue(ok)
            ok = tiP.reload()
            self.assertTrue(ok)
            #
        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def targetInteractionSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(TargetInteractionProviderTests("testTargetInteractionProviderBootstrap"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = targetInteractionSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
