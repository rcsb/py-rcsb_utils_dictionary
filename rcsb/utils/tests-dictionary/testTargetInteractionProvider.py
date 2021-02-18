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
import unittest

from rcsb.utils.dictionary.TargetInteractionProvider import TargetInteractionProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class TargetInteractionProviderTests(unittest.TestCase):
    skipFlag = True

    def setUp(self):
        mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        # configPath = os.path.join(mockTopPath, "config", "dbload-setup-example.yml")
        configPath = os.path.join(HERE, "test-data", "stash-config-example.yml")
        self.__configName = "site_info_configuration"
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=self.__configName, mockTopPath=mockTopPath)
        self.__dataPath = os.path.join(HERE, "test-data")
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        self.__stashUrl = None
        self.__stashRemotePath = os.path.join(self.__cachePath, "stash-remote")

    def tearDown(self):
        pass

    def testTargetInteractionProviderBootstrap(self):
        tiP = TargetInteractionProvider(self.__cfgOb, self.__configName, self.__cachePath, useCache=False)
        ok = tiP.generate(distLimit=5.0, fmt="json", indent=3)
        self.assertTrue(ok)
        ok = tiP.reload()
        self.assertTrue(ok)
        ok = tiP.testCache(minCount=85)
        self.assertTrue(ok)
        #
        tiP = TargetInteractionProvider(self.__cfgOb, self.__configName, self.__cachePath, useCache=True)
        ok = tiP.testCache(minCount=85)
        self.assertTrue(ok)

    @unittest.skipIf(skipFlag, "Private test")
    def testStashRemote(self):
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


def targetInteractionSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(TargetInteractionProviderTests("testTargetInteractionProviderBootstrap"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = targetInteractionSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
