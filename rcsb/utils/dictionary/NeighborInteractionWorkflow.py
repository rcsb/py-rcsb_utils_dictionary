##
# File:    NeighborInteractionWorkflow.py
# Author:  J. Westbrook
# Date:    18-Feb-2021
#
# Update:
#  16-Mar-2023 aae  Update configuration to use HERE and CACHE folder
#
##
"""
Worflow for generating and stashing ligand and target neighbor interactions

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os
import time

from rcsb.utils.dictionary.NeighborInteractionProvider import NeighborInteractionProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class NeighborInteractionWorkflow(object):
    def __init__(self, **kwargs):
        """Workflow --  Workflow to rebuild and stash ligand and target neighbor interactions

        Args:
            configPath (str, optional): path to configuration file (default: exdb-config-example.yml)
            configName (str, optional): configuration section name (default: site_info_remote_configuration)
            mockTopPath (str, optional):  mockTopPath is prepended to path configuration options if it specified (default=None)
            workPath (str, optional):  path to working directory (default: HERE)
            cachePath (str, optional):  path to cache directory (default: HERE/CACHE)
            numProc (int, optional):  number processors to include in multiprocessing mode (default: 10)
            chunkSize (int, optional):  incremental chunk size used for distribute work processes (default: 10)
            useCache (bool, optional):  use cached configuration assets (default: False)
            stashRemotePrefix (str, optional): file name prefix (channel) applied to remote stash file artifacts (default: None)
            debugFlag (bool, optional):  sets logger to debug mode (default: False)
        """
        configPath = kwargs.get("configPath", "exdb-config-example.yml")
        logger.info("Configuration file path %s", configPath)
        self.__configName = kwargs.get("configName", "site_info_remote_configuration")
        mockTopPath = kwargs.get("mockTopPath", None)
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=self.__configName, mockTopPath=mockTopPath)
        self.__workPath = kwargs.get("workPath", HERE)
        self.__cachePath = kwargs.get("cachePath", os.path.join(self.__workPath, "CACHE"))
        #
        self.__numProc = kwargs.get("numProc", 10)
        self.__chunkSize = kwargs.get("chunkSize", 10)
        self.__useCache = kwargs.get("useCache", False)
        #
        self.__stashRemotePrefix = kwargs.get("stashRemotePrefix", None)
        #
        self.__debugFlag = kwargs.get("debugFlag", False)
        if self.__debugFlag:
            logger.setLevel(logging.DEBUG)
            self.__startTime = time.time()
            logger.debug("Starting at %s", time.strftime("%Y %m %d %H:%M:%S", time.localtime()))
        #
        self.__tiP = NeighborInteractionProvider(
            self.__cachePath, useCache=self.__useCache, cfgOb=self.__cfgOb, configName=self.__configName, numProc=self.__numProc, chunkSize=self.__chunkSize
        )

    def update(self, incremental=True):
        ok = self.__tiP.generate(distLimit=5.0, updateOnly=incremental, fmt="pickle")
        return ok

    def backup(self):
        ok = self.__tiP.backup(self.__cfgOb, self.__configName, remotePrefix=self.__stashRemotePrefix, useStash=True, useGit=True)
        return ok

    def restore(self, minCount=0):
        ok1 = self.__tiP.restore(self.__cfgOb, self.__configName, remotePrefix=self.__stashRemotePrefix, useStash=True, useGit=True)
        ok2 = self.__tiP.reload()
        ok3 = self.__tiP.testCache(minCount=minCount)
        return ok1 and ok2 and ok3

    def convert(self, fmt1="json", fmt2="pickle"):
        ok = self.__tiP.convert(fmt1=fmt1, fmt2=fmt2)
        return ok


if __name__ == "__main__":
    tiWf = NeighborInteractionWorkflow()
    tiWf.update(incremental=False)
