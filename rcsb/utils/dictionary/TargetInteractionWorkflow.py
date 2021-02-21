##
# File:    TargetInteractionWorkflow.py
# Author:  J. Westbrook
# Date:    18-Feb-2021
#
# Update:
#
#
##
"""
Worflow for generating and stashing non-polymer instance target interactions

"""

__docformat__ = "restructuredtext en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os

from rcsb.utils.dictionary.TargetInteractionProvider import TargetInteractionProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class TargetInteractionWorkflow(object):
    def __init__(self, **kwargs):
        # - edit as needed -
        self.__configName = kwargs.get("configName", "site_info_remote_configuration")
        self.__configPath = kwargs.get("configPath", os.path.join(HERE, "exdb-config-example.yml"))
        self.__cachePath = kwargs.get("cachePath", os.path.join(HERE, "CACHE"))
        self.__mockTopPath = kwargs.get("mockTopPath", None)
        self.__numProc = kwargs.get("numProc", 10)
        self.__chunkSize = kwargs.get("chunkSize", 10)
        self.__useCache = kwargs.get("useCache", False)
        #
        self.__cfgOb = ConfigUtil(configPath=self.__configPath, defaultSectionName=self.__configName, mockTopPath=self.__mockTopPath)
        logger.info("Configuration file path %s", self.__configPath)
        self.__tiP = TargetInteractionProvider(self.__cfgOb, self.__configName, self.__cachePath, useCache=self.__useCache, numProc=self.__numProc, chunkSize=self.__chunkSize)

    def update(self, incremental=True):
        ok = self.__tiP.generate(distLimit=5.0, updateOnly=incremental, fmt="json", indent=0)
        return ok

    def backup(self):
        ok = self.__tiP.toStash()
        return ok

    def restore(self):
        ok1 = self.__tiP.fromStash()
        ok2 = self.__tiP.reload()
        ok3 = self.__tiP.testCache(minCount=80)
        return ok1 and ok2 and ok3


if __name__ == "__main__":
    tiWf = TargetInteractionWorkflow()
    tiWf.update(incremental=False)
