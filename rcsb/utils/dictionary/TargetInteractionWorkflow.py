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
        self.__numProc = kwargs.get("numProc", 10)
        self.__chunkSize = kwargs.get("chunkSize", 100)
        #
        self.__cfgOb = ConfigUtil(configPath=self.__configPath, defaultSectionName=self.__configName)
        logger.info("Configuration file path %s", self.__configPath)

    def update(self):
        tiP = TargetInteractionProvider(self.__cfgOb, self.__configName, self.__cachePath, useCache=False, numProc=self.__numProc, chunkSize=self.__chunkSize)
        ok1 = tiP.generate(distLimit=5.0, fmt="json", indent=3)
        ok2 = tiP.toStash()
        return ok1 & ok2

    def restore(self):
        tiP = TargetInteractionProvider(self.__cfgOb, self.__configName, self.__cachePath, useCache=False)
        ok1 = tiP.fromStash()
        ok2 = tiP.reload()
        ok3 = tiP.testCache(minCount=80)
        return ok1 and ok2 and ok3


if __name__ == "__main__":
    tiWf = TargetInteractionWorkflow()
    tiWf.update()
