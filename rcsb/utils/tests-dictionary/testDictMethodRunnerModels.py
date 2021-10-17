# File:    DictMethodRunnerTests.py
# Author:  J. Westbrook
# Date:    18-Aug-2018
# Version: 0.001
#
# Update:
#    12-Nov-2018 jdw add chemical component and bird chemical component tests
#     5-Jun-2019 jdw revise for new method runner api
#    16-Jul-2019 jdw remove schema processing.
##
"""
Tests for applying dictionary methods defined as references to helper plugin methods .

"""

__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import glob
import logging
import os
import platform
import resource
import time
import unittest

from mmcif.api.DictMethodRunner import DictMethodRunner
from rcsb.utils.dictionary.DictionaryApiProviderWrapper import DictionaryApiProviderWrapper
from rcsb.utils.dictionary.DictMethodResourceProvider import DictMethodResourceProvider
from rcsb.utils.repository.RepositoryProvider import RepositoryProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil

# from rcsb.utils.insilico3d.AlphaFoldModelProvider import AlphaFoldModelProvider
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


class DictMethodRunnerModelsTests(unittest.TestCase):
    #
    def setUp(self):

        self.__isMac = platform.system() == "Darwin"
        self.__excludeType = None if self.__isMac else "optional"
        self.__export = True
        self.__numProc = 2
        self.__fileLimit = 5
        self.__mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        self.__dataPath = os.path.join(HERE, "test-data")
        configPath = os.path.join(self.__mockTopPath, "config", "dbload-setup-example.yml")
        configName = "site_info_configuration"
        self.__configName = configName
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName, mockTopPath=self.__mockTopPath)
        self.__mU = MarshalUtil(workPath=self.__cachePath)
        self.__rpP = RepositoryProvider(cfgOb=self.__cfgOb, numProc=self.__numProc, fileLimit=self.__fileLimit, cachePath=self.__cachePath)
        #
        self.__testCaseList = [
            {"contentType": "pdbx_comp_model_core", "mockLength": self.__fileLimit, "mergeContent": None, "fileLimit": 5},
        ]
        #
        self.__modulePathMap = self.__cfgOb.get("DICT_METHOD_HELPER_MODULE_PATH_MAP", sectionName=configName)
        #
        self.__modelFixture()
        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 1.0e6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def __modelFixture(self):
        fU = FileUtil()
        modelSourcePath = os.path.join(self.__mockTopPath, "AF")
        for iPath in glob.iglob(os.path.join(modelSourcePath, "*.cif.gz")):
            fn = os.path.basename(iPath)
            uId = fn.split("-")[1]
            h3 = uId[-2:]
            h2 = uId[-4:-2]
            h1 = uId[-6:-4]
            oPath = os.path.join(self.__cachePath, "AlphaFold", h1, h2, h3, fn)
            fU.put(iPath, oPath)

    # def testFetchModels(self):
    #    aFMP = AlphaFoldModelProvider(cachePath=self.__cachePath, useCache=True, alphaFoldRequestedSpeciesList=["Staphylococcus aureus"])
    #    ok = aFMP.testCache()
    #    self.assertTrue(ok)

    def testMethodModelRunner(self):
        """Test method runner for multiple content types."""
        for tD in self.__testCaseList:
            self.__runContentType(tD["contentType"], tD["mergeContent"], tD["fileLimit"])

    def __runContentType(self, contentType, mergeContent, fileLimit):
        """Read and process test fixture data files from the input content type."""
        try:
            dP = DictionaryApiProviderWrapper(self.__cachePath, useCache=True, cfgOb=self.__cfgOb, configName=self.__configName)
            dictApi = dP.getApiByName(contentType)
            rP = DictMethodResourceProvider(
                self.__cfgOb,
                configName=self.__configName,
                cachePath=self.__cachePath,
                restoreUseStash=False,
                restoreUseGit=True,
                providerTypeExclude=self.__excludeType,
            )
            dmh = DictMethodRunner(dictApi, modulePathMap=self.__modulePathMap, resourceProvider=rP)
            locatorObjList = self.__rpP.getLocatorObjList(contentType=contentType, mergeContentTypes=mergeContent)
            logger.info("Length of locator list (%d)", len(locatorObjList))
            containerList = self.__rpP.getContainerList(locatorObjList[:fileLimit])

            logger.info("Processing container length (%d)", len(containerList))

            for container in containerList:
                cName = container.getName()
                logger.info("Processing container %s", cName)
                dmh.apply(container)
                if self.__export:
                    savePath = os.path.join(HERE, "test-output", "export-cif", cName + "-with-method.cif")
                    self.__mU.doExport(savePath, [container], fmt="mmcif")

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def dictMethodRunnerSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DictMethodRunnerModelsTests("testMethodModelRunner"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = dictMethodRunnerSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
