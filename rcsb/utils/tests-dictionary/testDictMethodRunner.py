# File:    DictMethodRunnerTests.py
# Author:  J. Westbrook
# Date:    18-Aug-2018
# Version: 0.001
#
# Update:
#    12-Nov-2018 jdw add chemical component and bird chemical component tests
#     5-Jun-2019 jdw revise for new method runner api
#    16-Jul-2019 jdw remove schema processing.
#    16-Oct-2026 dwp add testSubjectOfInvestigationNoLigandOfInterest (31vo) covering
#                    _pdbx_entry_details.has_ligand_of_interest = "N"
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

from mmcif.api.DictMethodRunner import DictMethodRunner
from rcsb.utils.dictionary.DictionaryApiProviderWrapper import DictionaryApiProviderWrapper
from rcsb.utils.dictionary.DictMethodResourceProvider import DictMethodResourceProvider
from rcsb.utils.repository.RepositoryProvider import RepositoryProvider
from rcsb.utils.config.ConfigUtil import ConfigUtil
from rcsb.utils.io.MarshalUtil import MarshalUtil

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()
logger.setLevel(logging.INFO)

HERE = os.path.abspath(os.path.dirname(__file__))
TOPDIR = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))


class DictMethodRunnerTests(unittest.TestCase):
    def setUp(self):
        # self.__isMac = platform.system() == "Darwin"
        # self.__excludeTypeL = None if self.__isMac else ["optional"]
        self.__isMac = False
        self.__excludeTypeL = ["optional"]
        self.__export = True
        self.__numProc = 2
        mockTopPath = os.path.join(TOPDIR, "rcsb", "mock-data")
        self.__cachePath = os.path.join(HERE, "test-output", "CACHE")
        configPath = os.path.join(mockTopPath, "config", "dbload-setup-example.yml")
        configName = "site_info_configuration"
        self.__configName = configName
        self.__cfgOb = ConfigUtil(configPath=configPath, defaultSectionName=configName, mockTopPath=mockTopPath)
        self.__mU = MarshalUtil(workPath=self.__cachePath)

        #
        pdbIdList = [
            "1ah1",
            "1b5f",
            "1bmv",
            "1c58",
            "1dsr",
            "1dul",
            "1kqe",
            "1o3q",
            "1sfo",
            "2hw3",
            "2hyv",
            "2osl",
            "2voo",
            "2wmg",
            "3ad7",
            "3hya",
            "3iyd",
            "3mbg",
            "3rer",
            "3vd8",
            "3vfj",
            "3x11",
            "3ztj",
            "4e2o",
            "4en8",
            "4mey",
            "5eu8",
            "5kds",
            "5tm0",
            "5vh4",
            "5vp2",
            "6fsz",
            "6lu7",
            "6nn7",
            "6q20",
            "6rfk",
            "6rku",
            "6yrq",
        ]
        ccIdList = [
            "0EG",
            "0F7",
            "0G7",
            "0PP",
            "1G1",
            "2RT",
            "2XL",
            "2XN",
            "ATP",
            "BJA",
            "BM3",
            "CNC",
            "DAL",
            "DDZ",
            "DHA",
            "DSN",
            "GTP",
            "HKL",
            "JS4",
            "NAC",
            "NAG",
            "NND",
            "PTR",
            "SEP",
            "SMJ",
            "STL",
            "UNK",
            "UNX",
            "UVL",
        ]
        prdIdList = [
            "PRD_000009",
            "PRD_000010",
            "PRD_000060",
            "PRD_000154",
            "PRD_000198",
            "PRD_000220",
            "PRD_000877",
            "PRD_000882",
            "PRD_000979",
        ]
        self.__testCaseList = [
            {
                "contentType": "pdbx_core",
                "mergeContent": ["vrpt"],
                "fileLimit": 10,
                "idCodeList": pdbIdList,
            },
            {
                "contentType": "chem_comp",
                "mergeContent": None,
                "fileLimit": None,
                "idCodeList": ccIdList,
            },
            {
                "contentType": "bird",
                "mergeContent": None,
                "fileLimit": None,
                "idCodeList": prdIdList,
            },
        ]
        #
        self.__modulePathMap = self.__cfgOb.get("DICT_METHOD_HELPER_MODULE_PATH_MAP", sectionName=configName)
        #
        self.__startTime = time.time()
        logger.debug("Starting %s at %s", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()))

    def tearDown(self):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Maximum resident memory size %.4f %s", rusageMax / 1.0e6, unitS)
        endTime = time.time()
        logger.info("Completed %s at %s (%.4f seconds)", self.id(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - self.__startTime)

    def testMethodRunner(self):
        """Test method runner for multiple content types."""
        for tD in self.__testCaseList:
            self.__runContentType(tD["contentType"], tD["mergeContent"], tD["fileLimit"], tD["idCodeList"])

    def __runContentType(self, contentType, mergeContent, fileLimit, idCodeList):
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
                providerTypeExcludeL=self.__excludeTypeL,
            )
            dmh = DictMethodRunner(dictApi, modulePathMap=self.__modulePathMap, resourceProvider=rP)
            rpP = RepositoryProvider(cfgOb=self.__cfgOb, numProc=self.__numProc, fileLimit=fileLimit, cachePath=self.__cachePath)
            locatorObjList = rpP.getLocatorObjList(contentType=contentType, mergeContentTypes=mergeContent, inputIdCodeList=idCodeList)
            logger.info("Length of locator list (%d)", len(locatorObjList))
            if self.__isMac:
                containerList = rpP.getContainerList(locatorObjList)
            else:
                # strip down tests for Azure linux low memory -
                tObjL = []
                for locatorObj in locatorObjList[:10]:
                    if "locator" in locatorObj and "5vp2" in locatorObj["locator"].lower():
                        continue
                    tObjL.append(locatorObj)
                containerList = rpP.getContainerList(tObjL)
            #
            logger.info("Processing container length (%d)", len(containerList))

            for container in containerList:
                cName = container.getName()
                #
                logger.info("Processing container %s", cName)
                dmh.apply(container)
                if self.__export:
                    savePath = os.path.join(HERE, "test-output", "export-cif", cName + "-with-method.cif")
                    self.__mU.doExport(savePath, [container], fmt="mmcif")

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    @unittest.skip("Redundant test")
    def testResourceCache(self):
        rP = DictMethodResourceProvider(
            self.__cfgOb,
            configName=self.__configName,
            cachePath=self.__cachePath,
            restoreUseStash=False,
            restoreUseGit=True,
            providerTypeExcludeL=self.__excludeTypeL,
        )
        ok = rP.cacheResources(useCache=True, doRestore=True, useStash=False, useGit=True)
        self.assertTrue(ok)

    def testSubjectOfInvestigationNoLigandOfInterest(self):
        """Entry 31vo carries _pdbx_entry_details.has_ligand_of_interest = "N".

        The author has stated that no ligand is the subject of the current study, so no
        instance in rcsb_nonpolymer_instance_validation_score may be flagged.  Without the
        flag, NAG (FW 221.21, absent from the ligand exclude list) would be designated
        "Y"/"RCSB" by the molecular-weight heuristic, so this entry exercises that branch.
        """
        try:
            cName, cObj = self.__applyMethodsForId("31vo")
            self.assertIsNotNone(cObj, "%s is missing rcsb_nonpolymer_instance_validation_score" % cName)
            #
            compIdL = [cObj.getValue("comp_id", ii) for ii in range(cObj.getRowCount())]
            logger.info("%s nonpolymer validation score comp_ids %r", cName, compIdL)
            # Guard against a vacuous pass if the entry yielded no scored nonpolymer instances
            self.assertIn("NAG", compIdL, "%s has no NAG row - was the validation report merged?" % cName)
            #
            for ii in range(cObj.getRowCount()):
                compId = cObj.getValue("comp_id", ii)
                isTarget = cObj.getValueOrDefault("is_subject_of_investigation", ii, defaultValue=None)
                isTargetProv = cObj.getValueOrDefault("is_subject_of_investigation_provenance", ii, defaultValue=None)
                self.assertEqual(isTarget, "N", "%s %s is_subject_of_investigation %r (expected 'N')" % (cName, compId, isTarget))
                self.assertIsNone(isTargetProv, "%s %s unexpected provenance %r" % (cName, compId, isTargetProv))

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()

    def __applyMethodsForId(self, idCode, catName="rcsb_nonpolymer_instance_validation_score"):
        """Apply the dictionary methods to a single pdbx_core entry and return (name, category object)."""
        dP = DictionaryApiProviderWrapper(self.__cachePath, useCache=True, cfgOb=self.__cfgOb, configName=self.__configName)
        dictApi = dP.getApiByName("pdbx_core")
        rP = DictMethodResourceProvider(
            self.__cfgOb,
            configName=self.__configName,
            cachePath=self.__cachePath,
            restoreUseStash=False,
            restoreUseGit=True,
            providerTypeExcludeL=self.__excludeTypeL,
        )
        dmh = DictMethodRunner(dictApi, modulePathMap=self.__modulePathMap, resourceProvider=rP)
        rpP = RepositoryProvider(cfgOb=self.__cfgOb, numProc=self.__numProc, fileLimit=None, cachePath=self.__cachePath)
        locatorObjList = rpP.getLocatorObjList(contentType="pdbx_core", mergeContentTypes=["vrpt"], inputIdCodeList=[idCode])
        self.assertEqual(len(locatorObjList), 1, "No locator for %s - is it in the current holdings list?" % idCode)
        containerList = rpP.getContainerList(locatorObjList)
        self.assertEqual(len(containerList), 1, "No container for %s" % idCode)
        #
        container = containerList[0]
        dmh.apply(container)
        return container.getName(), container.getObj(catName) if container.exists(catName) else None

    def testMethodRunnerSetup(self):
        """Test the setup methods for method runner class"""
        try:
            dP = DictionaryApiProviderWrapper(self.__cachePath, useCache=True, cfgOb=self.__cfgOb, configName=self.__configName)
            dictApi = dP.getApiByName("pdbx")
            rP = DictMethodResourceProvider(
                self.__cfgOb,
                configName=self.__configName,
                cachePath=self.__cachePath,
                restoreUseStash=False,
                restoreUseGit=True,
                providerTypeExcludeL=self.__excludeTypeL,
            )
            dmh = DictMethodRunner(dictApi, modulePathMap=self.__modulePathMap, resourceProvider=rP)
            ok = dmh is not None
            self.assertTrue(ok)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
            self.fail()


def dictMethodRunnerSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DictMethodRunnerTests("testMethodRunner"))
    suiteSelect.addTest(DictMethodRunnerTests("testSubjectOfInvestigationNoLigandOfInterest"))
    return suiteSelect


def dictMethodRunnerSetupSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DictMethodRunnerTests("testMethodRunnerSetup"))
    return suiteSelect


if __name__ == "__main__":

    mySuite = dictMethodRunnerSetupSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)

    mySuite = dictMethodRunnerSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
