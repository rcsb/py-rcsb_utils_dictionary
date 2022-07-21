##
# File:    DictMethodResourceProvider.py
# Author:  J. Westbrook
# Date:    25-Jul-2021
# Version: 0.002 Refactored version
#
#
# Updates:
#  17-Jul-2019 jdw add resource for common utilities and dictionary api
#   7-Aug-2019 jdw use dictionary locator map
#  13-Aug-2019 jdw return class instances in all cases. Add cache management support.
#   9-Sep-2019 jdw add AtcProvider() and SiftsSummaryProvider()
#  25-Nov-2019 jdw add CitationReferenceProvider(), ChemCompProvider() and  JournalTitleAbbreviationProvider()'s
#  16-Feb-2020 jdw add support for configuration of development resources
#  19-Mar-2020 jdw add ResidProvider() and send cachePath directly to all modules in rcsb.utils.chemref.
#  29-Jul-2020 jdw add PubChemProvider() from  rcsb.utils.chemref.
#  30-Jul-2020 jdw add PharosProvider() from  rcsb.utils.chemref.
#  29-Oct-2020 jdw add method getReferenceSequenceAlignmentOpt()
#  18-Feb-2021 jdw add TargerInteractionProvider()
#  25-Jul-2021 jdw refactored with common provider api
#  27-Jul-2021 jdw add exclusion filter for testing
#   3-Aug-2021 jdw add backup options for nonbuildable providers
#  29-Apr-2022 dwp add ModelCacheProvider()
#  21-Jul-2022  bv Update ModelCacheProvider to make providerType "core" and not stashable or buildable
##
##
"""
Resource provider for dictionary method runner and DictMethodHelper tools.

"""
__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import platform
import resource
import time

from rcsb.utils.chemref.AtcProvider import AtcProvider
from rcsb.utils.chemref.BirdProvider import BirdProvider
from rcsb.utils.chemref.ChemCompModelProvider import ChemCompModelProvider
from rcsb.utils.chemref.ChemCompProvider import ChemCompProvider
from rcsb.utils.chemref.DrugBankProvider import DrugBankProvider
from rcsb.utils.chemref.PharosProvider import PharosProvider
from rcsb.utils.chemref.PsiModProvider import PsiModProvider
from rcsb.utils.chemref.PubChemProvider import PubChemProvider
from rcsb.utils.chemref.RcsbLigandScoreProvider import RcsbLigandScoreProvider
from rcsb.utils.chemref.ResidProvider import ResidProvider

# ---
from rcsb.utils.citation.CitationReferenceProvider import CitationReferenceProvider
from rcsb.utils.citation.JournalTitleAbbreviationProvider import JournalTitleAbbreviationProvider

# ---
from rcsb.utils.dictionary.DictionaryApiProviderWrapper import DictionaryApiProviderWrapper
from rcsb.utils.dictionary.DictMethodCommonUtils import DictMethodCommonUtils
from rcsb.utils.dictionary.NeighborInteractionProvider import NeighborInteractionProvider
from rcsb.utils.ec.EnzymeDatabaseProvider import EnzymeDatabaseProvider
from rcsb.utils.io.FileUtil import FileUtil
from rcsb.utils.io.SingletonClass import SingletonClass

# ---
from rcsb.utils.seq.GlycanProvider import GlycanProvider
from rcsb.utils.seq.GlyGenProvider import GlyGenProvider
from rcsb.utils.seq.PfamProvider import PfamProvider
from rcsb.utils.seq.SiftsSummaryProvider import SiftsSummaryProvider

# ---
from rcsb.utils.struct.CathClassificationProvider import CathClassificationProvider
from rcsb.utils.struct.EcodClassificationProvider import EcodClassificationProvider
from rcsb.utils.struct.EntryInfoProvider import EntryInfoProvider
from rcsb.utils.struct.Scop2ClassificationProvider import Scop2ClassificationProvider
from rcsb.utils.struct.ScopClassificationProvider import ScopClassificationProvider

# ---
from rcsb.utils.targets.CARDTargetFeatureProvider import CARDTargetFeatureProvider
from rcsb.utils.targets.ChEMBLTargetCofactorProvider import ChEMBLTargetCofactorProvider
from rcsb.utils.targets.DrugBankTargetCofactorProvider import DrugBankTargetCofactorProvider
from rcsb.utils.targets.IMGTTargetFeatureProvider import IMGTTargetFeatureProvider
from rcsb.utils.targets.PharosTargetCofactorProvider import PharosTargetCofactorProvider
from rcsb.utils.targets.SAbDabTargetFeatureProvider import SAbDabTargetFeatureProvider

# --
from rcsb.utils.taxonomy.TaxonomyProvider import TaxonomyProvider

# --
from rcsb.utils.insilico3d.ModelCacheProvider import ModelCacheProvider

logger = logging.getLogger(__name__)


class DictMethodResourceProvider(SingletonClass):
    """Resource provider for dictionary method runner and DictMethodHelper tools."""

    def __init__(self, cfgOb, **kwargs):
        """Resource provider for dictionary method runner and DictMethodHelper tools.

        Arguments:
            cfgOb (object): instance ConfigUtils class
            configName (str, optional): configuration section name (default: default section name)
            cachePath (str, optional): path used for temporary file management (default: '.')
            restoreUseStash (bool, optional): use remote stash storage for restore operations
            restoreUseGit (bool, optional): use remote storage for restore operations
            providerTypeExclude (str, optional): exclude providers containing this name string. Defaults to None.
        """
        self.__cfgOb = cfgOb
        self.__configName = kwargs.get("configName", self.__cfgOb.getDefaultSectionName())
        self.__cachePath = kwargs.get("cachePath", ".")
        self.__restoreUseStash = kwargs.get("restoreUseStash", True)
        self.__restoreUseGit = kwargs.get("restoreUseGit", True)
        self.__providerTypeExclude = kwargs.get("providerTypeExclude", None)
        # --
        self.__providerInstanceD = {}
        self.__providerD = {
            "DictionaryAPIProviderWrapper instance": {
                "class": DictionaryApiProviderWrapper,
                "configArgMap": {
                    "cfgOb": (self.__cfgOb, "value"),
                    "configName": (self.__configName, "value"),
                },
                "stashable": False,
                "buildable": False,
                "providerType": "core",
            },
            "DictMethodCommonUtils instance": {
                "class": DictMethodCommonUtils,
                "configArgMap": {},
                "stashable": False,
                "buildable": False,
                "providerType": "core",
            },
            "NeighborInteractionProvider instance": {
                "class": NeighborInteractionProvider,
                "configArgMap": {
                    "cfgOb": (self.__cfgOb, "value"),
                    "configName": (self.__configName, "value"),
                },
                "stashable": True,
                "buildable": False,
                "providerType": "core",
            },
            "Scop2Provider instance": {
                "class": Scop2ClassificationProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "optional_1",
            },
            "EcodProvider instance": {
                "class": EcodClassificationProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "optional_1",
            },
            "ScopProvider instance": {
                "class": ScopClassificationProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "optional_1",
            },
            "CathProvider instance": {
                "class": CathClassificationProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "optional_1",
            },
            "DrugBankProvider instance": {
                "class": DrugBankProvider,
                "configArgMap": {
                    "username": ("_DRUGBANK_AUTH_USERNAME", "configItem"),
                    "password": ("_DRUGBANK_AUTH_PASSWORD", "configItem"),
                    # "urlTarget": ("DRUGBANK_MOCK_URL_TARGET", "configPath"),
                },
                "stashable": True,
                "buildable": True,
                "providerType": "core",
            },
            "AtcProvider instance": {
                "class": AtcProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "core",
            },
            "BirdProvider instance": {
                "class": BirdProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "core",
            },
            "ChemCompModelProvider instance": {
                "class": ChemCompModelProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "core",
            },
            "ChemCompProvider instance": {
                "class": ChemCompProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "core",
            },
            #
            "PharosProvider instance": {
                "class": PharosProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": False,
                "providerType": "optional_1",
            },
            "PsiModProvider instance": {
                "class": PsiModProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "optional_1",
            },
            "PubChemProvider instance": {
                "class": PubChemProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": False,
                "providerType": "core",
            },
            "RcsbLigandScoreProvider instance": {
                "class": RcsbLigandScoreProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "core",
            },
            "ResidProvider instance": {
                "class": ResidProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "core",
            },
            #
            "CitationReferenceProvider instance": {
                "class": CitationReferenceProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "core",
            },
            "JournalTitleAbbreviationProvider instance": {
                "class": JournalTitleAbbreviationProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "core",
            },
            "TaxonomyProvider instance": {
                "class": TaxonomyProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "core",
            },
            "EnzymeDatabaseProvider instance": {
                "class": EnzymeDatabaseProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "core",
            },
            #
            "GlyGenProvider instance": {
                "class": GlyGenProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "optional_1",
            },
            "GlycanProvider instance": {
                "class": GlycanProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": False,
                "providerType": "optional_1",
            },
            "SiftsSummaryProvider instance": {
                "class": SiftsSummaryProvider,
                "configArgMap": {
                    "abbreviated": ("TEST", "value"),
                    "srcDirPath": ("SIFTS_SUMMARY_DATA_PATH", "configPath"),
                },
                "stashable": True,
                "buildable": True,
                "providerType": "core",
            },
            "PfamProvider instance": {
                "class": PfamProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": True,
                "providerType": "optional_1",
            },
            "DrugBankTargetCofactorProvider instance": {
                "class": DrugBankTargetCofactorProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": False,
                "providerType": "optional_1",
            },
            "ChEMBLTargetCofactorProvider instance": {
                "class": ChEMBLTargetCofactorProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": False,
                "providerType": "optional_1",
            },
            "PharosTargetCofactorProvider instance": {
                "class": PharosTargetCofactorProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": False,
                "providerType": "optional_1",
            },
            "CARDTargetFeatureProvider instance": {
                "class": CARDTargetFeatureProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": False,
                "providerType": "optional_1",
            },
            "IMGTTargetFeatureProvider instance": {
                "class": IMGTTargetFeatureProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": False,
                "providerType": "optional_1",
            },
            "SAbDabTargetFeatureProvider instance": {
                "class": SAbDabTargetFeatureProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": False,
                "providerType": "optional_1",
            },
            "EntryInfoProvider instance": {
                "class": EntryInfoProvider,
                "configArgMap": {},
                "stashable": True,
                "buildable": False,
                "providerType": "optional_1",
            },
            "ModelCacheProvider instance": {
                "class": ModelCacheProvider,
                "configArgMap": {
                    "holdingsRemotePath": ("PDBX_COMP_MODEL_CACHE_LIST_PATH", "configPath"),
                },
                "stashable": False,
                "buildable": False,
                "providerType": "core",
            },
            # --
        }
        logger.info("Dictionary resource provider restoreUseGit %r restoreUseStash %r providerTypeExclude %r", self.__restoreUseGit, self.__restoreUseStash, self.__providerTypeExclude)
        self.__filterProviderWarnD = {}
        #

    def echo(self, msg):
        logger.info(msg)

    def getReferenceSequenceAlignmentOpt(self):
        return self.__cfgOb.get("REFERENCE_SEQUENCE_ALIGNMENTS", sectionName=self.__configName, default="SIFTS")

    def getResource(self, providerName, default=None, useCache=True, **kwargs):
        """Return the named input cached resource or the default value.

        Arguments:
            providerName (str): resource provider name
            useCache (bool, optional): use an existing provider instance in the current cache. Defaults to True.
            default (obj, optional): the default return value if the requested object cannot be returned. Defaults to None.
            doRestore (bool, optional): when useCache=True restore the cache from a prior backup. Defaults to False.
            doBackup (bool, optional): when building cache (useCache=False) backup the cached directory. Defaults to False.
            useGit (bool, optional): use git repository storage for backup operations. Defaults to False.
            useStash (bool, optional): use stash storage and backup operations. Defaults to True.
            remotePrefix (str, optional): remote prefix for a multi-channel stash rotation. Defaults to None.
            cacheInstance (bool, optional): hold a reference to the data for the cached provided. Defaults to True.

        Returns:
            (obj): instance of the resource object or default value
        """
        logger.debug("Requesting cached provider resource %r (useCache %r)", providerName, useCache)
        #
        # Return a cached instance
        if useCache and providerName in self.__providerInstanceD and self.__providerInstanceD[providerName]:
            return self.__providerInstanceD[providerName]
        #
        if providerName not in self.__providerD:
            logger.error("Request for unsupported provider resource %r returning %r", providerName, default)
            return default

        # Apply exclusions
        if self.__providerTypeExclude and self.__providerTypeExclude in self.__providerD[providerName]["providerType"]:
            if providerName not in self.__filterProviderWarnD:
                logger.info("Provider %r excluded by filter %r", providerName, self.__providerTypeExclude)
            self.__filterProviderWarnD[providerName] = True
            return default
        #
        # Reload the instance into the cache and optionally store the instance
        kwargs["cacheInstance"] = kwargs["cacheInstance"] if "cacheInstance" in kwargs else True
        ok = self.__cacheProvider(providerName, self.__cfgOb, self.__configName, self.__cachePath, useCache=useCache, **kwargs)
        if ok and kwargs["cacheInstance"]:
            return self.__providerInstanceD[providerName]
        #
        return default

    def cacheResources(self, useCache=False, **kwargs):
        """Update and optionally clear all resource caches.

        Args:
            useCache (bool, optional): use current cace. Defaults to False.
            doRestore (bool, optional): when useCache=True restore the cache from a prior backup. Defaults to False.
            doBackup (bool, optional): when building cache (useCache=False) backup the cached directory. Defaults to False.
            useGit (bool, optional): use git repository storage for  backup operations. Defaults to False.
            useStash (bool, optional): use stash storage for backup operations. Defaults to True.
            remotePrefix (str, optional): remote prefix for a multi-channel stash rotation. Defaults to None.
            providerSelect (str, optional): select buildable|nonbuildable|None. Defaults to None
            cacheInstance (bool, optional): hold a reference to the data for the cached provided. Defaults to False.
            clearCache (bool, optional): clear the cache directory before rebuilding.
        Returns:
            bool: True for success or False otherwise
        """
        ret = True
        tName = "CHECKING" if useCache else "REBUILDING/RELOADING"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("Begin %s cache for %d candidate resources", tName, len(self.__providerD))
        #
        kwargs["cacheInstance"] = kwargs["cacheInstance"] if "cacheInstance" in kwargs else True
        providerSelect = kwargs.get("providerSelect", None)
        clearCache = kwargs.get("clearCache", None)
        failList = []
        #
        if not useCache and clearCache:
            fU = FileUtil()
            fU.remove(self.__cachePath)
        #
        for providerName in sorted(self.__providerD):
            if self.__providerTypeExclude and self.__providerTypeExclude in self.__providerD[providerName]["providerType"]:
                logger.info("Provider excluded by filter %r", providerName)
                continue
            #
            rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            bFlag = self.__providerD[providerName]["buildable"]
            if providerSelect and (bFlag == (providerSelect != "buildable")):
                continue
            startTime = time.time()
            rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            #
            # Check the current instance cache --
            if useCache and providerName in self.__providerInstanceD:
                ok = self.__providerInstanceD[providerName].testCache()
                if ok:
                    continue
            #
            # Update the cache if necessary
            logger.debug("Updating cache resources for %r", providerName)
            ok = self.__cacheProvider(providerName, self.__cfgOb, self.__configName, self.__cachePath, useCache=useCache, **kwargs)
            if not ok:
                logger.error("%s %s fails", tName, providerName)
                failList.append(providerName)
            #
            ret = ret and ok
            if not ret:
                logger.info("%s resource %r step status %r cumulative status %r", tName, providerName, ok, ret)
            self.__resourceUsageReport(providerName, startTime, rusageMax)
        #
        logger.info("Completed %s %d resource instances status (%r) failures %r", tName, len(self.__providerInstanceD), ret, failList)
        return ret

    def __getClassArgs(self, providerName, cfgOb, configName):
        classArgs = {}
        for argName, configTup in self.__providerD[providerName]["configArgMap"].items():
            if configTup[1] == "configItem":
                classArgs[argName] = cfgOb.get(configTup[0], sectionName=configName)
            elif configTup[1] == "configPath":
                classArgs[argName] = cfgOb.getPath(configTup[0], sectionName=configName)
            elif configTup[1] == "value":
                classArgs[argName] = configTup[0]
        return classArgs

    def __cacheProvider(self, providerName, cfgOb, configName, cachePath, useCache=True, **kwargs):
        if not self.__providerD[providerName]["stashable"]:
            return self.__cacheProviderNonStashable(providerName, cfgOb, configName, cachePath, useCache=useCache, **kwargs)
        elif self.__providerD[providerName]["buildable"]:
            return self.__cacheProviderBuildable(providerName, cfgOb, configName, cachePath, useCache=useCache, **kwargs)
        elif not self.__providerD[providerName]["buildable"]:
            return self.__cacheProviderNonBuildable(providerName, cfgOb, configName, cachePath, useCache=useCache, **kwargs)
        else:
            return False

    def __cacheProviderNonStashable(self, providerName, cfgOb, configName, cachePath, useCache=True, **kwargs):
        """Instantiate a non-stashable resource provider.

        Args:
            providerName (str): provider name
            cfgOb (obj): instance of the configuration object ConfigUtil()
            configName (str): configuration section name
            cachePath (str): path to the directory containing the cached data
            useCache (bool, optional): use existing cache (with optional restore) otherwise rebuild the cache from scratch. Defaults to True.

        Returns:
            bool: True for success or False otherwise

         ---
        Provider types -
            stashable == False ->  class objects with no data payload

            stashable == True  ->  class with data payload

                buildable == True  both construct (useCache=False) and deliver (useCache=True) data payloads


                buildable == False  only deliver stashed payloads from precomputed workflows
                                  (useCache = True)   deliver locally cached data payload (with optional restore)
                                  (useCache = False)  restore stashed payload to local cache

        """
        logger.debug("providerName %r configName %s cachePath %s kwargs %r", providerName, configName, cachePath, kwargs)
        #
        cacheInstance = kwargs.get("cacheInstance", True)
        classArgs = self.__getClassArgs(providerName, cfgOb, configName)
        logger.debug("%r classArgs %r", providerName, classArgs)
        #
        try:
            prI = self.__providerD[providerName]["class"](cachePath=cachePath, useCache=useCache, **classArgs)
            ok = prI.testCache()
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        if ok and cacheInstance:
            self.__providerInstanceD[providerName] = prI
        else:
            del prI
        #
        return ok

    def __cacheProviderNonBuildable(self, providerName, cfgOb, configName, cachePath, useCache=True, **kwargs):
        """Load or restore/load the cached resources for a non-buildable resource provider.

        Args:
            providerName (str): provider name
            cfgOb (obj): instance of the configuration object ConfigUtil()
            configName (str): configuration section name
            cachePath (str): path to the directory containing the cached data
            useCache (bool, optional): use existing cache (with optional restore) otherwise rebuild the cache from scratch. Defaults to True.
            doRestore (bool, optional): when useCache=True restore the cache from a prior backup. Defaults to False.
                                        when useCache=False restore is done by default.
            remotePrefix (str, optional): remote prefix for a multi-channel stash rotation. Defaults to None.

        Returns:
            bool: True for success or False otherwise

         ---
        Provider types -
            stashable == False ->  class objects with no data payload

            stashable == True  ->  class with data payload

                buildable == True  both construct (useCache=False) and deliver (useCache=True) data payloads


                buildable == False  only deliver stashed payloads from precomputed workflows
                                  (useCache = True)   deliver locally cached data payload (with optional restore)
                                  (useCache = False)  restore stashed payload to local cache

        """
        logger.debug("providerName %r configName %s cachePath %s kwargs %r", providerName, configName, cachePath, kwargs)
        #
        classArgs = self.__getClassArgs(providerName, cfgOb, configName)

        isStashable = self.__providerD[providerName]["stashable"]
        logger.debug("%r classArgs %r", providerName, classArgs)
        #
        cacheInstance = kwargs.get("cacheInstance", True)
        remotePrefix = kwargs.get("remotePrefix", None)
        minCount = kwargs.get("remotePrefix", 5)
        if useCache:
            doRestore = kwargs.get("doRestore", True)
            try:
                prI = self.__providerD[providerName]["class"](cachePath=cachePath, useCache=useCache, **classArgs)
                ok = prI.testCache()
                if not ok and doRestore and isStashable:
                    prI.restore(cfgOb, configName, remotePrefix=remotePrefix, useStash=self.__restoreUseStash, useGit=self.__restoreUseGit)
                    prI = self.__providerD[providerName]["class"](cachePath=cachePath, useCache=True, **classArgs)
                    ok = prI.testCache(minCount=minCount)
            except Exception as e:
                logger.exception("Failing with %s", str(e))
        else:
            doBackup = kwargs.get("doBackup", False)
            try:
                prI = self.__providerD[providerName]["class"](cachePath=cachePath, useCache=useCache, **classArgs)
                prI.restore(cfgOb, configName, remotePrefix=remotePrefix, useStash=self.__restoreUseStash, useGit=self.__restoreUseGit)
                prI = self.__providerD[providerName]["class"](cachePath=cachePath, useCache=True, **classArgs)
                ok = prI.testCache(minCount=minCount)
                if ok and doBackup and isStashable:
                    useGit = kwargs.get("useGit", False)
                    useStash = kwargs.get("useStash", True)
                    okB = prI.backup(cfgOb, configName, remotePrefix=remotePrefix, useStash=useStash, useGit=useGit)
                    ok = ok and okB
            except Exception as e:
                logger.exception("Failing with %s", str(e))
        #
        if ok and cacheInstance:
            self.__providerInstanceD[providerName] = prI
        else:
            del prI
        #
        return ok

    def __cacheProviderBuildable(self, providerName, cfgOb, configName, cachePath, useCache=True, **kwargs):
        """Load or build the cached data for the input buildable resource provider.

        Args:
            providerName (str): provider name
            cfgOb (obj): instance of the configuration object ConfigUtil()
            configName (str): configuration section name
            cachePath (str): path to the directory containing the cached data
            useCache (bool, optional): use existing cache (with optional restore) otherwise rebuild the cache from scratch. Defaults to True.
            doRestore (bool, optional): when useCache=True restore the cache from a prior backup. Defaults to False.
            doBackup (bool, optional): when building cache (useCache=False) backup the cached directory. Defaults to False.
            useGit (bool, optional): use git repository storage for backup operations. Defaults to False.
            useStash (bool, optional): use stash storage for backup operations. Defaults to True.
            remotePrefix (str, optional): remote prefix for a multi-channel stash rotation. Defaults to None.

        Returns:
            bool: True for success or False otherwise

         ---
        Provider types -
            stashable == False ->  class objects with no data payload

            stashable == True  ->  class with data payload

                buildable == True  both construct (useCache=False) and deliver (useCache=True) data payloads


                buildable == False  only deliver stashed payloads from precomputed workflows
                                  (useCache = True)   deliver locally cached data payload
                                  (useCache = False)  restore stashed payload to local cache

        """
        logger.debug("providerName %r configName %s useCache %r cachePath %s kwargs %r", providerName, configName, useCache, cachePath, kwargs)
        #
        classArgs = self.__getClassArgs(providerName, cfgOb, configName)

        isStashable = self.__providerD[providerName]["stashable"]
        useGit = kwargs.get("useGit", False)
        useStash = kwargs.get("useStash", True)
        logger.debug("%r classArgs %r", providerName, classArgs)
        #
        ok = False
        cacheInstance = kwargs.get("cacheInstance", True)
        remotePrefix = kwargs.get("remotePrefix", None)
        if useCache:
            doRestore = kwargs.get("doRestore", True)
            try:
                prI = self.__providerD[providerName]["class"](cachePath=cachePath, useCache=useCache, **classArgs)
                ok = prI.testCache()
                if not ok and doRestore and isStashable:
                    prI.restore(cfgOb, configName, remotePrefix=remotePrefix, useStash=self.__restoreUseStash, useGit=self.__restoreUseGit)
                    prI = self.__providerD[providerName]["class"](cachePath=cachePath, useCache=useCache, **classArgs)
                    ok = prI.testCache()
            except Exception as e:
                logger.exception("Failing with %s", str(e))
        else:
            doBackup = kwargs.get("doBackup", False)
            try:
                prI = self.__providerD[providerName]["class"](cachePath=cachePath, useCache=useCache, **classArgs)
                ok = prI.testCache()
                if ok and doBackup and isStashable:
                    okB = prI.backup(cfgOb, configName, remotePrefix=remotePrefix, useStash=useStash, useGit=useGit)
                    ok = ok and okB
            except Exception as e:
                logger.exception("Failing with %s", str(e))
        #
        if ok and cacheInstance:
            self.__providerInstanceD[providerName] = prI
        else:
            del prI
        #
        return ok

    def syncCache(self, providerName, cfgOb, configName, cachePath, remotePrefix=None, sourceCache="stash"):
        """Synchronize cache data for the input provider from the input source cache to git stash storage

        Args:
            providerName (str): provider name
            cfgOb (obj): instance of the configuration object ConfigUtil()
            configName (str): configuration section name
            cachePath (str): path to the directory containing the cached data
            useCache (bool, optional): use existing cache (with optional restore) otherwise rebuild the cache from scratch. Defaults to True.
            remotePrefix (str, optional): remote prefix for a multi-channel stash rotation. Defaults to None.
            sourceCache (str, optional): source cache for the sync operation. Defaults to "stash".

        Returns:
            bool: True for success or False otherwise
        """
        logger.debug("providerName %r configName %s cachePath %s sourceCache %r", providerName, configName, cachePath, sourceCache)
        #
        classArgs = self.__getClassArgs(providerName, cfgOb, configName)
        logger.debug("%r classArgs %r", providerName, classArgs)
        #
        useCache = True
        try:
            if sourceCache == "stash":
                prI = self.__providerD[providerName]["class"](cachePath=cachePath, useCache=useCache, **classArgs)
                ok = prI.testCache()
                prI.restore(cfgOb, configName, remotePrefix=remotePrefix, useStash=True, useGit=False)
                prI = self.__providerD[providerName]["class"](cachePath=cachePath, useCache=useCache, **classArgs)
                ok = prI.testCache()
                okB = prI.backup(cfgOb, configName, remotePrefix=remotePrefix, useStash=False, useGit=True)
            elif sourceCache == "git":
                prI = self.__providerD[providerName]["class"](cachePath=cachePath, useCache=useCache, **classArgs)
                ok = prI.testCache()
                prI.restore(cfgOb, configName, remotePrefix=remotePrefix, useStash=False, useGit=True)
                prI = self.__providerD[providerName]["class"](cachePath=cachePath, useCache=useCache, **classArgs)
                ok = prI.testCache()
                okB = prI.backup(cfgOb, configName, remotePrefix=remotePrefix, useStash=True, useGit=False)
            else:
                logger.error("Unsupported source cache %r", sourceCache)

        except Exception as e:
            logger.exception("Failing with %s", str(e))
        #
        return ok and okB

    def __resourceUsageReport(self, providerName, startTime, startRusageMax):
        unitS = "MB" if platform.system() == "Darwin" else "GB"
        rusageMax = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        endTime = time.time()
        logger.info(
            "Step %s completed at %s (%.4f secs/Max %.3f %s/Delta %.3f %s)",
            providerName,
            time.strftime("%Y %m %d %H:%M:%S", time.localtime()),
            endTime - startTime,
            rusageMax / 1.0e6,
            unitS,
            (rusageMax - startRusageMax) / 1.0e6,
            unitS,
        )
