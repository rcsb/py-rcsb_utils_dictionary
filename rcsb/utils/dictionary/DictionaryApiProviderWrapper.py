##
# File:    DictionaryApiProviderWrapper.py
# Author:  J. Westbrook
# Date:   18-Aug-2019
# Version: 0.001 Initial version
#
# Updates:
#
##
"""
Wrapper for dictionary API provider.

"""
__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

import logging
import os.path

from rcsb.utils.dictionary.DictionaryApiProvider import DictionaryApiProvider
from rcsb.utils.io.SingletonClass import SingletonClass

logger = logging.getLogger(__name__)


class DictionaryApiProviderWrapper(SingletonClass):
    """Wrapper for dictionary API provider."""

    def __init__(self, cachePath, useCache=True, cfgOb=None, configName=None, **kwargs):
        """Wrapper for dictionary API provider.

        Args:
            cachePath (str): top path to contain the dictionary cache directory
            useCache (bool, optional): flag to use cached files. Defaults to True.
            cfgOb (object):  ConfigInfo() object instance
            configName (str): ConfigInfo() section

        """
        # self.__cfgOb = kwargs.get("cfgOb", None)
        # self.__configName = kwargs.get("configName", self.__cfgOb.getDefaultSectionName())
        self.__cfgOb = cfgOb
        self.__configName = configName
        self.__contentInfoConfigName = "content_info_helper_configuration"
        self.__dictLocatorMap = self.__cfgOb.get("DICT_LOCATOR_CONFIG_MAP", sectionName=self.__contentInfoConfigName)
        # self.__dictLocatorMap: dictionary of configured dictionary name mappings, e.g.:
        #   ordereddict([('pdbx_comp_model_core', ['PDBX_DICT_LOCATOR', 'RCSB_COMP_MODEL_DICT_LOCATOR', 'MA_DICT_LOCATOR']),
        #                ('pdbx', ['PDBX_DICT_LOCATOR', 'RCSB_DICT_LOCATOR']),
        #                ('pdbx_core', ['PDBX_DICT_LOCATOR', 'RCSB_DICT_LOCATOR', 'VRPT_DICT_LOCATOR']), ...
        dirPath = os.path.join(cachePath, self.__cfgOb.get("DICTIONARY_CACHE_DIR", sectionName=self.__configName))
        self.__dP = DictionaryApiProvider(dirPath, useCache=useCache, **kwargs)
        logger.debug("Leaving constructor")

    def testCache(self):
        return self.__cfgOb is not None

    def getApiByLocators(self, dictLocators, **kwargs):
        """Return a dictionary API object for the input dictionary locator list.

        Args:
            dictLocators (list str): list of dictionary locators

        Returns:
            (object): Instance of DictionaryApi()
        """
        return self.__dP.getApi(dictLocators, **kwargs)

    def getApiByName(self, databaseName, **kwargs):
        """Return a dictionary API object for the input schema name.

        Args:
            databaseName (str): database schema name

        Returns:
            (object): Instance of DictionaryApi()
        """
        if databaseName not in self.__dictLocatorMap:
            logger.error("Missing dictionary locator configuration for database schema %s", databaseName)
            dictLocators = []
        else:
            dictLocators = [self.__cfgOb.getPath(configLocator, sectionName=self.__configName) for configLocator in self.__dictLocatorMap[databaseName]]
            # Example dictLocators for databaseName 'pdbx_comp_model_core':
            #  ['https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets/development/dictionary_files/reference/mmcif_pdbx_v5_next.dic',
            #   'https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets/development/dictionary_files/dist/rcsb_mmcif_comp_model_ext.dic',
            #   'https://raw.githubusercontent.com/rcsb/py-rcsb_exdb_assets/development/dictionary_files/reference/mmcif_ma.dic']
        #
        logger.debug("Fetching dictionary API for %s using %r", databaseName, dictLocators)
        return self.__dP.getApi(dictLocators, **kwargs)
