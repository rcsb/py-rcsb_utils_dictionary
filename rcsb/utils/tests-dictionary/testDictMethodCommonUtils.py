##
# File:    testDictMethodCommonUtils.py
# Author:  D. Piehl
# Date:    28-Oct-2026
# Version: 0.001
#
# Update:
##
"""
Tests for selected DictMethodCommonUtils methods that don't need a data container.

"""

__docformat__ = "google en"
__author__ = "Dennis Piehl"
__email__ = "dennis.piehl@rcsb.org"
__license__ = "Apache 2.0"

import logging
import unittest

from rcsb.utils.dictionary.DictMethodCommonUtils import DictMethodCommonUtils

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s]-%(module)s.%(funcName)s: %(message)s")
logger = logging.getLogger()


class DictMethodCommonUtilsTests(unittest.TestCase):
    def setUp(self):
        self.__commonU = DictMethodCommonUtils()

    def testGetExtAndShortIds(self):
        """Test extended and short ID forms for short, extended, CSM and other entry IDs"""
        for entryId, expected in [
            ("4HHB", ("pdb_00004hhb", "4HHB")),
            ("4hhb", ("pdb_00004hhb", "4HHB")),
            ("pdb_00004hhb", ("pdb_00004hhb", "4HHB")),
            ("PDB_00004HHB", ("pdb_00004hhb", "4HHB")),
            ("pdb_00011abc", ("pdb_00011abc", None)),
            ("AF_AFP0CB50F1", ("AF_AFP0CB50F1", None)),
            ("ma_macoffeslacc102465g1i2", ("MA_MACOFFESLACC102465G1I2", None)),
            ("PDBDEV_00000020", ("PDBDEV_00000020", None)),
        ]:
            self.assertEqual(self.__commonU.getExtAndShortIds(entryId), expected, msg=entryId)


def commonUtilsSuite():
    suiteSelect = unittest.TestSuite()
    suiteSelect.addTest(DictMethodCommonUtilsTests("testGetExtAndShortIds"))
    return suiteSelect


if __name__ == "__main__":
    mySuite = commonUtilsSuite()
    unittest.TextTestRunner(verbosity=2).run(mySuite)
