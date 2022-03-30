##
# File:    DictMethodCompModelHelper.py
# Author:  J. Westbrook
# Date:    12-Oct-2021
# Version: 0.001 Initial version
#
#
# Updates:
#   07-Dec-2021  dwp Only add unassigned polymer entity-level taxonomy if _ma_target_ref_db_details.ncbi_taxonomy_id and .organism_scientific
#                    are both present, since these are not mandatory fields and thus may not be present in some cases
##
"""
Helper class implements computed model method references in the RCSB dictionary extension.

All data accessors and structures here refer to dictionary category and attribute names.

"""
__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

# pylint: disable=too-many-lines

import logging
from collections import defaultdict

from mmcif.api.DataCategory import DataCategory

logger = logging.getLogger(__name__)


class DictMethodCompModelHelper(object):
    """Helper class implements computed model method references in the RCSB dictionary extension."""

    aaFwDict3 = {
        "ALA": 89.093,
        "ARG": 175.209,
        "ASN": 132.118,
        "ASP": 133.103,
        "ASX": 100.096,
        "CYS": 121.158,
        "GLN": 146.144,
        "GLU": 147.129,
        "GLX": 114.123,
        "GLY": 75.067,
        "HIS": 156.162,
        "ILE": 131.173,
        "LEU": 131.173,
        "LYS": 147.195,
        "MET": 149.211,
        "PHE": 165.189,
        "PRO": 115.130,
        "SER": 105.093,
        "THR": 119.119,
        "TRP": 204.225,
        "TYR": 181.189,
        "VAL": 117.146,
        "PYL": 255.313,
        "SEC": 168.053,
    }

    def __init__(self, **kwargs):
        """
        Args:
            resourceProvider: (obj) instance of DictMethodResourceProvider()

        """
        #
        logger.debug("Dictionary entry method helper init with kwargs %r", kwargs)
        self._raiseExceptions = kwargs.get("raiseExceptions", False)
        #
        rP = kwargs.get("resourceProvider")
        self.__commonU = rP.getResource("DictMethodCommonUtils instance") if rP else None
        # dapw = rP.getResource("DictionaryAPIProviderWrapper instance") if rP else None
        # self.__dApi = dapw.getApiByName("pdbx_core") if dapw else None
        self.__dApi = kwargs.get("dictionaryApi", None)
        if self.__dApi:
            logger.debug("Loaded API for: %r", self.__dApi.getDictionaryTitle())
        else:
            logger.error("Missing dictionary API %r", kwargs)

        #
        # logger.debug("Dictionary entry method helper init")

    def echo(self, msg):
        logger.info(msg)

    def addPolymerEntityFormulaWeight(self, dataContainer, catName, **kwargs):
        """Add unassigned polymer entity-level formula weights.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name

        Returns:
            bool: True for success or False otherwise
        """
        logger.debug("Starting with %r %r %r", dataContainer.getName(), catName, kwargs)
        try:
            if not (dataContainer.exists("entity") and dataContainer.exists("entity_poly_seq") and dataContainer.exists("ma_data")):
                return False
            eObj = dataContainer.getObj("entity")
            epsObj = dataContainer.getObj("entity_poly_seq")
            fwD = defaultdict(float)
            for ii in range(epsObj.getRowCount()):
                monId = epsObj.getValue("mon_id", ii)
                entityId = epsObj.getValue("entity_id", ii)
                if monId not in DictMethodCompModelHelper.aaFwDict3:
                    logger.error("Unanticipated amino acid %r", monId)
                else:
                    fwD[entityId] += DictMethodCompModelHelper.aaFwDict3[monId]
            if not eObj.hasAttribute("formula_weight"):
                eObj.appendAttribute("formula_weight")
            for ii in range(eObj.getRowCount()):
                entityId = eObj.getValue("id", ii)
                if entityId in fwD:
                    eObj.setValue(round(fwD[entityId], 3), "formula_weight", ii)
            #
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def addPolymerEntityTaxonomy(self, dataContainer, catName, **kwargs):
        """Add unassigned polymer entity-level taxonomy (if both _ma_target_ref_db_details.ncbi_taxonomy_id and .organism_scientific are present).

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name

        Returns:
            bool: True for success or False otherwise

        Example:

            _ma_target_ref_db_details.db_accession        Q60392
            _ma_target_ref_db_details.db_code             Y084_METJA
            _ma_target_ref_db_details.db_name             UNP
            _ma_target_ref_db_details.ncbi_taxonomy_id    243232
            _ma_target_ref_db_details.organism_scientific "Methanocaldococcus jannaschii (strain ATCC 43067 / DSM 2661 / JAL-1 / JCM 10045 / NBRC 100440)"
            _ma_target_ref_db_details.seq_db_align_begin  1
            _ma_target_ref_db_details.seq_db_align_end    249
            _ma_target_ref_db_details.seq_db_isoform      ?
            _ma_target_ref_db_details.target_entity_id    1
        """
        atL = [
            "entity_id",
            "pdbx_organism_scientific",
            "nat_common_name",
            "pdbx_ncbi_taxonomy_id",
            "pdbx_src_id",
            "pdbx_beg_seq_num",
            "pdbx_end_seq_num",
            "rcsb_gene_name",
        ]
        logger.debug("Starting with %r %r %r", dataContainer.getName(), catName, kwargs)
        try:
            if not (dataContainer.exists("ma_data") and dataContainer.exists("ma_target_ref_db_details")):
                return False
            if not all([ai in dataContainer.getObj('ma_target_ref_db_details').getAttributeList() for ai in ["ncbi_taxonomy_id", "organism_scientific"]]):
                return False
            geneName = None
            if dataContainer.exists("af_target_ref_db_details"):
                tObj = dataContainer.getObj("af_target_ref_db_details")
                geneName = tObj.getValue("gene", 0)

            epLenD = self.__commonU.getPolymerEntityLengths(dataContainer)
            tObj = dataContainer.getObj("ma_target_ref_db_details")
            if not dataContainer.exists("entity_src_nat"):
                dataContainer.append(DataCategory("entity_src_nat", attributeNameList=atL))
            sObj = dataContainer.getObj("entity_src_nat")
            #
            for ii in range(tObj.getRowCount()):
                taxId = tObj.getValue("ncbi_taxonomy_id", ii)
                orgName = tObj.getValue("organism_scientific", ii)
                entityId = tObj.getValue("target_entity_id", ii)
                #
                sObj.setValue(entityId, "entity_id", ii)
                sObj.setValue(taxId, "pdbx_ncbi_taxonomy_id", ii)
                sObj.setValue(orgName, "pdbx_organism_scientific", ii)
                sObj.setValue("1", "pdbx_src_id", ii)
                sObj.setValue("1", "pdbx_beg_seq_num", ii)
                sObj.setValue(str(epLenD[entityId]), "pdbx_end_seq_num", ii)
                sObj.setValue(geneName, "rcsb_gene_name", ii)
            #
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False
