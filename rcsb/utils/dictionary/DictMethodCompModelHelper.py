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
#   02-Apr-2022   bv Add method 'consolidateGlobalQAScores'
#   29-Apr-2022  dwp Add method 'buildCompModelProvenance'
#   17-May-2022   bv Add method 'addStructInfo'
#   29-Jun-2022  dwp Update method 'buildCompModelProvenance' for populating rcsb_comp_model_provenance using mapping between internal and external IDs
#    6-Jul-2022   bv RO-3357: Fix taxonomy assignment in addPolymerEntityTaxonomy
#    6-Jul-2022  dwp Only populate rcsb_comp_model_provenance.source_url if it exists in CSM holdings file (may not exist for all AF fragments)
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
        self.__mcP = rP.getResource("ModelCacheProvider instance") if rP else None
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

            epLenD = self.__commonU.getPolymerEntityLengths(dataContainer)
            tObj = dataContainer.getObj("ma_target_ref_db_details")
            if not dataContainer.exists("entity_src_nat"):
                dataContainer.append(DataCategory("entity_src_nat", attributeNameList=atL))
            sObj = dataContainer.getObj("entity_src_nat")
            #
            jj = 0
            for ii in range(tObj.getRowCount()):
                taxId = tObj.getValue("ncbi_taxonomy_id", ii)
                orgName = tObj.getValue("organism_scientific", ii)
                entityId = tObj.getValue("target_entity_id", ii)
                geneName = tObj.getValueOrDefault("gene_name", ii, defaultValue=None)
                dbName = tObj.getValue("db_name", ii)
                #
                if dbName == "UNP":
                    sObj.setValue(entityId, "entity_id", jj)
                    sObj.setValue(taxId, "pdbx_ncbi_taxonomy_id", jj)
                    sObj.setValue(orgName, "pdbx_organism_scientific", jj)
                    sObj.setValue("1", "pdbx_src_id", jj)
                    sObj.setValue("1", "pdbx_beg_seq_num", jj)
                    sObj.setValue(str(epLenD[entityId]), "pdbx_end_seq_num", jj)
                    sObj.setValue(geneName, "rcsb_gene_name", jj)
                    jj += 1
            #
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def consolidateGlobalQAScores(self, dataContainer, catName, **kwargs):
        """Transform Global QA scores from ModelCIF to RCSB extension

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name

        Returns:
            bool: True for success or False otherwise

        Example:

            Data in ModelCIF:

            _entry.id       ma-bak-cepc-0001

            loop_
            _ma_qa_metric.id
            _ma_qa_metric.mode
            _ma_qa_metric.name
            _ma_qa_metric.type
            1   global  pLDDT   pLDDT
            2   global  ZDOPE   zscore

            loop_
            _ma_qa_metric_global.metric_id
            _ma_qa_metric_global.metric_value
            _ma_qa_metric_global.model_id
            _ma_qa_metric_global.ordinal_id
            1 81.23 1 1
            2 -2.0  1 2

            Transformed to:
            _rcsb_ma_qa_metric_global.entry_id                      ma-bak-cepc-0001
            _rcsb_ma_qa_metric_global.model_id                      1
            _rcsb_ma_qa_metric_global.ma_qa_global_metric_type      pLDDT,zscore
            _rcsb_ma_qa_metric_global.ma_qa_global_metric_name      pLDDT,ZDOPE
            _rcsb_ma_qa_metric_global.ma_qa_global_metric_value     81.23,-2.0

        """
        logger.debug("Starting with %r %r %r", dataContainer.getName(), catName, kwargs)
        try:
            # Check for all categories/data items required for the ETL (will return False for experimental models)
            if not dataContainer.exists("ma_qa_metric_global"):
                return False
            if not all([ai in dataContainer.getObj('ma_qa_metric_global').getAttributeList() for ai in ["model_id", "metric_id", "metric_value"]]):
                return False

            catName = "rcsb_ma_qa_metric_global"

            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=self.__dApi.getAttributeNameList(catName)))

            cObj = dataContainer.getObj(catName)

            if dataContainer.exists("ma_qa_metric_global"):
                bObj = dataContainer.getObj("ma_qa_metric_global")

            maQaMetricTypeD = self.__commonU.getMaQaMetricType(dataContainer)
            maQaMetricGlobalTypeD = maQaMetricTypeD["maQaMetricGlobalTypeD"]
            if not maQaMetricGlobalTypeD:
                return False

            modelIdList = bObj.getAttributeValueList("model_id")
            modelIdL = list(set(modelIdList))
            if not modelIdL:
                return False

            dObj = dataContainer.getObj("entry")
            entryId = dObj.getValue("id", 0)
            if not entryId:
                return False

            vD = {}
            mD = {}

            for ii, modelId in enumerate(modelIdL):
                cObj.setValue(entryId, "entry_id", ii)
                cObj.setValue(modelId, "model_id", ii)
                vD[modelId] = bObj.selectValuesWhere("metric_value", modelId, "model_id")
                mD[modelId] = bObj.selectValuesWhere("metric_id", modelId, "model_id")
            for ii in range(cObj.getRowCount()):
                modelId = cObj.getValue("model_id", ii)
                cObj.setValue(",".join(vD[modelId]), "ma_qa_metric_global_value", ii)
                cObj.setValue(",".join([str(maQaMetricGlobalTypeD[mId]["type"]) for mId in mD[modelId]]), "ma_qa_metric_global_type", ii)
                cObj.setValue(",".join([str(maQaMetricGlobalTypeD[mId]["name"]) for mId in mD[modelId]]), "ma_qa_metric_global_name", ii)

            return True
        except Exception as e:
            logger.exception("For %s populating rcsb_ma_qa_metric_global failing with %s", dataContainer.getName(), str(e))
        return False

    def buildCompModelProvenance(self, dataContainer, catName, **kwargs):
        """Load the input category with _rcsb_comp_model_provenance content.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name

        Returns:
            bool: True for success or False otherwise

        For example:

        _rcsb_comp_model_provenance.source_filename
        _rcsb_comp_model_provenance.source_url
        _rcsb_comp_model_provenance.source_db
        _rcsb_comp_model_provenance.entry_id
        ...

        """
        try:
            logger.debug("Starting with %r %r %r", dataContainer.getName(), catName, kwargs)
            if not (dataContainer.exists("entry") and dataContainer.exists("ma_data")):
                return False

            if catName != "rcsb_comp_model_provenance":
                logger.warning("input catName (%s) not 'rcsb_comp_model_provenance'")

            catName = "rcsb_comp_model_provenance"

            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=self.__dApi.getAttributeNameList(catName)))

            cObj = dataContainer.getObj(catName)

            tObj = dataContainer.getObj("entry")
            entryId = tObj.getValue("id", 0)

            # NOTE: Old approach - from when we were using the external/source ID as entry.id
            # compModelId = self.__mcP.getInternalCompModelId(entryId)
            # if not compModelId:
            #     logger.error("Unable to map computed-model entryId/sourceId (%s) to internal identifier - sourceId key not found.", entryId)
            #     return False

            compModelD = self.__mcP.getCompModelData(entryId)

            cObj.setValue(compModelD["sourceId"], "entry_id", 0)
            cObj.setValue(compModelD["sourceDb"], "source_db", 0)
            sourceModelUrl = compModelD.get("sourceModelUrl", None)
            if sourceModelUrl:
                cObj.setValue(compModelD["sourceModelUrl"], "source_url", 0)
            cObj.setValue(compModelD["sourceModelFileName"], "source_filename", 0)
            sourceModelPaeUrl = compModelD.get("sourceModelPaeUrl", None)
            if sourceModelPaeUrl:
                cObj.setValue(compModelD["sourceModelPaeUrl"], "source_pae_url", 0)

            return True

        except Exception as e:
            logger.exception("For %s failing with %s", catName, str(e))

        return False

    def addStructInfo(self, dataContainer, catName, **kwargs):
        """Add missing struct table

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name

        Returns:
            bool: True for success or False otherwise

        """
        logger.debug("Starting with %r %r %r", dataContainer.getName(), catName, kwargs)
        try:
            # Check if struct table is already present
            if dataContainer.exists("struct"):
                return False
            # Check for all categories/data items required for the ETL (will return False for experimental models)
            if not (dataContainer.exists("entry") and dataContainer.exists("entity") and dataContainer.exists("ma_data")):
                return False

            logger.debug("Adding struct table for %s", dataContainer.getName())
            catName = "struct"
            dataContainer.append(DataCategory(catName, attributeNameList=self.__dApi.getAttributeNameList(catName)))
            xObj = dataContainer.getObj(catName)

            tObj = dataContainer.getObj("entry")
            entryId = tObj.getValue("id", 0)

            methodType = "computational"
            tL = []
            qObj = dataContainer.getObj("entity")

            for ii in range(qObj.getRowCount()):
                pdbxDesc = qObj.getValue("pdbx_description", ii)
                eType = qObj.getValue("type", ii)
                if pdbxDesc not in [".", "?"] and eType.lower() == "polymer":
                    tL.append(str(pdbxDesc))

            if tL:
                strTitle = "Computed structure model of " + ", ".join([str(x) for x in tL])
            else:
                logger.debug("Missing entity description for %s. Using generic title", dataContainer.getName())
                strTitle = "Computed structure model"

            xObj.setValue(entryId, "entry_id", 0)
            xObj.setValue(strTitle, "title", 0)
            xObj.setValue(methodType, "pdbx_structure_determination_methodology", 0)

            return True
        except Exception as e:
            logger.exception("For %s populating struct failing with %s", dataContainer.getName(), str(e))
        return False
