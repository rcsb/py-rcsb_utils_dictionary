##
# File:    DictMethodCommonUtils.py
# Author:  J. Westbrook
# Date:    16-Jul-2019
# Version: 0.001 Initial version
#
# Updates:
# 26-Jul-2019 jdw Include struct_mon_prot_cis with secondary structure features
#                 Add general processing of intermolecular and other connections.
# 19-Sep-2019 jdw Add method getEntityReferenceAlignments()
# 13-Oct-2019 jdw add isoform support
# 17-Feb-2021 jdw add non-polymer neighbor calculation
# 19-Jan-2022 dwp add method 'filterStructureDeterminationMethodType'to return type of method used
#                 (experimental, computational, or integrative);
#                 remove 'THEORETICAL MODEL' from list of 'Other' experimental_method types, and
#                 add it to computational method list
# 28-Mar-2022 bv add method 'getRepresentativeModels' to get representative models for NMR ensembles
#                Fix pylint issues
#  2-Apr-2022 bv Add methods 'getCompModelDb2L', 'getMaQaMetricType', and 'getCompModelLocalScores'
# 20-Apr-2022 bv Update method 'getCompModelDb2L'
# 29-Apr-2022 dwp Use internal computed-model identifiers for 'rcsb_id'
#  3-May-2022 dwp Use internal computed-model identifiers for 'entry_id' in containter_identifiers
# 27-Jun-2022 bv  Update _rcsb_ma_qa_metric_global.ma_qa_metric_global_type to 'pLDDT' for AF models
# 29-Jun-2022 dwp Use internal computed-model identifiers everywhere (in same manner as experimental models)
#  3-Jul-2023 aae Update __getInstanceModelOutliers (old version is backed up as __getInstanceModelOutliersXML),
#                 add method __getValidationData and getLocalValidationData to get data from validation mmcif files
# 01-Feb-2024 bv  Add method 'getInstanceDeuWatMolCounts' to support deuterated water molecule count
#                 Update methods 'getDepositedAtomCounts' and '__getAtomSiteInfo'
# 18-Mar-2024 dwp Add method 'getPolymerEntityReferenceAlignments' to enable retrieval of all UniProt IDs
# 24-Jul-2024 dwp Adjust ligand interaction calculation and provider to not rely on struct_conn and instead
#                 base calculation on coordinates only
#  7-Jan-2025  bv Update '__getInstanceModelOutliers' to handle validation data
#                 Add getRepresentativeModelId and getMethodList
# 16-Jan-2025 dwp Consolidate getRepresentativeModelId and associated methods/calls;
#                 Skip non-representative models for calculating/building content;
#                 Fix modelId assignment in getNeighborInfo
# 13-Feb-2025  bv Add methods to support integrative structures
#
##
"""
Helper class implements common utility external method references supporting the RCSB dictionary extension.

"""
__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

# pylint: disable=too-many-lines

import datetime
import itertools
import logging
import re
import sys
import time
import copy
from collections import OrderedDict, namedtuple, defaultdict
from operator import itemgetter

import numpy as np
from scipy import spatial

from rcsb.utils.io.CacheUtils import CacheUtils
from rcsb.utils.seq.SeqAlign import SeqAlign

logger = logging.getLogger(__name__)

OutlierValueFields = ("compId", "seqId", "outlierType", "description", "reported", "reference", "uncertaintyValue", "uncertaintyType")
OutlierValue = namedtuple("OutlierValue", OutlierValueFields, defaults=(None,) * len(OutlierValueFields))

BoundEntityFields = ("targetCompId", "connectType", "partnerCompId", "partnerEntityId", "partnerEntityType")
NonpolymerBoundEntity = namedtuple("NonpolymerBoundEntity", BoundEntityFields, defaults=(None,) * len(BoundEntityFields))

BoundInstanceFields = (
    "targetCompId",
    "targetAtomId",
    "targetAltId",
    "connectType",
    "partnerEntityType",
    "partnerEntityId",
    "partnerCompId",
    "partnerAsymId",
    "partnerSeqId",
    "partnerAuthSeqId",
    "partnerAtomId",
    "partnerAltId",
    "bondDistance",
    "bondOrder",
    "role",
)
NonpolymerBoundInstance = namedtuple("NonpolymerBoundInstance", BoundInstanceFields, defaults=(None,) * len(BoundInstanceFields))

NonpolymerValidationFields = (
    "rsr",
    "rscc",
    "nAtomsEds",
    "mogul_bonds_rmsz",
    "mogul_angles_rmsz",
    "numAnglesRmsZ",
    "numBondsRmsZ",
    "avgOccupancy",
    "intermolecular_clashes",
    "mogul_bond_outliers",
    "mogul_angle_outliers",
    "stereo_outliers",
)
NonpolymerValidationInstance = namedtuple("NonpolymerValidationInstance", NonpolymerValidationFields, defaults=(None,) * len(NonpolymerValidationFields))

LigandTargetFields = (
    "ligandModelId",
    "ligandAsymId",
    "ligandCompId",
    "ligandAtomId",
    "ligandAltId",
    "connectType",
    "partnerModelId",
    "partnerEntityType",
    "partnerEntityId",
    "partnerCompId",
    "partnerAsymId",
    "partnerSeqId",
    "partnerAuthSeqId",
    "partnerAtomId",
    "partnerAltId",
    "distance",
)
LigandTargetInstance = namedtuple("LigandTargetInstance", LigandTargetFields, defaults=(None,) * len(LigandTargetFields))

ReferenceFields = ("entityId", "entityType", "asymId", "compId", "seqId", "authSeqId", "atomId", "altId", "modelId")
ReferenceInstance = namedtuple("ReferenceInstance", ReferenceFields, defaults=(None,) * len(ReferenceFields))


class DictMethodCommonUtils(object):
    """Helper class implements common utility external method references supporting the RCSB dictionary extension."""

    # Dictionary of current standard monomers -
    aaDict3 = {
        "ALA": "A",
        "ARG": "R",
        "ASN": "N",
        "ASP": "D",
        "ASX": "B",
        "CYS": "C",
        "GLN": "Q",
        "GLU": "E",
        "GLX": "Z",
        "GLY": "G",
        "HIS": "H",
        "ILE": "I",
        "LEU": "L",
        "LYS": "K",
        "MET": "M",
        "PHE": "F",
        "PRO": "P",
        "SER": "S",
        "THR": "T",
        "TRP": "W",
        "TYR": "Y",
        "VAL": "V",
        "PYL": "O",
        "SEC": "U",
    }
    dnaDict3 = {"DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U", "DI": "I"}
    rnaDict3 = {"A": "A", "C": "C", "G": "G", "I": "I", "N": "N", "T": "T", "U": "U"}
    # "UNK": "X",
    # "MSE":"M",
    # ".": "."
    monDict3 = {**aaDict3, **dnaDict3, **rnaDict3}

    def __init__(self, **kwargs):
        """
        Args:
            **kwargs: (dict)  Placeholder for future key-value arguments

        """
        #
        self._raiseExceptions = kwargs.get("raiseExceptions", False)
        self.__wsPattern = re.compile(r"\s+", flags=re.UNICODE | re.MULTILINE)
        self.__reNonDigit = re.compile(r"[^\d]+")
        #
        cacheSize = 8  # Set to max number of procs
        self.__entityAndInstanceMapCache = CacheUtils(size=cacheSize, label="instance mapping")
        self.__atomInfoCache = CacheUtils(size=cacheSize, label="atom site counts and mapping")
        self.__instanceConnectionCache = CacheUtils(size=cacheSize, label="instance connections")
        self.__entityReferenceSequenceDetailsCache = CacheUtils(size=cacheSize, label="entity reference sequence details")
        self.__entitySequenceFeatureCache = CacheUtils(size=cacheSize, label="entity sequence features")
        self.__instanceSiteInfoCache = CacheUtils(size=cacheSize, label="instance site details")
        self.__instanceUnobservedCache = CacheUtils(size=cacheSize, label="instance unobserved details")
        self.__modelOutliersCache = CacheUtils(size=cacheSize, label="model outlier details")
        self.__neighborInfoCache = CacheUtils(size=cacheSize, label="ligand and target nearest neighbors")
        self.__localValidationCache = CacheUtils(size=cacheSize, label="local validation data")
        #
        logger.debug("Dictionary common utilities init")

    def echo(self, msg):
        logger.info(msg)

    def testCache(self):
        return len(DictMethodCommonUtils.aaDict3) > 20

    def isFloat(self, val):
        try:
            float(val)
        except Exception:
            return False
        return True

    def __fetchEntityAndInstanceTypes(self, dataContainer):
        wD = self.__entityAndInstanceMapCache.get(dataContainer.getName())
        if not wD:
            wD = self.__getEntityAndInstanceTypes(dataContainer)
            self.__entityAndInstanceMapCache.set(dataContainer.getName(), wD)
        return wD

    def getFormulaWeightNonSolvent(self, dataContainer):
        """Return a formula weight of the non-solvent entities in the deposited entry.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            float: formula weight (kilodaltons)
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["fwNonSolvent"] if "fwNonSolvent" in wD else {}

    def getInstancePolymerTypes(self, dataContainer):
        """Return a dictionary of polymer types for each polymer instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {'asymId': <dictionary polymer type>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["instancePolymerTypeD"] if "instancePolymerTypeD" in wD else {}

    def getInstanceTypes(self, dataContainer):
        """Return a dictionary of entity types for each entity instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {'asymId': <entity type>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["instanceTypeD"] if "instanceTypeD" in wD else {}

    def getInstanceTypeCounts(self, dataContainer):
        """Return a dictionary of the counts entity types for each entity type.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {'entity type': <# of instances>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["instanceTypeCountD"] if "instanceTypeCountD" in wD else {}

    def getInstanceEntityMap(self, dataContainer):
        """Return a dictionary of entities corresponding to each entity instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {'asymId': <entity id>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["instEntityD"] if "instEntityD" in wD else {}

    def getEntityPolymerTypes(self, dataContainer):
        """Return a dictionary of polymer types for each polymer entity.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {'entityId': <dictionary polymer types>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["epTypeD"] if "epTypeD" in wD else {}

    def getEntityTypes(self, dataContainer):
        """Return a dictionary of entity types for each entity.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {'entityId': <dictionary entity types>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["eTypeD"] if "eTypeD" in wD else {}

    def getPolymerEntityFilteredTypes(self, dataContainer):
        """Return a dictionary of filtered entity polymer types for each polymer entity.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {'entityId': <filtered entity polymer types>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["epTypeFilteredD"] if "epTypeFilteredD" in wD else {}

    def getPolymerEntityLengths(self, dataContainer):
        """Return a dictionary of entity polymer lengths for each polymer entity.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {'entityId': <monomer length>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["epLengthD"] if "epLengthD" in wD else {}

    def getPolymerEntityLengthsEnumerated(self, dataContainer):
        """Return a dictionary of entity polymer lengths for each polymer entity.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {'entityId': <monomer length>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["entityPolymerLengthD"] if "entityPolymerLengthD" in wD else {}

    def getPolymerEntityMonomerCounts(self, dataContainer):
        """Return a dictionary of monomer counts for each polymer entity.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {'entityId': {'compId': <monomer count>, ... }}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["entityPolymerMonomerCountD"] if "entityPolymerMonomerCountD" in wD else {}

    def getPolymerEntityModifiedMonomers(self, dataContainer):
        """Return a dictionary of nonstandard monomers for each polymer entity.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {'entityId': [mod_comp_id, mod_comp_id,...]}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["entityPolymerModifiedMonomers"] if "entityPolymerModifiedMonomers" in wD else {}

    def getPolymerModifiedMonomerFeatures(self, dataContainer):
        """Return a dictionary of nonstandard monomer features.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: [(entityId, seqId, compId, 'modified_monomer')] = set(compId)

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["seqModMonomerFeatureD"] if "seqModMonomerFeatureD" in wD else {}

    def getEntityPolymerLengthBounds(self, dataContainer):
        """Return a dictionary of polymer length bounds by entity type.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            tuple: (minLen, maxLen)
        """
        if not dataContainer or not dataContainer.getName():
            return ()
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["entityPolymerLengthBounds"] if "entityPolymerLengthBounds" in wD else ()

    def getEntityFormulaWeightBounds(self, dataContainer):
        """Return a dictionary of formula weight bounds by entity type.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: [entityType] = (minFw, maxFw)
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["fwTypeBoundD"] if "fwTypeBoundD" in wD else {}

    def getTargetComponents(self, dataContainer):
        """Return the author identified components targeted in the current entry.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            list: [compId, compId,...]
        """
        if not dataContainer or not dataContainer.getName():
            return []
        wD = self.__fetchEntityAndInstanceTypes(dataContainer)
        return wD["ccTargets"] if "ccTargets" in wD else []

    def __getEntityAndInstanceTypes(self, dataContainer):
        """Internal method to collect and return entity/instance type, size and mapping information.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            (dict) : Return dictionary of entity types, type counts and polymer type (where applicable) for
                     each instance in the deposited unit.

            Type and count contents:

              instanceTypeD[asymId] = <entity_type>
              instanceTypeCountD[<entity_type>] = #
              instancePolymerTypeD[asymId] = <filtered polymer type>
              eTypeD[entityId] = <dictionary entity type>
              instEntityD[asymId] = entityId
              epTypeD[entityId] = <dictionary polymer type>
              epTypeFilteredD[entityId] = <dictionary polymer type>
              epLengthD[entityId] = polymer monomer length (from one-letter-code)
              entityPolymerLengthD[entityId] = polymer monomer length (from enumerated sequence)
              entityPolymerMonomerCountD[entityId][compId] = mononer count
              entityPolymerModifiedMonomers[entity]=[mod compId, mod compId]
              seqModMonomerFeatureD[(entityId, seqId, compId, 'modified_monomer')] = set(compId)
              fwNonSolvent = float value (kilodaltons)
              fwTypeBoundD[entityType] = (minFw, maxFw)
              entityPolymerLengthBounds = (minL, maxL)
              ccTargets = [compId, compId]
        """
        rD = {}
        #
        try:
            #
            if not dataContainer.exists("entity") or not dataContainer.exists("struct_asym"):
                return {}
            logger.debug("Starting for container %s", dataContainer.getName())
            eFwD = {}
            instanceTypeD = {}
            instancePolymerTypeD = {}
            instanceTypeCountD = {}
            #
            eObj = dataContainer.getObj("entity")
            eTypeD = {}
            for ii in range(eObj.getRowCount()):
                # logger.info("Attribute %r %r" % (ii, eObj.getAttributeList()))
                entityId = eObj.getValue("id", ii)
                eType = eObj.getValue("type", ii)
                eTypeD[entityId] = eType
                fw = eObj.getValue("formula_weight", ii)
                eFwD[entityId] = float(fw) if fw and fw not in [".", "?"] else 0.0
            #
            epTypeD = {}
            epLengthD = {}
            epTypeFilteredD = {}
            hasEntityPoly = False
            if dataContainer.exists("entity_poly"):
                hasEntityPoly = True
                epObj = dataContainer.getObj("entity_poly")
                for ii in range(epObj.getRowCount()):
                    entityId = epObj.getValue("entity_id", ii)
                    pType = epObj.getValue("type", ii)
                    epTypeFilteredD[entityId] = self.filterEntityPolyType(pType)
                    epTypeD[entityId] = pType
                    if epObj.hasAttribute("pdbx_seq_one_letter_code_can"):
                        sampleSeq = self.__stripWhiteSpace(epObj.getValue("pdbx_seq_one_letter_code_can", ii))
                        epLengthD[entityId] = len(sampleSeq) if sampleSeq and sampleSeq not in ["?", "."] else None

            #
            seqModMonomerFeatureD = {}
            entityPolymerMonomerCountD = {}
            entityPolymerLengthD = {}
            hasEntityPolySeq = False
            epsObj = None
            if dataContainer.exists("entity_poly_seq"):
                epsObj = dataContainer.getObj("entity_poly_seq")
                hasEntityPolySeq = True
                tSeqD = {}
                for ii in range(epsObj.getRowCount()):
                    entityId = epsObj.getValue("entity_id", ii)
                    seqNum = epsObj.getValue("num", ii)
                    compId = epsObj.getValue("mon_id", ii)
                    if compId not in DictMethodCommonUtils.monDict3:
                        seqModMonomerFeatureD.setdefault((entityId, seqNum, compId, "modified_monomer"), set()).add(compId)
                    # handle heterogeneity with the entityId,seqNum tuple
                    tSeqD.setdefault(entityId, set()).add((entityId, seqNum))
                    if entityId not in entityPolymerMonomerCountD:
                        entityPolymerMonomerCountD[entityId] = {}
                    entityPolymerMonomerCountD[entityId][compId] = entityPolymerMonomerCountD[entityId][compId] + 1 if compId in entityPolymerMonomerCountD[entityId] else 1
                #
                entityPolymerLengthD = {entityId: len(tSet) for entityId, tSet in tSeqD.items()}
            #
            if not hasEntityPoly and hasEntityPolySeq:
                for entityId, eType in eTypeD.items():
                    if eType in ["polymer"]:
                        monomerL = epsObj.selectValuesWhere("mon_id", entityId, "entity_id")
                        pType, fpType = self.guessEntityPolyTypes(monomerL)
                        epTypeFilteredD[entityId] = fpType
                        epTypeD[entityId] = pType
                        epLengthD[entityId] = len(monomerL)

            entityPolymerModifiedMonomers = {}
            for entityId, cD in entityPolymerMonomerCountD.items():
                tL = []
                for compId, _ in cD.items():
                    modFlag = "N" if compId in DictMethodCommonUtils.monDict3 else "Y"
                    if modFlag == "Y":
                        tL.append(compId)
                entityPolymerModifiedMonomers[entityId] = sorted(set(tL))
            #
            logger.debug("%s entityPolymerModifiedMonomers %r", dataContainer.getName(), entityPolymerModifiedMonomers)
            #  Add branched here
            #
            instEntityD = {}
            sObj = dataContainer.getObj("struct_asym")
            for ii in range(sObj.getRowCount()):
                entityId = sObj.getValue("entity_id", ii)
                asymId = sObj.getValue("id", ii)
                instEntityD[asymId] = entityId
                if entityId in eTypeD:
                    instanceTypeD[asymId] = eTypeD[entityId]
                else:
                    logger.warning("Missing entity id entry %r asymId %r entityId %r", dataContainer.getName(), entityId, asymId)
                if entityId in epTypeD:
                    instancePolymerTypeD[asymId] = epTypeFilteredD[entityId]
                #
            #
            # Count the instance by type - initialize all types
            #
            instanceTypeCountD = {k: 0 for k in ["polymer", "non-polymer", "branched", "macrolide", "water"]}
            for asymId, eType in instanceTypeD.items():
                instanceTypeCountD[eType] += 1
            #
            # Compute the total weight of polymer and non-polymer instances (full entities) - (kilodaltons)
            #
            fwNonSolvent = 0.0
            for asymId, eType in instanceTypeD.items():
                if eType not in ["water"]:
                    entityId = instEntityD[asymId]
                    fwNonSolvent += eFwD[entityId]
            fwNonSolvent = fwNonSolvent / 1000.0
            #
            # Get ligand of interest.
            #
            ccTargets = []
            if dataContainer.exists("pdbx_entity_instance_feature"):
                ifObj = dataContainer.getObj("pdbx_entity_instance_feature")
                for ii in range(ifObj.getRowCount()):
                    compId = ifObj.getValue("comp_id", ii)
                    ft = ifObj.getValue("feature_type", ii)
                    if ft.upper() in ["SUBJECT OF INVESTIGATION"]:
                        ccTargets.append(compId)
            #
            #
            fwTypeBoundD = {}
            tBoundD = {et: {"min": float("inf"), "max": -1.0} for eId, et in eTypeD.items()}
            for entityId, fw in eFwD.items():
                fw = fw / 1000.0
                eType = eTypeD[entityId]
                tBoundD[eType]["min"] = fw if fw < tBoundD[eType]["min"] else tBoundD[eType]["min"]
                tBoundD[eType]["max"] = fw if fw > tBoundD[eType]["max"] else tBoundD[eType]["max"]
            for eType in tBoundD:
                if tBoundD[eType]["min"] > 0.00000001:
                    fwTypeBoundD[eType] = tBoundD[eType]
            #

            entityPolymerLengthBounds = None
            maxL = -1
            minL = sys.maxsize
            if epLengthD:
                for entityId, pLen in epLengthD.items():
                    minL = pLen if pLen < minL else minL
                    maxL = pLen if pLen > maxL else maxL
                entityPolymerLengthBounds = (minL, maxL)
            #

            rD = {
                "instanceTypeD": instanceTypeD,
                "instancePolymerTypeD": instancePolymerTypeD,
                "instanceTypeCountD": instanceTypeCountD,
                "instEntityD": instEntityD,
                "eTypeD": eTypeD,
                "epLengthD": epLengthD,
                "epTypeD": epTypeD,
                "epTypeFilteredD": epTypeFilteredD,
                "entityPolymerMonomerCountD": entityPolymerMonomerCountD,
                "entityPolymerLengthD": entityPolymerLengthD,
                "entityPolymerModifiedMonomers": entityPolymerModifiedMonomers,
                "seqModMonomerFeatureD": seqModMonomerFeatureD,
                "fwNonSolvent": fwNonSolvent,
                "fwTypeBoundD": fwTypeBoundD,
                "entityPolymerLengthBounds": entityPolymerLengthBounds,
                "ccTargets": ccTargets,
            }
            logger.debug("%s length struct_asym %d (%d) instanceTypeD %r", dataContainer.getName(), sObj.getRowCount(), len(instanceTypeD), instanceTypeD)
        #
        except Exception as e:
            logger.exception("Failing with %r with %r", dataContainer.getName(), str(e))
        #
        return rD

    def getAsymAuthIdMap(self, dataContainer):
        """Return a dictionary of mapping between asymId and authAsymId.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {'asymId': authAsymId, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer)
        return wD["asymAuthIdD"] if "asymAuthIdD" in wD else {}

    def getInstanceHeavyAtomCounts(self, dataContainer, modelId="1"):
        """Return a dictionary of deposited heavy atom counts for each entity instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance
            modelId (str, optional): model index. Defaults to "1".


        Returns:
            dict: {'asymId': <# of deposited atoms>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer, modelId=modelId)
        return wD["instanceHeavyAtomCountD"] if "instanceHeavyAtomCountD" in wD else {}

    def getInstanceDeuWatMolCounts(self, dataContainer, modelId="1"):
        """Return a dictionary of deuterated water molecule counts for each entity instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance
            modelId (str, optional): model index. Defaults to "1".


        Returns:
            dict: {'asymId': <# of deuterated water molecules>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer, modelId=modelId)
        return wD["instanceDeuWatMolCountD"] if "instanceDeuWatMolCountD" in wD else {}

    def getInstanceHydrogenAtomCounts(self, dataContainer, modelId="1"):
        """Return a dictionary of deposited hydrogen atom counts for each entity instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance
            modelId (str, optional): model index. Defaults to "1".


        Returns:
            dict: {'asymId': <# of deposited atoms>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer, modelId=modelId)
        return wD["instanceHydrogenAtomCountD"] if "instanceHydrogenAtomCountD" in wD else {}

    def getModelIdList(self, dataContainer):
        """Return a list of model identifiers for the entry.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            list: [1,2,3]
        """
        if not dataContainer or not dataContainer.getName():
            return []
        wD = self.__fetchAtomSiteInfo(dataContainer)
        return wD["modelIdList"] if "modelIdList" in wD else []

    def getEntityTypeHeavyAtomCounts(self, dataContainer, modelId="1"):
        """Return a dictionary of deposited heavy atom counts for each entity type.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance
            modelId (str, optional): model index. Defaults to "1".

        Returns:
            dict: {'entity type': <# of deposited atoms>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer, modelId=modelId)
        return wD["typeHeavyAtomCountD"] if "typeHeavyAtomCountD" in wD else {}

    def getInstanceModeledMonomerCounts(self, dataContainer, modelId="1"):
        """Return a dictionary of deposited modeled monomer counts for each entity instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance
            modelId (str, optional): model index. Defaults to "1".

        Returns:
            dict: {'asymId': <# of deposited modeled monomers>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer, modelId=modelId)
        return wD["instancePolymerModeledMonomerCountD"] if "instancePolymerModeledMonomerCountD" in wD else {}

    def getInstanceUnModeledMonomerCounts(self, dataContainer, modelId="1"):
        """Return a dictionary of deposited unmodeled monomer counts for each entity instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance
            modelId (str, optional): model index. Defaults to "1".

        Returns:
            dict: {'asymId': <# of deposited unmodeled mononmers>, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer, modelId=modelId)
        return wD["instancePolymerUnmodeledMonomerCountD"] if "instancePolymerUnmodeledMonomerCountD" in wD else {}

    def getDepositedMonomerCounts(self, dataContainer, modelId="1"):
        """Return deposited modeled and unmodeled polymer monomer counts for the input modelid.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance
            modelId (str, optional): model index. Defaults to "1".


        Returns:
            (int,int):  modeled and unmodeled monomer counts
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer, modelId=modelId)
        modeledCount = sum(wD["instancePolymerModeledMonomerCountD"].values())
        unModeledCount = sum(wD["instancePolymerUnmodeledMonomerCountD"].values())
        return modeledCount, unModeledCount

    def getDepositedAtomCounts(self, dataContainer, modelId="1"):
        """Return the number of deposited heavy atoms in the input model, the total deposited atom
        and the total model count.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance
            modelId (str, optional): model index. Defaults to "1".

        Returns:
            (int, int, int, int)  deposited heavy atoms in input model, hydrogen atoms in input model, total deposited atom count, and total deposited model count
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer, modelId=modelId)
        numHeavyAtomsModel = wD["numHeavyAtomsModel"] if "numHeavyAtomsModel" in wD else 0
        numDeuWatMolModel = wD["numDeuWatMolModel"] if "numDeuWatMolModel" in wD else 0
        numHydrogenAtomsModel = wD["numHydrogenAtomsModel"] if "numHydrogenAtomsModel" in wD else 0
        numAtomsTotal = wD["numAtomsAll"] if "numAtomsAll" in wD else 0
        numModelsTotal = wD["numModels"] if "numModels" in wD else 0
        return numHeavyAtomsModel, numHydrogenAtomsModel, numAtomsTotal, numModelsTotal, numDeuWatMolModel

    def getInstancePolymerRanges(self, dataContainer):
        """Return a dictionary of polymer residue range and length for each entity instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {"asymId": , {"sampleSeqLen": sampleSeqLen,
                                "obsSeqLen": obsSeqLen,
                                "begSeqId": begSeqId,
                                "endSeqId": endSeqId,
                                "begAuthSeqId": begAuthSeqId,
                                "endAuthSeqId": endAuthSeqId,
                                "begInsCode": begAuthInsCode,
                                "endInsCode": endAuthInsCode,}...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer)
        return wD["asymIdPolymerRangesD"] if "asymIdPolymerRangesD" in wD else {}

    def getInstanceIdMap(self, dataContainer):
        """Return a dictionary of cardinal identifiers for each entity instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {"asymId":  {"entry_id": entryId,
                                "entity_id": entityId,
                                "entity_type": entityTypeD[entityId],
                                "asym_id": asymId,
                                "auth_asym_id": authAsymId,
                                "comp_id": monId,
                                "auth_seq_id": "?",}, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer)
        return wD["instanceIdMapD"] if "instanceIdMapD" in wD else {}

    def getNonPolymerIdMap(self, dataContainer):
        """Return a dictionary of cardinal identifiers for each non-polymer entity instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(authAsymId, resNum):   {"entry_id": entryId,
                                            "entity_id": entityId,
                                            "entity_type": entityTypeD[entityId],
                                            "asym_id": asymId,
                                            "auth_asym_id": authAsymId,
                                            "comp_id": monId,
                                            "auth_seq_id": resNum,
                                            }, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer)
        return wD["npAuthAsymIdMapD"] if "npAuthAsymIdMapD" in wD else {}

    def getPolymerIdMap(self, dataContainer):
        """Return a dictionary of cardinal identifiers for each polymer entity instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(authAsymId, authSeqId, insCode): {
                        "entry_id": entryId,
                        "entity_id": entityId,
                        "entity_type": entityTypeD[entityId],
                        "asym_id": asymId,
                        "comp_id": compId,
                        "seq_id": seqId,
                    }, ... }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer)
        return wD["pAuthAsymIdMapD"] if "pAuthAsymIdMapD" in wD else {}

    def getBranchedIdMap(self, dataContainer):
        """Return a dictionary of cardinal identifiers for each branched entity instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict:  {(authAsymId, authSeqNum): {
                        "entry_id": entryId,
                        "entity_id": entityId,
                        "entity_type": entityTypeD[entityId],
                        "asym_id": asymId,
                        "auth_asym_id": authAsymId,
                        "comp_id": monId,
                        "auth_seq_id": authSeqNum,
                        "seq_num": seqNum,
                    }, ...}

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer)
        return wD["brAuthAsymIdMapD"] if "brAuthAsymIdMapD" in wD else {}

    def getEntityTypeUniqueIds(self, dataContainer):
        """Return a nested dictionary of selected unique identifiers for entity types.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict:  [<entity_type>][<entity_id>] = {'asymIds': [...],'authAsymIds': [...], 'ccIds': [...]}


        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer)
        return wD["entityTypeUniqueIds"] if "entityTypeUniqueIds" in wD else {}

    def getAuthToSeqIdMap(self, dataContainer):
        """Return an instance (asymId) dictionary of auth to entity residue sequence mapping

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict:   seqIdMapAsymD[asymId] = [<authSeqId + insCode>, ... ]
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer)
        return wD["seqIdMapAsymD"] if "seqIdMapAsymD" in wD else {}

    def __fetchAtomSiteInfo(self, dataContainer, modelId="1"):
        wD = self.__atomInfoCache.get((dataContainer.getName(), modelId))
        if not wD:
            wD = self.__getAtomSiteInfo(dataContainer, modelId=modelId)
            self.__atomInfoCache.set((dataContainer.getName(), modelId), wD)
        return wD

    def __getAtomSiteInfo(self, dataContainer, modelId="1"):
        """Get counting information for each instance in the deposited coordinates for the input model.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance
            modelId (str, optional): model index. Defaults to "1".

        Returns:
            (dict): with atom site counting and instance mapping details.

            For instance, the following are calculated:

                instanceHeavyAtomCountD[asymId]:  number of deposited heavy atoms
                typeHeavyAtomCountD[entity type]: number of deposited heavy atoms

                numHeavyAtomsModel:  number of deposited heavy atoms in input model_id
                modelId: modelId

                instancePolymerModeledMonomerCountD[asymId]: number modeled polymer monomers in deposited coordinates
                instancePolymerUnmodeledMonomerCountD[asymId]: number of polymer unmodeled monomers in deposited coordinates

                numModels: total number of deposited models
                numAtomsAll: total number of deposited atoms

                asymAuthIdD = {asymId: authAsymId, ... }

                asymIdPolymerRangesD = {asymId: {"sampleSeqLen": sampleSeqLen,
                                                 "obsSeqLen": obsSeqLen,
                                                 "begSeqId": begSeqId,
                                                 "endSeqId": endSeqId,
                                                 "begAuthSeqId": begAuthSeqId,
                                                 "endAuthSeqId": endAuthSeqId,
                                                 "begInsCode": begAuthInsCode,
                                                 "endInsCode": endAuthInsCode,}, ...}
                instanceIdMapD = {asymId:  {"entry_id": entryId,
                                            "entity_id": entityId,
                                            "entity_type": entityTypeD[entityId],
                                            "asym_id": asymId,
                                            "auth_asym_id": authAsymId,
                                            "comp_id": monId,
                                            "auth_seq_id": "?",}, ...}

                 pAuthAsymIdMapD[(authAsymId, authSeqId, insCode)] = {
                        "entry_id": entryId,
                        "entity_id": entityId,
                        "entity_type": entityTypeD[entityId],
                        "asym_id": asymId,
                        "comp_id": compId,
                        "seq_id": seqId,
                    }

                npAuthAsymIdMapD[(authAsymId, resNum)] = {
                        "entry_id": entryId,
                        "entity_id": entityId,
                        "entity_type": entityTypeD[entityId],
                        "asym_id": asymId,
                        "auth_asym_id": authAsymId,
                        "comp_id": monId,
                        "auth_seq_id": resNum,
                    }

                brAuthAsymIdMapD[(authAsymId, authSeqNum)] = {
                        "entry_id": entryId,
                        "entity_id": entityId,
                        "entity_type": entityTypeD[entityId],
                        "asym_id": asymId,
                        "auth_asym_id": authAsymId,
                        "comp_id": monId,
                        "auth_seq_id": authSeqNum,
                        "seq_num": seqNum,
                    }
                entityTypeUniqueIds[<entity_type>][<entity_id>] = {'asymIds': [...],'authAsymIds': [...], 'ccIds': [...]}

                seqIdMapAsymD[asymId] = [<authSeqId + insCode>, ... ]

        """
        #
        numAtomsAll = 0
        numHeavyAtomsModel = 0
        numDeuWatMolModel = 0
        typeHeavyAtomCountD = {}
        instanceHeavyAtomCountD = {}
        instanceDeuWatMolCountD = {}
        #
        numHydrogenAtomsModel = 0
        typeHydrogenAtomCountD = {}
        instanceHydrogenAtomCountD = {}
        #
        instancePolymerModeledMonomerCountD = {}
        instancePolymerUnmodeledMonomerCountD = {}
        atomSiteInfoD = {}
        modelIdL = []
        asymAuthIdD = {}
        occupancySumD = {}
        instanceTypeD = self.getInstanceTypes(dataContainer)
        entityTypeD = self.getEntityTypes(dataContainer)
        #
        eObj = dataContainer.getObj("entity")
        entityIdL = eObj.getAttributeValueList("id")
        #
        try:
            if dataContainer.exists("ihm_model_list"):
                modelId = self.getIhmRepresentativeModelId(dataContainer)
            if dataContainer.exists("atom_site"):
                tObj = dataContainer.getObj("atom_site")
                # All atoms all types deposited -
                numAtomsAll = tObj.getRowCount()
                # Heavy atoms per model -
                cndL = [("type_symbol", "not in", ["H", "D", "T"]), ("pdbx_PDB_model_num", "eq", modelId)]
                numHeavyAtomsModel = tObj.countValuesWhereOpConditions(cndL)
                #
                # Deuterated water molecules per model -
                cndL2 = [("label_atom_id", "eq", "O"), ("label_comp_id", "eq", "DOD"), ("pdbx_PDB_model_num", "eq", modelId)]
                numDeuWatMolModel = tObj.countValuesWhereOpConditions(cndL2)
                #
                modelIdL = tObj.getAttributeUniqueValueList("pdbx_PDB_model_num")
                cD = tObj.getCombinationCountsWithConditions(["label_asym_id", "pdbx_PDB_model_num"], [("type_symbol", "not in", ["H", "D", "T"])])
                dwD = tObj.getCombinationCountsWithConditions(["label_asym_id", "pdbx_PDB_model_num"], [("label_atom_id", "eq", "O"), ("label_comp_id", "eq", "DOD")])
                #
                for asymId, _ in instanceTypeD.items():
                    instanceHeavyAtomCountD[asymId] = cD[(asymId, modelId)] if (asymId, modelId) in cD else 0
                    instanceDeuWatMolCountD[asymId] = dwD[(asymId, modelId)] if (asymId, modelId) in dwD else 0
                #
                # for eType in ['polymer', 'non-polymer', 'branched', 'macrolide', 'solvent']:
                typeHeavyAtomCountD = {k: 0 for k in ["polymer", "non-polymer", "branched", "macrolide", "water"]}
                for asymId, aCount in instanceHeavyAtomCountD.items():
                    tt = instanceTypeD[asymId]
                    typeHeavyAtomCountD[tt] += aCount

                # Hydrogen counts ...
                cndL = [("type_symbol", "in", ["H", "D", "T"]), ("pdbx_PDB_model_num", "eq", modelId)]
                numHydrogenAtomsModel = tObj.countValuesWhereOpConditions(cndL)
                #
                cD = tObj.getCombinationCountsWithConditions(["label_asym_id", "pdbx_PDB_model_num"], [("type_symbol", "in", ["H", "D", "T"])])
                for asymId, _ in instanceTypeD.items():
                    instanceHydrogenAtomCountD[asymId] = cD[(asymId, modelId)] if (asymId, modelId) in cD else 0
                #
                typeHydrogenAtomCountD = {k: 0 for k in ["polymer", "non-polymer", "branched", "macrolide", "water"]}
                for asymId, aCount in instanceHydrogenAtomCountD.items():
                    tt = instanceTypeD[asymId]
                    typeHydrogenAtomCountD[tt] += aCount
                # --
                for ii in range(tObj.getRowCount()):
                    tModelId = tObj.getValue("pdbx_PDB_model_num", ii)
                    if tModelId != modelId:
                        continue
                    asymId = tObj.getValue("label_asym_id", ii)
                    atomType = tObj.getValue("type_symbol", ii)
                    altId = tObj.getValueOrDefault("label_alt_id", ii, None)
                    occupancy = tObj.getValueOrDefault("occupancy", ii, "0.0")
                    if atomType not in ["H", "D", "T"]:
                        if not altId:
                            occupancySumD.setdefault(asymId, defaultdict(float))["FL"] += float(occupancy)
                        else:
                            occupancySumD.setdefault(asymId, defaultdict(float))[altId] += float(occupancy)
                # --
            elif dataContainer.exists("ihm_sphere_obj_site"):
                sObj = dataContainer.getObj("ihm_sphere_obj_site")
                numAtomsAll = numHeavyAtomsModel = numDeuWatMolModel = numHydrogenAtomsModel = 0
                modelIdL = sObj.getAttributeUniqueValueList("model_id")
                for asymId, _ in instanceTypeD.items():
                    instanceHeavyAtomCountD[asymId] = 0
                    instanceDeuWatMolCountD[asymId] = 0
                    instanceHydrogenAtomCountD[asymId] = 0
                    occupancySumD.setdefault(asymId, defaultdict(float))["FL"] = 0.0
                typeHeavyAtomCountD = {k: 0 for k in ["polymer", "non-polymer", "branched", "macrolide", "water"]}
                typeHydrogenAtomCountD = {k: 0 for k in ["polymer", "non-polymer", "branched", "macrolide", "water"]}
            else:
                logger.warning("Missing atom_site category for %s", dataContainer.getName())
            #
            numModels = len(modelIdL)
            if numModels < 1:
                logger.warning("Missing model details in atom_site category for %s", dataContainer.getName())
            #
            atomSiteInfoD = {
                "instanceHeavyAtomCountD": instanceHeavyAtomCountD,
                "typeHeavyAtomCountD": typeHeavyAtomCountD,
                "instanceDeuWatMolCountD": instanceDeuWatMolCountD,
                "numAtomsAll": numAtomsAll,
                "numHeavyAtomsModel": numHeavyAtomsModel,
                "numDeuWatMolModel": numDeuWatMolModel,
                "numModels": len(modelIdL),
                "modelId": modelId,
                "modelIdList": sorted(modelIdL),
                "instancePolymerModeledMonomerCountD": {},
                "instancePolymerUnmodeledMonomerCountD": {},
                "instanceHydrogenAtomCountD": instanceHydrogenAtomCountD,
                "typeHydrogenAtomCountD": typeHydrogenAtomCountD,
                "numHydrogenAtomsModel": numHydrogenAtomsModel,
                "occupancySumD": occupancySumD,
            }
        except Exception as e:
            logger.exception("Failing with %r with %r", dataContainer.getName(), str(e))

        #
        entityTypeUniqueIds = {}
        tAsymIdD = {}
        seqIdObsMapD = {}
        seqIdMapAsymD = {}
        epLengthD = self.getPolymerEntityLengths(dataContainer)
        asymIdPolymerRangesD = {}
        instanceIdMapD = {}
        npAuthAsymIdMapD = {}
        pAuthAsymIdMapD = {}
        brAuthAsymIdMapD = {}
        try:
            eObj = dataContainer.getObj("entry")
            entryId = eObj.getValue("id", 0)
            #
            psObj = dataContainer.getObj("pdbx_poly_seq_scheme")
            if psObj is not None:
                # --
                for eId in entityIdL:
                    if entityTypeD[eId] in ["polymer"]:
                        tAsymIdL = psObj.selectValuesWhere("asym_id", eId, "entity_id")
                        tAuthAsymIdL = psObj.selectValuesWhere("pdb_strand_id", eId, "entity_id")
                        tCcIdL = psObj.selectValuesWhere("mon_id", eId, "entity_id")
                        entityTypeUniqueIds.setdefault(entityTypeD[eId], {}).setdefault(eId, {"asymIds": tAsymIdL, "authAsymIds": tAuthAsymIdL, "ccIds": tCcIdL})
                # ---
                aSeqD = {}
                aOrgSeqD = {}
                hasPdbSeqNum = psObj.hasAttribute("pdb_seq_num")
                for ii in range(psObj.getRowCount()):
                    asymId = psObj.getValue("asym_id", ii)
                    #
                    # AF models do not provide - pdb_seq_num
                    if hasPdbSeqNum:
                        authSeqId = psObj.getValue("pdb_seq_num", ii)
                    else:
                        authSeqId = psObj.getValue("auth_seq_num", ii)

                    authOrgSeqId = psObj.getValue("auth_seq_num", ii)
                    seqId = psObj.getValue("seq_id", ii)
                    compId = psObj.getValue("mon_id", ii)
                    entityId = psObj.getValue("entity_id", ii)
                    authAsymId = psObj.getValue("pdb_strand_id", ii)
                    #
                    insCode = psObj.getValueOrDefault("pdb_ins_code", ii, defaultValue=None)
                    aSeqD.setdefault(asymId, []).append(authSeqId)
                    aOrgSeqD.setdefault(asymId, []).append(authOrgSeqId)
                    # ---
                    tC = authSeqId
                    if authSeqId not in [".", "?"]:
                        seqIdObsMapD.setdefault(asymId, {})[seqId] = (authSeqId, insCode)
                    else:
                        tC = "?"
                    if insCode and tC != "?":
                        tC += insCode
                    seqIdMapAsymD.setdefault(asymId, []).append(tC)
                    # ---
                    #
                    pAuthAsymIdMapD[(authAsymId, authSeqId, insCode)] = {
                        "entry_id": entryId,
                        "entity_id": entityId,
                        "entity_type": entityTypeD[entityId],
                        "asym_id": asymId,
                        "comp_id": compId,
                        "seq_id": seqId,
                    }
                    #
                    if asymId in tAsymIdD:
                        continue
                    tAsymIdD[asymId] = entityId
                    asymAuthIdD[asymId] = authAsymId
                    #
                    instanceIdMapD[asymId] = {
                        "entry_id": entryId,
                        "entity_id": entityId,
                        "entity_type": entityTypeD[entityId],
                        "asym_id": asymId,
                        "auth_asym_id": authAsymId,
                        "rcsb_id": entryId + "." + asymId,
                        "comp_id": "?",
                        "auth_seq_id": "?",
                    }
                    #

                #
                #  Get the modeled and unmodeled monomer counts by asymId
                #  JDW not use aOrgSeqD.items()
                for asymId, sL in aOrgSeqD.items():
                    instancePolymerModeledMonomerCountD[asymId] = len([t for t in sL if t not in ["?", "."]])
                    instancePolymerUnmodeledMonomerCountD[asymId] = len([t for t in sL if t in ["?", "."]])
                #  Get polymer range details for each polymer instance
                for asymId, entityId in tAsymIdD.items():
                    sampleSeqLen = epLengthD[entityId] if entityId in epLengthD else None
                    sL = list(seqIdObsMapD[asymId].items())
                    begSeqId, (begAuthSeqId, begAuthInsCode) = sL[0]
                    endSeqId, (endAuthSeqId, endAuthInsCode) = sL[-1]
                    obsSeqLen = len(sL)
                    #
                    asymIdPolymerRangesD[asymId] = {
                        "sampleSeqLen": sampleSeqLen,
                        "obsSeqLen": obsSeqLen,
                        "begSeqId": begSeqId,
                        "endSeqId": endSeqId,
                        "begAuthSeqId": begAuthSeqId,
                        "endAuthSeqId": endAuthSeqId,
                        "begInsCode": begAuthInsCode,
                        "endInsCode": endAuthInsCode,
                    }
            atomSiteInfoD["instancePolymerModeledMonomerCountD"] = instancePolymerModeledMonomerCountD
            atomSiteInfoD["instancePolymerUnmodeledMonomerCountD"] = instancePolymerUnmodeledMonomerCountD
            atomSiteInfoD["asymAuthIdD"] = asymAuthIdD
            atomSiteInfoD["asymIdPolymerRangesD"] = asymIdPolymerRangesD
            atomSiteInfoD["seqIdMapAsymD"] = seqIdMapAsymD
            # --------------
            logger.debug(
                "%s instancePolymerModeledMonomerCountD(%d) %r",
                dataContainer.getName(),
                sum(atomSiteInfoD["instancePolymerModeledMonomerCountD"].values()),
                atomSiteInfoD["instancePolymerModeledMonomerCountD"],
            )
            logger.debug("%s instancePolymerUnmodeledMonomerCountD %r", dataContainer.getName(), atomSiteInfoD["instancePolymerUnmodeledMonomerCountD"])
            #
            # -------------- -------------- -------------- -------------- -------------- -------------- -------------- --------------
            #  Add nonpolymer instance mapping
            #
            npsObj = dataContainer.getObj("pdbx_nonpoly_scheme")
            if npsObj is not None:
                # --
                for eId in entityIdL:
                    if entityTypeD[eId] in ["non-polymer", "water"]:
                        tAsymIdL = npsObj.selectValuesWhere("asym_id", eId, "entity_id")
                        tAuthAsymIdL = npsObj.selectValuesWhere("pdb_strand_id", eId, "entity_id")
                        tCcIdL = npsObj.selectValuesWhere("mon_id", eId, "entity_id")
                        entityTypeUniqueIds.setdefault(entityTypeD[eId], {}).setdefault(eId, {"asymIds": tAsymIdL, "authAsymIds": tAuthAsymIdL, "ccIds": tCcIdL})
                # ---
                for ii in range(npsObj.getRowCount()):
                    asymId = npsObj.getValue("asym_id", ii)
                    entityId = npsObj.getValue("entity_id", ii)
                    authAsymId = npsObj.getValue("pdb_strand_id", ii)
                    resNum = npsObj.getValue("pdb_seq_num", ii)
                    monId = npsObj.getValue("mon_id", ii)
                    asymAuthIdD[asymId] = authAsymId
                    if asymId not in instanceIdMapD:
                        instanceIdMapD[asymId] = {
                            "entry_id": entryId,
                            "entity_id": entityId,
                            "entity_type": entityTypeD[entityId],
                            "asym_id": asymId,
                            "auth_asym_id": authAsymId,
                            "rcsb_id": entryId + "." + asymId,
                            "comp_id": monId,
                            "auth_seq_id": resNum,
                        }
                    npAuthAsymIdMapD[(authAsymId, resNum)] = {
                        "entry_id": entryId,
                        "entity_id": entityId,
                        "entity_type": entityTypeD[entityId],
                        "asym_id": asymId,
                        "auth_asym_id": authAsymId,
                        "comp_id": monId,
                        "auth_seq_id": resNum,
                    }

            # ---------
            brsObj = dataContainer.getObj("pdbx_branch_scheme")
            if brsObj is not None:
                # --
                for eId in entityIdL:
                    if entityTypeD[eId] in ["branched"]:
                        tAsymIdL = brsObj.selectValuesWhere("asym_id", eId, "entity_id")
                        # changed to pdb_asym_id on 2020-07-29
                        tAuthAsymIdL = brsObj.selectValuesWhere("pdb_asym_id", eId, "entity_id")
                        tCcIdL = brsObj.selectValuesWhere("mon_id", eId, "entity_id")
                        entityTypeUniqueIds.setdefault(entityTypeD[eId], {}).setdefault(eId, {"asymIds": tAsymIdL, "authAsymIds": tAuthAsymIdL, "ccIds": tCcIdL})
                # ---
                for ii in range(brsObj.getRowCount()):
                    asymId = brsObj.getValue("asym_id", ii)
                    entityId = brsObj.getValue("entity_id", ii)
                    #
                    authAsymId = brsObj.getValue("pdb_asym_id", ii)
                    authSeqNum = brsObj.getValue("pdb_seq_num", ii)
                    monId = brsObj.getValue("mon_id", ii)
                    seqNum = brsObj.getValue("num", ii)
                    asymAuthIdD[asymId] = authAsymId
                    if asymId not in instanceIdMapD:
                        instanceIdMapD[asymId] = {
                            "entry_id": entryId,
                            "entity_id": entityId,
                            "entity_type": entityTypeD[entityId],
                            "asym_id": asymId,
                            "auth_asym_id": authAsymId,
                            "rcsb_id": entryId + "." + asymId,
                            "comp_id": monId,
                            "auth_seq_id": "?",
                        }
                    brAuthAsymIdMapD[(authAsymId, authSeqNum)] = {
                        "entry_id": entryId,
                        "entity_id": entityId,
                        "entity_type": entityTypeD[entityId],
                        "asym_id": asymId,
                        "auth_asym_id": authAsymId,
                        "comp_id": monId,
                        "auth_seq_id": authSeqNum,
                        "seq_num": seqNum,
                    }

            #
            atomSiteInfoD["instanceIdMapD"] = instanceIdMapD
            atomSiteInfoD["npAuthAsymIdMapD"] = npAuthAsymIdMapD
            atomSiteInfoD["pAuthAsymIdMapD"] = pAuthAsymIdMapD
            atomSiteInfoD["brAuthAsymIdMapD"] = brAuthAsymIdMapD
            atomSiteInfoD["entityTypeUniqueIds"] = entityTypeUniqueIds

        except Exception as e:
            logger.exception("Failing for %s with %s", dataContainer.getName(), str(e))

        #
        return atomSiteInfoD

    # Connection related
    def getInstanceConnectionCounts(self, dataContainer):
        """Return a dictionary instance connection counts.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {<connection type>: #count, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchInstanceConnections(dataContainer)
        return wD["instConnectCountD"] if "instConnectCountD" in wD else {}

    def getInstanceConnections(self, dataContainer):
        """Return a list of instance connections.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            list: [{"connect_type": <val>, "connect_target_label_comp_id": <val>, ... },...]

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchInstanceConnections(dataContainer)
        return wD["instConnectL"] if "instConnectL" in wD else {}

    def getBoundNonpolymersComponentIds(self, dataContainer):
        """Return a list of bound non-polymers in the entry.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {<entityId>: NonpolymerBoundEntity("targetCompId", "connectType", "partnerCompId", "partnerEntityId", "partnerEntityType"), }
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchInstanceConnections(dataContainer)
        return wD["boundNonpolymerComponentIdL"] if "boundNonpolymerComponentIdL" in wD else {}

    def getBoundNonpolymersByEntity(self, dataContainer):
        """Return a dictionary of bound non-polymers by entity.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {<entityId>: NonpolymerBoundEntity("targetCompId", "connectType", "partnerCompId", "partnerEntityId", "partnerEntityType"), }
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchInstanceConnections(dataContainer)
        return wD["boundNonpolymerEntityD"] if "boundNonpolymerEntityD" in wD else {}

    def getBoundNonpolymersByInstance(self, dataContainer):
        """Return a dictionary of bound non-polymers by instance.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {<asymId>: NonpolymerBoundInstance( "targetCompId", "targetAtomId", "targetAltId", "connectType", "partnerEntityType", "partnerEntityId",
                                                      "partnerCompId","partnerAsymId", "partnerSeqId", "partnerAuthSeqId", "partnerAtomId", "targetAltId",
                                                      "bondDistance", "bondOrder"), }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchInstanceConnections(dataContainer)
        return wD["boundNonpolymerInstanceD"] if "boundNonpolymerInstanceD" in wD else {}

    def __fetchInstanceConnections(self, dataContainer):
        wD = self.__instanceConnectionCache.get(dataContainer.getName())
        if not wD:
            wD = self.__getInstanceConnections(dataContainer)
            self.__instanceConnectionCache.set(dataContainer.getName(), wD)
        return wD

    def __getInstanceConnections(self, dataContainer):
        """Get instance connections (e.g., intermolecular bonds and non-primary connectivity)

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: instConnectCountD{<bond_type>: count, ... }

            For instance, the following are calculated:
                     {Get counting information about intermolecular linkages.
            covale  .
            disulf  .
            hydrog  .
            metalc

            loop_
            _struct_asym.id
            _struct_asym.pdbx_blank_PDB_chainid_flag
            _struct_asym.pdbx_modified
            _struct_asym.entity_id
            _struct_asym.details
            A N N 1 ?
            B N N 1 ?
            #
            _struct_biol.id   1
            #
            loop_
            _struct_conn.id
            _struct_conn.conn_type_id
            _struct_conn.pdbx_leaving_atom_flag
            _struct_conn.pdbx_PDB_id
            _struct_conn.ptnr1_label_asym_id
            _struct_conn.ptnr1_label_comp_id
            _struct_conn.ptnr1_label_seq_id
            _struct_conn.ptnr1_label_atom_id
            _struct_conn.pdbx_ptnr1_label_alt_id
            _struct_conn.pdbx_ptnr1_PDB_ins_code
            _struct_conn.pdbx_ptnr1_standard_comp_id
            _struct_conn.ptnr1_symmetry
            _struct_conn.ptnr2_label_asym_id
            _struct_conn.ptnr2_label_comp_id
            _struct_conn.ptnr2_label_seq_id
            _struct_conn.ptnr2_label_atom_id
            _struct_conn.pdbx_ptnr2_label_alt_id
            _struct_conn.pdbx_ptnr2_PDB_ins_code
            _struct_conn.ptnr1_auth_asym_id
            _struct_conn.ptnr1_auth_comp_id
            _struct_conn.ptnr1_auth_seq_id
            _struct_conn.ptnr2_auth_asym_id
            _struct_conn.ptnr2_auth_comp_id
            _struct_conn.ptnr2_auth_seq_id
            _struct_conn.ptnr2_symmetry
            _struct_conn.pdbx_ptnr3_label_atom_id
            _struct_conn.pdbx_ptnr3_label_seq_id
            _struct_conn.pdbx_ptnr3_label_comp_id
            _struct_conn.pdbx_ptnr3_label_asym_id
            _struct_conn.pdbx_ptnr3_label_alt_id
            _struct_conn.pdbx_ptnr3_PDB_ins_code
            _struct_conn.details
            _struct_conn.pdbx_dist_value
            _struct_conn.pdbx_value_order
            disulf1  disulf ? ? A CYS 31 SG ? ? ? 1_555 B CYS 31 SG ? ? A CYS 31 B CYS 31 1_555 ? ? ? ? ? ? ? 1.997 ?
            covale1  covale ? ? A VAL 8  C  ? ? ? 1_555 A DPR 9  N  ? ? A VAL 8  A DPR 9  1_555 ? ? ? ? ? ? ? 1.360 ?
            covale2  covale ? ? A DPR 9  C  ? ? ? 1_555 A GLY 10 N  ? ? A DPR 9  A GLY 10 1_555 ? ? ? ? ? ? ? 1.324 ?
            #
        """
        iAttMapD = {
            "id": "id",
            "connect_type": "conn_type_id",
            "connect_target_label_comp_id": "ptnr1_label_comp_id",
            "connect_target_label_asym_id": "ptnr1_label_asym_id",
            "connect_target_label_seq_id": "ptnr1_label_seq_id",
            "connect_target_auth_seq_id": "ptnr1_auth_seq_id",
            "connect_target_label_atom_id": "ptnr1_label_atom_id",
            "connect_target_label_alt_id": "pdbx_ptnr1_label_alt_id",
            "connect_target_symmetry": "ptnr1_symmetry",
            #
            "connect_partner_label_comp_id": "ptnr2_label_comp_id",
            "connect_partner_label_asym_id": "ptnr2_label_asym_id",
            "connect_partner_label_seq_id": "ptnr2_label_seq_id",
            "connect_partner_auth_seq_id": "ptnr2_auth_seq_id",
            "connect_partner_label_atom_id": "ptnr2_label_atom_id",
            "connect_partner_label_alt_id": "pdbx_ptnr2_label_alt_id",
            "connect_partner_symmetry": "ptnr2_symmetry",
            "value_order": "pdbx_value_order",
            "dist_value": "pdbx_dist_value",
            "description": "details",
            "role": "pdbx_role",
        }
        jAttMapD = {
            "id": "id",
            "connect_type": "conn_type_id",
            "connect_target_label_comp_id": "ptnr2_label_comp_id",
            "connect_target_label_asym_id": "ptnr2_label_asym_id",
            "connect_target_label_seq_id": "ptnr2_label_seq_id",
            "connect_target_auth_seq_id": "ptnr2_auth_seq_id",
            "connect_target_label_atom_id": "ptnr2_label_atom_id",
            "connect_target_label_alt_id": "pdbx_ptnr2_label_alt_id",
            "connect_target_symmetry": "ptnr2_symmetry",
            #
            "connect_partner_label_comp_id": "ptnr1_label_comp_id",
            "connect_partner_label_asym_id": "ptnr1_label_asym_id",
            "connect_partner_label_seq_id": "ptnr1_label_seq_id",
            "connect_partner_auth_seq_id": "ptnr1_auth_seq_id",
            "connect_partner_label_atom_id": "ptnr1_label_atom_id",
            "connect_partner_label_alt_id": "pdbx_ptnr1_label_alt_id",
            "connect_partner_symmetry": "ptnr1_symmetry",
            "value_order": "pdbx_value_order",
            "dist_value": "pdbx_dist_value",
            "description": "details",
            "role": "pdbx_role",
        }
        typeMapD = {
            "covale": "covalent bond",
            "disulf": "disulfide bridge",
            "hydrog": "hydrogen bond",
            "metalc": "metal coordination",
            "mismat": "mismatched base pairs",
            "saltbr": "ionic interaction",
            "modres": "covalent residue modification",
            "covale_base": "covalent modification of a nucleotide base",
            "covale_sugar": "covalent modification of a nucleotide sugar",
            "covale_phosphate": "covalent modification of a nucleotide phosphate",
        }
        #
        instConnectL = []
        instConnectCountD = {ky: 0 for ky in typeMapD}
        boundNonpolymerEntityD = {}
        boundNonpolymerInstanceD = {}
        boundNonpolymerComponentIdL = []
        #
        if dataContainer.exists("struct_conn"):
            tObj = dataContainer.getObj("struct_conn")
            for ii in range(tObj.getRowCount()):
                bt = str(tObj.getValue("conn_type_id", ii)).strip().lower()
                if bt not in instConnectCountD:
                    logger.error("Unsupported intermolecular bond type %r in %r", bt, dataContainer.getName())
                    continue
                instConnectCountD[bt] = instConnectCountD[bt] + 1 if bt in instConnectCountD else instConnectCountD[bt]
                #
                tD = OrderedDict()
                for ky, atName in iAttMapD.items():
                    if tObj.hasAttribute(atName):
                        val = tObj.getValue(atName, ii) if atName != "conn_type_id" else typeMapD[tObj.getValue(atName, ii).lower()]
                        tD[ky] = val
                instConnectL.append(tD)
                # Flip the bond sense so all target connections are accounted for
                tD = OrderedDict()
                for ky, atName in jAttMapD.items():
                    if tObj.hasAttribute(atName):
                        val = tObj.getValue(atName, ii) if atName != "conn_type_id" else typeMapD[tObj.getValue(atName, ii).lower()]
                        tD[ky] = val
                instConnectL.append(tD)

            boundNonpolymerEntityD, boundNonpolymerInstanceD, boundNonpolymerComponentIdL = self.__getBoundNonpolymers(dataContainer, instConnectL)

        return {
            "instConnectL": instConnectL,
            "instConnectCountD": instConnectCountD,
            "boundNonpolymerEntityD": boundNonpolymerEntityD,
            "boundNonpolymerInstanceD": boundNonpolymerInstanceD,
            "boundNonpolymerComponentIdL": boundNonpolymerComponentIdL,
        }

    def __getBoundNonpolymers(self, dataContainer, instConnectL):
        """Get nonpolymer bound

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            bool: True for success or False otherwise

        Example:
        """
        logger.debug("Starting with %r", dataContainer.getName())
        #
        boundNonpolymerEntityD = {}
        boundNonpolymerInstanceD = {}
        boundNonpolymerComponentIdL = []
        try:
            asymIdD = self.getInstanceEntityMap(dataContainer)
            # asymAuthIdD = self.getAsymAuthIdMap(dataContainer)
            eTypeD = self.getEntityTypes(dataContainer)
            #
            ts = set()
            for cD in instConnectL:
                tAsymId = cD["connect_target_label_asym_id"]
                tEntityId = asymIdD[tAsymId]
                if eTypeD[tEntityId] == "non-polymer" and cD["connect_type"] in ["covale", "covalent bond", "metalc", "metal coordination"]:
                    pAsymId = cD["connect_partner_label_asym_id"]
                    pEntityId = asymIdD[pAsymId]
                    pCompId = cD["connect_partner_label_comp_id"]
                    pSeqId = cD["connect_partner_label_seq_id"]
                    pAuthSeqId = cD["connect_partner_auth_seq_id"]
                    tCompId = cD["connect_target_label_comp_id"]
                    tAtomId = cD["connect_target_label_atom_id"]
                    pAtomId = cD["connect_partner_label_atom_id"]
                    tAltId = cD["connect_target_label_alt_id"]
                    pAltId = cD["connect_partner_label_alt_id"]
                    bondOrder = cD["value_order"]
                    bondDist = cD["dist_value"]
                    role = cD["role"] if "role" in cD else None
                    eType = eTypeD[pEntityId]
                    #
                    ts.add(tCompId)
                    boundNonpolymerInstanceD.setdefault(tAsymId, []).append(
                        NonpolymerBoundInstance(
                            tCompId,
                            tAtomId,
                            tAltId,
                            cD["connect_type"],
                            eType,
                            pEntityId,
                            pCompId,
                            pAsymId,
                            pSeqId,
                            pAuthSeqId,
                            pAtomId,
                            pAltId,
                            bondDist,
                            bondOrder,
                            role,
                        )
                    )
                    boundNonpolymerEntityD.setdefault(tEntityId, []).append(NonpolymerBoundEntity(tCompId, cD["connect_type"], pCompId, pEntityId, eType))
            #
            cloneD = copy.deepcopy(boundNonpolymerInstanceD)
            for asymId in cloneD:
                boundNonpolymerInstanceD[asymId] = sorted(set(cloneD[asymId]))
            cloneD = copy.deepcopy(boundNonpolymerEntityD)
            for entityId in cloneD:
                boundNonpolymerEntityD[entityId] = sorted(set(cloneD[entityId]))
            boundNonpolymerComponentIdL = sorted(ts)
        except Exception as e:
            logger.exception("%s failing with %s", dataContainer.getName(), str(e))
        return boundNonpolymerEntityD, boundNonpolymerInstanceD, boundNonpolymerComponentIdL

    def getEntitySequenceFeatureCounts(self, dataContainer):
        """Return a dictionary of sequence feature counts.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {<entity>: {'mutation': #, 'artifact': #, 'conflict': #, ...  }, }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchSequenceFeatures(dataContainer)
        return wD["seqFeatureCountsD"] if "seqFeatureCountsD" in wD else {}

    def getEntitySequenceMonomerFeatures(self, dataContainer):
        """Return a dictionary of sequence monomer features.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(entityId,seqId,compId,filteredFeature): {detail,detail},  .. }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchSequenceFeatures(dataContainer)
        return wD["seqMonomerFeatureD"] if "seqMonomerFeatureD" in wD else {}

    def getEntitySequenceRangeFeatures(self, dataContainer):
        """Return a dictionary of sequence range features.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(entityId,benSeqId,endSeqId,filteredFeature): {detail,detail},  .. }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchSequenceFeatures(dataContainer)
        return wD["seqRangeFeatureD"] if "seqRangeFeatureD" in wD else {}

    def getEntityReferenceAlignments(self, dataContainer):
        """Return a dictionary of reference sequence alignments for each entity.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {entityId: {'dbName': , 'dbAccession': , 'authAsymId': , 'entitySeqIdBeg':, 'dbSeqIdBeg':, ... },  .. }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchReferenceSequenceDetails(dataContainer)
        return wD["seqEntityAlignmentD"] if "seqEntityAlignmentD" in wD else {}

    def getEntityPolymerSequences(self, dataContainer):
        """Return a dictionary of the sequences (one-letter-codes) for each polymer entity.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {entityId: {'sequence': ..., 'polymerType': ... , 'polymerTypeFiltered': ... },  ... }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchReferenceSequenceDetails(dataContainer)
        return wD["entityPolymerSequenceD"] if "entityPolymerSequenceD" in wD else {}

    def getEntitySequenceReferenceCodes(self, dataContainer):
        """Return a dictionary of reference database accession codes.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {entityId: {'dbName': , 'dbAccession': },  ... }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchReferenceSequenceDetails(dataContainer)
        return wD["seqEntityRefDbD"] if "seqEntityRefDbD" in wD else {}

    def __fetchSequenceFeatures(self, dataContainer):
        wD = self.__entitySequenceFeatureCache.get(dataContainer.getName())
        if not wD:
            wD = self.__getSequenceFeatures(dataContainer)
            self.__entitySequenceFeatureCache.set(dataContainer.getName(), wD)
        return wD

    def __fetchReferenceSequenceDetails(self, dataContainer):
        wD = self.__entityReferenceSequenceDetailsCache.get(dataContainer.getName())
        if not wD:
            wD = self.__getReferenceSequenceDetails(dataContainer)
            self.__entityReferenceSequenceDetailsCache.set(dataContainer.getName(), wD)
        return wD

    def getDatabaseNameMap(self):
        dbNameMapD = {
            "UNP": "UniProt",
            "GB": "GenBank",
            "PDB": "PDB",
            # Can be deleted once all IHM entries are remediated
            # "PDB-DEV": "PDB-Dev",
            # "PDB-IHM": "PDB-IHM",
            "EMBL": "EMBL",
            "GENP": "GenBank",
            "NDB": "NDB",
            "NOR": "NORINE",
            "PIR": "PIR",
            "PRF": "PRF",
            "REF": "RefSeq",
            "TPG": "GenBank",
            "TREMBL": "UniProt",
            "SWS": "UniProt",
            "SWALL": "UniProt",
        }
        return dbNameMapD

    def getPolymerEntityReferenceAlignments(self, dataContainer, entityId=None, dbName=None):
        """Get list of polymer entity reference alignments from category 'rcsb_polymer_entity_align'

        Args:
            dataContainer (object): mmcif.api.DataContainer object instance
            entityId (optional): Only return reference alignments for a specfic entity.
            dbName (optional): _description_. Defaults to None.

        Returns:
            list: list of entity reference alignment dictionaries
                  e.g., [{
                    'ordinal': 1, 'entry_id': '7XIW', 'entity_id': '1', 'reference_database_name': 'UniProt',
                    'reference_database_accession': 'P0DTC2', 'reference_database_isoform': None, 'provenance_source': 'SIFTS',
                    'aligned_regions_entity_beg_seq_id': '1', 'aligned_regions_ref_beg_seq_id': '1', 'aligned_regions_length': '1270'
                  }]
        """
        pdbEntityAlignL = []

        if dataContainer.exists("rcsb_polymer_entity_align"):
            aObj = dataContainer.getObj("rcsb_polymer_entity_align")
            for idx in range(aObj.getRowCount()):
                aD = aObj.getRowAttributeDict(idx)
                # Example aD:
                # {
                #     "oridnal" : 1,
                #     "entry_id" : "1B5F",
                #     "entity_id" : 1,
                #     "reference_database_name" : "UniProt",
                #     "reference_database_accession" : "Q9XFX3",
                #     "provenance_source" : "SIFTS",
                #     "aligned_regions" : [
                #         {
                #             "ref_beg_seq_id" : 418,
                #             "entity_beg_seq_id" : 1,
                #             "length" : 87
                #         }
                #     ]
                # }
                if entityId and str(aD.get("entity_id", "")) != str(entityId):
                    continue
                if dbName and aD.get("reference_database_name") != dbName:
                    continue
                pdbEntityAlignL.append(aD)
        else:
            logger.warning("Missing rcsb_polymer_entity_align information for dataContainer %r", dataContainer.getName())

        return pdbEntityAlignL

    def __getReferenceSequenceDetails(self, dataContainer):
        """Get reference sequence and related alignment details.

        Args:
            dataContainer (object): mmcif.api.DataContainer object instance

        Returns:
            dict : {
                    "seqEntityAlignmentD" : {entityId: [{'dbName': 'UNP' , 'dbAccession': 'P000000', ... }]}
                    "seqEntityRefDbD":  {entityId: [{'dbName': 'UNP' , 'dbAccession': 'P000000'),  }]},
                    }

        Example source content:

            _struct_ref.id                         1
            _struct_ref.db_name                    UNP
            _struct_ref.db_code                    KSYK_HUMAN
            _struct_ref.pdbx_db_accession          P43405
            _struct_ref.entity_id                  1
            _struct_ref.pdbx_seq_one_letter_code
            ;ADPEEIRPKEVYLDRKLLTLEDKELGSGNFGTVKKGYYQMKKVVKTVAVKILKNEANDPALKDELLAEANVMQQLDNPYI
            VRMIGICEAESWMLVMEMAELGPLNKYLQQNRHVKDKNIIELVHQVSMGMKYLEESNFVHRDLAARNVLLVTQHYAKISD
            FGLSKALRADENYYKAQTHGKWPVKWYAPECINYYKFSSKSDVWSFGVLMWEAFSYGQKPYRGMKGSEVTAMLEKGERMG
            CPAGCPREMYDLMNLCWTYDVENRPGFAAVELRLRNYYYDVVN
            ;
            _struct_ref.pdbx_align_begin           353
            _struct_ref.pdbx_db_isoform            ?
            #
            _struct_ref_seq.align_id                      1
            _struct_ref_seq.ref_id                        1
            _struct_ref_seq.pdbx_PDB_id_code              1XBB
            _struct_ref_seq.pdbx_strand_id                A
            _struct_ref_seq.seq_align_beg                 1
            _struct_ref_seq.pdbx_seq_align_beg_ins_code   ?
            _struct_ref_seq.seq_align_end                 283
            _struct_ref_seq.pdbx_seq_align_end_ins_code   ?
            _struct_ref_seq.pdbx_db_accession             P43405
            _struct_ref_seq.db_align_beg                  353
            _struct_ref_seq.pdbx_db_align_beg_ins_code    ?
            _struct_ref_seq.db_align_end                  635
            _struct_ref_seq.pdbx_db_align_end_ins_code    ?
            _struct_ref_seq.pdbx_auth_seq_align_beg       353
            _struct_ref_seq.pdbx_auth_seq_align_end       635
            _struct_ref_seq.rcsb_entity_id                1
            #
            loop_
            _struct_ref_seq_dif.align_id
            _struct_ref_seq_dif.pdbx_pdb_id_code
            _struct_ref_seq_dif.mon_id
            _struct_ref_seq_dif.pdbx_pdb_strand_id
            _struct_ref_seq_dif.seq_num
            _struct_ref_seq_dif.pdbx_pdb_ins_code
            _struct_ref_seq_dif.pdbx_seq_db_name
            _struct_ref_seq_dif.pdbx_seq_db_accession_code
            _struct_ref_seq_dif.db_mon_id
            _struct_ref_seq_dif.pdbx_seq_db_seq_num
            _struct_ref_seq_dif.details
            _struct_ref_seq_dif.pdbx_auth_seq_num
            _struct_ref_seq_dif.pdbx_ordinal
            _struct_ref_seq_dif.rcsb_entity_id
            1 1XBB MET A 1   ? UNP P43405 ALA 353 'CLONING ARTIFACT' 353 1  1
            1 1XBB ALA A 2   ? UNP P43405 ASP 354 'CLONING ARTIFACT' 354 2  1
            1 1XBB LEU A 3   ? UNP P43405 PRO 355 'CLONING ARTIFACT' 355 3  1
            1 1XBB GLU A 284 ? UNP P43405 ?   ?   'CLONING ARTIFACT' 636 4  1
            1 1XBB GLY A 285 ? UNP P43405 ?   ?   'CLONING ARTIFACT' 637 5  1
            1 1XBB HIS A 286 ? UNP P43405 ?   ?   'EXPRESSION TAG'   638 6  1
            1 1XBB HIS A 287 ? UNP P43405 ?   ?   'EXPRESSION TAG'   639 7  1
            1 1XBB HIS A 288 ? UNP P43405 ?   ?   'EXPRESSION TAG'   640 8  1
            1 1XBB HIS A 289 ? UNP P43405 ?   ?   'EXPRESSION TAG'   641 9  1
            1 1XBB HIS A 290 ? UNP P43405 ?   ?   'EXPRESSION TAG'   642 10 1
            1 1XBB HIS A 291 ? UNP P43405 ?   ?   'EXPRESSION TAG'   643 11 1
            #
            #
            loop_
            _struct_ref_seq_dif.align_id
            _struct_ref_seq_dif.pdbx_pdb_id_code
            _struct_ref_seq_dif.mon_id
            _struct_ref_seq_dif.pdbx_pdb_strand_id
            _struct_ref_seq_dif.seq_num
            _struct_ref_seq_dif.pdbx_pdb_ins_code
            _struct_ref_seq_dif.pdbx_seq_db_name
            _struct_ref_seq_dif.pdbx_seq_db_accession_code
            _struct_ref_seq_dif.db_mon_id
            _struct_ref_seq_dif.pdbx_seq_db_seq_num
            _struct_ref_seq_dif.details
            _struct_ref_seq_dif.pdbx_auth_seq_num
            _struct_ref_seq_dif.pdbx_ordinal
            _struct_ref_seq_dif.rcsb_entity_id
            1 3RIJ TYR A 53  ? UNP Q5SHN1 PHE 54  'ENGINEERED MUTATION' 54  1  1
            1 3RIJ GLY A 54  ? UNP Q5SHN1 VAL 55  'ENGINEERED MUTATION' 55  2  1
            2 3RIJ ASP A 98  ? UNP Q5SHN1 ALA 99  'ENGINEERED MUTATION' 99  3  1
            2 3RIJ ALA A 99  ? UNP Q5SHN1 ILE 100 'ENGINEERED MUTATION' 100 4  1
            2 3RIJ LEU A 158 ? UNP Q5SHN1 ?   ?   INSERTION             159 5  1
            2 3RIJ GLU A 159 ? UNP Q5SHN1 ?   ?   INSERTION             160 6  1
            2 3RIJ HIS A 160 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      161 7  1
            2 3RIJ HIS A 161 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      162 8  1
            2 3RIJ HIS A 162 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      163 9  1
            2 3RIJ HIS A 163 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      164 10 1
            2 3RIJ HIS A 164 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      165 11 1
            2 3RIJ HIS A 165 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      166 12 1
            3 3RIJ TYR B 53  ? UNP Q5SHN1 PHE 54  'ENGINEERED MUTATION' 54  13 1
            3 3RIJ GLY B 54  ? UNP Q5SHN1 VAL 55  'ENGINEERED MUTATION' 55  14 1
            4 3RIJ ASP B 98  ? UNP Q5SHN1 ALA 99  'ENGINEERED MUTATION' 99  15 1
            4 3RIJ ALA B 99  ? UNP Q5SHN1 ILE 100 'ENGINEERED MUTATION' 100 16 1
            4 3RIJ LEU B 158 ? UNP Q5SHN1 ?   ?   INSERTION             159 17 1
            4 3RIJ GLU B 159 ? UNP Q5SHN1 ?   ?   INSERTION             160 18 1
            4 3RIJ HIS B 160 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      161 19 1
            4 3RIJ HIS B 161 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      162 20 1
            4 3RIJ HIS B 162 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      163 21 1
            4 3RIJ HIS B 163 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      164 22 1
            4 3RIJ HIS B 164 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      165 23 1
            4 3RIJ HIS B 165 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      166 24 1
            5 3RIJ TYR C 53  ? UNP Q5SHN1 PHE 54  'ENGINEERED MUTATION' 54  25 1
            5 3RIJ GLY C 54  ? UNP Q5SHN1 VAL 55  'ENGINEERED MUTATION' 55  26 1
            6 3RIJ ASP C 98  ? UNP Q5SHN1 ALA 99  'ENGINEERED MUTATION' 99  27 1
            6 3RIJ ALA C 99  ? UNP Q5SHN1 ILE 100 'ENGINEERED MUTATION' 100 28 1
            6 3RIJ LEU C 158 ? UNP Q5SHN1 ?   ?   INSERTION             159 29 1
            6 3RIJ GLU C 159 ? UNP Q5SHN1 ?   ?   INSERTION             160 30 1
            6 3RIJ HIS C 160 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      161 31 1
            6 3RIJ HIS C 161 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      162 32 1
            6 3RIJ HIS C 162 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      163 33 1
            6 3RIJ HIS C 163 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      164 34 1
            6 3RIJ HIS C 164 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      165 35 1
            6 3RIJ HIS C 165 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      166 36 1
            7 3RIJ TYR D 53  ? UNP Q5SHN1 PHE 54  'ENGINEERED MUTATION' 54  37 1
            7 3RIJ GLY D 54  ? UNP Q5SHN1 VAL 55  'ENGINEERED MUTATION' 55  38 1
            8 3RIJ ASP D 98  ? UNP Q5SHN1 ALA 99  'ENGINEERED MUTATION' 99  39 1
            8 3RIJ ALA D 99  ? UNP Q5SHN1 ILE 100 'ENGINEERED MUTATION' 100 40 1
            8 3RIJ LEU D 158 ? UNP Q5SHN1 ?   ?   INSERTION             159 41 1
            8 3RIJ GLU D 159 ? UNP Q5SHN1 ?   ?   INSERTION             160 42 1
            8 3RIJ HIS D 160 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      161 43 1
            8 3RIJ HIS D 161 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      162 44 1
            8 3RIJ HIS D 162 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      163 45 1
            8 3RIJ HIS D 163 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      164 46 1
            8 3RIJ HIS D 164 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      165 47 1
            8 3RIJ HIS D 165 ? UNP Q5SHN1 ?   ?   'EXPRESSION TAG'      166 48 1
            #
        """
        logger.debug("Starting with %r", dataContainer.getName())
        self.__addStructRefSeqEntityIds(dataContainer)
        #
        #  To exclude self references -
        excludeRefDbList = ["PDB"]
        rD = {"seqEntityAlignmentD": {}, "seqEntityRefDbD": {}, "entityPolymerSequenceD": {}}
        try:
            # Exit if source categories are missing
            if not (dataContainer.exists("struct_ref_seq") and dataContainer.exists("struct_ref") and dataContainer.exists("entity_poly")):
                return rD
            # ------- --------- ------- --------- ------- --------- ------- --------- ------- ---------
            entityPolymerSequenceD = {}
            if dataContainer.exists("entity_poly"):
                epObj = dataContainer.getObj("entity_poly")
                for ii in range(epObj.getRowCount()):
                    entityId = epObj.getValue("entity_id", ii)
                    pType = epObj.getValue("type", ii)
                    pTypeFiltered = self.filterEntityPolyType(pType)
                    if epObj.hasAttribute("pdbx_seq_one_letter_code_can"):
                        sampleSeq = self.__stripWhiteSpace(epObj.getValue("pdbx_seq_one_letter_code_can", ii))
                        if sampleSeq and sampleSeq not in ["?", "."]:
                            entityPolymerSequenceD[entityId] = {"sequence": sampleSeq, "polymerType": pType, "polymerTypeFiltered": pTypeFiltered}
            #
            srObj = None
            if dataContainer.exists("struct_ref"):
                srObj = dataContainer.getObj("struct_ref")
            #
            srsObj = None
            if dataContainer.exists("struct_ref_seq"):
                srsObj = dataContainer.getObj("struct_ref_seq")

            # srsdObj = None
            # if dataContainer.exists("struct_ref_seq_dif"):
            #    srsdObj = dataContainer.getObj("struct_ref_seq_dif")

            polymerEntityTypeD = self.getPolymerEntityFilteredTypes(dataContainer)
            # Map alignId -> entityId
            seqEntityRefDbD = {}
            tupSeqEntityRefDbD = {}
            alignEntityMapD = {}
            # entity alignment details
            seqEntityAlignmentD = {}
            for ii in range(srObj.getRowCount()):
                dbAccessionAlignS = set()
                entityId = srObj.getValue("entity_id", ii)
                refId = srObj.getValue("id", ii)
                dbName = str(srObj.getValue("db_name", ii)).strip().upper()
                #
                if dbName in excludeRefDbList:
                    continue
                #
                if entityId not in polymerEntityTypeD:
                    logger.debug("%s skipping non-polymer entity %r in sequence reference", dataContainer.getName(), entityId)
                    continue

                if dbName in ["UNP"] and polymerEntityTypeD[entityId] != "Protein":
                    logger.debug("%s skipping inconsistent reference assignment for %s polymer type %s", dataContainer.getName(), dbName, polymerEntityTypeD[entityId])
                    continue
                #
                tS = srObj.getValue("pdbx_db_accession", ii)
                dbAccession = tS if tS and tS not in [".", "?"] else None
                #
                tS = srObj.getValueOrDefault("pdbx_db_isoform", ii, defaultValue=None)
                dbIsoform = tS if tS and tS not in [".", "?"] else None
                # Look for a stray isoform
                if dbName in ["UNP"] and dbAccession and "-" in dbAccession:
                    if not dbIsoform:
                        dbIsoform = dbAccession
                    ff = dbAccession.split("-")
                    dbAccession = ff[0]

                #
                if dbIsoform and dbAccession not in dbIsoform:
                    logger.debug("entryId %r entityId %r accession %r isoform %r inconsistency", dataContainer.getName(), entityId, dbAccession, dbIsoform)
                # ---
                # Get indices for the target refId.
                iRowL = srsObj.selectIndices(refId, "ref_id")
                logger.debug("entryId %r entityId %r refId %r rowList %r", dataContainer.getName(), entityId, refId, iRowL)
                entitySeqIdBeg = entitySeqIdEnd = 0
                for iRow in iRowL:
                    try:
                        entitySeqIdBeg = srsObj.getValue("seq_align_beg", iRow)
                        entitySeqIdEnd = srsObj.getValue("seq_align_end", iRow)
                        entityAlignLength = int(entitySeqIdEnd) - int(entitySeqIdBeg) + 1
                    except Exception:
                        entityAlignLength = 0
                    #
                    if entityAlignLength <= 0:
                        logger.debug("%s entity %r skipping bad alignment seqBeg %r seqEnd %r", dataContainer.getName(), entityId, entitySeqIdBeg, entitySeqIdEnd)
                        continue

                    alignId = srsObj.getValue("align_id", iRow)
                    alignEntityMapD[alignId] = entityId
                    #
                    authAsymId = srsObj.getValueOrDefault("pdbx_strand_id", iRow, defaultValue=None)
                    dbSeqIdBeg = srsObj.getValue("db_align_beg", iRow)
                    dbSeqIdEnd = srsObj.getValue("db_align_end", iRow)
                    # ----
                    try:
                        idbSeqIdBeg = int(dbSeqIdBeg)
                        if idbSeqIdBeg == 0:
                            idbSeqIdBeg = 1
                            dbSeqIdBeg = str(idbSeqIdBeg)
                            idbSeqIdEnd = int(dbSeqIdEnd)
                            idbSeqIdEnd += 1
                            dbSeqIdEnd = str(idbSeqIdEnd)
                            logger.debug("%s offset reference sequence database position", dataContainer.getName())
                    except Exception:
                        pass
                    # ----
                    #
                    tS = srsObj.getValueOrDefault("pdbx_db_accession", iRow, defaultValue=None)
                    # use the parent pdbx_accession
                    dbAccessionAlign = tS if tS and tS not in [".", "?"] else dbAccession
                    # Look for a stray isoform
                    if dbName in ["UNP"] and dbAccessionAlign and "-" in dbAccessionAlign:
                        if not dbIsoform:
                            dbIsoform = dbAccessionAlign
                        ff = dbAccessionAlign.split("-")
                        dbAccessionAlign = ff[0]

                    dbAccessionAlignS.add(dbAccessionAlign)
                    #
                    #
                    seqEntityAlignmentD.setdefault(entityId, []).append(
                        SeqAlign(
                            "PDB",
                            **{
                                "authAsymId": authAsymId,
                                "entitySeqIdBeg": entitySeqIdBeg,
                                "entitySeqIdEnd": entitySeqIdEnd,
                                "dbSeqIdBeg": dbSeqIdBeg,
                                "dbSeqIdEnd": dbSeqIdEnd,
                                "dbName": dbName,
                                "dbAccession": dbAccessionAlign,
                                "dbIsoform": dbIsoform,
                                "entityAlignLength": entityAlignLength,
                            },
                        )
                    )
                # Check consistency
                try:
                    if len(dbAccessionAlignS) == 1 and list(dbAccessionAlignS)[0] == dbAccession:
                        tupSeqEntityRefDbD.setdefault(entityId, []).append((dbName, dbAccession, dbIsoform))
                    elif len(dbAccessionAlignS) == 1 and list(dbAccessionAlignS)[0]:
                        tupSeqEntityRefDbD.setdefault(entityId, []).append((dbName, list(dbAccessionAlignS)[0], None))
                    elif dbAccession:
                        tupSeqEntityRefDbD.setdefault(entityId, []).append((dbName, dbAccession, dbIsoform))
                    else:
                        logger.debug("%s entityId %r inconsistent reference sequence %r %r", dataContainer.getName(), entityId, dbAccession, dbAccessionAlignS)
                except Exception:
                    logger.exception("%s entityId %r inconsistent reference sequence %r %r", dataContainer.getName(), entityId, dbAccession, dbAccessionAlignS)

            # -----
            dbMapD = self.getDatabaseNameMap()
            for entityId, tupL in tupSeqEntityRefDbD.items():
                uTupL = list(OrderedDict({tup: True for tup in tupL}).keys())
                for tup in uTupL:
                    tS = dbMapD[tup[0]] if tup[0] in dbMapD else tup[0]
                    if tup[1]:
                        seqEntityRefDbD.setdefault(entityId, []).append({"dbName": tS, "dbAccession": tup[1], "dbIsoform": tup[2]})
                    else:
                        logger.debug("%s %s skipping incomplete sequence reference assignment %r", dataContainer.getName(), entityId, tup)

            return {
                "seqEntityAlignmentD": seqEntityAlignmentD,
                "seqEntityRefDbD": seqEntityRefDbD,
                "entityPolymerSequenceD": entityPolymerSequenceD,
            }
        except Exception as e:
            logger.exception("%s failing with %s", dataContainer.getName(), str(e))
        return rD

    def __getSequenceFeatures(self, dataContainer):
        """Get point and range sequence features.

        Args:
            dataContainer (object): mmcif.api.DataContainer object instance

        Returns:
            dict : {"seqFeatureCountsD": {entityId: {"mutation": #, "conflict": # ... }, }
                    "seqMonomerFeatureD": {(entityId, seqId, compId, filteredFeature): set(feature,...), ...}
                    "seqRangeFeatureD" : {(entityId, str(beg), str(end), "artifact"): set(details)}
                    }

        """
        logger.debug("Starting with %r", dataContainer.getName())
        self.__addStructRefSeqEntityIds(dataContainer)
        #
        #  To exclude self references -
        # excludeRefDbList = ["PDB"]
        rD = {"seqFeatureCountsD": {}, "seqMonomerFeatureD": {}, "seqRangeFeatureD": {}}
        try:
            # Exit if source categories are missing
            if not (dataContainer.exists("struct_ref_seq") and dataContainer.exists("struct_ref")):
                return rD
            # ------- --------- ------- --------- ------- --------- ------- --------- ------- ---------
            # srObj = None
            # if dataContainer.exists("struct_ref"):
            #    srObj = dataContainer.getObj("struct_ref")
            #
            # srsObj = None
            # if dataContainer.exists("struct_ref_seq"):
            #    srsObj = dataContainer.getObj("struct_ref_seq")

            srsdObj = None
            if dataContainer.exists("struct_ref_seq_dif"):
                srsdObj = dataContainer.getObj("struct_ref_seq_dif")

            # polymerEntityTypeD = self.getPolymerEntityFilteredTypes(dataContainer)
            #
            # ------- --------- ------- --------- ------- --------- ------- --------- ------- ---------
            #   (entityId, seqId, compId, filteredFeature) -> set{details, ...}
            #
            seqFeatureCountsD = {}
            seqMonomerFeatureD = {}
            seqRangeFeatureD = {}
            entityArtifactD = {}
            seqIdDetailsD = {}
            if srsdObj:
                for ii in range(srsdObj.getRowCount()):
                    # alignId = srsdObj.getValue("align_id", ii)
                    #
                    # entityId = alignEntityMapD[alignId]
                    entityId = srsdObj.getValueOrDefault("rcsb_entity_id", ii, defaultValue=None)
                    if not entityId:
                        continue
                    #
                    # authAsymId = srsdObj.getValue("pdbx_pdb_strand_id", ii)
                    # dbName = srsdObj.getValue("pdbx_seq_db_name", ii)
                    #
                    # Can't rely on alignId
                    # Keep difference records for self-referenced entity sequences.
                    # if alignId not in alignEntityMapD and dbName not in excludeRefDbList:
                    #    logger.warning("%s inconsistent alignment ID %r in difference record %d", dataContainer.getName(), alignId, ii + 1)
                    #    continue
                    #
                    seqId = srsdObj.getValueOrDefault("seq_num", ii, defaultValue=None)
                    if not seqId:
                        continue
                    compId = srsdObj.getValue("mon_id", ii)
                    #
                    details = srsdObj.getValue("details", ii)
                    filteredDetails = self.filterRefSequenceDif(details)
                    if filteredDetails == "artifact":
                        try:
                            entityArtifactD.setdefault(entityId, []).append(int(seqId))
                            seqIdDetailsD[int(seqId)] = details.lower()
                        except Exception:
                            logger.debug("Incomplete sequence difference for %r %r %r %r", dataContainer.getName(), entityId, seqId, details)
                    else:
                        seqMonomerFeatureD.setdefault((entityId, seqId, compId, filteredDetails), set()).add(details.lower())
                #
                # Consolidate the artifacts as ranges -
                for entityId, sL in entityArtifactD.items():
                    # logger.debug("%s artifact ranges SL %r ranges %r", dataContainer.getName(), sL, list(self.__toRangeList(sL)))
                    srL = self.__toRangeList(sL)
                    for sr in srL:
                        seqRangeFeatureD.setdefault((entityId, str(sr[0]), str(sr[1]), "artifact"), set()).update([seqIdDetailsD[sr[0]], seqIdDetailsD[sr[1]]])
                # JDW
                # logger.info("%s seqMonomerFeatureD %r ", dataContainer.getName(), seqMonomerFeatureD)
                #
                # Tabulate sequence monomer features by entity for the filtered cases -
                for (entityId, _, _, fDetails), _ in seqMonomerFeatureD.items():
                    if entityId not in seqFeatureCountsD:
                        seqFeatureCountsD[entityId] = {"mutation": 0, "artifact": 0, "insertion": 0, "deletion": 0, "conflict": 0, "other": 0}
                    seqFeatureCountsD[entityId][fDetails] += 1
                #
                #
                # Tabulate sequence range features by entity for the filtered cases -
                for (entityId, _, _, fDetails), _ in seqRangeFeatureD.items():
                    if entityId not in seqFeatureCountsD:
                        seqFeatureCountsD[entityId] = {"mutation": 0, "artifact": 0, "insertion": 0, "deletion": 0, "conflict": 0, "other": 0}
                    seqFeatureCountsD[entityId][fDetails] += 1

            return {
                "seqFeatureCountsD": seqFeatureCountsD,
                "seqMonomerFeatureD": seqMonomerFeatureD,
                "seqRangeFeatureD": seqRangeFeatureD,
            }
        except Exception as e:
            logger.exception("%s failing with %s", dataContainer.getName(), str(e))
        return rD

    def __addStructRefSeqEntityIds(self, dataContainer):
        """Add entity ids in categories struct_ref_seq and struct_ref_seq_dir instances.

        Args:
            dataContainer (object): mmcif.api.DataContainer object instance
            catName (str): Category name
            atName (str): Attribute name

        Returns:
            bool: True for success or False otherwise

        """
        catName = "struct_ref_seq"
        try:
            logger.debug("Starting with %r %r", dataContainer.getName(), catName)
            #
            if not (dataContainer.exists(catName) and dataContainer.exists("struct_ref")):
                return False
            #
            atName = "rcsb_entity_id"
            srsObj = dataContainer.getObj(catName)
            if not srsObj.hasAttribute(atName):
                srsObj.appendAttributeExtendRows(atName, defaultValue="?")
            else:
                # skip if attribute has already been added -
                return True
            #
            srObj = dataContainer.getObj("struct_ref")
            #
            srsdObj = None
            if dataContainer.exists("struct_ref_seq_dif"):
                srsdObj = dataContainer.getObj("struct_ref_seq_dif")
                if not srsdObj.hasAttribute(atName):
                    # srsdObj.appendAttribute(atName)
                    srsdObj.appendAttributeExtendRows(atName, defaultValue="?")

            for ii in range(srObj.getRowCount()):
                entityId = srObj.getValue("entity_id", ii)
                refId = srObj.getValue("id", ii)
                #
                # Get indices for the target refId.
                iRowL = srsObj.selectIndices(refId, "ref_id")
                for iRow in iRowL:
                    srsObj.setValue(entityId, "rcsb_entity_id", iRow)
                    alignId = srsObj.getValue("align_id", iRow)
                    #
                    if srsdObj:
                        jRowL = srsdObj.selectIndices(alignId, "align_id")
                        for jRow in jRowL:
                            srsdObj.setValue(entityId, "rcsb_entity_id", jRow)

            return True
        except Exception as e:
            logger.exception("%s %s failing with %s", dataContainer.getName(), catName, str(e))
        return False

    def filterRefSequenceDif(self, details):
        filteredDetails = details
        if details.upper() in [
            "ACETYLATION",
            "CHROMOPHORE",
            "VARIANT",
            "MODIFIED RESIDUE",
            "MODIFIED",
            "ENGINEERED",
            "ENGINEERED MUTATION",
            "AMIDATION",
            "FORMYLATION",
            "ALLELIC VARIANT",
            "AUTOPHOSPHORYLATION",
            "BENZOYLATION",
            "CHEMICAL MODIFICATION",
            "CHEMICALLY MODIFIED",
            "CHROMOPHOR, REM 999",
            "CHROMOPHORE, REM 999",
            "D-CONFIGURATION",
            "ENGINEERED AND OXIDIZED CYS",
            "ENGINEERED MUTANT",
            "ENGINERED MUTATION",
            "HYDROXYLATION",
            "METHYLATED ASN",
            "METHYLATION",
            "MICROHETEROGENEITY",
            "MODEIFED RESIDUE",
            "MODIFICATION",
            "MODIFIED AMINO ACID",
            "MODIFIED CHROMOPHORE",
            "MODIFIED GLN",
            "MODIFIED RESIDUES",
            "MUTATION",
            "MYC EPITOPE",
            "MYRISTOYLATED",
            "MYRISTOYLATION",
            "NATURAL VARIANT",
            "NATURAL VARIANTS",
            "OXIDIZED CY",
            "OXIDIZED CYS",
            "PHOSPHORYLATION",
            "POLYMORPHIC VARIANT",
            "PROPIONATION",
            "SOMATIC VARIANT",
            "SUBSTITUTION",
            "TRNA EDITING",
            "TRNA MODIFICATION",
            "TRNA",
            "VARIANT STRAIN",
            "VARIANTS",
        ]:
            filteredDetails = "mutation"
        elif details.upper() in [
            "LEADER SEQUENCE",
            "INITIATING METHIONINE",
            "INITIATOR METHIONINE",
            "LINKER",
            "EXPRESSION TAG",
            "CLONING",
            "CLONING ARTIFACT",
            "C-TERM CLONING ARTIFA",
            "C-TERMINAL HIS TAG",
            "C-TERMINLA HIS-TAG",
            "CLONING AETIFACT",
            "CLONING ARATIFACT",
            "CLONING ARTEFACT",
            "CLONING ARTFIACT",
            "CLONING ARTIACT",
            "CLONING ARTIFACTS",
            "CLONING ARTUFACT",
            "CLONING ATIFACT",
            "CLONING MUTATION",
            "CLONING REMNANT",
            "CLONING SITE RESIDUE",
            "CLONNG ARTIFACT",
            "CLONONG ARTIFACT",
            "DETECTION TAG",
            "ENGINEERED LINKER",
            "EXPRESSION ARTIFACT",
            "EXPRESSIOPN TAG",
            "EXPRSSION TAG",
            "FLAG TAG",
            "GCN4 TAG",
            "GPGS TAG",
            "GST TAG",
            "HIA TAG",
            "HIS TAG",
            "HIS-TAG",
            "INITIAL METHIONINE",
            "INITIATING MET",
            "INITIATING METHIONIE",
            "INITIATING MSE",
            "INITIATING RESIDUE",
            "INITIATOR N-FORMYL-MET",
            "INTIATING METHIONINE",
            "INTRACHAIN HIS TAG",
            "LINKER INSERTION",
            "LINKER PEPTIDE",
            "LINKER RESIDUE",
            "LINKER SEQUENCE",
            "LYS TAG",
            "MOD. RESIDUE/CLONING ARTIFACT",
            "MYC TAG",
            "N-TERMINAL EXTENSION",
            "N-TERMINAL HIS TAG",
            "PURIFICATION TAG",
            "RANDOM MUTAGENESIS",
            "RECOMBINANT HIS TAG",
            "RESIDUAL LINKER",
            "STREP-TAGII",
            "T7 EPITOPE TAG",
            "T7-TAG",
            "TAG",
        ]:
            filteredDetails = "artifact"
        elif details.upper() in ["INSERTION", "ENGINEERED INSERTION", "INSERTED", "INSERTION AT N-TERMINUS"]:
            filteredDetails = "insertion"
        elif details.upper() in ["DELETION", "CONFLICT/DELETION", "ENGINEERED DELETION"]:
            filteredDetails = "deletion"
        elif details.upper() in ["CONFLICT", "SEQUENCE CONFLICT", "SEQUENCE CONFLICT8"]:
            filteredDetails = "conflict"
        else:
            logger.debug("Unanticipated sequence difference details %r", details)
            filteredDetails = "other"
        #
        return filteredDetails

    def filterEntityPolyType(self, pType):
        """Map input dictionary polymer type to simplified molecular type.

        Args:
            pType (str): PDBx/mmCIF dictionary polymer type

        Returns:
            str: simplified mappings

        Returns mappings:
            'Protein'   'polypeptide(D) or polypeptide(L)'
            'DNA'       'polydeoxyribonucleotide'
            'RNA'       'polyribonucleotide'
            'NA-hybrid' 'polydeoxyribonucleotide/polyribonucleotide hybrid'
            'Other'      'polysaccharide(D), polysaccharide(L), cyclic-pseudo-peptide, peptide nucleic acid, or other'
        """
        polymerType = pType.lower()
        if polymerType in ["polypeptide(d)", "polypeptide(l)"]:
            rT = "Protein"
        elif polymerType in ["polydeoxyribonucleotide"]:
            rT = "DNA"
        elif polymerType in ["polyribonucleotide"]:
            rT = "RNA"
        elif polymerType in ["polydeoxyribonucleotide/polyribonucleotide hybrid"]:
            rT = "NA-hybrid"
        else:
            rT = "Other"
        return rT

    def guessEntityPolyTypes(self, monomerL):
        """Guess the polymer types to from the monomer list.

        Args:
            monomerL (list): list of monomers (chemical component ids)

        Returns:
            tuple: polymerType, filtered polymer Type.

        Returns mappings:
            'Protein'   'polypeptide(D) or polypeptide(L)'
            'DNA'       'polydeoxyribonucleotide'
            'RNA'       'polyribonucleotide'
            'NA-hybrid' 'polydeoxyribonucleotide/polyribonucleotide hybrid'
            'Other'      'polysaccharide(D), polysaccharide(L), cyclic-pseudo-peptide, peptide nucleic acid, or other'
        """
        hasAA = hasDNA = hasRNA = False
        pType = fpType = None
        for monomer in monomerL:
            if monomer in DictMethodCommonUtils.aaDict3:
                hasAA = True
            elif monomer in DictMethodCommonUtils.dnaDict3:
                hasDNA = True
            elif monomer in DictMethodCommonUtils.rnaDict3:
                hasRNA = True
        #
        if hasAA and not hasDNA and not hasRNA:
            pType = "polypeptide(d)"
        elif hasDNA and not hasAA and not hasRNA:
            pType = "polydeoxyribonucleotide"
        elif hasRNA and not hasAA and not hasDNA:
            pType = "polyribonucleotide"
        elif not hasAA and hasDNA and hasRNA:
            pType = "polydeoxyribonucleotide/polyribonucleotide hybrid"

        if pType:
            fpType = self.filterEntityPolyType(pType)
        else:
            pType = None
            fpType = "Other"
        #
        return pType, fpType

    def getPolymerComposition(self, polymerTypeList):
        """Map in list of dictionary entity polymer/branched types to a composition string.
            Input polymerTypeList contains entity_poly.type and pdbx_entity_branch.type values.

        Args:
            polymerTypeList (list): List of PDBx/mmCIF dictionary polymer/branched types

        Returns:
            tuple: compClass, ptClass, naClass, cD

                   compClass - simplified composition string
                   ptClass - subset class
                   naClass - nucleic acid subset class
                   cD (dict) - composition type counts

        Current polymer type list:
             'polypeptide(D)'
             'polypeptide(L)'
             'polydeoxyribonucleotide'
             'polyribonucleotide'
             'polysaccharide(D)'
             'polysaccharide(L)'
             'polydeoxyribonucleotide/polyribonucleotide hybrid'
             'cyclic-pseudo-peptide'
             'peptide nucleic acid'
             'other'
             "other type pair (polymer type count = 2)"
             "other composition (polymer type count >= 3)"

        Current branch type list:
             'oligosaccharide'

        Output composition classes:

            'homomeric protein' 'single protein entity'
            'heteromeric protein' 'multiple protein entities'
            'DNA' 'DNA entity/entities only'
            'RNA' 'RNA entity/entities only'
            'NA-hybrid' 'DNA/RNA hybrid entity/entities only'
            'protein/NA' 'Both protein and nucleic acid polymer entities'
            'DNA/RNA' 'Both DNA and RNA polymer entities'
            'oligosaccharide' 'One of more oligosaccharide entities'
            'protein/oligosaccharide' 'Both protein and oligosaccharide entities'
            'NA/oligosaccharide' 'Both NA and oligosaccharide entities'
            'other' 'Neither an individual protein, nucleic acid polymer nor oligosaccharide entity'
            'other type pair' 'Other combinations of 2 polymer types'
            'other type composition' 'Other combinations of 3 or more polymer types'

        And selected types (ptClass)-
            'Protein (only)' 'protein entity/entities only'
            'Nucleic acid (only)' 'DNA, RNA or NA-hybrid entity/entities only'
            'Protein/NA' 'Both protein and nucleic acid (DNA, RNA, or NA-hybrid) polymer entities'
            'Other' 'Another polymer type composition'

        And selected NA types (naClass) -
            'DNA (only)' 'DNA entity/entities only'
            'RNA (only)' 'RNA entity/entities only'
            'NA-hybrid (only)' 'NA-hybrid entity/entities only'
            'DNA/RNA (only)' 'Both DNA and RNA polymer entities only'
            'Other' 'Another polymer type composition'
        """

        compClass = "other"
        # get type counts
        cD = {}
        for polymerType in polymerTypeList:
            if polymerType in ["polypeptide(D)", "polypeptide(L)"]:
                cD["protein"] = cD["protein"] + 1 if "protein" in cD else 1
            elif polymerType in ["polydeoxyribonucleotide"]:
                cD["DNA"] = cD["DNA"] + 1 if "DNA" in cD else 1
            elif polymerType in ["polyribonucleotide"]:
                cD["RNA"] = cD["RNA"] + 1 if "RNA" in cD else 1
            elif polymerType in ["polydeoxyribonucleotide/polyribonucleotide hybrid"]:
                cD["NA-hybrid"] = cD["NA-hybrid"] + 1 if "NA-hybrid" in cD else 1
            elif polymerType in ["oligosaccharide"]:
                cD["oligosaccharide"] = cD["oligosaccharide"] + 1 if "oligosaccharide" in cD else 1
            else:
                cD["other"] = cD["other"] + 1 if "other" in cD else 1
        #
        if len(cD) == 1:
            ky = list(cD.keys())[0]
            if "protein" in cD:
                if cD["protein"] == 1:
                    compClass = "homomeric protein"
                else:
                    compClass = "heteromeric protein"
            elif ky in ["DNA", "RNA", "NA-hybrid", "oligosaccharide", "other"]:
                compClass = ky
        elif len(cD) == 2:
            if "protein" in cD:
                if ("DNA" in cD) or ("RNA" in cD) or ("NA-hybrid" in cD):
                    compClass = "protein/NA"
                elif "oligosaccharide" in cD:
                    compClass = "protein/oligosaccharide"
            elif "DNA" in cD and "RNA" in cD:
                compClass = "DNA/RNA"
            elif "oligosaccharide" in cD and ("RNA" in cD or "DNA" in cD):
                compClass = "NA/oligosaccharide"
            else:
                compClass = "other type pair"
        elif len(cD) == 3:
            if "DNA" in cD and "RNA" in cD and "NA-hybrid" in cD:
                compClass = "DNA/RNA"
            elif "oligosaccharide" in cD and all([j in ["oligosaccharide", "DNA", "RNA", "NA-hybrid"] for j in cD]):
                compClass = "NA/oligosaccharide"
            elif "protein" in cD and all([j in ["protein", "DNA", "RNA", "NA-hybrid"] for j in cD]):
                compClass = "protein/NA"
            elif "oligosaccharide" in cD and "protein" in cD and all([j in ["protein", "oligosaccharide", "DNA", "RNA", "NA-hybrid"] for j in cD]):
                compClass = "protein/NA/oligosaccharide"
            else:
                compClass = "other type composition"
        elif len(cD) >= 4:
            if "oligosaccharide" in cD and all([j in ["oligosaccharide", "DNA", "RNA", "NA-hybrid"] for j in cD]):
                compClass = "NA/oligosaccharide"
            elif "protein" in cD and all([j in ["protein", "DNA", "RNA", "NA-hybrid"] for j in cD]):
                compClass = "protein/NA"
            elif "oligosaccharide" in cD and "protein" in cD and all([j in ["protein", "oligosaccharide", "DNA", "RNA", "NA-hybrid"] for j in cD]):
                compClass = "protein/NA/oligosaccharide"
            else:
                compClass = "other type composition"
        else:
            compClass = "other type composition"

        # Subset type class --
        #
        if compClass in ["homomeric protein", "heteromeric protein"]:
            ptClass = "Protein (only)"
        elif compClass in ["DNA", "RNA", "NA-hybrid", "DNA/RNA"]:
            ptClass = "Nucleic acid (only)"
        elif compClass in ["protein/NA"]:
            ptClass = "Protein/NA"
        # JDW
        elif compClass in ["protein/oligosaccharide"]:
            ptClass = "Protein/Oligosaccharide"
        elif compClass in ["oligosaccharide"]:
            ptClass = "Oligosaccharide (only)"
        # elif compClass in ["protein/NA/oligosaccharide"]:
        #    ptClass = "Protein/NA/Oligosaccharide"
        # JDW
        else:
            ptClass = "Other"
        #
        # NA subtype class ---
        #
        if compClass in ["DNA"]:
            naClass = "DNA (only)"
        elif compClass in ["RNA"]:
            naClass = "RNA (only)"
        elif compClass in ["NA-hybrid"]:
            naClass = "NA-hybrid (only)"
        elif compClass in ["DNA/RNA"]:
            naClass = "DNA/RNA (only)"
        else:
            naClass = "Other"
        #
        return compClass, ptClass, naClass, cD

    def filterExperimentalMethod(self, dataContainer):
        """Apply a standard filter to the input experimental method list returning a method count and
            a simplified method name.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance

        Returns:
            tuple(int,str): methodCount, simpleMethodName

        For example:
        'X-ray'            'X-RAY DIFFRACTION, FIBER DIFFRACTION, or POWDER DIFFRACTION'
        'NMR'              'SOLUTION NMR or SOLID-STATE NMR'
        'EM'               'ELECTRON MICROSCOPY or ELECTRON CRYSTALLOGRAPHY or ELECTRON TOMOGRAPHY'
        'Neutron'          'NEUTRON DIFFRACTION'
        'Multiple methods' 'Multiple experimental methods'
        'Other'            'SOLUTION SCATTERING, EPR, INFRARED SPECTROSCOPY or FLUORESCENCE TRANSFER'
        """
        methodL = self.getMethodList(dataContainer)
        methodCount = len(methodL)
        if methodCount > 1:
            expMethod = "Multiple methods"
        else:
            #
            mS = methodL[0].upper()
            expMethod = "Other"
            if mS in ["X-RAY DIFFRACTION", "FIBER DIFFRACTION", "POWDER DIFFRACTION"]:
                expMethod = "X-ray"
            elif mS in ["SOLUTION NMR", "SOLID-STATE NMR"]:
                expMethod = "NMR"
            elif mS in ["ELECTRON MICROSCOPY", "ELECTRON CRYSTALLOGRAPHY", "ELECTRON DIFFRACTION", "CRYO-ELECTRON MICROSCOPY", "ELECTRON TOMOGRAPHY"]:
                expMethod = "EM"
            elif mS in ["NEUTRON DIFFRACTION"]:
                expMethod = "Neutron"
            elif mS in ["SOLUTION SCATTERING", "EPR", "INFRARED SPECTROSCOPY", "FLUORESCENCE TRANSFER"]:
                expMethod = "Other"
            elif mS in ["THEORETICAL MODEL", "AB INITIO MODEL", "HOMOLOGY MODEL"]:
                expMethod = None
            else:
                logger.error("Unexpected method ")

        return methodCount, expMethod

    def filterStructureDeterminationMethodType(self, dataContainer):
        """Apply a standard filter to the input experimental or computational method list to return the type
            of structure determination method used--experimental, computational, or integrative.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance

        Returns:
            str: methodType

        For example:
            experimental   "Experimentally based structure determination"
            integrative    "Integrative/Hybrid methods"
            computational  "Computational modeling"
        """
        methodL = self.getMethodList(dataContainer)
        #
        experimentalMethodList = [
            "X-RAY DIFFRACTION", "FIBER DIFFRACTION", "POWDER DIFFRACTION",
            "SOLUTION NMR", "SOLID-STATE NMR",
            "ELECTRON MICROSCOPY", "ELECTRON CRYSTALLOGRAPHY", "ELECTRON DIFFRACTION", "CRYO-ELECTRON MICROSCOPY", "ELECTRON TOMOGRAPHY",
            "NEUTRON DIFFRACTION",
            "SOLUTION SCATTERING", "EPR", "INFRARED SPECTROSCOPY", "FLUORESCENCE TRANSFER",
        ]
        #
        computationalMethodList = [
            "THEORETICAL MODEL", "AB INITIO MODEL", "HOMOLOGY MODEL",
        ]
        #
        methodType = None
        methodCount = len(methodL)
        if methodCount > 1:
            if any([mS.upper() in experimentalMethodList for mS in methodL]):
                methodType = "experimental"
            if all([mS.upper() in computationalMethodList for mS in methodL]):
                methodType = "computational"
            # if any([mS.upper() in experimentalMethodList for mS in methodL]) and any([mS.upper() in computationalMethodList for mS in methodL]):
            #     methodType = "integrative"
        else:
            #
            mS = methodL[0].upper()
            if mS in experimentalMethodList:
                methodType = "experimental"
            elif mS in computationalMethodList:
                methodType = "computational"
            else:
                logger.error("Unexpected method type")

        return methodType

    def hasMethodNMR(self, dataContainer):
        """Return if the input dictionary experimental method list contains an NMR experimental method.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance

        Returns:
            bool: True if the input contains NMR or False otherwise
        """
        ok = False
        #
        methodL = self.getMethodList(dataContainer)
        if methodL:
            for method in methodL:
                if method in ["SOLUTION NMR", "SOLID-STATE NMR"]:
                    return True
        return ok

    def __getTimeStamp(self):
        utcnow = datetime.datetime.utcnow()
        ts = utcnow.strftime("%Y-%m-%d:%H:%M:%S")
        return ts

    def __stripWhiteSpace(self, val):
        """Remove all white space from the input value."""
        if val is None:
            return val
        return self.__wsPattern.sub("", val)

    def __toRangeList(self, iterable):
        iterable = sorted(set(iterable))
        for _, group in itertools.groupby(enumerate(iterable), lambda t: t[1] - t[0]):
            group = list(group)
            yield group[0][1], group[-1][1]

    #
    def getTargetSiteInfo(self, dataContainer):
        """Return a dictionary of target site binding interactions using standard nomenclature.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {site_id: [{'asymId': , 'compId': , 'seqId': }, ...],  ... }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchInstanceSiteInfo(dataContainer)
        return wD["targetSiteD"] if "targetSiteD" in wD else {}

    def getLigandSiteInfo(self, dataContainer):
        """Return a dictionary of ligand site binding interactions.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {site_id: {"evCode": software|author,
                            "fromDetails": True|False,
                            "isRaw": True|False,
                            "entityType": polymer|non-polymer,
                            "polymerLigand": {"asymId": ., "entityId": ., "begSeqId": ., "endSeqId":. },
                            "nonPolymerLigands": [{"asymId": ., "entityId": ., "compId": .}, ...],
                            "description": raw or generated text,
                            }
                            }
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchInstanceSiteInfo(dataContainer)
        return wD["ligandSiteD"] if "ligandSiteD" in wD else {}

    def __fetchInstanceSiteInfo(self, dataContainer):
        wD = self.__instanceSiteInfoCache.get(dataContainer.getName())
        if not wD:
            wD = self.__getInstanceSiteInfo(dataContainer)
            self.__instanceSiteInfoCache.set(dataContainer.getName(), wD)
        return wD

    def __getInstanceSiteInfo(self, dataContainer):
        """[summary]

        Args:
            dataContainer (object): mmcif.api.DataContainer object instance

        Returns:
            dict : {"targetSiteD" = {<site_id>: {}}
                    "ligandSiteD": {<site_id>: {}}
                    }

        For example:

                loop_
                _struct_site.id
                _struct_site.pdbx_evidence_code
                _struct_site.pdbx_auth_asym_id
                _struct_site.pdbx_auth_comp_id
                _struct_site.pdbx_auth_seq_id
                _struct_site.pdbx_auth_ins_code # never used
                _struct_site.pdbx_num_residues
                _struct_site.details
                AC1 Software ? ? ? ? 7  'BINDING SITE FOR RESIDUE ADP A 105'
                AC2 Software ? ? ? ? 16 'BINDING SITE FOR RESIDUE ADP B 101'
                AC3 Software ? ? ? ? 6  'BINDING SITE FOR RESIDUE MG B 66'
                AC4 Software ? ? ? ? 13 'BINDING SITE FOR RESIDUE ADP C 102'
                AC5 Software ? ? ? ? 16 'BINDING SITE FOR RESIDUE ADP E 103'
                AC6 Software ? ? ? ? 10 'BINDING SITE FOR RESIDUE ADP F 104'
                AC7 Software ? ? ? ? 6  'BINDING SITE FOR RESIDUE MG K 9'
                #
                loop_
                _struct_site_gen.id
                _struct_site_gen.site_id
                _struct_site_gen.pdbx_num_res
                _struct_site_gen.label_comp_id
                _struct_site_gen.label_asym_id
                _struct_site_gen.label_seq_id
                _struct_site_gen.pdbx_auth_ins_code
                _struct_site_gen.auth_comp_id
                _struct_site_gen.auth_asym_id
                _struct_site_gen.auth_seq_id
                _struct_site_gen.label_atom_id
                _struct_site_gen.label_alt_id
                _struct_site_gen.symmetry
                _struct_site_gen.details
                1  AC1 7  TYR A 25 ? TYR A 25  . ? 1_555 ?
                2  AC1 7  GLY A 29 ? GLY A 29  . ? 1_555 ?
                3  AC1 7  THR A 61 ? THR A 61  . ? 1_555 ?
                4  AC1 7  VAL A 63 ? VAL A 63  . ? 1_555 ?
                5  AC1 7  ILE B 30 ? ILE B 30  . ? 1_555 ?
                6  AC1 7  LEU B 32 ? LEU B 32  . ? 1_555 ?
                7  AC1 7  GLN B 52 ? GLN B 52  . ? 1_555 ?
                8  AC2 16 TYR B 25 ? TYR B 25  . ? 1_555 ?
                9  AC2 16 LEU B 26 ? LEU B 26  . ? 1_555 ?
                10 AC2 16 GLY B 29 ? GLY B 29  . ? 1_555 ?
                11 AC2 16 LYS B 31 ? LYS B 31  . ? 1_555 ?
                12 AC2 16 SER B 60 ? SER B 60  . ? 1_555 ?
                13 AC2 16 THR B 61 ? THR B 61  . ? 1_555 ?
                14 AC2 16 HOH P .  ? HOH B 113 . ? 1_555 ?
                15 AC2 16 HOH P .  ? HOH B 116 . ? 1_555 ?
                16 AC2 16 HOH P .  ? HOH B 201 . ? 1_555 ?
                17 AC2 16 HOH P .  ? HOH B 241 . ? 1_555 ?
                18 AC2 16 LEU C 26 ? LEU C 26  . ? 1_555 ?
                19 AC2 16 ASN C 28 ? ASN C 28  . ? 1_555 ?
                20 AC2 16 ILE C 30 ? ILE C 30  . ? 1_555 ?
                21 AC2 16 LEU C 32 ? LEU C 32  . ? 1_555 ?
                22 AC2 16 ARG F 16 ? ARG F 16  . ? 1_565 ?
                23 AC2 16 ARG F 17 ? ARG F 17  . ? 1_565 ?
        """
        logger.debug("Starting with %r", dataContainer.getName())
        #
        rD = {"targetSiteD": {}, "ligandSiteD": {}}
        try:
            # Exit if source categories are missing
            if not (dataContainer.exists("struct_site") and dataContainer.exists("struct_site_gen")):
                return rD
            # ------- --------- ------- --------- ------- --------- ------- --------- ------- ---------
            ssObj = None
            if dataContainer.exists("struct_site"):
                ssObj = dataContainer.getObj("struct_site")
            #
            ssgObj = None
            if dataContainer.exists("struct_site_gen"):
                ssgObj = dataContainer.getObj("struct_site_gen")

            #
            ligandSiteD = {}
            for ii in range(ssObj.getRowCount()):
                ligL = []
                evCode = str(ssObj.getValue("pdbx_evidence_code", ii)).lower()
                if evCode not in ["software", "author"]:
                    continue
                sId = ssObj.getValue("id", ii)
                authAsymId = ssObj.getValueOrDefault("pdbx_auth_asym_id", ii, defaultValue=None)
                compId = ssObj.getValueOrDefault("pdbx_auth_comp_id", ii, defaultValue=None)
                authSeqId = ssObj.getValueOrDefault("pdbx_auth_seq_id", ii, defaultValue=None)
                ssDetails = ssObj.getValueOrDefault("details", ii, defaultValue=None)
                fromDetails = False
                if authAsymId:
                    ligL.append((authAsymId, compId, authSeqId, ssDetails))
                else:
                    fromDetails = True
                    if evCode == "software":
                        ligL = self.__parseStructSiteLigandDetails(ssDetails)
                    elif evCode == "author":
                        ligL.append((None, None, None, ssDetails))
                #
                ligandSiteD[sId] = self.__transStructSiteLigandDetails(dataContainer, ligL, evCode=evCode, fromDetails=fromDetails)
            #

            targetSiteD = {}
            instTypeD = self.getInstanceTypes(dataContainer)
            for ii in range(ssgObj.getRowCount()):
                sId = ssgObj.getValue("site_id", ii)
                asymId = ssgObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                compId = ssgObj.getValueOrDefault("label_comp_id", ii, defaultValue=None)
                seqId = ssgObj.getValueOrDefault("label_seq_id", ii, defaultValue=None)
                #
                if asymId and compId and seqId and asymId in instTypeD and instTypeD[asymId] == "polymer":
                    targetSiteD.setdefault(sId, []).append({"asymId": asymId, "compId": compId, "seqId": seqId})
            #
            return {"targetSiteD": targetSiteD, "ligandSiteD": ligandSiteD}
        except Exception as e:
            logger.exception("%s failing with %s", dataContainer.getName(), str(e))
        return rD

    def __transStructSiteLigandDetails(self, dataContainer, ligL, evCode="software", fromDetails=True):
        """Convert struct_site ligand details to standard nomenclature.

        Args:
            dataContainer (object): mmcif.api.DataContainer object instance
            ligL (list): list of raw ligand details in author nomenclature
            evCode (str):  string  (software|author)
            fromDetails (bool, optional): details parsed from descriptive text. Defaults to True.

        Returns:
            dict: {"evCode": software|author,
                   "fromDetails": True|False,
                   "isRaw": True|False,
                   "entityType": polymer|non-polymer,
                   "polymerLigand": {"asymId": ., "entityId": ., "begSeqId": ., "endSeqId":. },
                   "nonPolymerLigands": [{"asymId": ., "entityId": ., "compId": .}, ...],
                   "description": raw or generated text,
                   "siteLabel": replacement for data site id,
                   }

        """
        rD = {
            "evCode": evCode,
            "fromDetails": fromDetails,
            "isRaw": True,
            "entityType": None,
            "polymerLigand": None,
            "nonPolymerLigands": None,
            "description": None,
            "siteLabel": None,
        }
        npAuthAsymD = self.getNonPolymerIdMap(dataContainer)
        pAuthAsymD = self.getPolymerIdMap(dataContainer)
        asymAuthIdD = self.getAsymAuthIdMap(dataContainer)
        asymIdPolymerRangesD = self.getInstancePolymerRanges(dataContainer)
        iTypeD = self.getInstanceTypes(dataContainer)
        asymAuthIdD = self.getAsymAuthIdMap(dataContainer)
        # Note that this is a non-unique index inversion
        authAsymD = {v: k for k, v in asymAuthIdD.items()}
        instEntityD = self.getInstanceEntityMap(dataContainer)
        evS = "Software generated" if evCode == "software" else "Author provided"
        #
        if len(ligL) == 1:
            authAsymId, compId, authSeqId, ssDetails = ligL[0]
            #
            if not authAsymId:
                rD["description"] = ssDetails
                rD["isRaw"] = True
            elif not authSeqId:
                # An unqualified authAsymId -
                asymId = authAsymD[authAsymId] if authAsymId in authAsymD else None
                entityId = instEntityD[asymId] if asymId in instEntityD else None
                if entityId and asymId and asymId in iTypeD and iTypeD[asymId] == "polymer" and asymId in asymIdPolymerRangesD:
                    # insert the full residue range -
                    rD["entityType"] = iTypeD[asymId]
                    begSeqId = asymIdPolymerRangesD[asymId]["begSeqId"]
                    endSeqId = asymIdPolymerRangesD[asymId]["endSeqId"]
                    tD = {"asymId": asymId, "entityId": instEntityD[asymId], "begSeqId": begSeqId, "endSeqId": endSeqId}
                    rD["description"] = "%s binding site for entity %s (%s-%s) instance %s chain %s" % (evS, entityId, begSeqId, endSeqId, asymId, authAsymId)
                    rD["polymerLigand"] = tD
                    rD["siteLabel"] = "chain %s" % authAsymId
            elif (authAsymId, authSeqId) in npAuthAsymD:
                # single non-polymer-ligand -
                asymId = npAuthAsymD[(authAsymId, authSeqId)]["asym_id"]
                rD["entityType"] = iTypeD[asymId]
                entityId = instEntityD[asymId]
                tD = {"asymId": asymId, "entityId": instEntityD[asymId], "compId": compId}
                rD["nonPolymerLigands"] = [tD]
                rD["description"] = "%s binding site for ligand entity %s component %s instance %s chain %s" % (evS, entityId, compId, asymId, authAsymId)
                rD["siteLabel"] = "ligand %s" % compId
            elif (authAsymId, authSeqId, None) in pAuthAsymD:
                # single monomer ligand - an odd case
                asymId = pAuthAsymD[(authAsymId, authSeqId, None)]["asym_id"]
                entityId = pAuthAsymD[(authAsymId, authSeqId, None)]["entity_id"]
                seqId = pAuthAsymD[(authAsymId, authSeqId, None)]["seq_id"]
                rD["entityType"] = iTypeD[asymId]
                tD = {"asymId": asymId, "entityId": entityId, "begSeqId": seqId, "endSeqId": seqId}
                rD["description"] = "%s binding site for entity %s instance %s chainId %s (%s)" % (evS, entityId, asymId, authAsymId, authSeqId)
                rD["polymerLigand"] = tD
                rD["siteLabel"] = "chain %s" % authAsymId
            else:
                logger.debug("%s untranslated single ligand details %r", dataContainer.getName(), ligL)
                logger.debug("npAuthAsymD %r", npAuthAsymD)
                rD["description"] = ssDetails
                rD["isRaw"] = True
            #
        elif len(ligL) == 2:
            authAsymIdA, compIdA, authSeqIdA, ssDetailsA = ligL[0]
            authAsymIdB, compIdB, authSeqIdB, _ = ligL[1]
            #
            # is np
            if (authAsymIdA, authSeqIdA) in npAuthAsymD and (authAsymIdB, authSeqIdB) in npAuthAsymD:
                asymIdA = npAuthAsymD[(authAsymIdA, authSeqIdA)]["asym_id"]
                entityIdA = npAuthAsymD[(authAsymIdA, authSeqIdA)]["entity_id"]
                asymIdB = npAuthAsymD[(authAsymIdB, authSeqIdB)]["asym_id"]
                entityIdB = npAuthAsymD[(authAsymIdB, authSeqIdB)]["entity_id"]
                tDA = {"asymId": asymIdA, "entityId": entityIdA, "compId": compIdA}
                tDB = {"asymId": asymIdB, "entityId": entityIdB, "compId": compIdB}
                rD["nonPolymerLigands"] = [tDA, tDB]
                rD["entityType"] = iTypeD[asymIdA]
                rD["description"] = "%s binding site for ligands: entity %s component %s instance %s chain %s and entity %s component %s instance %s chain %s" % (
                    evS,
                    entityIdA,
                    compIdA,
                    asymIdA,
                    authAsymIdA,
                    entityIdB,
                    compIdB,
                    asymIdB,
                    authAsymIdB,
                )
                rD["siteLabel"] = "ligands %s/%s" % (compIdA, compIdB)
            elif (authAsymIdA, authSeqIdA, None) in pAuthAsymD and (authAsymIdB, authSeqIdB, None) in pAuthAsymD and authAsymIdA == authAsymIdB:
                asymIdA = pAuthAsymD[(authAsymIdA, authSeqIdA, None)]["asym_id"]
                entityIdA = pAuthAsymD[(authAsymIdA, authSeqIdA, None)]["entity_id"]
                asymIdB = pAuthAsymD[(authAsymIdB, authSeqIdB, None)]["asym_id"]
                entityIdB = pAuthAsymD[(authAsymIdB, authSeqIdB, None)]["entity_id"]
                begSeqId = pAuthAsymD[(authAsymIdA, authSeqIdA, None)]["seq_id"]
                endSeqId = pAuthAsymD[(authAsymIdB, authSeqIdB, None)]["seq_id"]
                tD = {"asymId": asymIdA, "entityId": instEntityD[asymIdA], "begSeqId": begSeqId, "endSeqId": endSeqId}
                rD["entityType"] = iTypeD[asymIdA]
                rD["description"] = "%s binding site for entity %s instance %s chain %s and entity %s instance %s chain %s" % (
                    evS,
                    entityIdA,
                    asymIdA,
                    authAsymIdA,
                    entityIdB,
                    asymIdB,
                    authAsymIdB,
                )
                rD["polymerLigand"] = tD
                rD["siteLabel"] = "chains %s/%s" % (authAsymIdA, authAsymIdB)
            else:
                logger.debug("%s untranslated ligand details %r", dataContainer.getName(), ligL)
                rD["description"] = ssDetailsA
                rD["isRaw"] = True
        else:
            logger.error("%s unexpected ligand expression %r", dataContainer.getName(), ligL)
        return rD

    def __parseStructSiteLigandDetails(self, ssDetails):
        """Parse the input site description text and returning structured details
        where possible.

        Args:
            ssDetails (str): struct_site.details text

        Returns:
            list: [(authAsymId, compId, authSeqId, ssDetails), ... ]

        """
        retL = []
        #
        try:
            if not ssDetails:
                retL.append((None, None, None, None))
                return retL
            prefixL = [
                "BINDING SITE FOR RESIDUE ",
                "binding site for residue ",
                "Binding site for Ligand ",
                "binding site for Ligand ",
                "Binding site for Mono-Saccharide ",
                "BINDING SITE FOR MONO-SACCHARIDE ",
                "binding site for Mono-Saccharide ",
                "binding site for Poly-Saccharide ",
                "binding site for nucleotide ",
            ]
            for prefix in prefixL:
                tup = ssDetails.partition(prefix)
                if tup[1] == prefix:
                    ff = tup[2].split(" ")
                    # binding site for Ligand residues POL d 4 through N7P d 1 bound to THR b 1
                    if ff[0] == "residues" and len(ff) > 8 and ff[4].lower() == "through":
                        compIdA = ff[1]
                        authAsymIdA = ff[2]
                        authSeqIdA = ff[3]
                        retL.append((authAsymIdA, compIdA, authSeqIdA, ssDetails))
                        #
                        compIdB = ff[5]
                        authAsymIdB = ff[6]
                        authSeqIdB = ff[7]
                        retL.append((authAsymIdB, compIdB, authSeqIdB, ssDetails))
                        return retL
                    elif len(ff) == 2:
                        compId = ff[0]
                        authAsymId = ff[1][0]
                        authSeqId = ff[1][1:]
                        retL.append((authAsymId, compId, authSeqId, ssDetails))
                        return retL
                    elif len(ff) == 3:
                        compId = ff[0]
                        authAsymId = ff[1]
                        authSeqId = ff[2]
                        retL.append((authAsymId, compId, authSeqId, ssDetails))
                        return retL

            #
            # Binding site for residues GCD A 900 and NGA A 901
            # Binding site for residues FUC A1118 and BGC A1119'
            prefixL = [
                "Binding site for residues ",
                "binding site for residues ",
                "BINDING SITE FOR DI-SACCHARIDE ",
                "Binding site for Di-Saccharide ",
                "binding site for Di-Saccharide ",
                "binding site for Di-peptide ",
                "Binding site for Di-peptide ",
                "binding site for Di-nucleotide ",
            ]
            for prefix in prefixL:
                tup = ssDetails.partition(prefix)
                if tup[1] == prefix:
                    ff = tup[2].split(" ")
                    if len(ff) == 5:
                        compIdA = ff[0]
                        authAsymIdA = ff[1][0]
                        authSeqIdA = ff[1][1:]
                        compIdB = ff[3]
                        authAsymIdB = ff[4][0]
                        authSeqIdB = ff[4][1:]
                    elif len(ff) == 7:
                        compIdA = ff[0]
                        authAsymIdA = ff[1]
                        authSeqIdA = ff[2]
                        compIdB = ff[4]
                        authAsymIdB = ff[5]
                        authSeqIdB = ff[6]
                    else:
                        compIdA = authAsymIdA = authSeqIdA = compIdB = authAsymIdB = authSeqIdB = None

                    retL.append((authAsymIdA, compIdA, authSeqIdA, ssDetails))
                    retL.append((authAsymIdB, compIdB, authSeqIdB, ssDetails))
                    return retL
            #
            # BINDING SITE FOR LINKED RESIDUES A 1519 A 1520 A 1521 A 1522 A 1523 A 1524 A 1525
            # BINDING SITE FOR LINKED RESIDUES A 801 to 802
            prefixL = ["BINDING SITE FOR LINKED RESIDUES "]
            for prefix in prefixL:
                tup = ssDetails.partition(prefix)
                if tup[1] == prefix:
                    ff = tup[2].split(" ")
                    if len(ff) == 2:
                        # BINDING SITE FOR LINKED RESIDUES A 502-507
                        try:
                            tff = ff[1].split("-")
                            authAsymIdA = ff[0]
                            authSeqIdA = tff[0]
                            authSeqIdB = tff[1]
                        except Exception:
                            continue
                    if len(ff) == 4 and ff[2].lower() == "to":
                        authAsymIdA = ff[0]
                        authSeqIdA = ff[1]
                        authSeqIdB = ff[3]
                    elif len(ff) == 4 and ff[2].lower() != "to":
                        authAsymIdA = ff[0]
                        authSeqIdA = ff[1]
                        authSeqIdB = ff[3]
                    elif len(ff) > 4:
                        authAsymIdA = ff[0]
                        authSeqIdA = ff[1]
                        authSeqIdB = ff[-1]
                    else:
                        continue
                    retL.append((authAsymIdA, None, authSeqIdA, ssDetails))
                    retL.append((authAsymIdA, None, authSeqIdB, ssDetails))
                    return retL

            #
            #
            prefixL = ["BINDING SITE FOR CHAIN ", "binding site for chain "]
            for prefix in prefixL:
                tup = ssDetails.partition(prefix)
                if tup[1] == prefix:
                    ff = tup[2].split(" ")
                    authAsymId = ff[0]
                    retL.append((authAsymId, None, None, ssDetails))
                    return retL
            # punt -
            retL.append((None, None, None, ssDetails))
            return retL
        except Exception as e:
            logger.exception("Failing with %s for %r", str(e), ssDetails)
        return [(None, None, None, ssDetails)]

    def getUnobservedPolymerResidueInfo(self, dataContainer):
        """Return a dictionary of unobserved regions of polymer instances.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(modelId, asymId, occFlag): [seqId range list], ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchUnobservedInfo(dataContainer)
        return wD["polyResRng"] if "polyResRng" in wD else {}

    def getUnobservedPolymerAtomInfo(self, dataContainer):
        """Return a dictionary of polymer regions containing unobserved atoms.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(modelId, asymId, occFlag): [seqId range list], ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchUnobservedInfo(dataContainer)
        return wD["polyAtomRng"] if "polyAtomRng" in wD else {}

    def getUnobservedNonPolymerAtomInfo(self, dataContainer):
        """Return a dictionary of nonpolymer instances containing unobserved atoms (std nomenclature).

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(modelId, compId, asymId, occFlag): [atomId, .. ], ...}

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchUnobservedInfo(dataContainer)
        return wD["nonPolyMissingAtomD"] if "nonPolyMissingAtomD" in wD else {}

    def getUnobservedNonPolymerAtomInfoAuth(self, dataContainer):
        """Return a dictionary of nonpolymer instances containing unobserved atoms (auth nomenclature)

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(modelId, compId, authtAsymId, authSeqIdm, occFlag): [atomId, .. ], ...}

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchUnobservedInfo(dataContainer)
        return wD["nonPolyMissingAtomAuthD"] if "nonPolyMissingAtomAuthD" in wD else {}

    def __fetchUnobservedInfo(self, dataContainer):
        wD = self.__instanceUnobservedCache.get(dataContainer.getName())
        if not wD:
            wD = self.__getUnobserved(dataContainer)
            self.__instanceUnobservedCache.set(dataContainer.getName(), wD)
        return wD

    def __getUnobserved(self, dataContainer):
        """Internal method to extract unobserved and zero occupancy features.

        Args:
            dataContainer ([type]): [description]

        Returns:
            {"polyResRng":  {(modelId, asymId, occFlag): [seqId range list], ...},
             "polyAtomRng": {(modelId, asymId, occFlag): [seqId range list], ...},
             "nonPolyMissingAtomD": {(modelId, compId, asymId, zeroOccFlag): [atomId,...], },
             "nonPolyMissingAtomAuthD": {(modelId, compId, authAsymId, authSeqId, zeroOccFlag): [atomId,...], },
             }

            occFlag = 0 -> zero occupancy
        Example:

                loop_
                _pdbx_unobs_or_zero_occ_atoms.id
                _pdbx_unobs_or_zero_occ_atoms.PDB_model_num
                _pdbx_unobs_or_zero_occ_atoms.polymer_flag
                _pdbx_unobs_or_zero_occ_atoms.occupancy_flag
                _pdbx_unobs_or_zero_occ_atoms.auth_asym_id
                _pdbx_unobs_or_zero_occ_atoms.auth_comp_id
                _pdbx_unobs_or_zero_occ_atoms.auth_seq_id
                _pdbx_unobs_or_zero_occ_atoms.PDB_ins_code
                _pdbx_unobs_or_zero_occ_atoms.auth_atom_id
                _pdbx_unobs_or_zero_occ_atoms.label_alt_id
                _pdbx_unobs_or_zero_occ_atoms.label_asym_id
                _pdbx_unobs_or_zero_occ_atoms.label_comp_id
                _pdbx_unobs_or_zero_occ_atoms.label_seq_id
                _pdbx_unobs_or_zero_occ_atoms.label_atom_id
                1  1 Y 1 B ARG 17  ? NE    ? B ARG 17 NE
                2  1 Y 1 B ARG 17  ? CZ    ? B ARG 17 CZ
                3  1 Y 1 B ARG 17  ? NH1   ? B ARG 17 NH1

                #
                loop_
                _pdbx_unobs_or_zero_occ_residues.id
                _pdbx_unobs_or_zero_occ_residues.PDB_model_num
                _pdbx_unobs_or_zero_occ_residues.polymer_flag
                _pdbx_unobs_or_zero_occ_residues.occupancy_flag
                _pdbx_unobs_or_zero_occ_residues.auth_asym_id
                _pdbx_unobs_or_zero_occ_residues.auth_comp_id
                _pdbx_unobs_or_zero_occ_residues.auth_seq_id
                _pdbx_unobs_or_zero_occ_residues.PDB_ins_code
                _pdbx_unobs_or_zero_occ_residues.label_asym_id
                _pdbx_unobs_or_zero_occ_residues.label_comp_id
                _pdbx_unobs_or_zero_occ_residues.label_seq_id
                1  1 Y 1 A MET 1 ? A MET 1
                2  1 Y 1 A ALA 2 ? A ALA 2
                3  1 Y 1 A LYS 3 ? A LYS 3
        """
        logger.debug("Starting with %r", dataContainer.getName())
        #
        rD = {}
        try:
            # Exit if source categories are missing
            if not (dataContainer.exists("pdbx_unobs_or_zero_occ_residues") or dataContainer.exists("pdbx_unobs_or_zero_occ_atoms")):
                return rD
            # ------- --------- ------- --------- ------- --------- ------- --------- ------- ---------
            resObj = None
            if dataContainer.exists("pdbx_unobs_or_zero_occ_residues"):
                resObj = dataContainer.getObj("pdbx_unobs_or_zero_occ_residues")
            #
            atomObj = None
            if dataContainer.exists("pdbx_unobs_or_zero_occ_atoms"):
                atomObj = dataContainer.getObj("pdbx_unobs_or_zero_occ_atoms")
            #
            # Get representative model
            repModelId = self.getRepresentativeModelId(dataContainer)
            #
            polyResRngD = {}
            if resObj:
                for ii in range(resObj.getRowCount()):
                    modelId = resObj.getValueOrDefault("PDB_model_num", ii, defaultValue=None)
                    if str(modelId) != repModelId:  # Skip non-representative models
                        continue
                    pFlag = resObj.getValueOrDefault("polymer_flag", ii, defaultValue=None)
                    if pFlag == "Y":
                        occFlag = resObj.getValueOrDefault("occupancy_flag", ii, defaultValue=None)
                        zeroOccFlag = int(occFlag) == 0
                        asymId = resObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                        # authAsymId = resObj.getValueOrDefault("auth_asym_id", ii, defaultValue=None)
                        seqId = resObj.getValueOrDefault("label_seq_id", ii, defaultValue=None)
                        if seqId:
                            polyResRngD.setdefault((modelId, asymId, zeroOccFlag), []).append(int(seqId))
                #
                cloneD = copy.deepcopy(polyResRngD)
                for tup in cloneD:
                    polyResRngD[tup] = list(self.__toRangeList(cloneD[tup]))
                logger.debug("polyResRngD %r", polyResRngD)
            #
            polyAtomRngD = {}
            nonPolyMissingAtomD = {}
            nonPolyMissingAtomAuthD = {}
            if atomObj:
                for ii in range(atomObj.getRowCount()):
                    modelId = atomObj.getValueOrDefault("PDB_model_num", ii, defaultValue=None)
                    if str(modelId) != repModelId:  # Skip non-representative models
                        continue
                    pFlag = atomObj.getValueOrDefault("polymer_flag", ii, defaultValue=None)
                    occFlag = atomObj.getValueOrDefault("occupancy_flag", ii, defaultValue=None)
                    zeroOccFlag = occFlag and int(occFlag) == 0
                    asymId = atomObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                    if pFlag == "Y":
                        # authAsymId = resObj.getValueOrDefault("auth_asym_id", ii, defaultValue=None)
                        seqId = atomObj.getValueOrDefault("label_seq_id", ii, defaultValue=None)
                        if seqId:
                            polyAtomRngD.setdefault((modelId, asymId, zeroOccFlag), []).append(int(seqId))
                    else:
                        authAsymId = atomObj.getValueOrDefault("auth_asym_id", ii, defaultValue=None)
                        authSeqId = atomObj.getValueOrDefault("auth_seq_id", ii, defaultValue=None)
                        atomId = atomObj.getValueOrDefault("label_atom_id", ii, defaultValue=None)
                        compId = atomObj.getValueOrDefault("label_comp_id", ii, defaultValue=None)
                        nonPolyMissingAtomD.setdefault((modelId, compId, asymId, zeroOccFlag), []).append(atomId)
                        nonPolyMissingAtomAuthD.setdefault((modelId, compId, authAsymId, authSeqId, zeroOccFlag), []).append(atomId)
                #
                cloneD = copy.deepcopy(polyAtomRngD)
                for tup in cloneD:
                    polyAtomRngD[tup] = list(self.__toRangeList(cloneD[tup]))
                logger.debug("polyAtomRngD %r", polyAtomRngD)
            #
            rD = {"polyResRng": polyResRngD, "polyAtomRng": polyAtomRngD, "nonPolyMissingAtomD": nonPolyMissingAtomD, "nonPolyMissingAtomAuthD": nonPolyMissingAtomAuthD}
        except Exception as e:
            logger.exception("%s failing with %s", dataContainer.getName(), str(e))
        return rD

    def getInstanceModelOutlierInfo(self, dataContainer):
        """Return a dictionary of polymer model outliers.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(modelId, asymId): (seqId,compId), ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchInstanceModelOutliers(dataContainer)
        return wD["instanceModelOutlierD"] if "instanceModelOutlierD" in wD else {}

    def getInstanceNonpolymerValidationInfo(self, dataContainer):
        """Return a dictionary of nonpolymer validation details.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(modelId, asymId): NonpolymerValidationInstance(rsr, rsrCc, bondsRmsZ, anglesRmsZ,
                                             intermolecular_clashes, mogul_bond_outliers, mogul_angle_outliers, stereo_outliers)}

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchInstanceModelOutliers(dataContainer)
        return wD["instanceModelValidationD"] if "instanceModelValidationD" in wD else {}

    def __fetchInstanceModelOutliers(self, dataContainer):
        wD = self.__modelOutliersCache.get(dataContainer.getName())
        if not wD:
            wD = self.__getInstanceModelOutliers(dataContainer)
            self.__modelOutliersCache.set(dataContainer.getName(), wD)
        return wD

    # TO BE DELETED
    def __getInstanceModelOutliersXML(self, dataContainer):
        """Internal method to assemble model outliers details.

        Args:
            dataContainer ([type]): [description]

        Returns:
            {"instanceModelOutlierD": {(modelId, asymId): [(compId, seqId, "BOND_OUTLIER", optional_description), ...}}
        #
            loop_
            _pdbx_vrpt_instance_results.ordinal
            _pdbx_vrpt_instance_results.entity_id
            _pdbx_vrpt_instance_results.auth_asym_id
            _pdbx_vrpt_instance_results.label_asym_id
            _pdbx_vrpt_instance_results.label_comp_id
            _pdbx_vrpt_instance_results.auth_seq_id
            _pdbx_vrpt_instance_results.label_seq_id
            _pdbx_vrpt_instance_results.PDB_ins_code
            _pdbx_vrpt_instance_results.label_alt_id
            _pdbx_vrpt_instance_results.PDB_model_num
            _pdbx_vrpt_instance_results.num_H_reduce
            _pdbx_vrpt_instance_results.cis_peptide
            _pdbx_vrpt_instance_results.natoms_eds
            _pdbx_vrpt_instance_results.RSR
            _pdbx_vrpt_instance_results.RSRCC
            _pdbx_vrpt_instance_results.RSRZ
            _pdbx_vrpt_instance_results.OWAB
            _pdbx_vrpt_instance_results.average_occupancy
            _pdbx_vrpt_instance_results.ramachandran_class
            _pdbx_vrpt_instance_results.rotamer_class
            _pdbx_vrpt_instance_results.phi
            _pdbx_vrpt_instance_results.psi
            _pdbx_vrpt_instance_results.mogul_angles_RMSZ
            _pdbx_vrpt_instance_results.mogul_bonds_RMSZ
            _pdbx_vrpt_instance_results.mogul_RMSZ_num_angles
            _pdbx_vrpt_instance_results.mogul_RMSZ_num_bonds
            # ...
            302 1 A A TYR 340 343 ? ? 1 9  ?   12 0.108 0.943 0.117  71.350  1.000 Favored m-85    -111.8 6.4    ?    ?    ?  ?
            303 1 A A LYS 341 344 ? ? 1 13 ?   9  0.120 0.955 -0.380 67.860  1.000 Favored mttt    -73.3  139.6  ?    ?    ?  ?
            304 1 A A ILE 342 345 ? ? 1 11 ?   8  0.147 0.964 0.799  76.030  1.000 Favored pt      -140.0 171.7  ?    ?    ?  ?
            305 1 A A ASN 343 346 ? ? 1 6  ?   8  0.182 0.948 1.114  82.730  1.000 Favored m-80    52.8   49.6   ?    ?    ?  ?
            306 1 A A GLN 344 347 ? ? 1 2  ?   5  0.193 0.807 1.002  97.730  1.000 ?       ?       ?      ?      ?    ?    ?  ?
            # ...
            307 2 A B PEG 401 .   ? A 1 10 ?   14 0.154 0.914 ?      36.150  1.000 ?       ?       ?      ?      0.76 0.64 5  6
            308 2 A B PEG 401 .   ? B 1 10 ?   14 0.154 0.914 ?      36.150  1.000 ?       ?       ?      ?      0.97 0.68 5  6
            309 3 A C HYO 402 .   ? ? 1 ?  ?   21 0.108 0.947 ?      35.530  1.000 ?       ?       ?      ?      2.18 4.96 32 23
            310 4 A D NI  403 .   ? ? 1 ?  ?   1  0.096 0.999 ?      28.080  1.000 ?       ?       ?      ?      ?    ?    ?  ?
            311 5 A E OGA 404 .   ? ? 1 3  ?   10 0.104 0.976 ?      30.510  1.000 ?       ?       ?      ?      1.87 3.23 4  3
            312 6 A F EDO 405 .   ? ? 1 6  ?   4  0.097 0.941 ?      42.000  1.000 ?       ?       ?      ?      0.32 0.80 2  3
            313 6 A G EDO 406 .   ? ? 1 6  ?   4  0.252 0.797 ?      57.320  1.000 ?       ?       ?      ?      0.73 0.61 2  3
            314 7 A H SR  407 .   ? ? 1 ?  ?   1  0.143 1.000 ?      30.560  0.840 ?       ?       ?      ?      ?    ?    ?  ?
            315 8 A I UNX 408 .   ? ? 1 ?  ?   1  0.321 0.940 ?      41.340  1.000 ?       ?       ?      ?      ?    ?    ?  ?
            316 8 A J UNX 409 .   ? ? 1 ?  ?   1  0.611 0.922 ?      61.040  1.000 ?       ?       ?      ?      ?    ?    ?  ?
            # ...
        """
        logger.debug("Starting with %r", dataContainer.getName())
        #
        rD = {}
        try:
            # Exit if no source categories are present
            if not (
                dataContainer.exists("pdbx_vrpt_instance_results")
                or dataContainer.exists("pdbx_vrpt_bond_outliers")
                or dataContainer.exists("pdbx_vrpt_angle_outliers")
                or dataContainer.exists("pdbx_vrpt_mogul_bond_outliers")
                or dataContainer.exists("pdbx_vrpt_mogul_angle_outliers")
                or dataContainer.exists("pdbx_vrpt_stereo_outliers")
                or dataContainer.exists("pdbx_vrpt_clashes")
            ):
                return rD
            # ------- --------- ------- --------- ------- --------- ------- --------- ------- ---------
            instanceTypeD = self.getInstanceTypes(dataContainer)
            #
            instanceModelOutlierD = {}
            instanceModelValidationD = {}
            #
            npMogulBondOutlierD = defaultdict(int)
            npMogulAngleOutlierD = defaultdict(int)
            npStereoOutlierD = defaultdict(int)
            #
            vObj = None
            if dataContainer.exists("pdbx_vrpt_bond_outliers"):
                vObj = dataContainer.getObj("pdbx_vrpt_bond_outliers")
            if vObj:
                for ii in range(vObj.getRowCount()):
                    seqId = vObj.getValueOrDefault("label_seq_id", ii, defaultValue=None)
                    if seqId:
                        modelId = vObj.getValueOrDefault("PDB_model_num", ii, defaultValue=None)
                        asymId = vObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                        compId = vObj.getValueOrDefault("label_comp_id", ii, defaultValue=None)
                        altId = vObj.getValueOrDefault("label_alt_id", ii, defaultValue=None)
                        #
                        atomI = vObj.getValueOrDefault("atom0", ii, defaultValue=None)
                        atomJ = vObj.getValueOrDefault("atom1", ii, defaultValue=None)
                        obsDist = vObj.getValueOrDefault("obs", ii, defaultValue=None)
                        zVal = vObj.getValueOrDefault("Z", ii, defaultValue=None)
                        tS = "%s-%s (altId=%s) dist=%s Z=%s" % (atomI, atomJ, altId, obsDist, zVal) if altId else "%s-%s dist=%s Z=%s" % (atomI, atomJ, obsDist, zVal)
                        #
                        instanceModelOutlierD.setdefault((modelId, asymId, altId, True), []).append(
                            OutlierValue(
                                compId,
                                int(seqId),
                                "BOND_OUTLIER",
                                tS,
                            )
                        )

                #
                logger.debug("length instanceModelOutlierD %d", len(instanceModelOutlierD))
            # ----
            vObj = None
            if dataContainer.exists("pdbx_vrpt_angle_outliers"):
                vObj = dataContainer.getObj("pdbx_vrpt_angle_outliers")
            if vObj:
                for ii in range(vObj.getRowCount()):
                    seqId = vObj.getValueOrDefault("label_seq_id", ii, defaultValue=None)
                    if seqId:
                        modelId = vObj.getValueOrDefault("PDB_model_num", ii, defaultValue=None)
                        asymId = vObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                        compId = vObj.getValueOrDefault("label_comp_id", ii, defaultValue=None)
                        altId = vObj.getValueOrDefault("label_alt_id", ii, defaultValue=None)
                        #
                        atomI = vObj.getValueOrDefault("atom0", ii, defaultValue=None)
                        atomJ = vObj.getValueOrDefault("atom1", ii, defaultValue=None)
                        atomK = vObj.getValueOrDefault("atom2", ii, defaultValue=None)
                        obsDist = vObj.getValueOrDefault("obs", ii, defaultValue=None)
                        zVal = vObj.getValueOrDefault("Z", ii, defaultValue=None)
                        tS = (
                            "%s-%s-%s (altId %s) angle=%s Z=%s" % (atomI, atomJ, atomK, altId, obsDist, zVal)
                            if altId
                            else "%s-%s-%s angle=%s Z=%s" % (atomI, atomJ, atomK, obsDist, zVal)
                        )
                        #
                        instanceModelOutlierD.setdefault((modelId, asymId, altId, True), []).append(
                            OutlierValue(
                                compId,
                                int(seqId),
                                "ANGLE_OUTLIER",
                                tS,
                            )
                        )

                #
                logger.debug("length instanceModelOutlierD %d", len(instanceModelOutlierD))
            # ----
            vObj = None
            if dataContainer.exists("pdbx_vrpt_mogul_bond_outliers"):
                vObj = dataContainer.getObj("pdbx_vrpt_mogul_bond_outliers")
            if vObj:
                for ii in range(vObj.getRowCount()):
                    seqId = vObj.getValueOrDefault("label_seq_id", ii, defaultValue=None)

                    modelId = vObj.getValueOrDefault("PDB_model_num", ii, defaultValue=None)
                    asymId = vObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                    compId = vObj.getValueOrDefault("label_comp_id", ii, defaultValue=None)
                    altId = vObj.getValueOrDefault("label_alt_id", ii, defaultValue=None)
                    #
                    atoms = vObj.getValueOrDefault("atoms", ii, defaultValue=None)
                    obsDist = vObj.getValueOrDefault("obsval", ii, defaultValue=None)
                    meanValue = vObj.getValueOrDefault("mean", ii, defaultValue=None)
                    zVal = vObj.getValueOrDefault("Zscore", ii, defaultValue=None)
                    tS = "%s (altIt %s) angle=%s Z=%s" % (atoms, altId, obsDist, zVal) if altId else "%s angle=%s Z=%s" % (atoms, obsDist, zVal)
                    # OutlierValue = collections.namedtuple("OutlierValue", "compId, seqId, outlierType, description, reported, reference, uncertaintyValue, uncertaintyType")
                    if seqId:
                        instanceModelOutlierD.setdefault((modelId, asymId, altId, True), []).append(
                            OutlierValue(
                                compId,
                                int(seqId),
                                "MOGUL_BOND_OUTLIER",
                                tS,
                            )
                        )
                    else:
                        instanceModelOutlierD.setdefault((modelId, asymId, altId, False), []).append(
                            OutlierValue(compId, None, "MOGUL_BOND_OUTLIER", tS, obsDist, meanValue, zVal, "Z-Score")
                        )
                        npMogulBondOutlierD[(modelId, asymId, altId, compId)] += 1
                #
                logger.debug("length instanceModelOutlierD %d", len(instanceModelOutlierD))

            vObj = None
            if dataContainer.exists("pdbx_vrpt_mogul_angle_outliers"):
                vObj = dataContainer.getObj("pdbx_vrpt_mogul_angle_outliers")
            if vObj:
                for ii in range(vObj.getRowCount()):
                    seqId = vObj.getValueOrDefault("label_seq_id", ii, defaultValue=None)

                    modelId = vObj.getValueOrDefault("PDB_model_num", ii, defaultValue=None)
                    asymId = vObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                    compId = vObj.getValueOrDefault("label_comp_id", ii, defaultValue=None)
                    altId = vObj.getValueOrDefault("label_alt_id", ii, defaultValue=None)
                    #
                    atoms = vObj.getValueOrDefault("atoms", ii, defaultValue=None)
                    obsDist = vObj.getValueOrDefault("obsval", ii, defaultValue=None)
                    meanValue = vObj.getValueOrDefault("mean", ii, defaultValue=None)
                    zVal = vObj.getValueOrDefault("Zscore", ii, defaultValue=None)
                    tS = "%s (altId %s) angle=%s Z=%s" % (atoms, altId, obsDist, zVal) if altId else "%s angle=%s Z=%s" % (atoms, obsDist, zVal)
                    if seqId:
                        instanceModelOutlierD.setdefault((modelId, asymId, altId, True), []).append(
                            OutlierValue(
                                compId,
                                int(seqId),
                                "MOGUL_ANGLE_OUTLIER",
                                tS,
                            )
                        )
                    else:
                        instanceModelOutlierD.setdefault((modelId, asymId, altId, False), []).append(
                            OutlierValue(compId, None, "MOGUL_ANGLE_OUTLIER", tS, obsDist, meanValue, zVal, "Z-Score")
                        )
                        npMogulAngleOutlierD[(modelId, asymId, altId, compId)] += 1
                logger.debug("length instanceModelOutlierD %d", len(instanceModelOutlierD))
                #
            # --
            vObj = None
            if dataContainer.exists("pdbx_vrpt_stereo_outliers"):
                vObj = dataContainer.getObj("pdbx_vrpt_stereo_outliers")
            if vObj:
                for ii in range(vObj.getRowCount()):
                    seqId = vObj.getValueOrDefault("label_seq_id", ii, defaultValue=None)
                    modelId = vObj.getValueOrDefault("PDB_model_num", ii, defaultValue=None)
                    asymId = vObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                    compId = vObj.getValueOrDefault("label_comp_id", ii, defaultValue=None)
                    altId = vObj.getValueOrDefault("label_alt_id", ii, defaultValue=None)
                    description = vObj.getValueOrDefault("problem", ii, defaultValue=None)
                    #
                    if seqId:
                        instanceModelOutlierD.setdefault((modelId, asymId, altId, True), []).append(
                            OutlierValue(
                                compId,
                                int(seqId),
                                "STEREO_OUTLIER",
                                description,
                            )
                        )
                    else:
                        instanceModelOutlierD.setdefault((modelId, asymId, altId, False), []).append(OutlierValue(compId, None, "STEREO_OUTLIER", description))
                        npStereoOutlierD[(modelId, asymId, altId, compId)] += 1
                logger.debug("length instanceModelOutlierD %d", len(instanceModelOutlierD))
                #
                #
            # ----  Capture/evaluate non-polymer intermolecular clashes ... here filter internal molecule clashes ...
            instanceTypeD = self.getInstanceTypes(dataContainer)
            npClashD = defaultdict(int)
            tClashD = {}
            vObj = None
            if dataContainer.exists("pdbx_vrpt_clashes"):
                vObj = dataContainer.getObj("pdbx_vrpt_clashes")
            if vObj:
                logger.debug("Row count for %s: %d", vObj.getName(), vObj.getRowCount())
                for ii in range(vObj.getRowCount()):
                    seqId = vObj.getValueOrDefault("label_seq_id", ii, defaultValue=None)
                    asymId = vObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                    if asymId in instanceTypeD and instanceTypeD[asymId] == "non-polymer":
                        clashId = vObj.getValueOrDefault("cid", ii, defaultValue=None)
                        modelId = vObj.getValueOrDefault("PDB_model_num", ii, defaultValue=None)
                        asymId = vObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                        compId = vObj.getValueOrDefault("label_comp_id", ii, defaultValue=None)
                        altId = vObj.getValueOrDefault("label_alt_id", ii, defaultValue=None)
                        tClashD.setdefault((modelId, asymId, altId, compId), []).append(clashId)
                    #
                for ky, clashIdL in tClashD.items():
                    cD = defaultdict(int)
                    for clashId in clashIdL:
                        cD[clashId] += 1
                    for clashId, clashCount in cD.items():
                        if clashCount == 1:
                            npClashD[ky] += 1
            #
            logger.debug("%s npClashD %r", dataContainer.getName(), npClashD)
            # ----
            vObj = None
            if dataContainer.exists("pdbx_vrpt_instance_results"):
                vObj = dataContainer.getObj("pdbx_vrpt_instance_results")

            if vObj:
                logger.debug("Row count for %s: %d", vObj.getName(), vObj.getRowCount())
                for ii in range(vObj.getRowCount()):
                    seqId = vObj.getValueOrDefault("label_seq_id", ii, defaultValue=None)
                    modelId = vObj.getValueOrDefault("PDB_model_num", ii, defaultValue=None)
                    asymId = vObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                    compId = vObj.getValueOrDefault("label_comp_id", ii, defaultValue=None)
                    altId = vObj.getValueOrDefault("label_alt_id", ii, defaultValue=None)
                    #
                    rotamerClass = vObj.getValueOrDefault("rotamer_class", ii, defaultValue=None)
                    ramaClass = vObj.getValueOrDefault("ramachandran_class", ii, defaultValue=None)
                    rsr = vObj.getValueOrDefault("RSR", ii, defaultValue=None)
                    rsrZ = vObj.getValueOrDefault("RSRZ", ii, defaultValue=None)
                    rsrCc = vObj.getValueOrDefault("RSRCC", ii, defaultValue=None)
                    #
                    anglesRmsZ = vObj.getValueOrDefault("mogul_angles_RMSZ", ii, defaultValue=None)
                    bondsRmsZ = vObj.getValueOrDefault("mogul_bonds_RMSZ", ii, defaultValue=None)
                    # ---

                    # ---
                    if seqId:
                        if rotamerClass and rotamerClass.upper() == "OUTLIER":
                            instanceModelOutlierD.setdefault((modelId, asymId, altId, True), []).append(
                                OutlierValue(
                                    compId,
                                    int(seqId),
                                    "ROTAMER_OUTLIER",
                                    None,
                                )
                            )
                        if ramaClass and ramaClass.upper() == "OUTLIER":
                            instanceModelOutlierD.setdefault((modelId, asymId, altId, True), []).append(
                                OutlierValue(
                                    compId,
                                    int(seqId),
                                    "RAMACHANDRAN_OUTLIER",
                                    None,
                                )
                            )
                        if rsrZ and float(rsrZ) > 2.0:
                            tS = "%s > 2.0 (altId %s)" % (rsrZ, altId) if altId else "%s > 2.0 " % rsrZ
                            instanceModelOutlierD.setdefault((modelId, asymId, altId, True), []).append(
                                OutlierValue(
                                    compId,
                                    int(seqId),
                                    "RSRZ_OUTLIER",
                                    tS,
                                )
                            )
                        if rsrCc and float(rsrCc) < 0.650:
                            tS = "RSCC < 0.65 (altId %s)" % altId if altId else "RSCC < 0.65"
                            instanceModelOutlierD.setdefault((modelId, asymId, altId, True), []).append(
                                OutlierValue(
                                    compId,
                                    int(seqId),
                                    "RSCC_OUTLIER",
                                    tS,
                                )
                            )
                    else:
                        if rsrZ and float(rsrZ) > 2.0:
                            tS = "%s > 2.0 (altId %s)" % (rsrZ, altId) if altId else "%s > 2.0" % rsrZ
                            instanceModelOutlierD.setdefault((modelId, asymId, altId, False), []).append(OutlierValue(compId, None, "RSRZ_OUTLIER", tS, rsr, None, rsrZ, "Z-Score"))
                        if rsrCc and float(rsrCc) < 0.650:
                            tS = "RSCC < 0.65 (altId %s)" % altId if altId else "RSCC < 0.65"
                            instanceModelOutlierD.setdefault((modelId, asymId, altId, False), []).append(OutlierValue(compId, None, "RSCC_OUTLIER", tS, rsrCc))
                        if asymId in instanceTypeD and instanceTypeD[asymId] == "non-polymer":
                            instanceModelValidationD[(modelId, asymId, altId, compId)] = NonpolymerValidationInstance(
                                float(rsr) if rsr else None,
                                float(rsrCc) if rsrCc else None,
                                float(bondsRmsZ) if bondsRmsZ else None,
                                float(anglesRmsZ) if anglesRmsZ else None,
                                npClashD[(modelId, asymId, altId, compId)] if (modelId, asymId, altId, compId) in npClashD else 0,
                                npMogulBondOutlierD[(modelId, asymId, altId, compId)] if (modelId, asymId, altId, compId) in npMogulBondOutlierD else 0,
                                npMogulAngleOutlierD[(modelId, asymId, altId, compId)] if (modelId, asymId, altId, compId) in npMogulAngleOutlierD else 0,
                                npStereoOutlierD[(modelId, asymId, altId, compId)] if (modelId, asymId, altId, compId) in npStereoOutlierD else 0,
                            )
                #
            logger.debug("instanceModelOutlierD %r", instanceModelOutlierD)
            logger.debug("instanceModelValidationD %r", instanceModelValidationD)

            rD = {"instanceModelOutlierD": instanceModelOutlierD, "instanceModelValidationD": instanceModelValidationD}
        except Exception as e:
            logger.exception("%s failing with %s", dataContainer.getName(), str(e))
        return rD

    # TO BE DELETED -- Needs additional testing before deletion to ensure all downstream dependencies are addressed
    def __getInstanceModelOutliers(self, dataContainer):
        """Internal method to assemble model outliers details.

        Args:
            dataContainer ([type]): [description]

        Returns:
            {"instanceModelOutlierD": {(modelId, asymId): [(compId, seqId, "BOND_OUTLIER", optional_description), ...}}

        """
        logger.debug("Starting with %r", dataContainer.getName())
        #
        rD = {}
        try:
            # Exit if no source categories are present
            if not (
                    dataContainer.exists("pdbx_vrpt_model_instance")
                    or dataContainer.exists("pdbx_vrpt_model_instance_geometry")
                    or dataContainer.exists("pdbx_vrpt_model_instance_density")
                    or dataContainer.exists("pdbx_vrpt_model_instance_map_fitting")
            ) or not dataContainer.exists("pdbx_vrpt_model_instance"):
                return rD
            # ------- --------- ------- --------- ------- --------- ------- --------- ------- ---------
            instanceModelOutlierD = {}
            instanceModelValidationD = {}
            #
            npMogulBondOutlierD = {}
            npMogulAngleOutlierD = {}
            npStereoOutlierD = {}
            npClashD = {}
            #
            iObj = dataContainer.getObj("pdbx_vrpt_model_instance")
            if not iObj:
                return rD

            repModelId = self.getRepresentativeModelId(dataContainer)

            instanceTypeD = self.getInstanceTypes(dataContainer)
            npAsymL = [k for k, v in instanceTypeD.items() if v == "non-polymer"]
            if not npAsymL:
                return rD
            cndL2 = [("label_asym_id", "in", npAsymL), ("PDB_model_num", "eq", repModelId)]
            kL = iObj.selectIndicesWhereOpConditions(cndL2)
            instL = []
            cD = {}
            for ii in kL:
                modelId = iObj.getValueOrDefault("PDB_model_num", ii, defaultValue=None)
                if str(modelId) != repModelId:  # Skip non-representative models
                    continue
                asymId = iObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                if asymId in npAsymL:
                    instId = iObj.getValueOrDefault("id", ii, defaultValue=None)
                    instL.append(instId)
                    altId = iObj.getValueOrDefault("label_alt_id", ii, defaultValue=None)
                    compId = iObj.getValueOrDefault("label_comp_id", ii, defaultValue=None)
                    seqId = iObj.getValueOrDefault("label_seq_id", ii, defaultValue=None)
                    cD[instId] = [modelId, asymId, compId, altId, seqId]
                    countClashes = iObj.getValueOrDefault("count_clashes", ii, defaultValue=None)
                    countMogulAngleOutliers = iObj.getValueOrDefault("count_mogul_angle_outliers", ii, defaultValue=None)
                    countMogulBondOutliers = iObj.getValueOrDefault("count_mogul_bond_outliers", ii, defaultValue=None)
                    countStereoOutliers = iObj.getValueOrDefault("count_chiral_outliers", ii, defaultValue=None)
                    npMogulBondOutlierD[(modelId, asymId, altId, compId)] = countMogulBondOutliers if countMogulBondOutliers else 0
                    npMogulAngleOutlierD[(modelId, asymId, altId, compId)] = countMogulAngleOutliers if countMogulAngleOutliers else 0
                    npClashD[(modelId, asymId, altId, compId)] = countClashes if countClashes else 0
                    npStereoOutlierD[(modelId, asymId, altId, compId)] = countStereoOutliers if countStereoOutliers else 0
            cndL3 = [("instance_id", "in", instL)]
            # ----
            vObj = None
            if dataContainer.exists("pdbx_vrpt_model_instance_map_fitting"):
                vObj = dataContainer.getObj("pdbx_vrpt_model_instance_map_fitting")
            elif dataContainer.exists("pdbx_vrpt_model_instance_density"):
                vObj = dataContainer.getObj("pdbx_vrpt_model_instance_density")

            vD = {}
            if vObj:
                kL = vObj.selectIndicesWhereOpConditions(cndL3)
                for ii in kL:
                    instId = vObj.getValueOrDefault("instance_id", ii, defaultValue=None)
                    rsrCc = vObj.getValueOrDefault("RSRCC", ii, defaultValue=None)
                    rsr = vObj.getValueOrDefault("RSR", ii, defaultValue=None)
                    rsrZ = vObj.getValueOrDefault("RSRZ", ii, defaultValue=None)
                    nAtomsEds = vObj.getValueOrDefault("natoms_eds", ii, defaultValue=None)
                    vD[instId] = [rsrCc, rsr, rsrZ, nAtomsEds]

            gObj = None
            if dataContainer.exists("pdbx_vrpt_model_instance_geometry"):
                gObj = dataContainer.getObj("pdbx_vrpt_model_instance_geometry")
            if gObj:
                logger.debug("Row count for %s: %d", gObj.getName(), gObj.getRowCount())
                kL = gObj.selectIndicesWhereOpConditions(cndL3)
                for ii in kL:
                    instId = gObj.getValueOrDefault("instance_id", ii, defaultValue=None)
                    # [fL] = iObj.selectValueListWhere(["PDB_model_num", "label_asym_id", "label_comp_id", "label_alt_id", "label_seq_id"], instId, "id")
                    # [modelId, asymId, compId, altId, seqId] = [x if not (x in [".", "?", "None"]) else None for x in fL]
                    [modelId, asymId, compId, altId, seqId] = cD[instId]
                    if str(modelId) != repModelId:  # Skip non-representative models
                        continue
                    # ---
                    if not seqId:
                        # Get the matching data for non-polymers from pdbx_vrpt_model_instance_density or pdbx_vrpt_model_instance_map_fitting
                        rsr = None
                        rsrZ = None
                        rsrCc = None
                        nAtomsEds = None
                        if vObj:
                            [rsrCc, rsr, rsrZ, nAtomsEds] = vD[instId]
                            # TO DELETE
                            # iiL = vObj.selectIndices(instId, "instance_id")
                            # if len(iiL) == 1:
                            #   rsr = vObj.getValueOrDefault("RSR", iiL[0], defaultValue=None)
                            #   rsrZ = vObj.getValueOrDefault("RSRZ", iiL[0], defaultValue=None)
                            #   rsrCc = vObj.getValueOrDefault("RSRCC", iiL[0], defaultValue=None)
                            #   nAtomsEds = vObj.getValueOrDefault("natoms_eds", iiL[0], defaultValue=None)
                        # Only need mogul values from pdbx_vrpt_model_instance_geometry for non-polymers
                        anglesRmsZ = None
                        bondsRmsZ = None
                        numAnglesRmsZ = None
                        numBondsRmsZ = None
                        avgOccupancy = None
                        software = gObj.getValueOrDefault("program_for_bond_angle_geometry", ii, defaultValue=None)
                        if str(software).lower() == "mogul":
                            anglesRmsZ = gObj.getValueOrDefault("angles_RMSZ", ii, defaultValue=None)
                            bondsRmsZ = gObj.getValueOrDefault("bonds_RMSZ", ii, defaultValue=None)
                            numAnglesRmsZ = gObj.getValueOrDefault("num_angles_RMSZ", ii, defaultValue=None)
                            numBondsRmsZ = gObj.getValueOrDefault("num_bonds_RMSZ", ii, defaultValue=None)
                            avgOccupancy = gObj.getValueOrDefault("average_occupancy", ii, defaultValue=None)
                        # Evaluate RSRZ_OUTLIER and RSCC_OUTLIER for non-polymer annotations
                        if rsrZ and float(rsrZ) > 2.0:
                            tS = "%s > 2.0 (altId %s)" % (rsrZ, altId) if altId else "%s > 2.0" % rsrZ
                            instanceModelOutlierD.setdefault((modelId, asymId, altId, False), []).append(OutlierValue(compId, None, "RSRZ_OUTLIER", tS, rsr, None, rsrZ, "Z-Score"))
                        if rsrCc and float(rsrCc) < 0.650:
                            tS = "RSCC < 0.65 (altId %s)" % altId if altId else "RSCC < 0.65"
                            instanceModelOutlierD.setdefault((modelId, asymId, altId, False), []).append(OutlierValue(compId, None, "RSCC_OUTLIER", tS, rsrCc))
                        # Set all values for non-polymer validation score
                        if asymId in instanceTypeD and instanceTypeD[asymId] == "non-polymer":
                            instanceModelValidationD[(modelId, asymId, altId, compId)] = NonpolymerValidationInstance(
                                float(rsr) if rsr else None,
                                float(rsrCc) if rsrCc else None,
                                int(nAtomsEds) if nAtomsEds else None,
                                float(bondsRmsZ) if bondsRmsZ else None,
                                float(anglesRmsZ) if anglesRmsZ else None,
                                int(numAnglesRmsZ) if numAnglesRmsZ else None,
                                int(numBondsRmsZ) if numBondsRmsZ else None,
                                float(avgOccupancy) if avgOccupancy else None,
                                npClashD[(modelId, asymId, altId, compId)] if (modelId, asymId, altId, compId) in npClashD else 0,
                                npMogulBondOutlierD[(modelId, asymId, altId, compId)] if (modelId, asymId, altId, compId) in npMogulBondOutlierD else 0,
                                npMogulAngleOutlierD[(modelId, asymId, altId, compId)] if (modelId, asymId, altId, compId) in npMogulAngleOutlierD else 0,
                                npStereoOutlierD[(modelId, asymId, altId, compId)] if (modelId, asymId, altId, compId) in npStereoOutlierD else 0,
                            )
            # --
            logger.debug("instanceModelOutlierD %r", instanceModelOutlierD)
            logger.debug("instanceModelValidationD %r", instanceModelValidationD)

            rD = {"instanceModelOutlierD": instanceModelOutlierD, "instanceModelValidationD": instanceModelValidationD}
        except Exception as e:
            logger.exception("%s failing with %s", dataContainer.getName(), str(e))
        return rD

    def getNearestNeighborList(self, dataContainer):
        """Return a list of nearest neighbors for ligands and targets referenced
        by the target and ligand indices.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            list: [LigandTargetInstance(), ...]

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchNeighborInfo(dataContainer)
        return wD["nearestNeighbors"] if "nearestNeighbors" in wD else {}

    def getInteractionIndex(self, dataContainer):
        """Return the index dictionary of polymer and branched entity targets interactions
        with ligand entity instances.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            (dict): {targetAsymId: {(ligandAsymId, ligandCompId): [nnIndex1, nnIndex2]}}

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchNeighborInfo(dataContainer)
        return wD["interactionIndexD"] if "interactionIndexD" in wD else {}

    def getLigandNeighborIndex(self, dataContainer):
        """Return the index of nearest neighbors for ligands interacting
        with targets polymer and branched entity instances.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            (dict): {ligandAsymId: {(targetAsymId, targetAuthSeqId): nnIndex1, (): nnIndex2}

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchNeighborInfo(dataContainer)
        return wD["ligandNeighborIndexD"] if "ligandNeighborIndexD" in wD else {}

    def getTargetNeighborIndex(self, dataContainer):
        """Return the index of nearest neighbors for polymer and branched entity targets interacting
        with ligand entity instances.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            (dict): {(targetAsymId, targetAuthSeqId): {(ligandAsymId): nnIndex1, (): nnIndex2}

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchNeighborInfo(dataContainer)
        return wD["targetNeighborIndexD"] if "targetNeighborIndexD" in wD else {}

    def getLigandAtomCountD(self, dataContainer):
        """Return the ligand atom counts for all observed alternate conformers
        of ligand instances. (all atom types)

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            (dict): {ligandAsymId: {'FL': ### 'A': ###, 'B': ###}, ...  }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchNeighborInfo(dataContainer)
        return wD["ligandAtomCountD"] if "ligandAtomCountD" in wD else {}

    def getInstanceOccupancySumD(self, dataContainer):
        """Return the instance occupancy sums for all observed instances. (heavy atom types)

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            (dict): {asymId: {'FL': ### 'A': ###, 'B': ###}, ...  }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchAtomSiteInfo(dataContainer)
        return wD["occupancySumD"] if "occupancySumD" in wD else {}

    def getLigandHydrogenAtomCountD(self, dataContainer):
        """Return the ligand hydrogen atom counts for all observed alternate conformers
        of ligand instances. (all atom types)

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            (dict): {ligandAsymId: {'FL': ### 'A': ###, 'B': ###}, ...  }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchNeighborInfo(dataContainer)
        return wD["ligandHydrogenAtomCountD"] if "ligandHydrogenAtomCountD" in wD else {}

    def getLigandNeighborBoundState(self, dataContainer):
        """Return the dicitonary of ligand instances with isBound boolean status.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            (dict): {ligandAsymId: True if isBound,  ...  }

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchNeighborInfo(dataContainer)
        return wD["ligandIsBoundD"] if "ligandIsBoundD" in wD else {}

    def __fetchNeighborInfo(self, dataContainer):
        wD = self.__neighborInfoCache.get(dataContainer.getName())
        if not wD:
            wD = self.getNeighborInfo(dataContainer)
            self.__neighborInfoCache.set(dataContainer.getName(), wD)
        return wD

    def getNeighborInfo(self, dataContainer, distLimit=5.0, targetModelId=None):
        """Get bound and unbound neighbors for each non-polymer instance and ligand atom counts.

        Args:
            dataContainer (obj): DataContainer object
            distLimit (float, optional): neighbor distance limit (Angstroms). Defaults to 5.0.
            targetModelId (str, optional):  select only for this model identifier. Will default to representative model.

        Returns:
              (dict, dict, list, dict): {asymId: [LigandTargetInstance()]}, {asymId: {altId: atomcount}, }
        """
        try:
            startTime = time.time()
            #
            nearestNeighbors = []
            ligandIndexD = {}
            targetIndexD = {}
            interactionIndexD = {}
            ligandIsBoundD = {}
            ligandAtomCountD = {}
            ligandHydrogenAtomCountD = {}
            # occupancySumD = {}
            rD = {
                "targetNeighborIndexD": targetIndexD,
                "ligandNeighborIndexD": ligandIndexD,
                "nearestNeighbors": nearestNeighbors,
                "interactionIndexD": interactionIndexD,
                "ligandIsBoundD": ligandIsBoundD,
                "ligandAtomCountD": ligandAtomCountD,
                "ligandHydrogenAtomCountD": ligandHydrogenAtomCountD,
                # "occupancySumD": occupancySumD,
            }
            #
            ligandTargetInstanceD = {}
            instanceTypeD = self.getInstanceTypes(dataContainer)
            if "non-polymer" not in instanceTypeD.values():
                return rD
            #
            targetModelId = targetModelId if targetModelId else self.getRepresentativeModelId(dataContainer)
            #
            entryId = dataContainer.getName()
            logger.debug("Starting with entry %s", entryId)
            nonPolymerBoundD = self.getBoundNonpolymersByInstance(dataContainer)
            #
            instancePolymerTypeD = self.getInstancePolymerTypes(dataContainer)
            instanceEntityD = self.getInstanceEntityMap(dataContainer)
            # -----
            targetXyzL = []
            targetRefL = []
            ligandXyzD = {}
            ligandRefD = {}

            # partition the cooordinates between ligands and candidate targets
            aObj = dataContainer.getObj("atom_site")
            for ii in range(aObj.getRowCount()):
                modelId = aObj.getValue("pdbx_PDB_model_num", ii)
                if str(modelId) != targetModelId:
                    continue

                asymId = aObj.getValue("label_asym_id", ii)
                instanceType = instanceTypeD[asymId]
                polymerType = instancePolymerTypeD[asymId] if asymId in instancePolymerTypeD else None
                selectType = None
                if (instanceType == "polymer" and polymerType in ["Protein", "DNA", "RNA", "NA-hybrid"]) or instanceType == "branched":
                    selectType = "target"
                elif instanceType == "non-polymer":
                    selectType = "ligand"
                if selectType not in ["target", "ligand"]:
                    continue
                #
                atomType = aObj.getValue("type_symbol", ii)
                atomId = aObj.getValue("label_atom_id", ii)
                seqId = aObj.getValue("label_seq_id", ii)
                authSeqId = aObj.getValueOrDefault("auth_seq_id", ii, seqId)
                compId = aObj.getValue("label_comp_id", ii)
                altId = aObj.getValueOrDefault("label_alt_id", ii, None)
                xC = aObj.getValue("Cartn_x", ii)
                yC = aObj.getValue("Cartn_y", ii)
                zC = aObj.getValue("Cartn_z", ii)
                # occupancy = aObj.getValueOrDefault("occupancy", ii, "0.0")
                entityId = instanceEntityD[asymId]
                # if atomType != "H" and atomType != "D" and atomType != "T":
                #     if not altId:
                #        occupancySumD.setdefault(asymId, defaultdict(float))["FL"] += float(occupancy)
                #    else:
                #        occupancySumD.setdefault(asymId, defaultdict(float))[altId] += float(occupancy)
                #
                if selectType == "target":
                    targetXyzL.append((float(xC), float(yC), float(zC)))
                    targetRefL.append(
                        ReferenceInstance(entityId, instanceType, asymId, compId, int(seqId) if seqId not in [".", "?"] else None, authSeqId, atomId, altId, targetModelId)
                    )
                elif selectType == "ligand":
                    ligandXyzD.setdefault(asymId, []).append((float(xC), float(yC), float(zC)))
                    ligandRefD.setdefault(asymId, []).append(ReferenceInstance(entityId, instanceType, asymId, compId, None, authSeqId, atomId, altId, modelId))

                    if not altId:
                        ligandAtomCountD.setdefault(asymId, defaultdict(int))["FL"] += 1
                        if atomType == "H":
                            ligandHydrogenAtomCountD.setdefault(asymId, defaultdict(int))["FL"] += 1
                    else:
                        ligandAtomCountD.setdefault(asymId, defaultdict(int))[altId] += 1
                        if atomType == "H":
                            ligandHydrogenAtomCountD.setdefault(asymId, defaultdict(int))[altId] += 1
            #
            # ------
            logger.debug("%s targetXyzL (%d) targetRef (%d) ligandXyzD (%d) ", entryId, len(targetXyzL), len(targetXyzL), len(ligandXyzD))
            if not targetXyzL:
                return rD
            tArr = np.array(targetXyzL, order="F")
            logger.debug("targetXyzL[0] %r tArr.shape %r tArr[0] %r", targetXyzL[0], tArr.shape, tArr[0])
            tree = spatial.cKDTree(tArr)
            #
            # Calculate ligand - target interactions
            for asymId, ligXyzL in ligandXyzD.items():
                # Set kn nearest neighbors to find PER ATOM (6 for single atom, e.g., a metal ion; else default to 3)
                kn = 6 if len(ligXyzL) == 1 else 3
                # logger.info("%r %r %r", entryId, asymId, ligXyzL)
                lArr = np.array(ligXyzL, order="F")
                distance, index = tree.query(lArr, k=kn, distance_upper_bound=distLimit)  # Find the first k neighbors for the given ligand's atom
                logger.debug("%s lig asymId %s distance %r  index %r", entryId, asymId, distance, index)
                for ligIndex, (distL, indL) in enumerate(zip(distance, index)):
                    for (dist, ind) in zip(distL, indL):
                        if dist == np.inf:
                            continue
                        # ----
                        connectType = "non-bonded"
                        bondDist = dist
                        # Check if the ligand-polymer interaction is defined in struct_conn, and if so then get the bond type (metal coord. or covalent) and distance
                        if asymId in nonPolymerBoundD:
                            for tup in nonPolymerBoundD[asymId]:
                                if tup.partnerEntityType not in ["non-polymer", "water"]:
                                    tmpCompareL = [
                                        (tup.targetCompId, ligandRefD[asymId][ligIndex].compId),
                                        (tup.targetAtomId, ligandRefD[asymId][ligIndex].atomId),
                                        (tup.partnerEntityId, targetRefL[ind].entityId),
                                        (tup.partnerCompId, targetRefL[ind].compId),
                                        (tup.partnerAsymId, targetRefL[ind].asymId),
                                        (tup.partnerSeqId, targetRefL[ind].seqId),
                                        (tup.partnerAuthSeqId, targetRefL[ind].authSeqId),
                                        (tup.partnerAtomId, targetRefL[ind].atomId),
                                    ]
                                    # logger.info(
                                    #     "entryId %r, asymId %r, tmpCompareL %r, tup.connectType %r, tup.bondDistance %r (vs dist %r)",
                                    #     entryId, asymId, tmpCompareL, tup.connectType, float(tup.bondDistance), round(dist, 3)
                                    # )
                                    if all([str(k) == str(v) for k, v in tmpCompareL]) and tup.connectType in ["metal coordination", "covalent bond"]:
                                        # logger.info("match for entryId %r, asymId %r", entryId, asymId)
                                        connectType = tup.connectType
                                        bondDist = float(tup.bondDistance) if tup.bondDistance else None
                                        break
                        # ----
                        ligandTargetInstanceD.setdefault(asymId, []).append(
                            LigandTargetInstance(
                                ligandRefD[asymId][ligIndex].modelId,
                                asymId,
                                ligandRefD[asymId][ligIndex].compId,
                                ligandRefD[asymId][ligIndex].atomId,
                                ligandRefD[asymId][ligIndex].altId,
                                connectType,
                                targetRefL[ind].modelId,
                                targetRefL[ind].entityType,
                                targetRefL[ind].entityId,
                                targetRefL[ind].compId,
                                targetRefL[ind].asymId,
                                targetRefL[ind].seqId,
                                targetRefL[ind].authSeqId,
                                targetRefL[ind].atomId,
                                targetRefL[ind].altId,
                                round(bondDist, 3),
                            )
                        )
                        # ----
                if not (asymId in ligandTargetInstanceD and len(ligandTargetInstanceD[asymId])):
                    logger.debug("%s no neighbors for ligand asymId %s within %.2f", entryId, asymId, distLimit)
                    continue
            #
            # re-sort by distance -
            cloneD = copy.deepcopy(ligandTargetInstanceD)
            for asymId in cloneD:
                ligandTargetInstanceD[asymId] = sorted(cloneD[asymId], key=itemgetter(-1))
            # logger.info("Resorted ligandTargetInstanceD %r", ligandTargetInstanceD)
            #
            # --- ----
            tnD = {}
            for asymId, neighborL in ligandTargetInstanceD.items():
                isBound = False
                for neighbor in neighborL:
                    tnD.setdefault(asymId, {}).setdefault((neighbor.partnerAsymId, neighbor.partnerAuthSeqId), []).append(neighbor)
                    if neighbor.connectType != "non-bonded":
                        isBound = True
                ligandIsBoundD[asymId] = isBound
            # --- ----
            # Example tnD (for 3FJQ):
            #     {'C': {('A', '184'): [LigandTargetInstance1, LigandTargetInstance2, LigandTargetInstance3, ...],
            #      'D': {('A', '171'): [LigandTargetInstance1, LigandTargetInstance2],
            #            ('A', '184'): [LigandTargetInstance1]}}
            #
            jj = 0
            for asymId, nD in tnD.items():
                for (pAsymId, pAuthSeqId), nL in nD.items():
                    #
                    ligandIndexD.setdefault(asymId, {})[(pAsymId, pAuthSeqId)] = jj  # This is used for building rcsb_target_neighbors (see buildInstanceTargetNeighbors)
                    # E.g.: {'C': {('A', '184'): 0}, 'D': {('A', '171'): 1, ('A', '184'): 2}}
                    # (Organized by ligand asymId first, then polymer chain and residue. Indices correspond to the index in nearestNeighbors)
                    #
                    targetIndexD.setdefault((pAsymId, pAuthSeqId), {})[asymId] = jj  # This is used for building rcsb_ligand_neighbors (see buildInstanceLigandNeighbors)
                    # E.g.,: {('A', '184'): {'C': 0, 'D': 2}, ('A', '171'): {'D': 1}}
                    # (Organized by polymer chain and residue first, then by ligand)
                    #
                    # Below, pick only the "first" (i.e., closest) interaction of the list (for a given "pAsymId, pAuthSeqId" pair).
                    # So, for a given ligand that has multiple interactions with the same residue of the polymer, it picks the shortest one.
                    # Note that this will lead to loss of information for cases where different atoms of a given ligand have interactions with the same polymer residue
                    nearestNeighborD = nL[0]
                    nearestNeighbors.append(nearestNeighborD)
                    #
                    ligandCompId = nearestNeighborD.ligandCompId if nearestNeighborD.ligandCompId else None
                    if not ligandCompId:
                        logger.error("Missing ligand compId for asymId %r", asymId)
                    interactionIndexD.setdefault(pAsymId, {}).setdefault((asymId, ligandCompId), []).append(jj)
                    #
                    jj += 1
            #
            # --- ----
        except Exception as e:
            logger.exception("Failing for %r with %r", dataContainer.getName() if dataContainer else None, str(e))
        #
        logger.debug("Completed %s at %s (%.4f seconds)", dataContainer.getName(), time.strftime("%Y %m %d %H:%M:%S", time.localtime()), time.time() - startTime)
        return rD

    def getCompModelDb2L(self, dataContainer):
        """Get list of comp model database_ids from database_2

        Args:
            dataContainer (object): mmif.api.DataContainer object instance

        Returns:
            compModelDb2L = ["AlphaFoldDB", ...]

        """

        compModelDb2L = []

        try:
            db2EnumL = ["AlphaFoldDB", "MODBASE", "ModelArchive", "SWISS-MODEL_REPOSITORY", "AF", "MA", "SMR"]
            if dataContainer.exists("database_2"):
                eObj = dataContainer.getObj("database_2")
                dbL = eObj.getAttributeValueList("database_2")
                compModelDb2L = [db for db in dbL if db in db2EnumL]
            else:
                dObj = dataContainer.getObj("entry")
                entryId = dObj.getValue("id", 0)
                if entryId.upper().startswith("MA"):
                    compModelDb2L = ["MA"]
                elif entryId.upper().startswith("AF"):
                    compModelDb2L = ["AF"]

        except Exception as e:
            logger.exception("Missing database_2 information. %r failing with %s", dataContainer.getName(), str(e))

        if not compModelDb2L:
            logger.debug("Missing database_2 information for %r", dataContainer.getName())

        return compModelDb2L

    def getMaQaMetricType(self, dataContainer):
        """ Get mapping of metric ids to metric names and types for computed models
            from ma_qa_metric

            Args:
                dataContainer (object): mmif.api.DataContainer object instance

        """

        maQaMetricTypeD = {}
        maQaMetricLocalTypeD = {}
        maQaMetricGlobalTypeD = {}

        compModelScoreTypeEnumD = {
            "zscore": "MA_QA_METRIC_LOCAL_TYPE_ZSCORE",
            "energy": "MA_QA_METRIC_LOCAL_TYPE_ENERGY",
            "distance": "MA_QA_METRIC_LOCAL_TYPE_DISTANCE",
            "normalized score": "MA_QA_METRIC_LOCAL_TYPE_NORMALIZED_SCORE",
            "pLDDT": "MA_QA_METRIC_LOCAL_TYPE_PLDDT",
            "pLDDT in [0,1]": "MA_QA_METRIC_LOCAL_TYPE_PLDDT_[0,1]",
            "pLDDT all-atom": "MA_QA_METRIC_LOCAL_TYPE_PLDDT_ALL-ATOM",
            "pLDDT all-atom in [0,1]": "MA_QA_METRIC_LOCAL_TYPE_PLDDT_ALL-ATOM_[0,1]",
            "PAE": "MA_QA_METRIC_LOCAL_TYPE_PAE",
            "pTM": "MA_QA_METRIC_LOCAL_TYPE_PTM",
            "ipTM": "MA_QA_METRIC_LOCAL_TYPE_IPTM",
            "contact probability": "MA_QA_METRIC_LOCAL_TYPE_CONTACT_PROBABILITY",
            "other": "MA_QA_METRIC_LOCAL_TYPE_OTHER"
        }

        try:
            if dataContainer.exists("ma_qa_metric"):
                aObj = dataContainer.getObj("ma_qa_metric")
                for ii in range(aObj.getRowCount()):
                    mId = aObj.getValue("id", ii)
                    mMode = aObj.getValue("mode", ii)
                    mType = aObj.getValue("type", ii)
                    mName = aObj.getValue("name", ii)
                    if mName == "pLDDT" and mType.lower() == "other":
                        mType = "pLDDT"
                    if mMode == "local":
                        maQaMetricLocalTypeD[mId] = {"type": compModelScoreTypeEnumD[mType], "name": mName}
                    if mMode == "global":
                        maQaMetricGlobalTypeD[mId] = {"type": mType, "name": mName}

                maQaMetricTypeD["maQaMetricLocalTypeD"] = maQaMetricLocalTypeD
                maQaMetricTypeD["maQaMetricGlobalTypeD"] = maQaMetricGlobalTypeD

        except Exception as e:
            logger.exception("Failing for %s with %s", dataContainer.getName(), str(e))

        return maQaMetricTypeD

    def getCompModelLocalScores(self, dataContainer):
        """ Get Local QA Scores for computed models from the ModelCIF file
            (ma_qa_metric_local) and convert to objects corresponding to
            rcsb_entity_instance_feature in the RCSB extension dictionary

            Args:
                dataContainer (object): mmif.api.DataContainer object instance

        """
        metricValD = OrderedDict()
        compModelLocalScoresD = OrderedDict()

        try:
            if dataContainer.exists("ma_qa_metric_local"):
                tObj = dataContainer.getObj("ma_qa_metric_local")
                repModelId = self.getRepresentativeModelId(dataContainer)
                dL = []
                for ii in range(tObj.getRowCount()):
                    modelId = tObj.getValue("model_id", ii)
                    if str(modelId) != repModelId:  # Skip non-representative models
                        continue
                    seqId = tObj.getValue("label_seq_id", ii)
                    asymId = tObj.getValue("label_asym_id", ii)
                    metricId = tObj.getValue("metric_id", ii)
                    metricV = tObj.getValue("metric_value", ii)
                    tId = str(modelId) + "_" + str(asymId) + "_" + str(metricId) + "_" + str(seqId)
                    if seqId and seqId not in [".", "?"]:   # Eliminates non-polymers and branched
                        if tId not in dL:
                            metricValD.setdefault((modelId, asymId, metricId), []).append((seqId, metricV))
                            dL.append(tId)

                for (modelId, asymId, metricId), aL in metricValD.items():
                    tD = {}
                    sL = sorted(aL, key=lambda i: int(i[0]))
                    mL = [int(s[0]) for s in sL]
                    for ii in range(len(sL)):
                        seqId = sL[ii][0]
                        metricV = sL[ii][1]
                        for tup in list(self.__toRangeList(mL)):
                            beg = tup[0]
                            end = tup[1]
                            if int(beg) <= int(seqId) <= int(end):
                                tD.setdefault(int(beg), []).append(float(metricV))
                    compModelLocalScoresD.setdefault((modelId, asymId, metricId), []).append(tD)

        except Exception as e:
            logger.exception("Failing for %s with %s", dataContainer.getName(), str(e))

        return compModelLocalScoresD

    def getRepresentativeModels(self, dataContainer):
        """Return the list of representative models

        Example:
            #
            _pdbx_nmr_ensemble.entry_id                                      5TM0
            _pdbx_nmr_ensemble.conformers_calculated_total_number            15
            _pdbx_nmr_ensemble.conformers_submitted_total_number             15
            _pdbx_nmr_ensemble.conformer_selection_criteria                  'all calculated structures submitted'
            _pdbx_nmr_ensemble.representative_conformer                      ?
            _pdbx_nmr_ensemble.average_constraints_per_residue               ?
            _pdbx_nmr_ensemble.average_constraint_violations_per_residue     ?
            _pdbx_nmr_ensemble.maximum_distance_constraint_violation         ?
            _pdbx_nmr_ensemble.average_distance_constraint_violation         ?
            _pdbx_nmr_ensemble.maximum_upper_distance_constraint_violation   ?
            _pdbx_nmr_ensemble.maximum_lower_distance_constraint_violation   ?
            _pdbx_nmr_ensemble.distance_constraint_violation_method          ?
            _pdbx_nmr_ensemble.maximum_torsion_angle_constraint_violation    ?
            _pdbx_nmr_ensemble.average_torsion_angle_constraint_violation    ?
            _pdbx_nmr_ensemble.torsion_angle_constraint_violation_method     ?
            #
            _pdbx_nmr_representative.entry_id             5TM0
            _pdbx_nmr_representative.conformer_id         1
            _pdbx_nmr_representative.selection_criteria   'fewest violations'
        """
        eObj = dataContainer.getObj("entry")
        entryId = eObj.getValue("id", 0)
        repModelL = []
        mIdL = self.getModelIdList(dataContainer)
        #
        if mIdL:
            # If NMR structure...
            if self.hasMethodNMR(dataContainer):
                if dataContainer.exists("pdbx_nmr_representative"):
                    tObj = dataContainer.getObj("pdbx_nmr_representative")
                    if tObj.hasAttribute("conformer_id"):
                        for ii in range(tObj.getRowCount()):
                            nn = tObj.getValue("conformer_id", ii)
                            if nn is not None and nn.isdigit() and nn in mIdL and nn not in repModelL:
                                repModelL.append(nn)
                if dataContainer.exists("pdbx_nmr_ensemble"):
                    tObj = dataContainer.getObj("pdbx_nmr_ensemble")
                    if tObj.hasAttribute("representative_conformer"):
                        nn = tObj.getValue("representative_conformer", 0)
                        if nn is not None and nn and nn.isdigit() and nn in mIdL and nn not in repModelL:
                            repModelL.append(nn)
            #
            if not repModelL:
                logger.debug("Using the first model as representative model for %s.", dataContainer.getName())
                repModelL = ["1"] if "1" in mIdL else [mIdL[0]]

        else:
            logger.error("Missing model data for %s. Check file for data issue (model ID should exist in atom_site.pdbx_PDB_model_num).", dataContainer.getName())
            raise ValueError(
                "Missing model data for %s. Check file for data issue (model ID should exist in atom_site.pdbx_PDB_model_num)."
                % (dataContainer.getName())
            )
        #
        logger.debug("Representative model list for entryId %r %r", entryId, repModelL)
        #
        return repModelL

    def getMethodList(self, dataContainer):
        """Return experimental or computational method list.
        Args:
            dataContainer (object): mmif.api.DataContainer object instance
        Returns:
            methodL (list): List of dictionary experimental method names
        """
        #
        methodL = []
        try:
            if dataContainer.exists("exptl"):
                xObj = dataContainer.getObj("exptl")
                methodL = xObj.getAttributeValueList("method")
            elif dataContainer.exists("ma_model_list"):
                mObj = dataContainer.getObj("ma_model_list")
                methodL = mObj.getAttributeUniqueValueList("model_type")
        except Exception as e:
            logger.debug("Failed to get method list with %s", str(e))
        #
        return methodL

    def getRepresentativeModelId(self, dataContainer):
        """Return the first representative model ID.
        """
        repModelL = self.getRepresentativeModels(dataContainer)
        repModelId = repModelL[0]
        return str(repModelId)

    def getIhmRepresentativeModelId(self, dataContainer):
        """Return the first representative model ID for integrative structures
        """
        modelIdL = []
        repModelIdL = []
        repModelId = None
        try:
            if dataContainer.exists("ihm_model_list"):
                mObj = dataContainer.getObj("ihm_model_list")
                modelIdL = mObj.getAttributeUniqueValueList("model_id")
                modelIdL.sort(key=int)
            if dataContainer.exists("ihm_model_representative"):
                mrObj = dataContainer.getObj("ihm_model_representative")
                repModelIdL = mrObj.getAttributeUniqueValueList("model_id")
                repModelIdL.sort(key=int)
                if repModelIdL:
                    repModelId = self.__getIhmLargestModel(dataContainer, repModelIdL)
            if not repModelId:
                if modelIdL:
                    repModelId = self.__getIhmLargestModel(dataContainer, modelIdL)
            if not repModelId:
                repModelId = modelIdL[0]
        except Exception as e:
            logger.exception("Failing for %s with %s", dataContainer.getName(), str(e))
        return str(repModelId)

    def __getIhmLargestModel(self, dataContainer, modelIdL):
        """Return the largest model in assembly size
        """
        sD = {}
        mObj = dataContainer.getObj("ihm_model_list")
        for modelId in modelIdL:
            repAssemblyId = mObj.selectValuesWhere("assembly_id", modelId, "model_id")[0]
            if dataContainer.exists("ihm_struct_assembly_details"):
                sdObj = dataContainer.getObj("ihm_struct_assembly_details")
                sdL = sdObj.selectValuesWhere("asym_id", repAssemblyId, "assembly_id")
                sD[modelId] = len(list(set(sdL)))
        if sD:
            rModelId = max(sD, key=sD.get)
        else:
            rModelId = modelIdL[0]
        return str(rModelId)

    def filterIntegrativeMethod(self, dataContainer):
        """Apply a standard filter to the input integrative method list returning a method count and
            a simplified method name.

        """
        methodCount = 1
        expMethod = None
        mc = 0
        dsNameTypeMapD = self.getIhmDsNameTypeMapD()
        if dataContainer.exists("ihm_model_list"):
            expMethod = "Integrative"
        if dataContainer.exists("ihm_dataset_list"):
            dObj = dataContainer.getObj("ihm_dataset_list")
            dtL = dObj.getAttributeUniqueValueList("data_type")
            for dt in dtL:
                if dt in dsNameTypeMapD:
                    if dsNameTypeMapD[dt] == "Experimental data":
                        mc += 1
        if mc:
            methodCount = mc
        return methodCount, expMethod

    def getIhmDsNameTypeMapD(self):
        dsNameTypeMapD = {
            "NMR data": "Experimental data",
            "3DEM volume": "Experimental data",
            "2DEM class average": "Experimental data",
            "EM raw micrographs": "Experimental data",
            "X-ray diffraction data": "Experimental data",
            "SAS data": "Experimental data",
            "CX-MS data": "Experimental data",
            "Crosslinking-MS data": "Experimental data",
            "Mass Spectrometry data": "Experimental data",
            "EPR data": "Experimental data",
            "H/D exchange data": "Experimental data",
            "Single molecule FRET data": "Experimental data",
            "Ensemble FRET data": "Experimental data",
            "Experimental model": "Starting model",
            "Comparative model": "Starting model",
            "Integrative model": "Starting model",
            "De Novo model": "Starting model",
            "Predicted contacts": "Computed restraints",
            "Mutagenesis data": "Experimental data",
            "DNA footprinting data": "Experimental data",
            "Hydroxyl radical footprinting data": "Experimental data",
            "Yeast two-hybrid screening data": "Experimental data",
            "Quantitative measurements of genetic interactions": "Experimental data",
            "Other": "Other"
        }
        return dsNameTypeMapD

    # METHOD TO BE DELETED
    def __getValidationData(self, tObj, iObj, fields, idField, metricValD, dL, instL):
        cndL4 = [(idField, "in", instL)]
        kL = tObj.selectIndicesWhereOpConditions(cndL4)
        for ii in kL:
            instId = tObj.getValue(idField, ii)
            [[entityId, asymId, compId, authAsymId, seqId, modelNum]] = iObj.selectValueListWhere(
                [
                    "entity_id",
                    "label_asym_id",
                    "label_comp_id",
                    "auth_asym_id",
                    "label_seq_id",
                    "PDB_model_num"
                ], instId, "id"
            )
            for iFd, iFdv in fields.items():
                value = tObj.getValueOrDefault(iFdv, ii, None)
                if value is not None:
                    tId = iFd + "_" + entityId + "_" + asymId + "_" + modelNum + "_" + seqId
                    if tId not in dL:
                        metricValD.setdefault((entityId, asymId, authAsymId, modelNum, iFd, seqId and seqId not in [".", "?"]), []).append((compId, seqId, value))
                        dL.append(tId)

    # METHOD TO BE DELETED
    def __getLocalValidationPrev(self, dataContainer):
        """ Get Local validation data from the Validation report
            (e.g., pdbx_vrpt_model_instance_map_fitting.Q_score) and convert to objects corresponding to
            rcsb_entity_instance_feature in the RCSB extension dictionary

            Args:
                dataContainer (object): mmif.api.DataContainer object instance

        """
        metricValD = OrderedDict()
        localDataD = OrderedDict()

        # Exit if no source categories are present
        if not (
                dataContainer.exists("pdbx_vrpt_model_instance_map_fitting")
                or dataContainer.exists("pdbx_vrpt_model_instance_density")
                or dataContainer.exists("pdbx_vrpt_model_instance_geometry")
                or dataContainer.exists("pdbx_vrpt_model_instance")
        ):
            return localDataD

        ValidationFields = {
            "RSCC": "RSRCC",
            "RSR": "RSR",
            "RSRZ": "RSRZ",
            "Q_SCORE": "Q_score",
            "NATOMS_EDS": "natoms_eds"
        }
        OutlierCountFields = {
            "BOND_OUTLIERS": "count_bond_outliers",
            "ANGLE_OUTLIERS": "count_angle_outliers",
            "CLASHES": "count_clashes",
            "SYMM_CLASHES": "count_symm_clashes",
            "CHIRAL_OUTLIERS": "count_chiral_outliers",
            "PLANE_OUTLIERS": "count_plane_outliers",
            "MOGUL_BOND_OUTLIERS": "count_mogul_bond_outliers",
            "MOGUL_ANGLE_OUTLIERS": "count_mogul_angle_outliers",
            "MOGUL_TORSION_OUTLIERS": "count_mogul_torsion_outliers",
            "MOGUL_RING_OUTLIERS": "count_mogul_ring_outliers"
        }

        try:
            # Get representative model
            repModelId = self.getRepresentativeModelId(dataContainer)

            iObj = dataContainer.getObj("pdbx_vrpt_model_instance")

            instanceTypeD = self.getInstanceTypes(dataContainer)
            pAsymL = [k for k, v in instanceTypeD.items() if v == "polymer"]
            if not pAsymL:
                return localDataD
            cndL2 = [("label_asym_id", "in", pAsymL), ("PDB_model_num", "eq", repModelId)]
            instL = iObj.selectValuesWhereOpConditions("id", cndL2)

            eObj = dataContainer.getObj("entry")
            entryId = eObj.getValue("id", 0)
            logger.debug("Starting validation report feature for %s.", entryId)
            dL = []

            # pdbx_vrpt_model_instance_map_fitting
            if dataContainer.exists("pdbx_vrpt_model_instance_map_fitting"):
                tObj = dataContainer.getObj("pdbx_vrpt_model_instance_map_fitting")
                self.__getValidationData(tObj, iObj, ValidationFields, "instance_id", metricValD, dL, instL)

            # pdbx_vrpt_model_instance_density
            if dataContainer.exists("pdbx_vrpt_model_instance_density"):
                tObj = dataContainer.getObj("pdbx_vrpt_model_instance_density")
                self.__getValidationData(tObj, iObj, ValidationFields, "instance_id", metricValD, dL, instL)

            # pdbx_vrpt_model_instance_geometry
            if dataContainer.exists("pdbx_vrpt_model_instance_geometry"):
                tObj = dataContainer.getObj("pdbx_vrpt_model_instance_geometry")
                self.__getValidationData(tObj, iObj, {"OWAB": "OWAB", "AVERAGE_OCCUPANCY": "average_occupancy"}, "instance_id", metricValD, dL, instL)

            # pdbx_vrpt_model_instance
            if dataContainer.exists("pdbx_vrpt_model_instance"):
                tObj = dataContainer.getObj("pdbx_vrpt_model_instance")
                self.__getValidationData(tObj, iObj, OutlierCountFields, "id", metricValD, dL, instL)

            for (entityId, asymId, authAsymId, modelId, attrId, hasSeq), aL in metricValD.items():
                tD = {}
                if hasSeq:
                    sL = sorted(aL, key=lambda i: int(i[1]))
                    mL = [int(s[1]) for s in sL]
                    begCompId = ""
                    for ii in range(len(sL)):
                        compId = sL[ii][0]
                        seqId = sL[ii][1]
                        metricV = sL[ii][2]
                        for tup in list(self.__toRangeList(mL)):
                            beg = tup[0]
                            end = tup[1]
                            if int(beg) == int(seqId):
                                begCompId = compId
                            if int(beg) <= int(seqId) <= int(end) and metricV is not None:
                                tD.setdefault((begCompId, int(beg)), []).append(float(metricV))
                elif len(aL) > 0:
                    ii = 0
                    compId = aL[ii][0]
                    metricV = aL[ii][2]
                    tD.setdefault((compId, 0), []).append(float(metricV))
                localDataD.setdefault((entityId, asymId, authAsymId, modelId, attrId, hasSeq), []).append(tD)

        except Exception as e:
            logger.exception("Failing for %s with %s", dataContainer.getName(), str(e))

        return localDataD

    def __processValidationData(self, ii, cL, vObj, vFields, metricValD):
        [modelId, asymId, compId, _altId, seqId, entityId, authAsymId] = cL
        for k, v in vFields.items():
            value = vObj.getValueOrDefault(k, ii, defaultValue=None)
            if value is not None:
                metricValD.setdefault((entityId, asymId, authAsymId, modelId, v, seqId and seqId not in [".", "?"]), {}).setdefault((seqId, compId), (compId, seqId, value))

    def __getLocalValidation(self, dataContainer):
        """ Get Local validation data from the Validation report
            (e.g., pdbx_vrpt_model_instance_map_fitting.Q_score) and convert to objects corresponding to
            rcsb_entity_instance_feature in the RCSB extension dictionary

            Args:
                dataContainer (object): mmif.api.DataContainer object instance

        """
        metricValD = OrderedDict()
        localDataD = OrderedDict()
        instanceModelOutlierD = {}
        instanceModelValidationD = {}
        rD = {}

        # Exit if no source category is present
        if not dataContainer.exists("pdbx_vrpt_model_instance"):
            return rD

        v1Fields = {
            "RSRCC": "RSCC",
            "RSR": "RSR",
            "RSRZ": "RSRZ",
            "Q_score": "Q_SCORE",
            "natoms_eds": "NATOMS_EDS"
        }
        v2Fields = {
            "RSRCC": "RSCC",
            "RSR": "RSR",
            "RSRZ": "RSRZ",
            "natoms_eds": "NATOMS_EDS"
        }
        gFields = {
            "OWAB": "OWAB",
            "average_occupancy": "AVERAGE_OCCUPANCY",
        }
        iFields = {
            "count_bond_outliers": "BOND_OUTLIERS",
            "count_angle_outliers": "ANGLE_OUTLIERS",
            "count_clashes": "CLASHES",
            "count_symm_clashes": "SYMM_CLASHES",
            "count_chiral_outliers": "CHIRAL_OUTLIERS",
            "count_plane_outliers": "PLANE_OUTLIERS",
            "count_mogul_bond_outliers": "MOGUL_BOND_OUTLIERS",
            "count_mogul_angle_outliers": "MOGUL_ANGLE_OUTLIERS",
            "count_mogul_torsion_outliers": "MOGUL_TORSION_OUTLIERS",
            "count_mogul_ring_outliers": "MOGUL_RING_OUTLIERS"
        }
        startTime = time.time()

        try:
            logger.debug("Starting validation report feature for %s", dataContainer.getName())
            npMogulBondOutlierD = {}
            npMogulAngleOutlierD = {}
            npStereoOutlierD = {}
            npClashD = {}
            #
            iObj = dataContainer.getObj("pdbx_vrpt_model_instance")
            if not iObj:
                return rD

            # Get representative model
            repModelId = self.getRepresentativeModelId(dataContainer)

            # Get number of residues - Leaving it here for debugging purposes
            # modeledCount, unModeledCount = self.getDepositedMonomerCounts(dataContainer, modelId=repModelId)
            # depPolyMonomerCount = modeledCount + unModeledCount

            instanceTypeD = self.getInstanceTypes(dataContainer)
            npAsymL = [k for k, v in instanceTypeD.items() if v in ["polymer", "non-polymer"]]
            # Skip polymer validation features for large entries - Leaving it here for debugging purposes
            # if modeledCount > 25000:
            #     npAsymL = [k for k, v in instanceTypeD.items() if v in ["non-polymer"]]
            #     logger.info("Skipping polymer validation features for large entry PDBID %s", dataContainer.getName())
            if not npAsymL:
                return rD
            cD = {}
            nD = {}
            # cndL2 = [("label_asym_id", "in", npAsymL), ("PDB_model_num", "eq", repModelId)]
            # kL = iObj.selectIndicesWhereOpConditions(cndL2)
            # for ii in kL:
            for ii in range(iObj.getRowCount()):
                modelId = iObj.getValueOrDefault("PDB_model_num", ii, defaultValue=None)
                if str(modelId) != repModelId:  # Skip non-representative models
                    continue
                asymId = iObj.getValueOrDefault("label_asym_id", ii, defaultValue=None)
                if asymId in npAsymL:
                    instId = iObj.getValueOrDefault("id", ii, defaultValue=None)
                    entityId = iObj.getValueOrDefault("entity_id", ii, defaultValue=None)
                    authAsymId = iObj.getValueOrDefault("auth_asym_id", ii, defaultValue=None)
                    altId = iObj.getValueOrDefault("label_alt_id", ii, defaultValue=None)
                    compId = iObj.getValueOrDefault("label_comp_id", ii, defaultValue=None)
                    seqId = iObj.getValueOrDefault("label_seq_id", ii, defaultValue=None)
                    cD[instId] = [modelId, asymId, compId, altId, seqId, entityId, authAsymId]
                    if seqId:
                        self.__processValidationData(ii, cD[instId], iObj, iFields, metricValD)
                    else:
                        countClashes = iObj.getValueOrDefault("count_clashes", ii, defaultValue=None)
                        countMogulAngleOutliers = iObj.getValueOrDefault("count_mogul_angle_outliers", ii, defaultValue=None)
                        countMogulBondOutliers = iObj.getValueOrDefault("count_mogul_bond_outliers", ii, defaultValue=None)
                        countStereoOutliers = iObj.getValueOrDefault("count_chiral_outliers", ii, defaultValue=None)
                        npMogulBondOutlierD[(modelId, asymId, altId, compId)] = countMogulBondOutliers if countMogulBondOutliers else 0
                        npMogulAngleOutlierD[(modelId, asymId, altId, compId)] = countMogulAngleOutliers if countMogulAngleOutliers else 0
                        npClashD[(modelId, asymId, altId, compId)] = countClashes if countClashes else 0
                        npStereoOutlierD[(modelId, asymId, altId, compId)] = countStereoOutliers if countStereoOutliers else 0
                        nD[instId] = [modelId, asymId, compId, altId, seqId, entityId, authAsymId]
            # cndL3 = [("instance_id", "in", cD)]

            vObj = None
            vD = {}
            if dataContainer.exists("pdbx_vrpt_model_instance_map_fitting"):
                vObj = dataContainer.getObj("pdbx_vrpt_model_instance_map_fitting")
                # kL = vObj.selectIndicesWhereOpConditions(cndL3)
                # for ii in kL:
                for ii in range(vObj.getRowCount()):
                    instId = vObj.getValueOrDefault("instance_id", ii, defaultValue=None)
                    if instId not in cD:
                        continue
                    [modelId, asymId, compId, altId, seqId, entityId, authAsymId] = cD[instId]
                    if seqId:
                        self.__processValidationData(ii, cD[instId], vObj, v1Fields, metricValD)
                    else:
                        rsrCc = vObj.getValueOrDefault("RSRCC", ii, defaultValue=None)
                        rsr = vObj.getValueOrDefault("RSR", ii, defaultValue=None)
                        rsrZ = vObj.getValueOrDefault("RSRZ", ii, defaultValue=None)
                        nAtomsEds = vObj.getValueOrDefault("natoms_eds", ii, defaultValue=None)
                        vD[instId] = [rsrCc, rsr, rsrZ, nAtomsEds]
            elif dataContainer.exists("pdbx_vrpt_model_instance_density"):
                vObj = dataContainer.getObj("pdbx_vrpt_model_instance_density")
                # kL = vObj.selectIndicesWhereOpConditions(cndL3)
                # for ii in kL:
                for ii in range(vObj.getRowCount()):
                    instId = vObj.getValueOrDefault("instance_id", ii, defaultValue=None)
                    if instId not in cD:
                        continue
                    [modelId, asymId, compId, altId, seqId, entityId, authAsymId] = cD[instId]
                    if seqId:
                        self.__processValidationData(ii, cD[instId], vObj, v2Fields, metricValD)
                    else:
                        rsrCc = vObj.getValueOrDefault("RSRCC", ii, defaultValue=None)
                        rsr = vObj.getValueOrDefault("RSR", ii, defaultValue=None)
                        rsrZ = vObj.getValueOrDefault("RSRZ", ii, defaultValue=None)
                        nAtomsEds = vObj.getValueOrDefault("natoms_eds", ii, defaultValue=None)
                        vD[instId] = [rsrCc, rsr, rsrZ, nAtomsEds]

            gObj = None
            gD = {}
            if dataContainer.exists("pdbx_vrpt_model_instance_geometry"):
                gObj = dataContainer.getObj("pdbx_vrpt_model_instance_geometry")
                # kL = gObj.selectIndicesWhereOpConditions(cndL3)
                # for ii in kL:
                for ii in range(gObj.getRowCount()):
                    instId = gObj.getValueOrDefault("instance_id", ii, defaultValue=None)
                    if instId not in cD:
                        continue
                    [modelId, asymId, compId, altId, seqId, entityId, authAsymId] = cD[instId]
                    if seqId:
                        self.__processValidationData(ii, cD[instId], gObj, gFields, metricValD)
                    else:
                        anglesRmsZ = gObj.getValueOrDefault("angles_RMSZ", ii, defaultValue=None)
                        bondsRmsZ = gObj.getValueOrDefault("bonds_RMSZ", ii, defaultValue=None)
                        numAnglesRmsZ = gObj.getValueOrDefault("num_angles_RMSZ", ii, defaultValue=None)
                        numBondsRmsZ = gObj.getValueOrDefault("num_bonds_RMSZ", ii, defaultValue=None)
                        avgOccupancy = gObj.getValueOrDefault("average_occupancy", ii, defaultValue=None)
                        software = gObj.getValueOrDefault("program_for_bond_angle_geometry", ii, defaultValue=None)
                        gD[instId] = [anglesRmsZ, bondsRmsZ, numAnglesRmsZ, numBondsRmsZ, avgOccupancy, software]

            for instId, nL in nD.items():
                [modelId, asymId, compId, altId, seqId, entityId, authAsymId] = nL
                if str(modelId) != repModelId:  # Skip non-representative models
                    continue
                if not seqId:
                    rsr = rsrZ = rsrCc = nAtomsEds = None
                    anglesRmsZ = bondsRmsZ = numAnglesRmsZ = numBondsRmsZ = avgOccupancy = software = None
                    # Get the matching data from pdbx_vrpt_model_instance_density or pdbx_vrpt_model_instance_map_fitting
                    if instId in vD:
                        [rsrCc, rsr, rsrZ, nAtomsEds] = vD[instId]
                    # Get the matching data from pdbx_vrpt_model_instance_geometry
                    if instId in gD:
                        [anglesRmsZ, bondsRmsZ, numAnglesRmsZ, numBondsRmsZ, avgOccupancy, software] = gD[instId]
                    # Evaluate RSRZ_OUTLIER and RSCC_OUTLIER for non-polymer annotations
                    if rsrZ and float(rsrZ) > 2.0:
                        tS = "%s > 2.0 (altId %s)" % (rsrZ, altId) if altId else "%s > 2.0" % rsrZ
                        instanceModelOutlierD.setdefault((modelId, asymId, altId, False), []).append(OutlierValue(compId, None, "RSRZ_OUTLIER", tS, rsr, None, rsrZ, "Z-Score"))
                    if rsrCc and float(rsrCc) < 0.650:
                        tS = "RSCC < 0.65 (altId %s)" % altId if altId else "RSCC < 0.65"
                        instanceModelOutlierD.setdefault((modelId, asymId, altId, False), []).append(OutlierValue(compId, None, "RSCC_OUTLIER", tS, rsrCc))
                    # Set all values for non-polymer validation score
                    if asymId in instanceTypeD and instanceTypeD[asymId] == "non-polymer":
                        instanceModelValidationD[(modelId, asymId, altId, compId)] = NonpolymerValidationInstance(
                            float(rsr) if rsr else None,
                            float(rsrCc) if rsrCc else None,
                            int(nAtomsEds) if nAtomsEds else None,
                            float(bondsRmsZ) if bondsRmsZ else None,
                            float(anglesRmsZ) if anglesRmsZ else None,
                            int(numAnglesRmsZ) if numAnglesRmsZ else None,
                            int(numBondsRmsZ) if numBondsRmsZ else None,
                            float(avgOccupancy) if avgOccupancy else None,
                            npClashD[(modelId, asymId, altId, compId)] if (modelId, asymId, altId, compId) in npClashD else 0,
                            npMogulBondOutlierD[(modelId, asymId, altId, compId)] if (modelId, asymId, altId, compId) in npMogulBondOutlierD else 0,
                            npMogulAngleOutlierD[(modelId, asymId, altId, compId)] if (modelId, asymId, altId, compId) in npMogulAngleOutlierD else 0,
                            npStereoOutlierD[(modelId, asymId, altId, compId)] if (modelId, asymId, altId, compId) in npStereoOutlierD else 0,
                        )
            # Set validation data for polymer features
            for (entityId, asymId, authAsymId, modelId, attrId, hasSeq), aD in metricValD.items():
                tD = {}
                if hasSeq:
                    sL = sorted(aD.values(), key=lambda item: int(item[1]))
                    mL = [int(s[1]) for s in sL]
                    begCompId = ""
                    for ii in range(len(sL)):
                        compId = sL[ii][0]
                        seqId = sL[ii][1]
                        metricV = sL[ii][2]
                        for tup in list(self.__toRangeList(mL)):
                            beg = tup[0]
                            end = tup[1]
                            if int(beg) == int(seqId):
                                begCompId = compId
                            if int(beg) <= int(seqId) <= int(end) and metricV is not None:
                                tD.setdefault((begCompId, int(beg)), []).append(float(metricV))
                    localDataD.setdefault((entityId, asymId, authAsymId, modelId, attrId, hasSeq), []).append(tD)

            rD = {"localDataD": localDataD, "instanceModelValidationD": instanceModelValidationD, "instanceModelOutlierD": instanceModelOutlierD}
        except Exception as e:
            logger.exception("Failing for %s with %s", dataContainer.getName(), str(e))

        endTime = time.time()
        logger.debug("Completed at %s (%.4f seconds) PDBID %s", time.strftime("%Y %m %d %H:%M:%S", time.localtime()), endTime - startTime, dataContainer.getName())
        return rD

    def __fetchLocalValidationData(self, dataContainer):
        wD = self.__localValidationCache.get(dataContainer.getName())
        if not wD:
            # wD = self.__getLocalValidationPrev(dataContainer)
            wD = self.__getLocalValidation(dataContainer)
            self.__localValidationCache.set(dataContainer.getName(), wD)
        return wD

    def getLocalValidationData(self, dataContainer):
        """Return a dictionary of polymer model outliers.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(modelId, asymId): (seqId,compId), ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchLocalValidationData(dataContainer)
        return wD["localDataD"] if "localDataD" in wD else {}

    def getNonpolyValidationData(self, dataContainer):
        """Return a dictionary of nonpolymer validation details.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(modelId, asymId): NonpolymerValidationInstance(rsr, rsrCc, bondsRmsZ, anglesRmsZ,
                                             intermolecular_clashes, mogul_bond_outliers, mogul_angle_outliers, stereo_outliers)}

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchLocalValidationData(dataContainer)
        return wD["instanceModelValidationD"] if "instanceModelValidationD" in wD else {}

    def getNonpolyOutlierData(self, dataContainer):
        """Return a dictionary of nonpolymer validation details.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {(modelId, asymId): NonpolymerValidationInstance(rsr, rsrCc, bondsRmsZ, anglesRmsZ,
                                             intermolecular_clashes, mogul_bond_outliers, mogul_angle_outliers, stereo_outliers)}

        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchLocalValidationData(dataContainer)
        return wD["instanceModelOutlierD"] if "instanceModelOutlierD" in wD else {}
