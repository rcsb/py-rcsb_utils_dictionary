##
# File:    DictMethodSecStructUtils.py
# Author:  J. Westbrook
# Date:    30-Sep-2021
# Version: 0.001 Initial version
#
# Updates:
# 08-Aug-2022  bv RO-3382 Map DSSP SS types to PROMOTIF SS types (temporary fix until DSSP can be run for all experimental structures)
#
##
"""
Helper class implements utilities to process secondary structure records.

"""
__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

# pylint: disable=too-many-lines

import itertools
import logging
from collections import OrderedDict

from rcsb.utils.io.CacheUtils import CacheUtils

logger = logging.getLogger(__name__)


class DictMethodSecStructUtils(object):
    """Helper class implements utilities to process secondary structure records struct_conf* categories."""

    dsspTypeNames = ["HELX_RH_AL_P", "STRN", "HELX_RH_3T_P", "HELX_RH_PI_P", "HELX_LH_PP_P", "TURN_TY1_P", "BEND"]

    def __init__(self, rpObj, **kwargs):
        """
        Args:
            resourceProvider: (obj) instance of DictMethodResourceProvider()
            raiseExceptions: (bool, optional) flag to raise rather than handle exceptions

        """
        #
        self.__commonU = rpObj.getResource("DictMethodCommonUtils instance") if rpObj else None
        self.__raiseExceptions = kwargs.get("raiseExceptions", False)
        #

        self.__promotifTypes = ["HELX_P"]
        cacheSize = 2
        self.__protSSCache = CacheUtils(size=cacheSize, label="protein secondary structure")
        self.__cisPeptideCache = CacheUtils(size=cacheSize, label="cis-peptide instances")
        #
        logger.debug("Dictionary secondary structure utilities")

    def testCache(self):
        return len(DictMethodSecStructUtils.dsspTypeNames) > 8

    def getCisPeptides(self, dataContainer):
        """Return a dictionary cis-peptides linkages using standard nomenclature.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {<id>: (begAsymId, begSeqId, endSeqId, modelId, omegaAngle), ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchCisPeptideFeatures(dataContainer)
        return wD["cisPeptideD"] if "cisPeptideD" in wD else {}

    def getProtSecStructFeatures(self, dataContainer, ssType):
        """Return a dictionary protein secondary structure features (entity/label sequence coordinates).

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance
            ssType (str): one of helix, sheet|strn, bend, turn

        Returns:
            dict: {<ss_id>: (asymId, begSeqId, endSeqId, confType, provCode), ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchProtSecStructFeatures(dataContainer)
        if ssType.lower() in ["helix"]:
            return wD["helixRangeD"] if "helixRangeD" in wD else {}
        elif ssType.lower() in ["sheet", "strn"]:
            # return wD["sheetRangeD"] if "sheetRangeD" in wD else {}
            return wD["instSheetRangeD"] if "instSheetRangeD" in wD else {}
        elif ssType.lower() in ["bend"]:
            return wD["bendRangeD"] if "bendRangeD" in wD else {}
        elif ssType.lower() in ["turn"]:
            return wD["turnRangeD"] if "turnRangeD" in wD else {}
        else:
            logger.error("Unknown secondary structure type %r", ssType)
        return {}

    def getProtHelixFeatures(self, dataContainer):
        """Return a dictionary protein helical features (entity/label sequence coordinates).

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {<helix_id>: (asymId, begSeqId, endSeqId, confType, provCode), ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchProtSecStructFeatures(dataContainer)
        return wD["helixRangeD"] if "helixRangeD" in wD else {}

    def getProtUnassignedSecStructFeatures(self, dataContainer):
        """Return a dictionary protein regions lacking SS feature assignments (entity/label sequence coordinates).

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {<id>: (asymId, begSeqId, endSeqId, confType, provCode), ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchProtSecStructFeatures(dataContainer)
        return wD["unassignedRangeD"] if "unassignedRangeD" in wD else {}

    def getProtUnassignedSecStructProvenance(self, dataContainer):
        """Return a dictionary of provenance and version details for unassigned secondary structure features.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {provenance: <provenanceString>, "version": <versionString>}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchProtSecStructFeatures(dataContainer)
        return wD["unassignedProvenanceD"] if "unassignedProvenanceD" in wD else {}

    def getProtSheetFeatures(self, dataContainer):
        """Return a dictionary protein beta strand features (entity/label sequence coordinates).

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {<sheet_id>: {asymId: [(begSeqId, endSeqId), ...], }
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchProtSecStructFeatures(dataContainer)
        return wD["instSheetRangeD"] if "instSheetRangeD" in wD else {}

    def getProtSheetSense(self, dataContainer):
        """Return a dictionary protein beta strand sense .

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {<sheet_id>: mixed|parallel|anti-parallel, ...}
        """
        if not dataContainer or not dataContainer.getName():
            return {}
        wD = self.__fetchProtSecStructFeatures(dataContainer)
        return wD["senseTypeD"] if "senseTypeD" in wD else {}

    def getProtSecStructFeaturesAll(self, dataContainer):
        """Return a dictionary of protein secondary structure features, counts, and coverage statistics.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            dict: {}
        """
        wD = self.__fetchProtSecStructFeatures(dataContainer)
        return wD

    def __fetchProtSecStructFeatures(self, dataContainer):
        wD = self.__protSSCache.get(dataContainer.getName())
        if not wD:
            flavor = self.__getFlavor(dataContainer)
            if flavor == "PROMOTIF":
                wD = self.__assembleProtSecStructFeaturesProMotif(dataContainer)
            elif flavor == "DSSP":
                wD = self.__assembleProtSecStructFeaturesDssp(dataContainer)
            else:
                wD = self.__initSecStructFeatures()
            #
            self.__protSSCache.set(dataContainer.getName(), wD)
        return wD

    def __getFlavor(self, dataContainer):
        if dataContainer.exists("struct_conf"):
            tObj = dataContainer.getObj("struct_conf")
            for ii in range(tObj.getRowCount()):
                confType = str(tObj.getValue("conf_type_id", ii)).strip().upper()
                if confType not in self.__promotifTypes:
                    return "DSSP"
            return "PROMOTIF"
        if dataContainer.exists("struct_sheet_range"):
            return "PROMOTIF"
        return None

    def __initSecStructFeatures(self):
        return {
            "helixCountD": {},
            "sheetStrandCountD": {},
            "unassignedCountD": {},
            "helixLengthD": {},
            "sheetStrandLengthD": {},
            "unassignedLengthD": {},
            "helixFracD": {},
            "sheetStrandFracD": {},
            "unassignedFracD": {},
            "sheetSenseD": {},
            "sheetFullStrandCountD": {},
            "featureMonomerSequenceD": {},
            "featureSequenceD": {},
            #
            "unassignedRangeD": {},
            "helixRangeD": {},
            "instHelixD": {},
            "sheetRangeD": {},
            "instSheetD": {},
            "senseTypeD": {},
        }

    def __assembleProtSecStructFeaturesProMotif(self, dataContainer):
        """Get secondary structure features using standard nomenclature .

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            (dict): with secondary structure details

            For instance, the following are calculated:
                     {
                        "helixCountD": {},
                        "sheetStrandCountD": {},
                        "unassignedCountD": {},
                        "helixLengthD": {},
                        "sheetStrandLengthD": {},
                        "unassignedLengthD": {},
                        "helixFracD": {},
                        "sheetStrandFracD": {},
                        "unassignedFracD": {},
                        "sheetSenseD": {},
                        "sheetFullStrandCountD": {},
                        "featureMonomerSequenceD": {},
                        "featureSequenceD": {},
                        #
                        "unassignedRangeD": {},
                        "helixRangeD": {},
                        "instHelixD": {},
                        "sheetRangeD": {},
                        "instSheetD": {},
                         "senseTypeD": {}
                    }

            # -- Target data categories ---
            loop_
            _struct_conf.conf_type_id
            _struct_conf.id
            _struct_conf.pdbx_PDB_helix_id
            _struct_conf.beg_label_comp_id
            _struct_conf.beg_label_asym_id
            _struct_conf.beg_label_seq_id
            _struct_conf.pdbx_beg_PDB_ins_code
            _struct_conf.end_label_comp_id
            _struct_conf.end_label_asym_id
            _struct_conf.end_label_seq_id
            _struct_conf.pdbx_end_PDB_ins_code

            _struct_conf.beg_auth_comp_id
            _struct_conf.beg_auth_asym_id
            _struct_conf.beg_auth_seq_id
            _struct_conf.end_auth_comp_id
            _struct_conf.end_auth_asym_id
            _struct_conf.end_auth_seq_id
            _struct_conf.pdbx_PDB_helix_class
            _struct_conf.details
            _struct_conf.pdbx_PDB_helix_length
            HELX_P HELX_P1 AA1 SER A 5   ? LYS A 19  ? SER A 2   LYS A 16  1 ? 15
            HELX_P HELX_P2 AA2 GLU A 26  ? LYS A 30  ? GLU A 23  LYS A 27  5 ? 5
            HELX_P HELX_P3 AA3 GLY A 47  ? LYS A 60  ? GLY A 44  LYS A 57  1 ? 14
            HELX_P HELX_P4 AA4 ASP A 111 ? LEU A 125 ? ASP A 108 LEU A 122 1 ? 15
            #
            _struct_conf_type.id          HELX_P
            _struct_conf_type.criteria    ?
            _struct_conf_type.reference   ?
            # -------------------------------------------------------------------

            loop_
            _struct_asym.id
            _struct_asym.pdbx_blank_PDB_chainid_flag
            _struct_asym.pdbx_modified
            _struct_asym.entity_id
            _struct_asym.details
            A N N 1 ?
            B N N 1 ?
            #
            _struct_sheet.id               A
            _struct_sheet.type             ?
            _struct_sheet.number_strands   8
            _struct_sheet.details          ?
            #
            loop_
            _struct_sheet_order.sheet_id
            _struct_sheet_order.range_id_1
            _struct_sheet_order.range_id_2
            _struct_sheet_order.offset
            _struct_sheet_order.sense
            A 1 2 ? anti-parallel
            A 2 3 ? anti-parallel
            A 3 4 ? anti-parallel
            A 4 5 ? anti-parallel
            A 5 6 ? anti-parallel
            A 6 7 ? anti-parallel
            A 7 8 ? anti-parallel
            #
            loop_
            _struct_sheet_range.sheet_id
            _struct_sheet_range.id
            _struct_sheet_range.beg_label_comp_id
            _struct_sheet_range.beg_label_asym_id
            _struct_sheet_range.beg_label_seq_id
            _struct_sheet_range.pdbx_beg_PDB_ins_code
            _struct_sheet_range.end_label_comp_id
            _struct_sheet_range.end_label_asym_id
            _struct_sheet_range.end_label_seq_id
            _struct_sheet_range.pdbx_end_PDB_ins_code

            _struct_sheet_range.beg_auth_comp_id
            _struct_sheet_range.beg_auth_asym_id
            _struct_sheet_range.beg_auth_seq_id
            _struct_sheet_range.end_auth_comp_id
            _struct_sheet_range.end_auth_asym_id
            _struct_sheet_range.end_auth_seq_id
            A 1 LYS A 5  ? VAL A 8  ? LYS A 5  VAL A 8
            A 2 ARG A 11 ? THR A 16 ? ARG A 11 THR A 16
            A 3 VAL A 19 ? LEU A 26 ? VAL A 19 LEU A 26
            A 4 TYR A 29 ? ALA A 35 ? TYR A 29 ALA A 35
            A 5 TYR B 29 ? ALA B 35 ? TYR B 29 ALA B 35
            A 6 VAL B 19 ? LEU B 26 ? VAL B 19 LEU B 26
            A 7 ARG B 11 ? THR B 16 ? ARG B 11 THR B 16
            A 8 LYS B 5  ? VAL B 8  ? LYS B 5  VAL B 8
            #
            _struct_mon_prot_cis.pdbx_id                1
            _struct_mon_prot_cis.label_comp_id          ASN
            _struct_mon_prot_cis.label_seq_id           189
            _struct_mon_prot_cis.label_asym_id          C
            _struct_mon_prot_cis.label_alt_id           .
            _struct_mon_prot_cis.pdbx_PDB_ins_code      ?
            _struct_mon_prot_cis.auth_comp_id           ASN
            _struct_mon_prot_cis.auth_seq_id            2007
            _struct_mon_prot_cis.auth_asym_id           2

            _struct_mon_prot_cis.pdbx_label_comp_id_2   PRO
            _struct_mon_prot_cis.pdbx_label_seq_id_2    190
            _struct_mon_prot_cis.pdbx_label_asym_id_2   C
            _struct_mon_prot_cis.pdbx_PDB_ins_code_2    ?
            _struct_mon_prot_cis.pdbx_auth_comp_id_2    PRO
            _struct_mon_prot_cis.pdbx_auth_seq_id_2     2008
            _struct_mon_prot_cis.pdbx_auth_asym_id_2    2

            _struct_mon_prot_cis.pdbx_PDB_model_num     1
            _struct_mon_prot_cis.pdbx_omega_angle       -6.45
        """

        rD = {
            "helixCountD": {},
            "sheetStrandCountD": {},
            "unassignedCountD": {},
            "helixLengthD": {},
            "sheetStrandLengthD": {},
            "unassignedLengthD": {},
            "helixFracD": {},
            "sheetStrandFracD": {},
            "unassignedFracD": {},
            "sheetSenseD": {},
            "sheetFullStrandCountD": {},
            "featureMonomerSequenceD": {},
            "featureSequenceD": {},
            #
            "unassignedRangeD": {},
            "helixRangeD": {},
            "instHelixD": {},
            "sheetRangeD": {},
            "instSheetD": {},
            "senseTypeD": {},
            "unassignedProvenanceD": {},
            "turnRangeD": {},
            "bendRangeD": {},
        }
        try:
            instancePolymerTypeD = self.__commonU.getInstancePolymerTypes(dataContainer)
            instEntityD = self.__commonU.getInstanceEntityMap(dataContainer)
            epLengthD = self.__commonU.getPolymerEntityLengths(dataContainer)
            #
            helixRangeD = {}
            sheetRangeD = {}
            sheetSenseD = {}
            unassignedRangeD = {}
            bendRangeD = {}
            turnRangeD = {}

            if dataContainer.exists("struct_conf"):
                tObj = dataContainer.getObj("struct_conf")
                helixRangeD = OrderedDict()
                for ii in range(tObj.getRowCount()):
                    confType = str(tObj.getValue("conf_type_id", ii)).strip().upper()
                    if confType in ["HELX_P"]:
                        hId = tObj.getValue("id", ii)
                        begAsymId = tObj.getValue("beg_label_asym_id", ii)
                        endAsymId = tObj.getValue("end_label_asym_id", ii)
                        try:
                            tbegSeqId = int(tObj.getValue("beg_label_seq_id", ii))
                            tendSeqId = int(tObj.getValue("end_label_seq_id", ii))
                            begSeqId = min(tbegSeqId, tendSeqId)
                            endSeqId = max(tbegSeqId, tendSeqId)
                        except Exception:
                            continue
                        #
                        if (begAsymId == endAsymId) and (begSeqId <= endSeqId):
                            helixRangeD.setdefault(hId, []).append((begAsymId, begSeqId, endSeqId, "HELIX_P", "PROMOTIF", "V1.0"))
                        else:
                            logger.debug("%s inconsistent struct_conf description id = %s", dataContainer.getName(), hId)

            if dataContainer.exists("struct_sheet_range"):
                tObj = dataContainer.getObj("struct_sheet_range")
                sheetRangeD = OrderedDict()
                for ii in range(tObj.getRowCount()):
                    sId = tObj.getValue("sheet_id", ii)
                    begAsymId = tObj.getValue("beg_label_asym_id", ii)
                    endAsymId = tObj.getValue("end_label_asym_id", ii)
                    # Most obsolete entries do no define this
                    try:
                        tbegSeqId = int(tObj.getValue("beg_label_seq_id", ii))
                        tendSeqId = int(tObj.getValue("end_label_seq_id", ii))
                        begSeqId = min(tbegSeqId, tendSeqId)
                        endSeqId = max(tbegSeqId, tendSeqId)
                    except Exception:
                        continue
                    if (begAsymId == endAsymId) and (begSeqId <= endSeqId):
                        sheetRangeD.setdefault(sId, []).append((begAsymId, begSeqId, endSeqId, "SHEET", "PROMOTIF", "V1.0"))
                    else:
                        logger.debug("%s inconsistent struct_sheet_range description id = %s", dataContainer.getName(), sId)

            logger.debug("%s sheetRangeD %r", dataContainer.getName(), sheetRangeD.items())
            #
            if dataContainer.exists("struct_sheet_order"):
                tObj = dataContainer.getObj("struct_sheet_order")
                #
                sheetSenseD = OrderedDict()
                for ii in range(tObj.getRowCount()):
                    sId = tObj.getValue("sheet_id", ii)
                    sense = str(tObj.getValue("sense", ii)).strip().lower()
                    sheetSenseD.setdefault(sId, []).append(sense)
            #
            logger.debug("%s sheetSenseD %r", dataContainer.getName(), sheetSenseD.items())
            # --------

            unassignedCoverageD = {}
            unassignedCountD = {}
            unassignedLengthD = {}
            unassignedFracD = {}

            helixCoverageD = {}
            helixCountD = {}
            helixLengthD = {}
            helixFracD = {}
            instHelixD = {}

            sheetCoverageD = {}
            sheetStrandCountD = {}
            sheetStrandLengthD = {}
            strandsPerBetaSheetD = {}
            sheetFullStrandCountD = {}
            sheetStrandFracD = {}
            instSheetD = {}
            instSheetSenseD = {}
            #
            featureMonomerSequenceD = {}
            featureSequenceD = {}
            #
            # ------------
            # Initialize over all protein instances
            for asymId, filteredType in instancePolymerTypeD.items():
                if filteredType != "Protein":
                    continue
                helixCoverageD[asymId] = []
                helixLengthD[asymId] = []
                helixCountD[asymId] = 0
                helixFracD[asymId] = 0.0
                instHelixD[asymId] = []
                #
                sheetCoverageD[asymId] = []
                sheetStrandCountD[asymId] = 0
                sheetStrandLengthD[asymId] = []
                sheetFullStrandCountD[asymId] = []
                sheetStrandFracD[asymId] = 0.0
                instSheetD[asymId] = []
                instSheetSenseD[asymId] = []
                #
                unassignedCountD[asymId] = 0
                unassignedLengthD[asymId] = []
                unassignedFracD[asymId] = 0.0
                #
                featureMonomerSequenceD[asymId] = None
                featureSequenceD[asymId] = None
            # -------------
            #
            for hId, hL in helixRangeD.items():
                for (asymId, begSeqId, endSeqId, _, _, _) in hL:
                    helixCoverageD.setdefault(asymId, []).extend(range(begSeqId, endSeqId + 1))
                    helixLengthD.setdefault(asymId, []).append(abs(begSeqId - endSeqId) + 1)
                    helixCountD[asymId] = helixCountD[asymId] + 1 if asymId in helixCountD else 0
                    instHelixD.setdefault(asymId, []).append(hId)
            #
            # ---------
            # betaSheetCount = len(sheetRangeD)
            #
            for sId, sL in sheetRangeD.items():
                strandsPerBetaSheetD[sId] = len(sL)
                for (asymId, begSeqId, endSeqId, _, _, _) in sL:
                    sheetCoverageD.setdefault(asymId, []).extend(range(begSeqId, endSeqId + 1))
                    sheetStrandLengthD.setdefault(asymId, []).append(abs(begSeqId - endSeqId) + 1)
                    sheetStrandCountD[asymId] = sheetStrandCountD[asymId] + 1 if asymId in sheetStrandCountD else 0
                    instSheetD.setdefault(asymId, []).append(sId)
            #
            instSheetRangeD = {}
            for sId, sL in sheetRangeD.items():
                aD = {}
                for (asymId, begSeqId, endSeqId, confType, provCode, provVer) in sL:
                    aD.setdefault(asymId, []).append((begSeqId, endSeqId, confType, provCode, provVer))
                instSheetRangeD[sId] = aD
            #
            # ---------
            senseTypeD = {}
            for sheetId, sL in sheetSenseD.items():
                if not sL:
                    continue
                usL = list(set(sL))
                if len(usL) == 1:
                    senseTypeD[sheetId] = usL[0]
                else:
                    senseTypeD[sheetId] = "mixed"
            # ---------
            #
            for asymId, filteredType in instancePolymerTypeD.items():
                logger.debug("%s processing %s type %r", dataContainer.getName(), asymId, filteredType)
                if filteredType != "Protein":
                    continue
                entityId = instEntityD[asymId]
                entityLen = epLengthD[entityId]
                entityS = set(range(1, entityLen + 1))
                eLen = len(entityS)
                #
                helixS = set(helixCoverageD[asymId])
                sheetS = set(sheetCoverageD[asymId])
                commonS = helixS & sheetS
                if commonS:
                    logger.debug("%s asymId %s overlapping secondary structure assignments for monomers %r", dataContainer.getName(), asymId, commonS)
                    # continue

                hLen = len(helixS) if asymId in helixCoverageD else 0
                sLen = len(sheetS) if asymId in sheetCoverageD else 0
                unassignedS = entityS - helixS if hLen else entityS
                unassignedS = unassignedS - sheetS if sLen else unassignedS
                tLen = len(unassignedS)
                #
                logger.debug("%s (%s) helix (%d) sheet (%d) unassigned (%d)", dataContainer.getName(), asymId, hLen, sLen, tLen)
                #
                # if eLen != hLen + sLen + tLen:
                #    logger.warning("%s overlapping secondary structure assignments for asymId %s", dataContainer.getName(), asymId)
                #    continue
                #
                unassignedCoverageD[asymId] = list(unassignedS)
                helixFracD[asymId] = float(hLen) / float(eLen)
                sheetStrandFracD[asymId] = float(sLen) / float(eLen)
                unassignedFracD[asymId] = float(tLen) / float(eLen)
                #
                unassignedRangeD[asymId] = list(self.__toRangeList(unassignedS))
                unassignedCountD[asymId] = len(unassignedRangeD[asymId])
                unassignedLengthD[asymId] = [abs(i - j) + 1 for (i, j) in unassignedRangeD[asymId]]
                #
                # ------
                sIdL = instSheetD[asymId]
                #
                instSheetSenseD[asymId] = [senseTypeD[sId] for sId in sIdL if sId in senseTypeD]
                sheetFullStrandCountD[asymId] = [strandsPerBetaSheetD[sId] for sId in sIdL if sId in strandsPerBetaSheetD]
                #

                # ------
                ssTypeL = ["_"] * eLen
                if hLen:
                    for idx in helixCoverageD[asymId]:
                        ssTypeL[idx - 1] = "H"
                if sLen:
                    for idx in sheetCoverageD[asymId]:
                        ssTypeL[idx - 1] = "S"
                if tLen:
                    for idx in unassignedCoverageD[asymId]:
                        ssTypeL[idx - 1] = "_"
                #
                featureMonomerSequenceD[asymId] = "".join(ssTypeL)
                featureSequenceD[asymId] = "".join([t[0] for t in itertools.groupby(ssTypeL)])
            # ---------
            unassignedProvenanceD = {"provenance": "PROMOTIF", "version": "V1.0"}
            rD = {
                "helixCountD": helixCountD,
                "sheetStrandCountD": sheetStrandCountD,
                "unassignedCountD": unassignedCountD,
                "helixLengthD": helixLengthD,
                "sheetStrandLengthD": sheetStrandLengthD,
                "unassignedLengthD": unassignedLengthD,
                "helixFracD": helixFracD,
                "sheetStrandFracD": sheetStrandFracD,
                "unassignedFracD": unassignedFracD,
                "sheetSenseD": instSheetSenseD,
                "sheetFullStrandCountD": sheetFullStrandCountD,
                "featureMonomerSequenceD": featureMonomerSequenceD,
                "featureSequenceD": featureSequenceD,
                #
                "unassignedRangeD": unassignedRangeD,
                "helixRangeD": helixRangeD,
                "instHelixD": instHelixD,
                # "sheetRangeD": sheetRangeD,
                "instSheetRangeD": instSheetRangeD,
                "instSheetD": instSheetD,
                "senseTypeD": senseTypeD,
                "unassignedProvenanceD": unassignedProvenanceD,
                "turnRangeD": turnRangeD,
                "bendRangeD": bendRangeD,
            }
        except Exception as e:
            logger.exception("Failing for %s with %s", dataContainer.getName(), str(e))
        #
        return rD

    def __toRangeList(self, iterable):
        iterable = sorted(set(iterable))
        for _, group in itertools.groupby(enumerate(iterable), lambda t: t[1] - t[0]):
            group = list(group)
            yield group[0][1], group[-1][1]

    def __fetchCisPeptideFeatures(self, dataContainer):
        wD = self.__cisPeptideCache.get(dataContainer.getName())
        if not wD:
            wD = self.__assembleCisPeptideFeatures(dataContainer)
            self.__cisPeptideCache.set(dataContainer.getName(), wD)
        return wD

    def __assembleCisPeptideFeatures(self, dataContainer):
        """Assemble cis-peptide instances from from the struct_mon_prot_cis category.

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            (dict): cis-peptide instances

            For instance, the following are calculated:
                     {"cisPeptideD": {} }

            # Example category:

            _struct_mon_prot_cis.pdbx_id                1
            _struct_mon_prot_cis.label_comp_id          ASN
            _struct_mon_prot_cis.label_seq_id           189
            _struct_mon_prot_cis.label_asym_id          C
            _struct_mon_prot_cis.label_alt_id           .
            _struct_mon_prot_cis.pdbx_PDB_ins_code      ?
            _struct_mon_prot_cis.auth_comp_id           ASN
            _struct_mon_prot_cis.auth_seq_id            2007
            _struct_mon_prot_cis.auth_asym_id           2

            _struct_mon_prot_cis.pdbx_label_comp_id_2   PRO
            _struct_mon_prot_cis.pdbx_label_seq_id_2    190
            _struct_mon_prot_cis.pdbx_label_asym_id_2   C
            _struct_mon_prot_cis.pdbx_PDB_ins_code_2    ?
            _struct_mon_prot_cis.pdbx_auth_comp_id_2    PRO
            _struct_mon_prot_cis.pdbx_auth_seq_id_2     2008
            _struct_mon_prot_cis.pdbx_auth_asym_id_2    2

            _struct_mon_prot_cis.pdbx_PDB_model_num     1
            _struct_mon_prot_cis.pdbx_omega_angle       -6.45
        """
        #
        rD = {"cisPeptideD": {}}
        try:
            cisPeptideD = OrderedDict()
            #
            if dataContainer.exists("struct_mon_prot_cis"):
                tObj = dataContainer.getObj("struct_mon_prot_cis")
                for ii in range(tObj.getRowCount()):
                    cId = tObj.getValue("pdbx_id", ii)
                    begAsymId = tObj.getValue("label_asym_id", ii)
                    # begCompId = tObj.getValue("label_comp_id", ii)
                    begSeqId = int(tObj.getValue("label_seq_id", ii))
                    endAsymId = tObj.getValue("pdbx_label_asym_id_2", ii)
                    # endCompId = int(tObj.getValue("pdbx_label_comp_id_2", ii))
                    endSeqId = int(tObj.getValue("pdbx_label_seq_id_2", ii))
                    modelId = int(tObj.getValue("pdbx_PDB_model_num", ii))
                    omegaAngle = float(tObj.getValue("pdbx_omega_angle", ii))
                    #
                    if (begAsymId == endAsymId) and (begSeqId <= endSeqId):
                        cisPeptideD.setdefault(cId, []).append((begAsymId, begSeqId, endSeqId, modelId, omegaAngle))
                    else:
                        logger.debug("%s inconsistent cis peptide description id = %s", dataContainer.getName(), cId)

            rD = {"cisPeptideD": cisPeptideD}
        except Exception as e:
            logger.exception("Failing for %s with %s", dataContainer.getName(), str(e))
        #
        return rD

    def __assembleProtSecStructFeaturesDssp(self, dataContainer):
        """Get secondary structure features using standard nomenclature (DSSP).

        Args:
            dataContainer (object):  mmcif.api.DataContainer object instance

        Returns:
            (dict): with secondary structure details

            For instance, the following are calculated:
                     {
                        "helixCountD": {},
                        "sheetStrandCountD": {},
                        "unassignedCountD": {},
                        "helixLengthD": {},
                        "sheetStrandLengthD": {},
                        "unassignedLengthD": {},
                        "helixFracD": {},
                        "sheetStrandFracD": {},
                        "unassignedFracD": {},
                        "sheetSenseD": {},
                        "sheetFullStrandCountD": {},
                        "featureMonomerSequenceD": {},
                        "featureSequenceD": {},
                        #
                        "unassignedRangeD": {},
                        "helixRangeD": {},
                        "instHelixD": {},
                        "sheetRangeD": {},
                        "instSheetD": {},
                         "senseTypeD": {}
                    }

            # -- Target data categories for DSSP ---
            #
            loop_
            _struct_conf.beg_auth_asym_id
            _struct_conf.beg_auth_comp_id
            _struct_conf.beg_auth_seq_id
            _struct_conf.beg_label_asym_id
            _struct_conf.beg_label_comp_id
            _struct_conf.beg_label_seq_id
            _struct_conf.conf_type_id
            _struct_conf.end_auth_asym_id
            _struct_conf.end_auth_comp_id
            _struct_conf.end_auth_seq_id
            _struct_conf.end_label_asym_id
            _struct_conf.end_label_comp_id
            _struct_conf.end_label_seq_id
            _struct_conf.id
            _struct_conf.pdbx_beg_PDB_ins_code
            _struct_conf.pdbx_end_PDB_ins_code
            A LYS 0 A LYS 2   STRN         A GLY 0 A GLY 7   STRN1          ? ?
            A LYS 0 A LYS 8   BEND         A LYS 0 A LYS 8   BEND1          ? ?
            A GLY 0 A GLY 9   TURN_TY1_P   A GLY 0 A GLY 10  TURN_TY1_P1    ? ?
            A CYS 0 A CYS 11  BEND         A GLY 0 A GLY 12  BEND2          ? ?
            A LYS 0 A LYS 13  HELX_RH_AL_P A LYS 0 A LYS 26  HELX_RH_AL_P1  ? ?
            A LYS 0 A LYS 27  TURN_TY1_P   A GLY 0 A GLY 28  TURN_TY1_P2    ? ?
            A ASN 0 A ASN 30  STRN         A GLY 0 A GLY 36  STRN2          ? ?
            # ....
            loop_
            _struct_conf_type.criteria
            _struct_conf_type.id
            DSSP STRN
            DSSP BEND
            DSSP TURN_TY1_P
            DSSP HELX_RH_AL_P
            DSSP HELX_LH_PP_P
            DSSP HELX_RH_3T_P
            #
            # -------------------------------------------------------------------


        """

        dsspTypeMapD = {
            "HELX_RH_AL_P": "HELIX_P",
            "STRN": "SHEET",
            "HELX_RH_3T_P": "HELIX_P",
            "HELX_RH_PI_P": "HELIX_P",
            "HELX_LH_PP_P": "HELIX_P",
            "TURN_TY1_P": "TURN_TY1_P",
            "BEND": "BEND"
        }

        rD = {
            "helixCountD": {},
            "sheetStrandCountD": {},
            "unassignedCountD": {},
            "helixLengthD": {},
            "sheetStrandLengthD": {},
            "unassignedLengthD": {},
            "helixFracD": {},
            "sheetStrandFracD": {},
            "unassignedFracD": {},
            "sheetSenseD": {},
            "sheetFullStrandCountD": {},
            "featureMonomerSequenceD": {},
            "featureSequenceD": {},
            #
            "unassignedRangeD": {},
            "helixRangeD": {},
            "instHelixD": {},
            "sheetRangeD": {},
            "instSheetD": {},
            "senseTypeD": {},
            "unassignedProvenanceD": {},
            "turnRangeD": {},
            "bendRangeD": {},
        }
        try:
            instancePolymerTypeD = self.__commonU.getInstancePolymerTypes(dataContainer)
            instEntityD = self.__commonU.getInstanceEntityMap(dataContainer)
            epLengthD = self.__commonU.getPolymerEntityLengths(dataContainer)
            #
            helixRangeD = OrderedDict()
            sheetRangeD = OrderedDict()
            turnRangeD = OrderedDict()
            bendRangeD = OrderedDict()
            # sheetSenseD = {}
            unassignedRangeD = {}

            if dataContainer.exists("struct_conf"):
                tObj = dataContainer.getObj("struct_conf")
                numH = 0
                numS = 0
                numB = 0
                numT = 0
                for ii in range(tObj.getRowCount()):
                    confType = str(tObj.getValue("conf_type_id", ii)).strip().upper()
                    if confType in DictMethodSecStructUtils.dsspTypeNames:
                        ssId = tObj.getValue("id", ii)
                        begAsymId = tObj.getValue("beg_label_asym_id", ii)
                        endAsymId = tObj.getValue("end_label_asym_id", ii)
                        try:
                            tbegSeqId = int(tObj.getValue("beg_label_seq_id", ii))
                            tendSeqId = int(tObj.getValue("end_label_seq_id", ii))
                            begSeqId = min(tbegSeqId, tendSeqId)
                            endSeqId = max(tbegSeqId, tendSeqId)
                        except Exception:
                            continue
                        #

                        if confType.startswith("HELX") and (begAsymId == endAsymId) and (begSeqId <= endSeqId):
                            numH += 1
                            sId = dsspTypeMapD[confType] + str(numH)
                            helixRangeD.setdefault(sId, []).append((begAsymId, begSeqId, endSeqId, dsspTypeMapD[confType], "DSSP", "4"))
                        else:
                            logger.debug("%s inconsistent struct_conf description id = %s", dataContainer.getName(), ssId)
                        #
                        if confType.startswith("STRN") and (begAsymId == endAsymId) and (begSeqId <= endSeqId):
                            numS += 1
                            sId = dsspTypeMapD[confType] + str(numS)
                            sheetRangeD.setdefault(sId, []).append((begAsymId, begSeqId, endSeqId, dsspTypeMapD[confType], "DSSP", "4"))
                        else:
                            logger.debug("%s inconsistent struct_conf description id = %s", dataContainer.getName(), ssId)

                        if confType.startswith("BEND") and (begAsymId == endAsymId) and (begSeqId <= endSeqId):
                            numB += 1
                            sId = dsspTypeMapD[confType] + str(numB)
                            bendRangeD.setdefault(sId, []).append((begAsymId, begSeqId, endSeqId, dsspTypeMapD[confType], "DSSP", "4"))
                        else:
                            logger.debug("%s inconsistent struct_conf description id = %s", dataContainer.getName(), ssId)

                        if confType.startswith("TURN") and (begAsymId == endAsymId) and (begSeqId <= endSeqId):
                            numT += 1
                            sId = dsspTypeMapD[confType] + str(numT)
                            turnRangeD.setdefault(sId, []).append((begAsymId, begSeqId, endSeqId, dsspTypeMapD[confType], "DSSP", "4"))
                        else:
                            logger.debug("%s inconsistent struct_conf description id = %s", dataContainer.getName(), ssId)

            # --------
            unassignedCoverageD = {}
            unassignedCountD = {}
            unassignedLengthD = {}
            unassignedFracD = {}

            helixCoverageD = {}
            helixCountD = {}
            helixLengthD = {}
            helixFracD = {}
            instHelixD = {}

            bendCoverageD = {}
            bendCountD = {}
            bendLengthD = {}
            bendFracD = {}
            instBendD = {}

            turnCoverageD = {}
            turnCountD = {}
            turnLengthD = {}
            turnFracD = {}
            instTurnD = {}

            sheetCoverageD = {}
            sheetStrandCountD = {}
            sheetStrandLengthD = {}
            strandsPerBetaSheetD = {}
            sheetFullStrandCountD = {}
            sheetStrandFracD = {}
            instSheetD = {}
            # instSheetSenseD = {}
            #
            featureMonomerSequenceD = {}
            featureSequenceD = {}
            #
            # ------------
            # Initialize over all protein instances
            for asymId, filteredType in instancePolymerTypeD.items():
                if filteredType != "Protein":
                    continue
                helixCoverageD[asymId] = []
                helixLengthD[asymId] = []
                helixCountD[asymId] = 0
                helixFracD[asymId] = 0.0
                instHelixD[asymId] = []
                #
                bendCoverageD[asymId] = []
                bendLengthD[asymId] = []
                bendCountD[asymId] = 0
                bendFracD[asymId] = 0.0
                instBendD[asymId] = []
                #
                turnCoverageD[asymId] = []
                turnLengthD[asymId] = []
                turnCountD[asymId] = 0
                turnFracD[asymId] = 0.0
                instTurnD[asymId] = []
                #
                sheetCoverageD[asymId] = []
                sheetStrandCountD[asymId] = 0
                sheetStrandLengthD[asymId] = []
                sheetFullStrandCountD[asymId] = []
                sheetStrandFracD[asymId] = 0.0
                instSheetD[asymId] = []
                # instSheetSenseD[asymId] = []
                #
                unassignedCountD[asymId] = 0
                unassignedLengthD[asymId] = []
                unassignedFracD[asymId] = 0.0
                #
                featureMonomerSequenceD[asymId] = None
                featureSequenceD[asymId] = None
            # -------------
            #
            for hId, hL in helixRangeD.items():
                for (asymId, begSeqId, endSeqId, _, _, _) in hL:
                    helixCoverageD.setdefault(asymId, []).extend(range(begSeqId, endSeqId + 1))
                    helixLengthD.setdefault(asymId, []).append(abs(begSeqId - endSeqId) + 1)
                    helixCountD[asymId] = helixCountD[asymId] + 1 if asymId in helixCountD else 0
                    instHelixD.setdefault(asymId, []).append(hId)

            for bId, bL in bendRangeD.items():
                for (asymId, begSeqId, endSeqId, _, _, _) in bL:
                    bendCoverageD.setdefault(asymId, []).extend(range(begSeqId, endSeqId + 1))
                    bendLengthD.setdefault(asymId, []).append(abs(begSeqId - endSeqId) + 1)
                    bendCountD[asymId] = bendCountD[asymId] + 1 if asymId in bendCountD else 0
                    instBendD.setdefault(asymId, []).append(bId)

            for tId, tL in turnRangeD.items():
                for (asymId, begSeqId, endSeqId, _, _, _) in tL:
                    turnCoverageD.setdefault(asymId, []).extend(range(begSeqId, endSeqId + 1))
                    turnLengthD.setdefault(asymId, []).append(abs(begSeqId - endSeqId) + 1)
                    turnCountD[asymId] = turnCountD[asymId] + 1 if asymId in turnCountD else 0
                    instTurnD.setdefault(asymId, []).append(tId)

            for sId, sL in sheetRangeD.items():
                strandsPerBetaSheetD[sId] = len(sL)
                for (asymId, begSeqId, endSeqId, _, _, _) in sL:
                    sheetCoverageD.setdefault(asymId, []).extend(range(begSeqId, endSeqId + 1))
                    sheetStrandLengthD.setdefault(asymId, []).append(abs(begSeqId - endSeqId) + 1)
                    sheetStrandCountD[asymId] = sheetStrandCountD[asymId] + 1 if asymId in sheetStrandCountD else 0
                    instSheetD.setdefault(asymId, []).append(sId)
            #
            instSheetRangeD = {}
            for sId, sL in sheetRangeD.items():
                aD = {}
                for (asymId, begSeqId, endSeqId, confType, provCode, provVer) in sL:
                    aD.setdefault(asymId, []).append((begSeqId, endSeqId, confType, provCode, provVer))
                instSheetRangeD[sId] = aD
            #
            # ---------
            senseTypeD = {}
            # ---------
            #
            for asymId, filteredType in instancePolymerTypeD.items():
                logger.debug("%s processing %s type %r", dataContainer.getName(), asymId, filteredType)
                if filteredType != "Protein":
                    continue
                entityId = instEntityD[asymId]
                entityLen = epLengthD[entityId]
                entityS = set(range(1, entityLen + 1))
                eLen = len(entityS)
                #
                helixS = set(helixCoverageD[asymId])
                sheetS = set(sheetCoverageD[asymId])
                bendS = set(bendCoverageD[asymId])
                turnS = set(turnCoverageD[asymId])
                commonS = helixS & sheetS & bendS & turnS
                if commonS:
                    logger.info("%s asymId %s overlapping secondary structure assignments for monomers %r", dataContainer.getName(), asymId, commonS)
                    # continue

                hLen = len(helixS) if asymId in helixCoverageD else 0
                sLen = len(sheetS) if asymId in sheetCoverageD else 0
                turnLen = len(turnS) if asymId in turnCoverageD else 0
                bendLen = len(bendS) if asymId in bendCoverageD else 0
                #
                unassignedS = entityS - helixS if hLen else entityS
                unassignedS = unassignedS - sheetS if sLen else unassignedS
                unassignedS = unassignedS - turnS if turnLen else unassignedS
                unassignedS = unassignedS - bendS if bendLen else unassignedS
                uLen = len(unassignedS)
                #
                # if eLen != hLen + sLen + turnLen + bendLen + uLen:
                #    logger.warning("%s overlapping secondary structure assignments for asymId %s", dataContainer.getName(), asymId)
                #    continue
                #
                unassignedCoverageD[asymId] = list(unassignedS)
                helixFracD[asymId] = float(hLen) / float(eLen)
                sheetStrandFracD[asymId] = float(sLen) / float(eLen)
                unassignedFracD[asymId] = float(uLen) / float(eLen)
                #
                unassignedRangeD[asymId] = list(self.__toRangeList(unassignedS))
                unassignedCountD[asymId] = len(unassignedRangeD[asymId])
                unassignedLengthD[asymId] = [abs(i - j) + 1 for (i, j) in unassignedRangeD[asymId]]

                #
                # ------
                sIdL = instSheetD[asymId]
                #
                # instSheetSenseD[asymId] = [senseTypeD[sId] for sId in sIdL if sId in senseTypeD]
                sheetFullStrandCountD[asymId] = [strandsPerBetaSheetD[sId] for sId in sIdL if sId in strandsPerBetaSheetD]
                # ------
                ssTypeL = ["_"] * eLen
                if hLen:
                    for idx in helixCoverageD[asymId]:
                        ssTypeL[idx - 1] = "H"
                if sLen:
                    for idx in sheetCoverageD[asymId]:
                        ssTypeL[idx - 1] = "S"
                if bendLen:
                    for idx in bendCoverageD[asymId]:
                        ssTypeL[idx - 1] = "B"
                if turnLen:
                    for idx in turnCoverageD[asymId]:
                        ssTypeL[idx - 1] = "T"
                if uLen:
                    for idx in unassignedCoverageD[asymId]:
                        ssTypeL[idx - 1] = "_"
                #
                featureMonomerSequenceD[asymId] = "".join(ssTypeL)
                featureSequenceD[asymId] = "".join([t[0] for t in itertools.groupby(ssTypeL)])
            # ---------
            unassignedProvenanceD = {"provenance": "DSSP", "version": "V4"}
            rD = {
                "helixCountD": helixCountD,
                "sheetStrandCountD": sheetStrandCountD,
                "unassignedCountD": unassignedCountD,
                "helixLengthD": helixLengthD,
                "sheetStrandLengthD": sheetStrandLengthD,
                "unassignedLengthD": unassignedLengthD,
                "helixFracD": helixFracD,
                "sheetStrandFracD": sheetStrandFracD,
                "unassignedFracD": unassignedFracD,
                # "sheetSenseD": instSheetSenseD,
                "sheetFullStrandCountD": sheetFullStrandCountD,
                "featureMonomerSequenceD": featureMonomerSequenceD,
                "featureSequenceD": featureSequenceD,
                #
                "unassignedRangeD": unassignedRangeD,
                "helixRangeD": helixRangeD,
                "instHelixD": instHelixD,
                # "sheetRangeD": sheetRangeD,
                "instSheetRangeD": instSheetRangeD,
                "instSheetD": instSheetD,
                "senseTypeD": senseTypeD,
                "unassignedProvenanceD": unassignedProvenanceD,
                "turnRangeD": turnRangeD,
                "bendRangeD": bendRangeD,
            }

        except Exception as e:
            logger.exception("Failing for %s with %s", dataContainer.getName(), str(e))
        #
        return rD
