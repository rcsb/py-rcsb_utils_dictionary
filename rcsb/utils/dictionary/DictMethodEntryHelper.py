##
# File:    DictMethodEntryHelper.py (DictMethodRunnerHelper.py)
# Author:  J. Westbrook
# Date:    18-Aug-2018
# Version: 0.001 Initial version
#
#
# Updates:
#  4-Sep-2018 jdw add methods to construct entry and entity identier categories.
# 10-Sep-2018 jdw add method for citation author aggregation
# 22-Sep-2018 jdw add method assignAssemblyCandidates()
# 27-Oct-2018 jdw add method consolidateAccessionDetails()
# 30-Oct-2018 jdw add category methods addChemCompRelated(), addChemCompInfo(),
#                 addChemCompDescriptor()
# 10-Nov-2018 jdw add addChemCompSynonyms(), addChemCompTargets(), filterBlockByMethod()
# 12-Nov-2018 jdw add InChIKey matching in addChemCompRelated()
# 15-Nov-2018 jdw add handling for antibody misrepresentation of multisource organisms
# 28-Nov-2018 jdw relax constraints on the production of rcsb_entry_info
#  1-Dec-2018 jdw add ncbi source and host organism info
# 11-Dec-2018 jdw add addStructRefSeqEntityIds and buildEntityPolySeq
# 10-Jan-2019 jdw better handle initialization in filterBlockByMethod()
# 11-Jan-2019 jdw revise classification in assignAssemblyCandidates()
# 16-Feb-2019 jdw add buildContainerEntityInstanceIds()
# 19-Feb-2019 jdw add internal method __addPdbxValidateAsymIds() to add cardinal identifiers to
#                 pdbx_validate_* categories
# 28-Feb-2019 jdw change criteria for adding rcsb_chem_comp_container_identifiers to work with ion definitions
# 11-Mar-2019 jdw replace taxonomy file handling with calls to TaxonomyUtils()
# 11-Mar-2019 jdw add EC lineage using EnzymeDatabaseUtils()
# 17-Mar-2019 jdw add support for entity subcategory rcsb_macromolecular_names_combined
# 23-Mar-2019 jdw change criteria chem_comp collection criteria to _chem_comp.pdbx_release_status
# 25-Mar-2019 jdw remap merged taxons and adjust exception handling for taxonomy lineage generation
#  7-Apr-2019 jdw add CathClassificationUtils and CathClassificationUtils and sequence difference type counts
# 25-Apr-2019 jdw For source and host organism add ncbi_parent_scientific_name
#                 add rcsb_entry_info.deposited_modeled_polymer_monomer_count and
#                     rcsb_entry_info.deposited_unmodeled_polymer_monomer_count,
#  1-May-2019 jdw add support for _rcsb_entry_info.deposited_polymer_monomer_count,
#                   _rcsb_entry_info.polymer_entity_count_protein,
#                   _rcsb_entry_info.polymer_entity_count_nucleic_acid,
#                   _rcsb_entry_info.polymer_entity_count_nucleic_acid_hybrid,
#                   _rcsb_entry_info.polymer_entity_count_DNA,
#                   _rcsb_entry_info.polymer_entity_count_RNA,
#                   _rcsb_entry_info.nonpolymer_ligand_entity_count
#                   _rcsb_entry_info.selected_polymer_entity_types
#                   _rcsb_entry_info.polymer_entity_taxonomy_count
#                   _rcsb_entry_info.assembly_count
#                    add categories rcsb_entity_instance_domain_scop and rcsb_entity_instance_domain_cath
#  4-May-2019 jdw extend content in categories rcsb_entity_instance_domain_scop and rcsb_entity_instance_domain_cath
# 13-May-2019 jdw add rcsb_entry_info.deposited_polymer_entity_instance_count and deposited_nonpolymer_entity_instance_count
#                 add entity_poly.rcsb_non_std_monomer_count and rcsb_non_std_monomers
# 15-May-2019 jdw add _rcsb_entry_info.na_polymer_entity_types update enumerations for _rcsb_entry_info.selected_polymer_entity_types
# 19-May-2019 jdw add method __getStructConfInfo()
# 21-May-2019 jdw handle odd ordering of records in struct_ref_seq_dif.
# 25-Nov-2019 jdw add method normalizeCitationJournalAbbrev() and dependencies
# 11-Mar-2022 bv Fix _rcsb_entry_info.deposited_model_count not being populated for certain NMR entries
# 28-Mar-2022 bv Move 'getRepresentativeModels' method to DictMethodCommonUtils
# 26-Apr-2022 bv Add missing pdbx_database_status for MA models and _rcsb_entry_info.structure_determination_methodology
# 29-Apr-2022 dwp Use internal computed-model identifiers for 'rcsb_id'
# 30 Apr-2022 bv Update consolidateAccessionDetails
#  3-May-2022 dwp Use internal computed-model identifiers for 'entry_id' in containter_identifiers
# 29-Jun-2022 dwp Use internal computed-model identifiers everywhere (in same manner as experimental models)
# 06-Jul-2022 dwp Add addtional filters for populating _rcsb_accession_info
# 01-Aug-2022 dwp Override rcsb_entry_info.structure_determination_methodology if struct.pdbx_structure_determination_methodology is '?' or '.' in the CIF file
# 08-Aug-2022  bv Set values for rcsb_entry_info.structure_determination_methodology_priority
# 03-Oct-2022  bv Set values for rcsb_entry_info.ndb_struct_conf_na_feature_combined
# 03-Jan-2023  bv Include _pdbx_database_status.status_code_nmr_data for experimental data availability
# 26-Jan-2023 dwp Populate or update pdbx_database_status attributes for CSMs to make ready for RELease
# 21-Feb-2023  bv Update '__filterExperimentalResolution' method to handle experimental resolutions properly (see RO-3559)
# 01-Feb-2024  bv Update method 'addEntryInfo' to support deuterated water molecule count
# 16-Jan-2025 dwp Use simplified method call for getting representative model ID
# 03-Feb-2025  bv Add method 'filterRevisionHistory' to remove data not relevant to structure model
#
##
"""
Helper class implements entry-level method references in the RCSB dictionary extension.

All data accessors and structures here refer to dictionary category and attribute names.

"""
__docformat__ = "google en"
__author__ = "John Westbrook"
__email__ = "jwest@rcsb.rutgers.edu"
__license__ = "Apache 2.0"

# pylint: disable=too-many-lines

import logging
from string import capwords

from mmcif.api.DataCategory import DataCategory
from rcsb.utils.dictionary.DictMethodSecStructUtils import DictMethodSecStructUtils

logger = logging.getLogger(__name__)


def cmpElements(lhs, rhs):
    return 0 if (lhs[-1].isdigit() or lhs[-1] in ["R", "S"]) and rhs[0].isdigit() else -1


class DictMethodEntryHelper(object):
    """Helper class implements entry-level method references in the RCSB dictionary extension."""

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
        self.__crP = rP.getResource("CitationReferenceProvider instance") if rP else None
        self.__jtaP = rP.getResource("JournalTitleAbbreviationProvider instance") if rP else None
        #
        self.__ssU = DictMethodSecStructUtils(rP, raiseExceptions=self._raiseExceptions)
        # logger.debug("Dictionary entry method helper init")

    def echo(self, msg):
        logger.info(msg)

    def deferredItemMethod(self, dataContainer, catName, atName, **kwargs):
        """Placeholder for an item method."""
        _ = kwargs
        logger.debug("Called deferred item method %r %r for %r", catName, atName, dataContainer.getName())
        return True

    def deferredCategoryMethod(self, dataContainer, catName, **kwargs):
        """Placeholder for a category method."""
        _ = kwargs
        logger.debug("Called deferred category method %r for %r", catName, dataContainer.getName())
        return True

    def setDatablockId(self, dataContainer, catName, atName, **kwargs):
        """Item-level method to set the value of the input item to the current container name.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name
            atName (str): Attribute name

        Returns:
            bool: True for success or False otherwise
        """
        logger.debug("Starting catName %s atName %s kwargs %r", catName, atName, kwargs)
        try:
            val = dataContainer.getName()
            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=[atName]))
            #
            cObj = dataContainer.getObj(catName)
            if not cObj.hasAttribute(atName):
                cObj.appendAttribute(atName)
            #
            rc = cObj.getRowCount()
            numRows = rc if rc else 1
            for ii in range(numRows):
                cObj.setValue(val, atName, ii)
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def setLoadDateTime(self, dataContainer, catName, atName, **kwargs):
        """Set the value of the input data item with container load date.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name
            atName (str): Attribute name

        Returns:
            bool: True for success or False otherwise
        """
        logger.debug("Starting catName %s atName %s kwargs %r", catName, atName, kwargs)
        try:
            val = dataContainer.getProp("load_date")
            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=[atName]))
            #
            cObj = dataContainer.getObj(catName)
            if not cObj.hasAttribute(atName):
                cObj.appendAttribute(atName)
            #
            rc = cObj.getRowCount()
            numRows = rc if rc else 1
            for ii in range(numRows):
                cObj.setValue(val, atName, ii)
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def setLocator(self, dataContainer, catName, atName, **kwargs):
        """Set the value of the input data item with container locator path.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name
            atName (str): Attribute name

        Returns:
            bool: True for success or False otherwise
        """
        logger.debug("Starting catName %s atName %s kwargs %r", catName, atName, kwargs)
        try:
            val = dataContainer.getProp("locator")
            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=[atName]))
            #
            cObj = dataContainer.getObj(catName)
            if not cObj.hasAttribute(atName):
                cObj.appendAttribute(atName)
            #
            rc = cObj.getRowCount()
            numRows = rc if rc else 1
            for ii in range(numRows):
                cObj.setValue(val, atName, ii)
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def setRowIndex(self, dataContainer, catName, atName, **kwargs):
        """Set the values of the input data item with the category row index.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name
            atName (str): Attribute name

        Returns:
            bool: True for success or False otherwise
        """
        logger.debug("Starting catName %s atName %s kwargs %r", catName, atName, kwargs)
        try:
            if not dataContainer.exists(catName):
                # exit if there is no category to index
                return False
            #
            cObj = dataContainer.getObj(catName)
            if not cObj.hasAttribute(atName):
                cObj.appendAttribute(atName)
            #
            rc = cObj.getRowCount()
            numRows = rc if rc else 1
            for ii, iRow in enumerate(range(numRows), 1):
                # Note - we set the integer value as a string  -
                cObj.setValue(str(ii), atName, iRow)
            return True
        except Exception as e:
            logger.exception("Failing with %s", str(e))
        return False

    def aggregateCitationOrcidIdentifiers(self, dataContainer, catName, atName, **kwargs):
        """Set the value of the input data item with list of citation authors.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name
            atName (str): Attribute name

        Returns:
            bool: True for success or False otherwise
        """
        logger.debug("Starting catName %s atName %s kwargs %r", catName, atName, kwargs)
        try:
            if not dataContainer.exists(catName) or not dataContainer.exists("citation_author"):
                return False
            #
            cObj = dataContainer.getObj(catName)
            if not cObj.hasAttribute(atName):
                cObj.appendAttribute(atName)
            citIdL = cObj.getAttributeValueList("id")
            #
            tObj = dataContainer.getObj("citation_author")
            #

            citIdL = list(set(citIdL))
            tD = {}
            for ii, citId in enumerate(citIdL):
                if tObj.hasAttribute("identifier_ORCID"):
                    tD[citId] = tObj.selectValuesWhere("identifier_ORCID", citId, "citation_id")
                else:
                    tD[citId] = []
            for ii in range(cObj.getRowCount()):
                citId = cObj.getValue("id", ii)
                if tD[citId]:
                    cObj.setValue(",".join(tD[citId]), atName, ii)
                else:
                    cObj.setValue("?", atName, ii)
            return True
        except Exception as e:
            logger.exception("Failing for %r with %s", dataContainer.getName(), str(e))
        return False

    def aggregateCitationAuthors(self, dataContainer, catName, atName, **kwargs):
        """Set the value of the input data item with list of citation authors.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name
            atName (str): Attribute name

        Returns:
            bool: True for success or False otherwise
        """
        logger.debug("Starting catName %s atName %s kwargs %r", catName, atName, kwargs)
        try:
            if not dataContainer.exists(catName) or not dataContainer.exists("citation_author"):
                return False
            #
            cObj = dataContainer.getObj(catName)
            if not cObj.hasAttribute(atName):
                cObj.appendAttribute(atName)
            citIdL = cObj.getAttributeValueList("id")
            #
            tObj = dataContainer.getObj("citation_author")
            #
            citIdL = list(set(citIdL))
            tD = {}
            for ii, citId in enumerate(citIdL):
                tD[citId] = tObj.selectValuesWhere("name", citId, "citation_id")
            for ii in range(cObj.getRowCount()):
                citId = cObj.getValue("id", ii)
                cObj.setValue("|".join(tD[citId]), atName, ii)
            return True
        except Exception as e:
            logger.exception("Failing for %r with %s", dataContainer.getName(), str(e))
        return False

    def normalizeCitationJournalAbbrev(self, dataContainer, catName, atName, **kwargs):
        """Normalize citation journal abbrev.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name
            atName (str): Attribute name

        Returns:
            bool: True for success or False otherwise
        """
        logger.debug("Starting catName %s atName %s kwargs %r", catName, atName, kwargs)
        revAbbrev = None
        try:
            if not dataContainer.exists(catName):
                return False
            #
            cObj = dataContainer.getObj(catName)
            if not cObj.hasAttribute("journal_abbrev") or not cObj.hasAttribute("id") or not cObj.hasAttribute("journal_id_ISSN"):
                return False
            #
            if not cObj.hasAttribute(atName):
                cObj.appendAttribute(atName)
            #
            rcsbId = dataContainer.getName()
            for ii in range(cObj.getRowCount()):
                # citId = cObj.getValue("id", ii)
                issn = cObj.getValueOrDefault("journal_id_ISSN", ii, defaultValue=None)
                curAbbrev = cObj.getValueOrDefault("journal_abbrev", ii, defaultValue=None)
                if curAbbrev:
                    revAbbrev = self.__updateJournalAbbreviation(rcsbId, issn, curAbbrev)
                revAbbrev = revAbbrev if revAbbrev else curAbbrev
                #
                logger.debug("%s journal abbreviation issn %r current %r normalized %r", rcsbId, issn, curAbbrev, revAbbrev)
                cObj.setValue(revAbbrev, atName, ii)
            return True
        except Exception as e:
            logger.exception("Failing for %r with %s", dataContainer.getName(), str(e))
        return False

    def __updateJournalAbbreviation(self, rcsbId, issn, curAbbrev):
        revAbbrev = None
        try:
            if issn:
                medlineAbbrev = self.__crP.getMedlineJournalAbbreviation(issn)
                # medlineIsoAbbrev = self.__crP.getMedlineJournalIsoAbbreviation(issn)
                crIssn = issn.replace("-", "")
                crTitle = self.__crP.getCrossRefJournalTitle(crIssn)
                #
                revAbbrev = medlineAbbrev
                if not medlineAbbrev and not crTitle:
                    logger.debug("%s: missing information for issn %r curAbbrev %r", rcsbId, issn, curAbbrev)
                    revAbbrev = capwords(curAbbrev.replace(".", " "))
                elif not medlineAbbrev:
                    revAbbrev = self.__jtaP.getJournalAbbreviation(crTitle, usePunctuation=False)
            else:
                if curAbbrev.upper() in ["TO BE PUBLISHED", "IN PREPARATION"]:
                    revAbbrev = "To be published"
                elif curAbbrev.upper().startswith("THESIS"):
                    revAbbrev = "Thesis"
                else:
                    revAbbrev = capwords(curAbbrev.replace(".", " "))
                    logger.debug("%r: missing issn and non-standard abbrev for %r", rcsbId, curAbbrev)

                if not curAbbrev:
                    logger.info("%r: missing issn and journal abbrev", rcsbId)
                #
            logger.debug("%s: revised: %r current: %r", rcsbId, revAbbrev, curAbbrev)
        except Exception as e:
            logger.exception("Failing on %r %r %r with %r", rcsbId, issn, curAbbrev, str(e))

        return revAbbrev

    def assignPrimaryCitation(self, dataContainer, catName, atName, **kwargs):
        """Normalize citation journal abbrev.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name
            atName (str): Attribute name

        Returns:
            bool: True for success or False otherwise
        """
        logger.debug("Starting catName %s atName %s kwargs %r", catName, atName, kwargs)
        try:
            if not dataContainer.exists(catName):
                return False
            #
            cObj = dataContainer.getObj(catName)
            if not cObj.hasAttribute(atName):
                cObj.appendAttribute(atName)
            #
            for ii in range(cObj.getRowCount()):
                citId = cObj.getValue("id", ii)
                if citId.upper() == "PRIMARY":
                    cObj.setValue("Y", atName, ii)
                else:
                    cObj.setValue("N", atName, ii)
            return True
        except Exception as e:
            logger.exception("Failing for %r with %s", dataContainer.getName(), str(e))
        return False

    def __getEmdbIdentifiers(self, dataContainer):
        """[summary]

        Args:
            dataContainer ([type]): [description]

        Returns:
            [type]: [description]

            #
            loop_
            _database_2.database_id
            _database_2.database_code
            PDB   6QUY
            WWPDB D_1292100913
            EMDB  EMD-4644
            #
            loop_
            _pdbx_database_related.db_name
            _pdbx_database_related.details
            _pdbx_database_related.db_id
            _pdbx_database_related.content_type
            EMDB 'HsCKK (human CAMSAP1) decorated 13pf taxol-GDP microtubule (asymmetric unit)' EMD-4643 'other EM volume'
            PDB  'HsCKK (human CAMSAP1) decorated 13pf taxol-GDP microtubule (asymmetric unit)' 6QUS     unspecified
            EMDB 'NgCKK (N.Gruberi CKK) decorated 13pf taxol-GDP microtubule'                   EMD-4644 'associated EM volume'
            #
        """
        emdbIdD = {}
        emdbIdAltD = {}
        if dataContainer.exists("database_2"):
            dbObj = dataContainer.getObj("database_2")
            for ii in range(dbObj.getRowCount()):
                dbId = dbObj.getValue("database_id", ii)
                dbCode = dbObj.getValue("database_code", ii)
                if dbId.upper() == "EMDB":
                    emdbIdD[dbCode] = "associated EM volume"

        if dataContainer.exists("pdbx_database_related"):
            drObj = dataContainer.getObj("pdbx_database_related")
            for ii in range(drObj.getRowCount()):
                dbCode = drObj.getValue("db_id", ii)
                dbName = drObj.getValue("db_name", ii)
                contentType = drObj.getValue("content_type", ii)
                if dbName.upper() == "EMDB" and contentType.upper() == "ASSOCIATED EM VOLUME" and dbCode not in emdbIdD:
                    emdbIdD[dbCode] = "associated EM volume"
                elif dbName.upper() == "EMDB" and contentType.upper() != "ASSOCIATED EM VOLUME" and dbCode not in emdbIdAltD:
                    emdbIdAltD[dbCode] = contentType
        return emdbIdD, emdbIdAltD

    def buildContainerEntryIds(self, dataContainer, catName, **kwargs):
        """Load the input category with rcsb_entry_container_identifiers content.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name

        Returns:
            bool: True for success or False otherwise

        For example:

        loop_
        _rcsb_entry_container_identifiers.entry_id
        _rcsb_entry_container_identifiers.entity_ids
        _rcsb_entry_container_identifiers.polymer_entity_ids_polymer
        _rcsb_entry_container_identifiers.non-polymer_entity_ids
        _rcsb_entry_container_identifiers.assembly_ids
        _rcsb_entry_container_identifiers.rcsb_id
        ...

        """
        logger.debug("Starting catName  %s kwargs %r", catName, kwargs)
        try:
            if not dataContainer.exists("entry"):
                return False
            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=self.__dApi.getAttributeNameList(catName)))
            #
            cObj = dataContainer.getObj(catName)

            tObj = dataContainer.getObj("entry")
            entryId = tObj.getValue("id", 0)

            cObj.setValue(entryId, "entry_id", 0)
            cObj.setValue(entryId, "rcsb_id", 0)

            #
            tObj = dataContainer.getObj("entity")
            entityIdL = tObj.getAttributeValueList("id")
            cObj.setValue(",".join(entityIdL), "entity_ids", 0)
            #
            #
            tIdL = tObj.selectValuesWhere("id", "polymer", "type")
            tV = ",".join(tIdL) if tIdL else "?"
            cObj.setValue(tV, "polymer_entity_ids", 0)

            tIdL = tObj.selectValuesWhere("id", "non-polymer", "type")
            tV = ",".join(tIdL) if tIdL else "?"
            cObj.setValue(tV, "non-polymer_entity_ids", 0)
            #
            tIdL = tObj.selectValuesWhere("id", "branched", "type")
            tV = ",".join(tIdL) if tIdL else "?"
            cObj.setValue(tV, "branched_entity_ids", 0)
            #
            # tIdL = tObj.selectValuesWhere("id", "water", "type")
            # tV = ",".join(tIdL) if tIdL else "?"
            # cObj.setValue(tV, "water_entity_ids", 0)
            #
            tObj = dataContainer.getObj("pdbx_struct_assembly")
            assemblyIdL = tObj.getAttributeValueList("id") if tObj else []
            tV = ",".join(assemblyIdL) if assemblyIdL else "?"
            cObj.setValue(tV, "assembly_ids", 0)
            #
            #
            emdbIdD, emdbIdAltD = self.__getEmdbIdentifiers(dataContainer)
            tV = ",".join([tId for tId in emdbIdD]) if emdbIdD else "?"
            cObj.setValue(tV, "emdb_ids", 0)
            tV = ",".join([tId for tId in emdbIdAltD]) if emdbIdAltD else "?"
            cObj.setValue(tV, "related_emdb_ids", 0)
            #
            modelIdList = self.__commonU.getModelIdList(dataContainer)
            tV = ",".join([str(tId) for tId in modelIdList]) if modelIdList else "?"
            cObj.setValue(tV, "model_ids", 0)
            #
            return True
        except Exception as e:
            logger.exception("For %s failing with %s", catName, str(e))
        return False

    def consolidateAccessionDetails(self, dataContainer, catName, **kwargs):
        """Consolidate accession details into the rcsb_accession_info category. Also include
        a flag for the availability of any supporting experimental data.

        Args:
            dataContainer (object): mmif.api.DataContainer object instance
            catName (str): Category name

        Returns:
            bool: True for success or False otherwise

        For example:
            For example -
             _rcsb_accession_info.entry_id                1ABC
             _rcsb_accession_info.status_code             REL
             _rcsb_accession_info.deposit_date            2018-01-11
             _rcsb_accession_info.initial_release_date    2018-03-23
             _rcsb_accession_info.major_revision          1
             _rcsb_accession_info.minor_revision          2
             _rcsb_accession_info.revision_date           2018-10-25


            Taking data values from:

            _pdbx_database_status.entry_id                        3OQP
            _pdbx_database_status.deposit_site                    RCSB
            _pdbx_database_status.process_site                    RCSB
            _pdbx_database_status.recvd_initial_deposition_date   2010-09-03
            _pdbx_database_status.status_code                     REL
            _pdbx_database_status.status_code_sf                  REL
            _pdbx_database_status.status_code_mr                  ?
            _pdbx_database_status.status_code_cs                  ?
            _pdbx_database_status.pdb_format_compatible           Y
            _pdbx_database_status.methods_development_category    ?
            _pdbx_database_status.SG_entry                        Y
            #
            loop_
            _pdbx_audit_revision_history.ordinal
            _pdbx_audit_revision_history.data_content_type
            _pdbx_audit_revision_history.major_revision
            _pdbx_audit_revision_history.minor_revision
            _pdbx_audit_revision_history.revision_date
            1 'Structure model' 1 0 2010-10-13
            2 'Structure model' 1 1 2011-07-13
            3 'Structure model' 1 2 2011-07-20
            4 'Structure model' 1 3 2014-11-12
            5 'Structure model' 1 4 2017-10-25
            #

            #  - For EM and SAS -
            _pdbx_database_related.db_name        EMDB
            _pdbx_database_related.details
            'pseudo-atomic model of the RNA polymerase lambda-based antitermination complex solved by cryo-EM'
            _pdbx_database_related.db_id          EMD-3561
            _pdbx_database_related.content_type   'associated EM volume'
        """
        ##
        try:
            logger.debug("Starting with  %r %r %r", dataContainer.getName(), catName, kwargs)
            #
            entryId = None
            # Add missing pdbx_database_status for MA or AF models (if absent in mmCIF file)
            cName = "pdbx_database_status"
            if dataContainer.exists("entry") and dataContainer.exists("ma_data"):
                if not dataContainer.exists(cName):
                    dObj = dataContainer.getObj("entry")
                    entryId = dObj.getValue("id", 0)
                    dataContainer.append(DataCategory(cName, attributeNameList=self.__dApi.getAttributeNameList(cName)))
                    eObj = dataContainer.getObj(cName)
                    eObj.setValue(entryId, "entry_id", 0)
                    eObj.setValue("REL", "status_code", 0)
                    eObj.setValue("?", "recvd_initial_deposition_date", 0)
                else:
                    # if it does exist but is missing one or more attributes, fill them in (for CSMs only!)
                    pdsAttrL = dataContainer.getObj(cName).getAttributeList()
                    eObj = dataContainer.getObj(cName)
                    if "entry_id" not in pdsAttrL:
                        dObj = dataContainer.getObj("entry")
                        entryId = dObj.getValue("id", 0)
                        eObj.appendAttribute("entry_id")
                        eObj.setValue(entryId, "entry_id", 0)
                    if "status_code" not in pdsAttrL:
                        eObj.appendAttribute("status_code")
                        eObj.setValue("REL", "status_code", 0)
                    if "recvd_initial_deposition_date" not in pdsAttrL:
                        eObj.appendAttribute("recvd_initial_deposition_date")
                        eObj.setValue("?", "recvd_initial_deposition_date", 0)
                # Make sure status_code is set to "REL" (and not "HPUB" or something else)
                if eObj.getValue("status_code", 0) != "REL":
                    eObj.setValue("REL", "status_code", 0)

            # if there is incomplete accessioninformation then exit
            if not dataContainer.exists("pdbx_database_status"):
                return False
            # Create the new target category
            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=self.__dApi.getAttributeNameList(catName)))

            cObj = dataContainer.getObj(catName)
            #
            tObj = dataContainer.getObj("pdbx_database_status")
            entryId = tObj.getValue("entry_id", 0)
            statusCode = tObj.getValue("status_code", 0)
            depositDate = tObj.getValueOrDefault("recvd_initial_deposition_date", 0, "?")
            #
            cObj.setValue(entryId, "entry_id", 0)
            cObj.setValue(statusCode, "status_code", 0)
            cObj.setValue(depositDate, "deposit_date", 0)
            # cObj.setValue(depositDate[:4], "deposit_year", 0)
            #
            # -- Experimental data availability --
            #
            expDataRelFlag = "N"
            statusSf = tObj.getValueOrDefault("status_code_sf", 0, defaultValue=None)
            statusMr = tObj.getValueOrDefault("status_code_mr", 0, defaultValue=None)
            statusCs = tObj.getValueOrDefault("status_code_cs", 0, defaultValue=None)
            statusNmrData = tObj.getValueOrDefault("status_code_nmr_data", 0, defaultValue=None)
            #
            if statusSf == "REL" or statusMr == "REL" or statusCs == "REL" or statusNmrData == "REL":
                expDataRelFlag = "Y"
            else:
                if dataContainer.exists("pdbx_database_related"):
                    rObj = dataContainer.getObj("pdbx_database_related")
                    ctL = rObj.getAttributeValueList("content_type")
                    if "associated EM volume" in ctL or "associated SAS data" in ctL:
                        expDataRelFlag = "Y"
            #
            cObj.setValue(expDataRelFlag, "has_released_experimental_data", 0)
            #
            if dataContainer.exists("pdbx_audit_revision_history"):
                tObj = dataContainer.getObj("pdbx_audit_revision_history")
                nRows = tObj.getRowCount()
                # Assuming the default sorting order from the release module -
                releaseDate = tObj.getValue("revision_date", 0)
                minorRevision = tObj.getValue("minor_revision", nRows - 1)
                majorRevision = tObj.getValue("major_revision", nRows - 1)
                revisionDate = tObj.getValue("revision_date", nRows - 1)
                cObj.setValue(releaseDate, "initial_release_date", 0)
                # cObj.setValue(releaseDate[:4], "initial_release_year", 0)
                cObj.setValue(minorRevision, "minor_revision", 0)
                cObj.setValue(majorRevision, "major_revision", 0)
                cObj.setValue(revisionDate, "revision_date", 0)
            #
            return True
        except Exception as e:
            logger.exception("In %s for %s failing with %s", dataContainer.getName(), catName, str(e))
        return False

    def filterRedundantRecords(self, dataContainer, catName, **kwargs):
        """Filter redundant records from input category subject to excluded/included attributes."""
        try:
            logger.debug("Starting with %r %r %r", dataContainer.getName(), catName, kwargs)
            # Exit if source categories are missing
            if not dataContainer.exists(catName):
                return False
            #
            cObj = dataContainer.getObj(catName)
            if cObj.getRowCount() < 2:
                return False
            #
            if catName == "pdbx_related_exp_data_set":
                logger.debug("Filtering %r %r", dataContainer.getName(), catName)
                try:
                    cObj.removeAttribute("ordinal")
                    cObj.removeDuplicateRows()
                    cObj.appendAttribute("ordinal")
                    for ii in range(cObj.getRowCount()):
                        cObj.setValue(ii + 1, "ordinal", ii)
                except Exception as e:
                    logger.exception("%s failing with %s", dataContainer.getName(), str(e))
                    return False

                return True
        except Exception as e:
            logger.exception("For %s %r failing with %s", dataContainer.getName(), catName, str(e))
        #
        return False

    def addEntryInfo(self, dataContainer, catName, **kwargs):
        """
        Add  _rcsb_entry_info, for example:
             _rcsb_entry_info.entry_id                              1ABC
             _rcsb_entry_info.polymer_composition                   'heteromeric protein'
             _rcsb_entry_info.structure_determination_methodology   'experimental'
             _rcsb_entry_info.experimental_method                   'multiple methods'
             _rcsb_entry_info.experimental_method_count             2
             _rcsb_entry_info.polymer_entity_count                  2
             _rcsb_entry_info.entity_count                          2
             _rcsb_entry_info.nonpolymer_entity_count               2
             _rcsb_entry_info.branched_entity_count                 0
             _rcsb_entry_info.software_programs_combined            'Phenix;RefMac'
             ....

        Also add the related field:

        _entity_poly.rcsb_entity_polymer_type

            'Protein'   'polypeptide(D) or polypeptide(L)'
            'DNA'       'polydeoxyribonucleotide'
            'RNA'       'polyribonucleotide'
            'NA-hybrid' 'polydeoxyribonucleotide/polyribonucleotide hybrid'
            'Other'      'polysaccharide(D), polysaccharide(L), cyclic-pseudo-peptide, peptide nucleic acid, or other'
            #
          _rcsb_entry_info.deposited_polymer_monomer_count
          'polymer_entity_count_protein',
          'polymer_entity_count_nucleic_acid',
          'polymer_entity_count_nucleic_acid_hybrid',
          'polymer_entity_count_DNA',
          'polymer_entity_count_RNA',

        """
        try:
            logger.debug("Starting with %r %r %r", dataContainer.getName(), catName, kwargs)
            # Exit if source categories are missing
            if not (dataContainer.exists("entity") and dataContainer.exists("entry")):
                return False
            if not (dataContainer.exists("exptl") or dataContainer.exists("ma_model_list")):
                return False
            #
            # Create the new target category rcsb_entry_info
            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=self.__dApi.getAttributeNameList(catName)))
            # --------------------------------------------------------------------------------------------------------
            # catName = rcsb_entry_info
            cObj = dataContainer.getObj(catName)
            #
            # --------------------------------------------------------------------------------------------------------
            #  Filter experimental methods
            #
            methodCount = 0
            expMethod = None
            methodType = None
            entryId = None
            methodPriority = None
            #
            if dataContainer.exists("struct"):
                xObj = dataContainer.getObj("struct")
                if xObj.hasAttribute("pdbx_structure_determination_methodology"):
                    methodType = xObj.getValue("pdbx_structure_determination_methodology", 0)
            if not methodType or methodType == "?" or methodType == ".":
                if dataContainer.exists("exptl"):
                    methodType = "experimental"
                if dataContainer.exists("ma_data"):
                    methodType = "computational"
            if methodType == "experimental":
                methodPriority = 10
            elif methodType == "computational":
                methodPriority = 100
            #
            if dataContainer.exists("exptl"):
                xObj = dataContainer.getObj("exptl")
                entryId = xObj.getValue("entry_id", 0)
                methodCount, expMethod = self.__commonU.filterExperimentalMethod(dataContainer)
                cObj.setValue(expMethod, "experimental_method", 0)
            elif dataContainer.exists("ma_model_list"):
                tObj = dataContainer.getObj("entry")
                entryId = tObj.getValue("id", 0)
                methodCount, expMethod = self.__commonU.filterExperimentalMethod(dataContainer)
            #
            if methodType not in ["experimental", "computational"]:
                logger.error("Unexpected methodType %r found for entry %r", methodType, entryId)
            #
            cObj.setValue(entryId, "entry_id", 0)
            cObj.setValue(methodCount, "experimental_method_count", 0)
            cObj.setValue(methodType, "structure_determination_methodology", 0)
            cObj.setValue(methodPriority, "structure_determination_methodology_priority", 0)
            #
            # --------------------------------------------------------------------------------------------------------
            #  Experimental resolution -
            #
            resL = self.__filterExperimentalResolution(dataContainer)
            if resL:
                cObj.setValue(",".join(resL), "resolution_combined", 0)
            #
            # ---------------------------------------------------------------------------------------------------------
            # Consolidate software details -
            #
            swNameL = []
            if dataContainer.exists("software"):
                swObj = dataContainer.getObj("software")
                swNameL.extend(swObj.getAttributeUniqueValueList("name"))
            if dataContainer.exists("pdbx_nmr_software"):
                swObj = dataContainer.getObj("pdbx_nmr_software")
                swNameL.extend(swObj.getAttributeUniqueValueList("name"))
            if dataContainer.exists("em_software"):
                swObj = dataContainer.getObj("em_software")
                swNameL.extend(swObj.getAttributeUniqueValueList("name"))
            if swNameL:
                swNameD = {swName.upper().strip(): True for swName in swNameL if swName not in [".", "?"]}
                swNameL = sorted(swNameD.keys())
                cObj.setValue(";".join(swNameL), "software_programs_combined", 0)
            # ---------------------------------------------------------------------------------------------------------
            #  ENTITY FEATURES
            #
            # Nucleic acid secondary structure features
            naFeatureL = []
            if dataContainer.exists("ndb_struct_conf_na"):
                naObj = dataContainer.getObj("ndb_struct_conf_na")
                naFeatureL.extend(naObj.getAttributeUniqueValueList("feature"))
            if naFeatureL:
                nL = [naFeature.strip() for naFeature in naFeatureL if naFeature not in [".", "?"]]
                cObj.setValue(",".join(nL), "ndb_struct_conf_na_feature_combined", 0)
            #  entity and polymer entity counts -
            ##
            eObj = dataContainer.getObj("entity")
            eTypeL = eObj.getAttributeValueList("type")
            #
            numPolymers = 0
            numNonPolymers = 0
            numBranched = 0
            numSolvent = 0
            for eType in eTypeL:
                if eType == "polymer":
                    numPolymers += 1
                elif eType == "non-polymer":
                    numNonPolymers += 1
                elif eType == "branched":
                    numBranched += 1
                elif eType == "water":
                    numSolvent += 1
                else:
                    logger.error("Unexpected entity type for %s %s", dataContainer.getName(), eType)
            totalEntities = numPolymers + numNonPolymers + numBranched + numSolvent
            #
            # Simplified entity polymer type: 'Protein', 'DNA', 'RNA', 'NA-hybrid', or 'Other'
            pTypeL = []
            if dataContainer.exists("entity_poly"):
                epObj = dataContainer.getObj("entity_poly")
                pTypeL = epObj.getAttributeValueList("type")
                #
                atName = "rcsb_entity_polymer_type"
                if not epObj.hasAttribute(atName):
                    epObj.appendAttribute(atName)
                for ii in range(epObj.getRowCount()):
                    epObj.setValue(self.__commonU.filterEntityPolyType(pTypeL[ii]), atName, ii)
            #
            # Add any branched entity types to the type list -
            if dataContainer.exists("pdbx_entity_branch"):
                ebObj = dataContainer.getObj("pdbx_entity_branch")
                pTypeL.extend(ebObj.getAttributeValueList("type"))
            #
            polymerCompClass, ptClass, naClass, eptD = self.__commonU.getPolymerComposition(pTypeL)
            if eptD and len(eptD) > 2:
                logger.debug("%s entity type count=%d class=%s typeD %r", dataContainer.getName(), len(eptD), polymerCompClass, eptD)
            #
            cObj.setValue(polymerCompClass, "polymer_composition", 0)
            cObj.setValue(ptClass, "selected_polymer_entity_types", 0)
            cObj.setValue(naClass, "na_polymer_entity_types", 0)
            cObj.setValue(numPolymers, "polymer_entity_count", 0)
            cObj.setValue(numNonPolymers, "nonpolymer_entity_count", 0)
            cObj.setValue(numBranched, "branched_entity_count", 0)
            cObj.setValue(numSolvent, "solvent_entity_count", 0)
            cObj.setValue(totalEntities, "entity_count", 0)
            #
            num = eptD["protein"] if "protein" in eptD else 0
            cObj.setValue(num, "polymer_entity_count_protein", 0)
            #
            num = eptD["NA-hybrid"] if "NA-hybrid" in eptD else 0
            cObj.setValue(num, "polymer_entity_count_nucleic_acid_hybrid", 0)
            #
            numDNA = eptD["DNA"] if "DNA" in eptD else 0
            cObj.setValue(numDNA, "polymer_entity_count_DNA", 0)
            #
            numRNA = eptD["RNA"] if "RNA" in eptD else 0
            cObj.setValue(numRNA, "polymer_entity_count_RNA", 0)
            cObj.setValue(numDNA + numRNA, "polymer_entity_count_nucleic_acid", 0)
            #
            # ---------------------------------------------------------------------------------------------------------
            # INSTANCE FEATURES
            #
            repModelId = self.__commonU.getRepresentativeModelId(dataContainer)
            #
            instanceTypeCountD = self.__commonU.getInstanceTypeCounts(dataContainer)
            cObj.setValue(instanceTypeCountD["polymer"], "deposited_polymer_entity_instance_count", 0)
            cObj.setValue(instanceTypeCountD["non-polymer"], "deposited_nonpolymer_entity_instance_count", 0)

            #
            # Various atom counts -
            #
            numHeavyAtomsModel, numHydrogenAtomsModel, numAtomsTotal, numModelsTotal, numDeuWatMolModel = self.__commonU.getDepositedAtomCounts(dataContainer, modelId=repModelId)
            #
            logger.debug("numAtomsTotal %d numHeavyAtomsModel %d numModelsTotal %d", numAtomsTotal, numHeavyAtomsModel, numModelsTotal)
            logger.debug("entity type atom counts %r", self.__commonU.getEntityTypeHeavyAtomCounts(dataContainer, modelId=repModelId))
            logger.debug("instance atom counts %r", self.__commonU.getEntityTypeHeavyAtomCounts(dataContainer, modelId=repModelId))
            #

            if numHeavyAtomsModel > 0:
                cObj.setValue(numHeavyAtomsModel, "deposited_atom_count", 0)
                cObj.setValue(numModelsTotal, "deposited_model_count", 0)
                cObj.setValue(numHydrogenAtomsModel, "deposited_hydrogen_atom_count", 0)
                cObj.setValue(numDeuWatMolModel, "deposited_deuterated_water_count", 0)
                tCD = self.__commonU.getEntityTypeHeavyAtomCounts(dataContainer, modelId=repModelId)
                wCount = tCD["water"] if tCD and "water" in tCD else 0
                cObj.setValue(wCount, "deposited_solvent_atom_count", 0)
            #
            # ---------------------------------------------------------------------------------------------------------
            #  Deposited monomer/residue instance counts
            #
            #  Get modeled and unmodeled residue counts
            #
            modeledCount, unModeledCount = self.__commonU.getDepositedMonomerCounts(dataContainer, modelId=repModelId)
            cObj.setValue(modeledCount, "deposited_modeled_polymer_monomer_count", 0)
            cObj.setValue(unModeledCount, "deposited_unmodeled_polymer_monomer_count", 0)
            cObj.setValue(modeledCount + unModeledCount, "deposited_polymer_monomer_count", 0)
            #
            # ---------------------------------------------------------------------------------------------------------
            #  Counts of intermolecular bonds/linkages
            #
            #
            bCountsD = self.__commonU.getInstanceConnectionCounts(dataContainer)
            cObj.setValue(bCountsD["disulf"], "disulfide_bond_count", 0)
            cObj.setValue(bCountsD["metalc"], "inter_mol_metalic_bond_count", 0)
            cObj.setValue(bCountsD["covale"], "inter_mol_covalent_bond_count", 0)
            #
            cisPeptideD = self.__ssU.getCisPeptides(dataContainer)
            cObj.setValue(len(cisPeptideD), "cis_peptide_count", 0)
            #
            # This is reset in anothor method - filterSourceOrganismDetails()
            cObj.setValue(None, "polymer_entity_taxonomy_count", 0)
            #
            fw = self.__commonU.getFormulaWeightNonSolvent(dataContainer)
            cObj.setValue(str(round(fw, 2)), "molecular_weight", 0)
            #
            # nonpolymer_bound_components
            #
            bcL = self.__commonU.getBoundNonpolymersComponentIds(dataContainer)
            if bcL:
                cObj.setValue(";".join(bcL), "nonpolymer_bound_components", 0)
            #
            # polymer_molecular_weight_minimum
            # polymer_molecular_weight_maximum
            # nonpolymer_molecular_weight_minimum
            # nonpolymer_molecular_weight_maximum
            # branched_molecular_weight_minimum
            # branched_molecular_weight_maximum
            #
            fwBoundD = self.__commonU.getEntityFormulaWeightBounds(dataContainer)
            if "polymer" in fwBoundD and fwBoundD["polymer"]["min"] and fwBoundD["polymer"]["max"]:
                cObj.setValue(str(round(fwBoundD["polymer"]["min"], 2)), "polymer_molecular_weight_minimum", 0)
                cObj.setValue(str(round(fwBoundD["polymer"]["max"], 2)), "polymer_molecular_weight_maximum", 0)
            if "non-polymer" in fwBoundD and fwBoundD["non-polymer"]["min"] and fwBoundD["non-polymer"]["max"]:
                cObj.setValue(str(round(fwBoundD["non-polymer"]["min"], 2)), "nonpolymer_molecular_weight_minimum", 0)
                cObj.setValue(str(round(fwBoundD["non-polymer"]["max"], 2)), "nonpolymer_molecular_weight_maximum", 0)
            if "branched" in fwBoundD and fwBoundD["branched"]["min"] and fwBoundD["branched"]["max"]:
                cObj.setValue(str(round(fwBoundD["branched"]["min"], 2)), "branched_molecular_weight_minimum", 0)
                cObj.setValue(str(round(fwBoundD["branched"]["max"], 2)), "branched_molecular_weight_maximum", 0)
            #
            # polymer_monomer_count_maximum
            # polymer_monomer_count_minimum
            #
            polymerLengthBounds = self.__commonU.getEntityPolymerLengthBounds(dataContainer)
            if polymerLengthBounds:
                cObj.setValue(str(polymerLengthBounds[0]), "polymer_monomer_count_minimum", 0)
                cObj.setValue(str(polymerLengthBounds[1]), "polymer_monomer_count_maximum", 0)
            #
            # ---------------------------------------------------------------------------------------------------------
            # Consolidate diffraction wavelength details -
            wL = []
            try:
                if dataContainer.exists("diffrn_radiation_wavelength"):
                    swObj = dataContainer.getObj("diffrn_radiation_wavelength")
                    wL.extend(swObj.getAttributeUniqueValueList("wavelength"))
                if dataContainer.exists("diffrn_radiation"):
                    swObj = dataContainer.getObj("diffrn_radiation")
                    if swObj.hasAttribute("pdbx_wavelength"):
                        wL.extend(swObj.getAttributeUniqueValueList("pdbx_wavelength"))
                    if swObj.hasAttribute("pdbx_wavelength_list"):
                        tL = []
                        for tS in swObj.getAttributeUniqueValueList("pdbx_wavelength_list"):
                            tL.extend(tS.split(","))
                        if tL:
                            wL.extend(tL)
                if dataContainer.exists("diffrn_source"):
                    swObj = dataContainer.getObj("diffrn_source")
                    if swObj.hasAttribute("pdbx_wavelength"):
                        wL.extend(swObj.getAttributeUniqueValueList("pdbx_wavelength"))
                    if swObj.hasAttribute("pdbx_wavelength_list"):
                        tL = []
                        for tS in swObj.getAttributeUniqueValueList("pdbx_wavelength_list"):
                            tL.extend(tS.split(","))
                        if tL:
                            wL.extend(tL)
                fL = []
                for wS in wL:
                    try:
                        fL.append(float(wS))
                    except Exception:
                        pass
                if fL:
                    cObj.setValue("%.4f" % min(fL), "diffrn_radiation_wavelength_minimum", 0)
                    cObj.setValue("%.4f" % max(fL), "diffrn_radiation_wavelength_maximum", 0)

            except Exception as e:
                logger.exception("%s failing wavelength processing with %s", entryId, str(e))
            #
            # JDW
            self.__updateReflnsResolution(dataContainer)
            return True
        except Exception as e:
            logger.exception("For %s %r failing with %s", dataContainer.getName(), catName, str(e))
        #
        return False

    def filterBlockByMethod(self, dataContainer, blockName, **kwargs):
        """Filter empty placeholder data categories by experimental method."""
        logger.debug("Starting with %r blockName %r kwargs %r", dataContainer.getName(), blockName, kwargs)
        try:
            if not dataContainer.exists("exptl"):
                return False
            #
            xObj = dataContainer.getObj("exptl")
            methodL = xObj.getAttributeValueList("method")
            objNameL = []
            # Test for a diffraction method in the case of multiple methods
            if len(methodL) > 1:
                isXtal = False
                for method in methodL:
                    if method in ["X-RAY DIFFRACTION", "FIBER DIFFRACTION", "POWDER DIFFRACTION", "ELECTRON CRYSTALLOGRAPHY", "NEUTRON DIFFRACTION", "ELECTRON DIFFRACTION"]:
                        isXtal = True
                        break
                if not isXtal:
                    objNameL = ["cell", "symmetry", "refine", "refine_hist", "software", "diffrn", "diffrn_radiation"]
            else:
                #
                mS = methodL[0].upper()
                if mS in ["X-RAY DIFFRACTION", "FIBER DIFFRACTION", "POWDER DIFFRACTION", "ELECTRON CRYSTALLOGRAPHY", "NEUTRON DIFFRACTION", "ELECTRON DIFFRACTION"]:
                    objNameL = []
                elif mS in ["SOLUTION NMR", "SOLID-STATE NMR"]:
                    objNameL = ["cell", "symmetry", "refine", "refine_hist", "software", "diffrn", "diffrn_radiation"]
                elif mS in ["ELECTRON MICROSCOPY", "CRYO-ELECTRON MICROSCOPY"]:
                    objNameL = ["cell", "symmetry", "refine", "refine_hist", "software", "diffrn", "diffrn_radiation"]
                elif mS in ["SOLUTION SCATTERING", "EPR", "THEORETICAL MODEL", "INFRARED SPECTROSCOPY", "FLUORESCENCE TRANSFER"]:
                    objNameL = ["cell", "symmetry", "refine", "refine_hist", "software", "diffrn", "diffrn_radiation"]
                else:
                    logger.error("%s Unexpected method %r", dataContainer.getName(), mS)
            #
            for objName in objNameL:
                dataContainer.remove(objName)
            return True
        except Exception as e:
            logger.exception("For %s failing with %s", dataContainer.getName(), str(e))
        return False

    def filterEnumerations(self, dataContainer, catName, atName, **kwargs):
        """Standardize the item value to conform to enumeration specifications."""
        logger.debug("Starting with %r %r %r %r", dataContainer.getName(), atName, catName, kwargs)
        subD = {("pdbx_reference_molecule", "class"): [("Anti-tumor", "Antitumor")]}
        try:
            if not dataContainer.exists(catName):
                return False
            #
            cObj = dataContainer.getObj(catName)
            if not cObj.hasAttribute(atName):
                return False
            #
            subL = subD[(catName, atName)] if (catName, atName) in subD else []
            #
            for ii in range(cObj.getRowCount()):
                tV = cObj.getValue(atName, ii)
                if tV and tV not in [".", "?"]:
                    for sub in subL:
                        if sub[0] in tV:
                            tV = tV.replace(sub[0], sub[1])
                            cObj.setValue(tV, atName, ii)
            return True
        except Exception as e:
            logger.exception("%s %s %s failing with %s", dataContainer.getName(), catName, atName, str(e))
        return False

    def __filterExperimentalResolution(self, dataContainer):
        """Collect resolution estimates from method specific sources."""

        rL = []
        fL = []
        eL = []

        # _refine.pdbx_refine_id and refine.entry_id are joint composite keys for the refine category
        # So, it is not possible to have multiple resolutions values for the same method in this category
        # But _refine.pdbx_refine_id does not have a controlled vocabulary, and this can lead to accidental errors
        if dataContainer.exists("refine"):
            tObj = dataContainer.getObj("refine")
            if tObj.hasAttribute("ls_d_res_high") and tObj.hasAttribute("pdbx_refine_id"):
                for ii in range(tObj.getRowCount()):
                    rv = tObj.getValue("ls_d_res_high", ii)
                    rid = tObj.getValue("pdbx_refine_id", ii)
                    rM = rid.upper()
                    if self.__commonU.isFloat(rv):
                        if rM in ["X-RAY DIFFRACTION", "FIBER DIFFRACTION", "POWDER DIFFRACTION", "ELECTRON CRYSTALLOGRAPHY", "NEUTRON DIFFRACTION", "ELECTRON DIFFRACTION"]:
                            rL.append(rv)

        if dataContainer.exists("em_3d_reconstruction"):
            tObj = dataContainer.getObj("em_3d_reconstruction")
            if tObj.hasAttribute("resolution") and tObj.hasAttribute("resolution_method"):
                for ii in range(tObj.getRowCount()):
                    rv = tObj.getValue("resolution", ii)
                    rM = tObj.getValue("resolution_method", ii)
                    if self.__commonU.isFloat(rv):
                        if rM.upper() in ["FSC 0.143 CUT-OFF"]:
                            fL.append(rv)
                        else:
                            eL.append(rv)
                if fL:
                    rL.append(min(fL))
                elif eL:
                    rL.append(min(eL))
                else:
                    pass

        return rL

    def addCategoryPrimaryCitation(self, dataContainer, blockName, **kwargs):
        """
        Add  rcsb_primary_citation category as a copy or the citation category
        with rcsb extensions.
        """
        catName = None
        try:
            logger.debug("Starting with %r %r %r", dataContainer.getName(), blockName, kwargs)
            # Exit if source categories are missing
            if not dataContainer.exists("citation"):
                return False
            cObj = dataContainer.getObj("citation")
            catName = "rcsb_primary_citation"
            #
            if not dataContainer.exists(catName):
                dataContainer.append(DataCategory(catName, attributeNameList=self.__dApi.getAttributeNameList(catName)))
            # --------------------------------------------------------------------------------------------------------
            rObj = dataContainer.getObj(catName)
            atNameList = self.__dApi.getAttributeNameList(catName)
            logger.debug("Category %s dict attributes %r", catName, atNameList)
            #
            for ii in range(cObj.getRowCount()):
                pv = cObj.getValue("id", ii)
                if pv.upper() == "PRIMARY":
                    for atName in atNameList:
                        if cObj.hasAttribute(atName):
                            rObj.setValue(cObj.getValue(atName, ii), atName, 0)

            return True
        except Exception as e:
            logger.exception("%s %s failing with %s", dataContainer.getName(), catName, str(e))
        return False

    def __updateReflnsResolution(self, dataContainer):
        """Find a plausable data collection diffraction high resolution limit from one of the following sources.
        #
        _rcsb_entry_info.diffrn_resolution_high_value
        _rcsb_entry_info.diffrn_resolution_high_provenance_source

        Update category 'reflns' with any missing resolution extrema data using limits in category reflns_shell.

            _reflns.entry_id                     2DCG
            _reflns.d_resolution_high            0.900
            _reflns.pdbx_diffrn_id               1
            _reflns.pdbx_ordinal                 1


            _refine.entry_id                                 2DCG
            _refine.ls_number_reflns_obs                     15000
            _refine.ls_number_reflns_all                     ?
            _refine.pdbx_ls_sigma_I                          2.000
            _refine.ls_d_res_low                             ?
            _refine.ls_d_res_high                            0.900
            _refine.pdbx_refine_id                           'X-RAY DIFFRACTION'
            _refine.pdbx_diffrn_id                           1

            _reflns_shell.d_res_high             1.18
            _reflns_shell.d_res_low              1.25
            _reflns_shell.pdbx_ordinal           1
            _reflns_shell.pdbx_diffrn_id         1
            #

        """
        try:
            logger.debug("Starting with %r", dataContainer.getName())
            #
            if not dataContainer.exists("exptl") or not dataContainer.exists("rcsb_entry_info"):
                return False
            # --------------------------------------------------------------------------------------------------------
            #  Only applicable to X-ray
            #
            _, expMethod = self.__commonU.filterExperimentalMethod(dataContainer)
            if expMethod not in ["X-ray", "Neutron", "Multiple methods"]:
                return False
            #
            resValue = resProvSource = None
            #
            # Here are the various cases -
            if dataContainer.exists("reflns"):
                rObj = dataContainer.getObj("reflns")
                if rObj.hasAttribute("d_resolution_high"):
                    rvL = rObj.getAttributeValueList("d_resolution_high")
                    fvL = [float(rv) for rv in rvL if self.__commonU.isFloat(rv)]
                    if fvL:
                        resValue = round(min(fvL), 2)
                        resProvSource = "Depositor assigned"

            if not resValue and dataContainer.exists("reflns_shell"):
                rObj = dataContainer.getObj("reflns_shell")
                if rObj.hasAttribute("d_res_high"):
                    rvL = rObj.getAttributeValueList("d_res_high")
                    fvL = [float(rv) for rv in rvL if self.__commonU.isFloat(rv)]
                    if fvL:
                        resValue = round(min(fvL), 2)
                        resProvSource = "From the high resolution shell"

            if not resValue and dataContainer.exists("refine"):

                rObj = dataContainer.getObj("refine")
                if rObj.hasAttribute("ls_d_res_high"):
                    fvL = []
                    for ii in range(rObj.getRowCount()):
                        rId = rObj.getValue("pdbx_refine_id", ii)
                        if rId in ["X-RAY DIFFRACTION", "NEUTRON DIFFRACTION", "FIBER DIFFRACTION"]:
                            rv = rObj.getValue("ls_d_res_high", ii)
                            if self.__commonU.isFloat(rv):
                                fvL.append(float(rv))
                    if fvL:
                        resValue = round(min(fvL), 2)
                        resProvSource = "From refinement resolution cutoff"
            #
            if not resValue:
                logger.debug("No source of data collection resolution available for %r", dataContainer.getName())
            else:
                logger.debug("Data collection diffraction limit %r PS %r", resValue, resProvSource)

            if resValue:
                eObj = dataContainer.getObj("rcsb_entry_info")
                for atName in ["diffrn_resolution_high_value", "diffrn_resolution_high_provenance_source"]:
                    if not eObj.hasAttribute(atName):
                        eObj.appendAttribute(atName)
                eObj.setValue(resValue, "diffrn_resolution_high_value", 0)
                eObj.setValue(resProvSource, "diffrn_resolution_high_provenance_source", 0)
                # --------------------------------------------------------------------------------------------------------
                return True
        except Exception as e:
            logger.exception("%s failing with %s", dataContainer.getName(), str(e))
        return False

    def filterRevisionHistory(self, dataContainer, catName, **kwargs):
        """Remove rows that don't belong to "data_content_type" == "Structure model"
           in revision history categories

        Example:
        loop_
        _pdbx_audit_revision_history.ordinal
        _pdbx_audit_revision_history.data_content_type
        _pdbx_audit_revision_history.major_revision
        _pdbx_audit_revision_history.minor_revision
        _pdbx_audit_revision_history.revision_date
        _pdbx_audit_revision_history.part_number
        1 'Structure model' 1 0 2025-01-22 ?
        2 'EM metadata' 1 0 2025-01-22 ?
        3 'Structure model' 1 1 2025-01-29 ?
        4 'EM metadata' 1 1 2025-01-29 ?
        #
        loop_
        _pdbx_audit_revision_details.ordinal
        _pdbx_audit_revision_details.revision_ordinal
        _pdbx_audit_revision_details.data_content_type
        _pdbx_audit_revision_details.provider
        _pdbx_audit_revision_details.type
        _pdbx_audit_revision_details.description
        _pdbx_audit_revision_details.details
        1 1 'Structure model' repository 'Initial release' ? ?
        2 2 'EM metadata' repository 'Initial release' ? ?
        3 4 'EM metadata' repository 'Data updated' ? ?
        #
        loop_
        _pdbx_audit_revision_group.ordinal
        _pdbx_audit_revision_group.revision_ordinal
        _pdbx_audit_revision_group.data_content_type
        _pdbx_audit_revision_group.group
        1 3 'Structure model' 'Data collection'
        2 3 'Structure model' Other
        3 3 'Structure model' 'Structure summary'
        4 4 'EM metadata' 'Experimental summary'
        5 4 'EM metadata' 'Structure summary'
        #
        loop_
        _pdbx_audit_revision_category.ordinal
        _pdbx_audit_revision_category.revision_ordinal
        _pdbx_audit_revision_category.data_content_type
        _pdbx_audit_revision_category.category
        1 3 'Structure model' em_admin
        2 3 'Structure model' pdbx_database_status
        3 3 'Structure model' pdbx_prerelease_seq
        4 3 'Structure model' struct_keywords
        5 4 'EM metadata' em_admin
        6 4 'EM metadata' struct_keywords
        """
        logger.debug("Starting with %s %r %r", dataContainer.getName(), catName, kwargs)
        try:
            if not dataContainer.exists("pdbx_audit_revision_history"):
                return False

            cndL = [("data_content_type", "not in", "Structure model")]
            cNameL = ["pdbx_audit_revision_history", "pdbx_audit_revision_details", "pdbx_audit_revision_group", "pdbx_audit_revision_category", "pdbx_audit_revision_item"]

            for catName in cNameL:
                if dataContainer.exists(catName):
                    cObj = dataContainer.getObj(catName)
                    rL = cObj.selectIndicesWhereOpConditions(cndL)
                    if rL:
                        logger.debug("For %s removing %s rows that don't correspond to structure model in %s", dataContainer.getName(), rL, catName)
                        cObj.removeRows(list(set(rL)))

            return True
        except Exception as e:
            logger.exception("For %s removing rows in revision history categories failing with %s", dataContainer.getName(), str(e))
        return False
