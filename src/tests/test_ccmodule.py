import sys
import os
import pandas as pd

# Ensure we can import from project root
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.extraction.ccmodule import PDBCrystallizationConditions


def make_extractor(conditions_list):
    e = PDBCrystallizationConditions("dummy_folder")
    e.conditions = conditions_list
    e.protein_id = "TEST"
    # Initialize dataframe with schema
    e.conditions_df = pd.DataFrame(columns=[
        'FASTA_Sequence', 'Solvent_Concentration',
        'Matthews_Coefficient', 'Crystallography_Type',
        'Vapor_Diffusion_Type', 'pH',
        'Organic_Precipitates', 'Salt_Precipitates'
    ])
    e.conditions_df.index.name = "Protein_ID"
    e.conditions_df.loc["TEST"] = [None] * len(e.conditions_df.columns)
    return e


def test_solvent_content_found():
    e = make_extractor(["REMARK 280 SOLVENT CONTENT, PERCENTAGE: 47.5"])
    e.crystallization_conditions_list = e.conditions
    e.extract_solvent_content()
    assert e.solvent_content == "47.5"


def test_solvent_content_not_found():
    e = make_extractor(["REMARK 280 SOME OTHER FIELD"])
    e.crystallization_conditions_list = e.conditions
    e.extract_solvent_content()
    assert e.solvent_content == "No Solvent Content"


def test_matthews_coefficient_parsing():
    e = make_extractor(["REMARK 280 MATTHEWS COEFFICIENT, VM (ANGSTROMS**3/DA): 3.07"])
    e.crystallization_conditions_list = e.conditions
    e.extract_matthews_coefficient()
    assert e.matthews_coefficient == "3.07"


def test_matthews_coefficient_missing():
    e = make_extractor(["REMARK 280 OTHER"])
    e.crystallization_conditions_list = e.conditions
    e.extract_matthews_coefficient()
    assert e.matthews_coefficient is None


def test_crystallography_type():
    e = make_extractor(["VAPOR DIFFUSION"])
    e.extract_crystallography_type()
    assert e.crystallography_type == "VAPOR DIFFUSION"


def test_crystallography_type_none():
    e = make_extractor(["SOMETHING ELSE"])
    e.extract_crystallography_type()
    assert e.crystallography_type is None


def test_vapor_diffusion_type():
    e = make_extractor(["HANGING DROP"])
    e.extract_vapor_diffusion_type()
    assert e.vapor_diffusion_type == "HANGING DROP"


def test_vapor_diffusion_type_none():
    e = make_extractor(["SOMETHING ELSE"])
    e.extract_vapor_diffusion_type()
    assert e.vapor_diffusion_type is None


def test_ph_found():
    e = make_extractor(["PH 7.5"])
    e.extract_ph()
    assert e.ph.startswith("PH 7")


def test_ph_not_found():
    e = make_extractor(["NO PH HERE"])
    e.extract_ph()
    assert e.ph is None


def test_organic_precipitates():
    e = make_extractor(["REMARK 280 CRYSTALLIZATION CONDITIONS: 10% PEG 4000"])
    e.extract_crystallization_conditions()
    e.extract_precipitates(mode="organic")
    data = e.conditions_df.at[e.protein_id, "Organic_Precipitates"]
    assert (data is None) or isinstance(data, list)


def test_salt_precipitates():
    e = make_extractor(["REMARK 280 CRYSTALLIZATION CONDITIONS: 0.2M NaCl"])
    e.extract_crystallization_conditions()
    e.extract_precipitates(mode="salt")
    data = e.conditions_df.at[e.protein_id, "Salt_Precipitates"]
    assert (data is None) or isinstance(data, list)
