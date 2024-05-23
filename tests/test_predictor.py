import pytest
import sys
from predictor import prediction


def test_prediction_returns_float():
    result = prediction("CCCCO")
    assert isinstance(result, float)

def test_prediction_handles_invalid_input():
    with pytest.raises(Exception):
        prediction("C1CC==C")

def test_prediction_handles_empty_input():
    with pytest.raises(Exception):
        prediction("")

def test_prediction_handles_whitespace_input():
    with pytest.raises(Exception):
        prediction("   ")

def test_prediction_complex_isomeric_smiles():
    # Check that the predictor can handle complex isomeric SMILES
    isomeric_smile = "CC[C@H](C)[C@H]1C(=O)N[C@@H](C(=O)N[C@@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)NCC(=O)N[C@@H](C(=O)N[C@H]2CSSC[C@@H]3C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@@H](C(=O)N[C@@H](CSSC[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N4CCC[C@@H]4C(=O)N5CCC[C@H]5C(=O)N[C@@H](C(=O)N[C@@H](C(=O)NCC(=O)N6CCC[C@H]6C(=O)N[C@H](CSSC[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N3)CC(=O)O)CCC(=O)O)C)CO)CCCCN)CC7=CC=CC=C7)CC(=O)N)CC(=O)N)CCCNC(=N)N)CCCCN)C)CCCNC(=N)N)NC(=O)CNC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](NC2=O)CCC(=O)N)[C@@H](C)O)CC8=CC=CC=C8)C(C)C)CC9=CC=C(C=C9)O)C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N1)CCCNC(=N)N)C)CCCCN)[C@H](C)O)CC1=CC=C(C=C1)O)CCC(=O)O)CC(C)C)NC(=O)[C@@H](CC1=CC=CC=C1)NC(=O)[C@H](CC(=O)O)NC(=O)[C@H]1CCCN1C(=O)[C@H](CCCNC(=N)N)N)C(=O)NCC(=O)NCC(=O)N[C@H](C)C(=O)O)[C@@H](C)O)CCCNC(=N)N)CCSC)CC(C)C)C)CCCCN)C)CC(=O)N)CC1=CC=C(C=C1)O)CC1=CC=CC=C1)CC1=CC=C(C=C1)O)CCCNC(=N)N)[C@H](C)CC"
    result = prediction(isomeric_smile)
    assert isinstance(result, float)

def test_prediction_isomeric_canonical_specific():
    # The predictor is isomer-specific, the prediction for the 2 following SMILES should be different
    isomeric_smile = "CC[C@H](C)[C@@H]1C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)NCC(=O)N[C@@H](C(=O)N[C@@H]2CSSC[C@@H]3C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@H](CSSC[C@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N4CCC[C@@H]4C(=O)N5CCC[C@H]5C(=O)N[C@H](C(=O)N[C@H](C(=O)NCC(=O)N6CCC[C@@H]6C(=O)N[C@@H](CSSC[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N[C@H](C(=O)N3)CC(=O)O)CCC(=O)O)C)CO)CCCCN)CC7=CC=CC=C7)CC(=O)N)CC(=O)N)CCCNC(=N)N)CCCCN)C)CCCNC(=N)N)NC(=O)CNC(=O)CNC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@H](NC2=O)CCC(=O)N)[C@H](C)O)CC8=CC=CC=C8)C(C)C)CC9=CC=C(C=C9)O)C(=O)N[C@@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N1)CCCNC(=N)N)C)CCCCN)[C@H](C)O)CC1=CC=C(C=C1)O)CCC(=O)O)CC(C)C)NC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@@H](CC(=O)O)NC(=O)[C@@H]1CCCN1C(=O)[C@@H](CCCNC(=N)N)N)C(=O)NCC(=O)NCC(=O)N[C@@H](C)C(=O)O)[C@H](C)O)CCCNC(=N)N)CCSC)CC(C)C)C)CCCCN)C)CC(=O)N)CC1=CC=C(C=C1)O)CC1=CC=CC=C1)CC1=CC=C(C=C1)O)CCCNC(=N)N)[C@H](C)CC"
    canonical_smile = "CCC(C)C1C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NCC(=O)NC(C(=O)NC2CSSCC3C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CSSCC(C(=O)NC(C(=O)NC(C(=O)N4CCCC4C(=O)N5CCCC5C(=O)NC(C(=O)NC(C(=O)NCC(=O)N6CCCC6C(=O)NC(CSSCC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)N3)CC(=O)O)CCC(=O)O)C)CO)CCCCN)CC7=CC=CC=C7)CC(=O)N)CC(=O)N)CCCNC(=N)N)CCCCN)C)CCCNC(=N)N)NC(=O)CNC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC2=O)CCC(=O)N)C(C)O)CC8=CC=CC=C8)C(C)C)CC9=CC=C(C=C9)O)C(=O)NC(C(=O)NC(C(=O)NC(C(=O)N1)CCCNC(=N)N)C)CCCCN)C(C)O)CC1=CC=C(C=C1)O)CCC(=O)O)CC(C)C)NC(=O)C(CC1=CC=CC=C1)NC(=O)C(CC(=O)O)NC(=O)C1CCCN1C(=O)C(CCCNC(=N)N)N)C(=O)NCC(=O)NCC(=O)NC(C)C(=O)O)C(C)O)CCCNC(=N)N)CCSC)CC(C)C)C)CCCCN)C)CC(=O)N)CC1=CC=C(C=C1)O)CC1=CC=CC=C1)CC1=CC=C(C=C1)O)CCCNC(=N)N)C(C)CC"
    result1 = prediction(isomeric_smile)
    result2 = prediction(canonical_smile)
    assert result1 != result2

def test_prediction_consistent_results():
    """
    During the development of the predictor, the model predicted sometimes different values for the same input, this has been fixed but the source of the issue is still unknown
    """
    result1 = prediction("CCCCO")
    result2 = prediction("CCCCO")
    assert result1 == result2

