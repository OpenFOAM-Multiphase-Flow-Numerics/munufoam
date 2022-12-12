import os
from numpy import logspace
import pytest
import oftest
from oftest import run_reset_case


def simMod(interfaceCapturingMethod, interfaceRepresentation):
    dir_name = os.path.dirname(os.path.abspath(__file__))
    filemod = {
        "system/simulationParameter": [
            ("interfaceCapturingMethod", interfaceCapturingMethod),
            ("interfaceRepresentation", interfaceRepresentation),
        ]
    }
    # case_mod = oftest.Case_modifiers(filemod, dir_name)
    meta_data = {
        "interfaceCapturingMethod": interfaceCapturingMethod,
        "interfaceRepresentation": interfaceRepresentation,
    }
    case_mod = oftest.Case_modifiers(filemod, dir_name, meta_data)
    return case_mod


isoAdvection_plicRDF = simMod("isoAdvection", "plicRDF")
isoAdvection_isoAlpha = simMod("isoAdvection", "isoAlpha")
isoAdvection_gradAlpha = simMod("isoAdvection", "gradAlpha")

parameters = [isoAdvection_plicRDF, isoAdvection_isoAlpha, isoAdvection_gradAlpha]


@pytest.mark.parametrize(
    "run_reset_case",
    parameters,
    indirect=["run_reset_case"],
    ids=["isoAdvection_plicRDF", "isoAdvection_isoAlpha", "isoAdvection_gradAlpha"],
)
def test_damBreak(run_reset_case):
    logs = oftest.path_logs()
    assert len(logs) > 0
    for log in logs:
        assert oftest.case_status(log) == "completed"  # checks if run completes
