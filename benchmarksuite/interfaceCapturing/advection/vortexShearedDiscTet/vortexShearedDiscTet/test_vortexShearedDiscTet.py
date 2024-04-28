import os
from numpy import logspace
import pytest
import oftest
from oftest import run_reset_case

def test_vortexShearedDiscTet(run_reset_case):
    logs = oftest.path_logs()
    assert len(logs) > 0
    for log in logs:
        print(f"current log file:", log)
        if "gmshv493" in log:
            lines =  oftest.tail_log(log)
            foundKeyWord = False
            for line in lines:
                if "Done writing" in line:
                    foundKeyWord = True
            assert foundKeyWord
        else:
            assert oftest.case_status(log) == "completed"  # checks if run completes

