from enum import Enum
from .GENIEWrapper import GENIEWrapper
from .ARACNEWrapper import ARACNEWrapper
from .LSCONWrapper import LSCONWrapper

class AlgorithmSelector(Enum):
    """Enum specifying which network inference algorithm should be used."""
    GENIE3 = 'GENIE3'
    ARACNe = 'ARACNe'
    LSCON = 'LSCON'

    def __str__(self):
        return self.value

def get_algorithm_wrapper(algorithm_selector):
    """Returns the appropriate algorithm based on the selection.
    Parameters
    ----------
    algorithm_selector : AlgorithmSelector
        Specifies which algorithm should be used.
    """
    if algorithm_selector == AlgorithmSelector.GENIE3:
        return GENIEWrapper()
    elif algorithm_selector == AlgorithmSelector.ARACNe:
        return ARACNe()
    elif algorithm_selector == AlgorithmSelector.LSCON:
        return LSCON()