def pressureOverUnderEstimate(pressureTest, etaStage):
    """ Describe function """
    # Import 
    import numpy as np
    if pressureTest == np.nan:
        debug = 0                   # breakpoint for debugging
    if  pressureTest > 0:
        etaST = etaStage - 0.003    # todo: need better way of iterating etaStage
    elif  pressureTest < 0:
        etaST = etaStage + 0.003

    return etaST