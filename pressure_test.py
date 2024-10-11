def pressureOverUnderEstimate(pressureTest, etaStage):
    """ Describe function """
    if  pressureTest > 0:
        etaST = etaStage - 0.0001    # todo: need better way of iterating etaStage
    elif  pressureTest < 0:
        etaST = etaStage + 0.0001

    return etaST