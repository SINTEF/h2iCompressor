def pressureOverUnderEstimate(pressureTest, etaStage):
    """ Describe function """
    if  pressureTest > 0:
        etaST = etaStage - 0.00005    # todo: need better way of iterating etaStage
    elif  pressureTest < 0:
        etaST = etaStage + 0.00005
    else:
        print(pressureTest)
        etaST = etaStage
    return etaST