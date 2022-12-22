from matpy.matrix import Matrix

def names(input):

    if type(input) is Matrix:
        if input.get_datanam() is not None:
            return input.get_datanam()
        elif input.get_axisnam(1) is not None:
            return input.get_axisnam(1)
        elif input.get_axisnam(0) is not None:
            return input.get_axisnam(0)
    elif type(input) is dict:
        return input.keys()
    return None