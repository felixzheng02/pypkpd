from matpy.matrix import Matrix

def names(input):

    if type(input) is not Matrix:
        return None
    if type(input) is Matrix:
        return input.get_datanam()