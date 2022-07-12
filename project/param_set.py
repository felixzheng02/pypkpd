def param_set(param, pypkpdInput, key_1, key_2=None, key_3=None):
    if param is not None:
        return param
    if pypkpdInput is not None:
        if key_1 is not None:
            if key_1 in pypkpdInput.keys():
                if key_2 is not None:
                    if key_2 in pypkpdInput[key_1].keys():
                        if key_3 is not None:
                            if key_3 in pypkpdInput[key_1][key_2].keys():
                                return pypkpdInput[key_1][key_2][key_3]
                            else: raise Exception("key_3 not in the dict keys.")
                        else: return pypkpdInput[key_1][key_2]
                    else: raise Exception("key_2 not in the dict keys.")
                else: return pypkpdInput[key_1]
            else: raise Exception("key_1 not in the dict keys.")
    raise Exception("Must provide pypkpdInput dict.")