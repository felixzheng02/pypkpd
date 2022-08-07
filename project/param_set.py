from project.util import get_dict_value


def param_set(param, pypkpdInput, *argv):
    if param is not None:
        return param
    if pypkpdInput is not None:
        return get_dict_value(pypkpdInput, *argv)
    raise Exception("Must provide pypkpdInput dict.")