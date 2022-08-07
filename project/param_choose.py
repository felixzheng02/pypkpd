"""
For parameter choose in create_poped_database()

Author: Caiya Zhang, Yuchen Zheng
"""


from project.util import get_dict_value


def param_choose(param, pypkpdInput, replace, *argv):
    if param is not None:
        return param
    if pypkpdInput is not None:
        output = get_dict_value(pypkpdInput, argv)
        if output is None:
            return replace
        return output
    raise Exception("Must provide pypkpdInput dict.")