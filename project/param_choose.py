"""
For parameter choose in create_poped_database()

Author: Caiya Zhang, Yuchen Zheng
"""


from project.pypkpd_choose import pypkpd_choose


def param_choose(param, default, alternative):
    """
    if param is defined (param is not None), return param
    otherwise, pypkpd_choose(default, alternative)
    """
    if param is None:
        return pypkpd_choose(default, alternative)
    return param