"""
## @param file.name A function or a string that is the name of a function. 
## @param ... Arguments for the function.  Multiple arguments separated by a comma.
## 
## @return Output from the defined function.
## 
## @family MATLAB
## @export
## @keywords internal
## 
## Function written to match MATLAB function feval()
## Author: Caiya Zhang, Yuchen Zheng
"""


def feval (func_name, *argv):
    if callable(func_name):
        return eval(func_name.__name__ + str(argv))
    elif type(func_name) is str:
        return eval(func_name + str(argv))




    