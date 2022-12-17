import path
from project.models import Models
from project.sfg import Sfg

def call_func(func, *args, **kwargs):
    """
    func could be a function object or a string
    """
    arguments = list(locals().values())[1:]
    if callable(func):
        return func(*args, **kwargs)
    elif type(func) is str:
        try:
            method = getattr(Models, func)
            return method(args, kwargs)
        except AttributeError:
            method = getattr(Sfg, func)
            return method(args, kwargs)
        except:
            raise Exception("No function found.")