import path
from matpy.matrix import Matrix
from matpy.num import Num

class Visualization:

    def dict_visualization(self, input: dict):
        for key, value in input.items():
            print(key + ":")
            if type(value) is Matrix:
                print(value.get_data())
            elif type(value) is Num:
                print(value.get_value())
            elif type(value) is dict:
                self.dict_visualization(value)
            else:
                print(value)