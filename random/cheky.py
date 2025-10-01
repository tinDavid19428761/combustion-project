# cheky.py
from pyomo.environ import *

class MyBlock(Block):
    def build(self):
        self.x = Var(initialize=2)
        
        @Expression(self)
        def doubled_x(b):
            return 2 * b.x

m = ConcreteModel()
m.block = MyBlock()
value(m.block.doubled_x)