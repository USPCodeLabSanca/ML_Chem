from pyCheckmol import CheckMol
import numpy as np




smi = 'C(C(F)(F)F)(F)(F)F'
cm = CheckMol()
res = cm.functionalGroupSmiles(smiles=smi, isString=True, generate3D=False, justFGcode=False, returnDataframe=False,deleteTMP=False)
print(res['Functional Group'])

