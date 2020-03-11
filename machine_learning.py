from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors


file = r"C:\Users\blaz-\PycharmProjects\lek_ml\molekule.sdf"
"""
mol = pybel.readfile("fa", file)
print(mol)
for i in mol:
    print(i)

suppl = Chem.SDMolSupplier(file)

for i in range(10):
    print(suppl[i].GetNumAtoms())
    ms = [x for x in suppl if x is not None]
    for m in ms: tmp=AllChem.Compute2DCoords(m)
"""
#Draw.MolToFile(ms[0], r"C:\Users\blaz-\PycharmProjects\lek_ml\molekul1.png")
#Draw.MolToFile(ms[1], r"C:\Users\blaz-\PycharmProjects\lek_ml\molekula2.png")


BRDLigs = PandasTools.LoadSDF(file)
#print(BRDLigs.columns)
BRDLigs = BRDLigs.drop(["ROMol", "ID"], axis=1)

j = 0
for i in BRDLigs["smiles"]:
    molekula_smiles =  Chem.MolFromSmiles(i)
    BRDLigs.at[j, "TPSA"] = Descriptors.TPSA(molekula_smiles)
    BRDLigs.at[j, "logP"] = Descriptors.MolLogP(molekula_smiles)
    BRDLigs.at[j, "MW"] = Descriptors.MolWt(molekula_smiles)
    BRDLigs.at[j, "Hdonors"] = Descriptors.NumHDonors(molekula_smiles)
    BRDLigs.at[j, "Hacceptors"] = Descriptors.NumHAcceptors(molekula_smiles)
    BRDLigs.at[j, "Rotbonds"] = Descriptors.NumRotatableBonds(molekula_smiles)

    j += 1

print(BRDLigs.columns)
print(BRDLigs)
