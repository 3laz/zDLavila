from __future__ import print_function
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors

#Lokacija .sdf datoteke z zapisi molekul
file = r"C:\*****\******\*******\*****\molekule.sdf"

#Pandas Dataframe object iz .sdf datoteke
BRDLigs = PandasTools.LoadSDF(file)
BRDLigs = BRDLigs.drop(["ROMol", "ID"], axis=1)


#Loopanje čez SMILES zapise molekul ter izračun deskriptorjev za posamezne molekule
j = 0
for i in BRDLigs["smiles"]:
    
    #2D struktura
    molekula_smiles =  Chem.MolFromSmiles(i)
    BRDLigs.at[j, "TPSA"] = Descriptors.TPSA(molekula_smiles) #Topological surface area
    BRDLigs.at[j, "logP"] = Descriptors.MolLogP(molekula_smiles) #LogP
    BRDLigs.at[j, "MW"] = Descriptors.MolWt(molekula_smiles) #Molecular Weight
    BRDLigs.at[j, "Hdonors"] = Descriptors.NumHDonors(molekula_smiles) #H donors
    BRDLigs.at[j, "Hacceptors"] = Descriptors.NumHAcceptors(molekula_smiles) #H acceptors
    BRDLigs.at[j, "Rotbonds"] = Descriptors.NumRotatableBonds(molekula_smiles) #Rotating bonds
    
    #3D struktura    
    molekula_smiles2 = Chem.AddHs(molekula_smiles) #Doda H atome na molekulo
    Chem.AllChem.EmbedMolecule(molekula_smiles2) #Nevem kaj naredi ampak rabiš
    Chem.AllChem.MMFFOptimizeMolecule(molekula_smiles2) #Optimizacija nardi 2-3% razlike v vrednosti če jo uporabis
    BRDLigs.at[j, "Molvolume"] = Chem.AllChem.ComputeMolVolume(molekula_smiles2) #molekulski volumen
    
    j += 1

print(BRDLigs.columns)
print(BRDLigs)
