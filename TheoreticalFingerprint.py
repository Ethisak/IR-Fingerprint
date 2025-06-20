import pandas as pd
from rdkit import Chem
import numpy as np
from array import array
import json

d = pd.read_excel("")

def fpteor(smi):
  fp = [0] * 101
  mol = Chem.MolFromSmiles(smi)

#fp1
  fp[0] = 1

#fp2
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H]')):
    fp[1] = 1

#fp3
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[Si][OX2H]')):
    fp[2] = 1

#fp4+5
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H][CX4,c;!$(C([OX2H])[O,S,#7,#15])]')):
    fp[3] = 1
    fp[4] = 1
  elif mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H1][C]=[C]')):
    fp[3] = 1
    fp[4] = 1
  elif mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H1][C]#[C]')):
    fp[3] = 1
    fp[4] = 1

#fp6
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]=[N][OX2H1]')):
    fp[5] = 1

#fp7
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3](=O)[OX2H1]')):
    fp[6] = 1

#fp8
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[Nh]')):
    fp[7] = 1

#fp9
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3H2][CX3](=[OX1])')):
    fp[8] = 1

#fp10
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6][NX3H1][CX3](=[OX1])')):
    fp[9] = 1

#fp11
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]')):
    fp[10] = 1

#fp12
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3;H1;!$(NC=[!#6]);!$(NC#[!#6])][#6]')):
    fp[11] = 1

#fp13
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3H2][CX3](=[OX1])')):
    fp[12] = 1

#fp14
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6][NX3H1][CX3](=[OX1])')):
    fp[13] = 1

#fp15
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6h]')):
    fp[14] = 1

#fp16
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[Ch]#[C]')):
    fp[15] = 1

#fp17
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[Ch]=[C]')):
    fp[16] = 1

#fp18
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[Ch;R]')):
    fp[17] = 1

#fp19
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[ch;R]')):
    fp[18] = 1

#fp20
  if mol.HasSubstructMatch(Chem.MolFromSmarts('n1[ch]cccc1')):
    fp[19] = 1

#fp21
  if mol.HasSubstructMatch(Chem.MolFromSmarts('n1[ch]cncc1')):
    fp[20] = 1

#fp22
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[nh1]1[ch]ccc1')):
    fp[21] = 1

#fp23
  if mol.HasSubstructMatch(Chem.MolFromSmarts('o1[ch]ccc1')):
    fp[22] = 1

#fp24
  if mol.HasSubstructMatch(Chem.MolFromSmarts('s1[ch]ccc1')):
    fp[23] = 1

#fp25
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6][CH1]([#6])[#6]')):
    fp[24] = 1

#fp26
  if mol.HasSubstructMatch(Chem.MolFromSmarts('*[Ch]*')):
    fp[25] = 1

#fp27
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3H1](=O)[#6]')):
    fp[26] = 1
  elif mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3H2]=O')):
    fp[26] = 1

#fp28
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[SX2H1]')):
    fp[27] = 1

#fp29
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]C#N')):
    fp[28] = 1

#fp30
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]C#N')):
    fp[29] = 1

#fp31
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#7]#[#16]')):
    fp[30] = 1

#fp32
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]#[#16]')):
    fp[31] = 1

#fp33
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]=[#16]=[#8]')):
    fp[32] = 1

#fp34
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#7]=[#6]=[#8]')):
    fp[33] = 1

#fp35
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]=[#16]=[#7]')):
    fp[34] = 1

#fp36
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#7]=[#6]=[#16]')):
    fp[35] = 1

#fp37
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]=[#7]=[#8]')):
    fp[36] = 1

#fp38
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#7]=[#16]=[#8]')):
    fp[37] = 1

#fp39
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[C]#[C]')):
    fp[38] = 1

#fp40
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[C]=[C]=[C]')):
    fp[39] = 1

#fp41
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]=[#8]')):
    fp[40] = 1

#fp42
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3](=[OX1])[Cl,F,I,Br]')):
    fp[41] = 1

#fp43
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6X3H0;R](=O)[#8X2H0;R]')):
    fp[42] = 1

#fp44
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3](=O)[OX2H0][#6]')):
    fp[43] = 1

#fp45
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3H1](=O)[#6]')):
    fp[44] = 1
  elif mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3H2]=O')):
    fp[44] = 1

#fp46
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6][#6X3](=O)[#6]')):
    fp[45] = 1

#fp47
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3](=O)[OX2H1]')):
    fp[46] = 1

#fp48
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6][NX3H1][CX3](=[OX1])')):
    fp[47] = 1

#fp49
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3H2][CX3](=[OX1])')):
    fp[48] = 1

#fp50
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]~[#6]~[#9]')):
    fp[49] = 1

#fp51
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[C]=[C]')):
    fp[50] = 1

#fp52
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[C]=[C]')):
    fp[51] = 1

#fp53
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[C&R1]=[C&R1]')):
    fp[52] = 1

#fp54
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]O[CH1]=[CH2]')):
    fp[53] = 1

#fp55
  if mol.HasSubstructMatch(Chem.MolFromSmarts('*~@:~*(@:*)C=C')):
    fp[54] = 1

#fp56
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]')):
    fp[55] = 1

#fp57
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]')):
    fp[56] = 1

#fp58
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]')):
    fp[57] = 1

#fp59
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[Nh]')):
    fp[58] = 1

#fp60
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6]')):
    fp[59] = 1

#fp61
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3;H1;!$(NC=[!#6]);!$(NC#[!#6])][#6]')):
    fp[60] = 1

#fp62
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3H2][CX3](=[OX1])')):
    fp[61] = 1

#fp63
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6][NX3H1][CX3](=[OX1])')):
    fp[62] = 1

#fp64
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6X3;R1](=O)[#7X3H1;R1]')):
    fp[63] = 1

#fp65
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3;$([C]([#6])[#6]),$([CH][#6])]=[N;$([NX2H1]),$([NX2][#6])]')):
    fp[64] = 1

#fp66
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]=[N][OX2H1]')):
    fp[65] = 1

#fp67
  if mol.HasSubstructMatch(Chem.MolFromSmarts('c1ccccc1')):
    fp[66] = 1

#fp68
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[Ch]')):
    fp[67] = 1

#fp69
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[Ch;R]')):
    fp[68] = 1

#fp70
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]-*(-[#6])(~*)(~*)')):
    fp[69] = 1

#fp71
  if mol.HasSubstructMatch(Chem.MolFromSmarts('*[Ch]*')):
    fp[70] = 1

#fp72
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[Ch]#C')):
    fp[71] = 1

#fp73
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H]')):
    fp[72] = 1

#fp74
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#16]=[#8]')):
    fp[73] = 1

#fp75
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[$([#16X4](=[OX1])=[OX1]),$([#16X4+2]([OX1-])[OX1-])]')):
    fp[74] = 1

#fp76
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[$([#16X4](=[OX1])=[OX1]),$([#16X4+2]([OX1-])[OX1-])]')):
    fp[75] = 1

#fp77
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]')):
    fp[76] = 1

#fp78
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[$([C][#16X3]=[OX1]),$([C][#16X3+][OX1-])]')):
    fp[77] = 1

#fp79
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[$([c][#16X3]=[OX1]),$([c][#16X3+][OX1-])]')):
    fp[78] = 1

#fp80
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]-[#8]')):
    fp[79] = 1

#fp81
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[CX3](=O)[OX2H1]')):
    fp[80] = 1

#fp82
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[c][CX3](=O)[OX2H0][#6]')):
    fp[81] = 1

#fp83
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H][CX4,c;!$(C([OX2H])[O,S,#7,#15])]')):
    fp[82] = 1
  elif mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H1][C]=[C]')):
    fp[82] = 1
  elif mol.HasSubstructMatch(Chem.MolFromSmarts('[OX2H1][C]#[C]')):
    fp[82] = 1

#fp84
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[C]O[C]')):
    fp[83] = 1

#fp85
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]-[#7]')):
    fp[84] = 1

#fp86
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3H2][CX3](=[OX1])')):
    fp[85] = 1

#fp87
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6][NX3H1][CX3](=[OX1])')):
    fp[86] = 1

#fp88
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3;v3;H2,H1,H0;!$(NC=[!#6]);!$(NC#[!#6])][c]')):
    fp[87] = 1

#fp89
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3;v3;H2,H1,H0;!$(NC=[!#6]);!$(NC#[!#6])][C]')):
    fp[88] = 1

#fp90
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[NX3][CX3]=[SX1]')):
    fp[89] = 1

#fp91
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]-[#9]')):
    fp[90] = 1

#fp92
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#14]-[#8]-[#14]')):
    fp[91] = 1

#fp93
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#15]-[#8]-[C]')):
    fp[92] = 1

#fp94
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]=[N][OX2H1]')):
    fp[93] = 1

#fp95
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[C]1-O-[C]1')):
    fp[94] = 1

#fp96
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#14]-[#6]')):
    fp[95] = 1

#fp97
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]-[#17,#53,#16,#35]')):
    fp[96] = 1

#fp98
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]-[#17]')):
    fp[97] = 1

#fp99
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]-[#16]')):
    fp[98] = 1

#fp100
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]-[#35]')):
    fp[99] = 1

#fp101
  if mol.HasSubstructMatch(Chem.MolFromSmarts('[#6]-[#53]')):
    fp[100] = 1

  globals()['FP'] = fp

dd = range(1278)

for ii in dd:
  fpteor(d['SMILES'][ii])
  fk = [i for i, x in enumerate(FP) if x == 1]
  fk = [x + 1 for x in fk]
  np.savetxt(str(int(ii + 1)) + ".fp", fk)

