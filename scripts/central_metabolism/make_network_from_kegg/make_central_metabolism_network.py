# %%
# make a model network from kgml files
# merge mmm00010, mmu00020, mmu00030
# replace pyruvate->acCoA, 2-oxoglu->sucCoA with one reaction
# delete small compounds;H2O, ATP, NAD+, NADH, NADPH, NADP+, ADP, orthophasphate, CoA, CO2, AMP, GDP, GTP, H+, ITP, and IDP
# delete protein histidine, protein p-his, quinone, and hydroquinone
import pickle
import csv

from src import ReactionNetwork
from src import read_kgml
# %%
pathnames = ['mmu00010', 'mmu00020', 'mmu00030']
reaction_list_merge = []
for pathname in pathnames:
    reaction_list = read_kgml.get_reactionData(
        f'../../../data/central_metabolism/kgml/{pathname}.xml')
    reaction_list_merge += reaction_list

# remove duplicates
picked = []
for reac in reaction_list_merge:
    if reac not in picked:
        picked.append(reac)
reaction_list_merge = picked
print(len(reaction_list_merge))
# %%
'''
remove intermediate reactions

pyruvate->acCoA
 rn:R00014, rn:R03270, rn:07618, rn:R07618_2, rn:R02569, rn:R02569_2
 -> rn:R00209 C00022 + C00010 + C00003 <=> C00024 + C00011 + C00004 + C00080

2-oxo -> SucCoA
 rn:R00621, rn:R03316, rn:R07618, rn:R007618_2, rn:R02570, rn:R02570_2,
 -> rn:R08549 C00026 + C00010 + C00003 <=> C00091 + C00011 + C00004 + C00080
'''

reactions_removed = ['rn:R00014', 'rn:R03270', 'rn:R07618', 'rn:R07618_2', 'rn:R02569', 'rn:R02569_2',
                     'rn:R00621', 'rn:R03316', 'rn:R02570', 'rn:R02570_2']

reaction_list_noInterReac = []
for reac in reaction_list_merge:
    if reac[0] not in reactions_removed:
        # print(reac[0])
        reaction_list_noInterReac.append(reac)

newreac1 = ['rn:R00209', ['C00022', 'C00010', 'C00003'],
            ['C00024', 'C00011', 'C00004', 'C00080']]
newreac2 = ['rn:R08549', ['C00026', 'C00010', 'C00003'],
            ['C00091', 'C00011', 'C00004', 'C00080']]
reaction_list_noInterReac += [newreac1, newreac2]
print(len(reaction_list_noInterReac))
# %%
# some reactions are imcomplete
reaction_list_completed = []
from tqdm import tqdm
for reac in tqdm(reaction_list_noInterReac):
    reac_comp = read_kgml.complete_reaction(reac)
    reaction_list_completed.append(reac_comp)
print(len(reaction_list_completed))
# %%
# remove small molecules and proteins
cpds_removed = ['C00001', 'C00002', 'C00003', 'C00004', 'C00005', 'C00006', 'C00008', 'C00009', 'C00010', 'C00011', 'C00020',
                'C00035', 'C00044', 'C00080', 'C00081', 'C00104',
                'C15602', 'C15603', 'C00615', 'C04262']
reaction_list_noSmallCpds = []
for reac in reaction_list_completed:
    lhs = [c for c in reac[1] if c not in cpds_removed]
    rhs = [c for c in reac[2] if c not in cpds_removed]
    reaction_list_noSmallCpds.append([reac[0], lhs, rhs])
picked = []
for reac in reaction_list_noSmallCpds:
    if reac not in picked:
        picked.append(reac)
reaction_list_noSmallCpds = picked
print(len(reaction_list_noSmallCpds))
# %%

# %%
# delete duplicated reactions
reaction_list_nodup = []
exist = []
for reac in reaction_list_noSmallCpds:
    if [reac[1], reac[2]] not in exist:
        reaction_list_nodup.append(reac)
        exist.append([reac[1], reac[2]])
print(len(reaction_list_nodup))
# %%
'''
# There are two types of rn:R01070 (F1,6P <-> GlyP + G3P and F1,6P <-> G3P) 
# Only the former is used.
# '''
# reaction_list_identical=reaction_list_nodup.copy()
# reaction_list_identical.remove(['rn:R01070',['cpd:C05378'],['cpd:C00118']])
# reaction_list_identical.remove(['rn:R01070_2',['cpd:C00118'],['cpd:C05378']])
# %%
# sort by reaction name
reaction_list_sorted = []
reactionNames = sorted([reac[0] for reac in reaction_list_nodup])
for name in reactionNames:
    for reac in reaction_list_nodup:
        if reac[0] == name:
            reaction_list_sorted.append(reac)
            break


# %%
# rename cpd ids to readable names
with open('../../../data/central_metabolism/kgml/cpddata_dict3.pkl', 'rb') as f:
    cpddata_dict = pickle.load(f)

cpddata_dict['out'] = 'out'
reaction_list_renamed = []
for reac in reaction_list_sorted:
    lhs = [cpddata_dict['cpd:'+cpdid] for cpdid in reac[1]]
    rhs = [cpddata_dict['cpd:'+cpdid] for cpdid in reac[2]]
    reac_renamed = [reac[0], lhs, rhs]
    reaction_list_renamed.append(reac_renamed)

# %%
# add outflows add dead-end metabolites
reaction_list_renamed.append(['in_G1_5L', ['out'], ['D-Glucono-1,5-lactone']])
reaction_list_renamed.append(['in_Glycerate', ['out'], ['D-Glycerate']])
reaction_list_renamed.append(['out_Diphosphate', ['Diphosphate'], ['out']])
reaction_list_renamed.append(['in_HCO3-', ['out'], ['HCO3-']])
# add reactions for the convergence of the dynamics
reaction_list_renamed.append(['out_AcCoA', ['Acetyl-CoA'], ['out']])
reaction_list_renamed.append(
    ['in_DR5P', ['out'], ['2-Deoxy-D-ribose 5-phosphate']])
reaction_list_renamed.append(
    ['out_G3P', ['D-Glyceraldehyde 3-phosphate'], ['out']])
reaction_list_renamed.append(['in_Glu', ['out'], ['alpha-D-Glucose']])
reaction_list_renamed.append(['out_Pyruvate', ['Pyruvate'], ['out']])
# %%
reaction_list_nospace = []
for reac in reaction_list_renamed:
    reacName = reac[0].replace(' ', '_').replace(',', '_')
    lhs = [cpd.replace(' ', '_').replace(',', '_') for cpd in reac[1]]
    rhs = [cpd.replace(' ', '_').replace(',', '_') for cpd in reac[2]]
    reaction_list_nospace.append([reacName, lhs, rhs])

# %%
# shorten cpd names
dict_shortname = dict()
with open('../../../data/central_metabolism/kgml/cpd_shortname.csv', 'r') as f:
    reader = csv.reader(f)
    for row in reader:
        dict_shortname[row[0]] = row[1]
# %%
reaction_list_shortname = []
for reac in reaction_list_nospace:
    print(reac[0])
    lhs = [dict_shortname[cpd] for cpd in reac[1]]
    rhs = [dict_shortname[cpd] for cpd in reac[2]]
    reaction_list_shortname.append([reac[0], lhs, rhs])
# %%
with open('../../../data/central_metabolism/kgml/reaction_list_renamed.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    for reac in reaction_list_shortname:
        lhs = ' '.join(reac[1])
        rhs = ' '.join(reac[2])
        writer.writerow([reac[0], lhs, rhs])

# %%
