"""Representative mammalian glycans; diagrams describe references, not MS assignments.

Masses: underivatized monoisotopic residue formulas. Sources and conventions are
exposed by the catalog API. A free reducing glycan adds one water to the sum.
"""
from collections import Counter

from masses import WATER, formula_mass
SUGARS = {
    'Gal': dict(name='Galactose', formula={'C':6,'H':10,'O':5}, color='#f3cf37', shape='circle'),
    'Man': dict(name='Mannose', formula={'C':6,'H':10,'O':5}, color='#55a660', shape='circle'),
    'GlcNAc': dict(name='N-acetylglucosamine', formula={'C':8,'H':13,'N':1,'O':5}, color='#4299d4', shape='square'),
    'GalNAc': dict(name='N-acetylgalactosamine', formula={'C':8,'H':13,'N':1,'O':5}, color='#f3cf37', shape='square'),
    'Fuc': dict(name='Fucose', formula={'C':6,'H':10,'O':4}, color='#e66a5d', shape='triangle'),
    'Neu5Ac': dict(name='N-acetylneuraminic acid (sialic acid)', formula={'C':11,'H':17,'N':1,'O':8}, color='#b780be', shape='diamond'),
}
for sugar in SUGARS.values():
    sugar['residue_mass'] = float(formula_mass(sugar['formula']))

O_SOURCE = 'https://www.ncbi.nlm.nih.gov/books/NBK579921/'
N_SOURCE = 'https://www.ncbi.nlm.nih.gov/books/NBK579964/'

def node(sugar, parent=None, linkage=''):
    return dict(sugar=sugar,parent=parent,linkage=linkage)

def o_core(core=1, sialic=0):
    nodes = [node('GalNAc')]
    if core in (1,2):
        nodes.append(node('Gal',0,'β1–3'))
    if core == 2:
        nodes.append(node('GlcNAc',0,'β1–6'))
    if core == 3:
        nodes.append(node('GlcNAc',0,'β1–3'))
    if sialic:
        nodes.append(node('Neu5Ac',1,'α2–3'))
    if sialic == 2:
        if core == 2:
            nodes.append(node('Gal',2,'β1–4'))
            nodes.append(node('Neu5Ac',len(nodes)-1,'α2–3'))
        else:
            nodes.append(node('Neu5Ac',0,'α2–6'))
    return nodes

def n_core():
    return [node('GlcNAc'),node('GlcNAc',0,'β1–4'),node('Man',1,'β1–4'),node('Man',2,'α1–3'),node('Man',2,'α1–6')]

def mannose(count):
    nodes = n_core()
    # Man5: the α1–6 arm branches into α1–3 and α1–6 terminal mannoses.
    nodes.extend([node('Man',4,'α1–3'),node('Man',4,'α1–6')])
    if count == 9:
        nodes.extend([node('Man',3,'α1–2'),node('Man',7,'α1–2'),node('Man',5,'α1–2'),node('Man',6,'α1–2')])
    return nodes

def complex_glycan(gal=0, sialic=0, fucose=False, antennae=2):
    nodes = n_core()
    tips = []
    for parent,link in [(3,'β1–2'),(4,'β1–2')]+([(3,'β1–4')] if antennae == 3 else []):
        tips.append(len(nodes))
        nodes.append(node('GlcNAc',parent,link))
    gal_tips = []
    for parent in tips[:gal]:
        gal_tips.append(len(nodes))
        nodes.append(node('Gal',parent,'β1–4'))
    for parent in gal_tips[:sialic]:
        nodes.append(node('Neu5Ac',parent,'α2–6'))
    if fucose:
        nodes.append(node('Fuc',0,'α1–6'))
    return nodes

def entry(identifier, kind, name, nodes, note=''):
    counts = Counter(n['sugar'] for n in nodes)
    mass = sum(SUGARS[s]['residue_mass']*n for s,n in counts.items())
    return dict(id=identifier,type=kind,name=name,nodes=nodes,composition=dict(counts),
                attached_mass=mass,free_mass=mass+WATER,note=note,
                source=O_SOURCE if kind=='O' else N_SOURCE)

CATALOG = [
    entry('o_tn','O','Tn antigen', [node('GalNAc')]),
    entry('o_core1','O','Core 1 / T antigen',o_core()),
    entry('o_sialyl_tn','O','Sialyl-Tn', [node('GalNAc'),node('Neu5Ac',0,'α2–6')]),
    entry('o_sialyl_core1','O','Sialylated core 1',o_core(sialic=1)),
    entry('o_disialyl_core1','O','Disialylated core 1',o_core(sialic=2)),
    entry('o_core2','O','Core 2',o_core(2)),
    entry('o_disialyl_core2','O','Disialylated extended core 2',o_core(2,2)),
    entry('o_core3','O','Core 3',o_core(3)),
    entry('o_glcnac','O','O-GlcNAc', [node('GlcNAc')], 'Isobaric with Tn; different sugar and protein linkage.'),
    entry('n_core','N','Trimannosyl N-glycan core',n_core()),
    entry('n_man5','N','High mannose / Man5',mannose(5)),
    entry('n_man9','N','High mannose / Man9',mannose(9)),
    entry('n_g0','N','Biantennary / G0',complex_glycan()),
    entry('n_g0f','N','Core-fucosylated / G0F',complex_glycan(fucose=True)),
    entry('n_g1f','N','Core-fucosylated / G1F',complex_glycan(gal=1,fucose=True),'One representative arm is shown; positional isomers have the same mass.'),
    entry('n_g2','N','Digalactosylated / G2',complex_glycan(gal=2)),
    entry('n_g2f','N','Core-fucosylated / G2F',complex_glycan(gal=2,fucose=True)),
    entry('n_a2g2s1','N','Monosialylated / A2G2S1',complex_glycan(gal=2,sialic=1),'Representative α2–6 sialylation; linkage and arm isomers share this mass.'),
    entry('n_a2g2s2','N','Disialylated / A2G2S2',complex_glycan(gal=2,sialic=2),'Representative α2–6 sialylation; linkage isomers share this mass.'),
    entry('n_fa2g2s2','N','Fucosylated, disialylated / FA2G2S2',complex_glycan(gal=2,sialic=2,fucose=True),'Representative α2–6 sialylation; linkage isomers share this mass.'),
    entry('n_a3g3s3','N','Triantennary, trisialylated / A3G3S3',complex_glycan(gal=3,sialic=3,antennae=3),'Representative β1–4 third antenna and α2–6 sialylation; other isomers share this mass.'),
]
BY_ID = {g['id']:g for g in CATALOG}
BY_ID['o_glcnac']['source'] = 'https://www.ncbi.nlm.nih.gov/books/NBK1954/'

def model_mass(glycan, mode):
    return glycan['attached_mass'] if mode == 'mono_attached' else glycan['free_mass']

def catalog_payload():
    return dict(glycans=CATALOG,sugars=SUGARS,water=WATER,
                mass_source='https://web.expasy.org/glycomod/glycomod_masses.html',
                symbol_source='https://www.ncbi.nlm.nih.gov/glycans/snfg.html')
