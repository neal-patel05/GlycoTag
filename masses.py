"""Monoisotopic masses derived from molecular formulas without rounded AA tables.

Atomic isotope masses: NIST Atomic Weights and Isotopic Compositions v4.1,
https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?isotype=all
All published digits of 12C, 1H, 14N, 16O and 32S are retained. Decimal
arithmetic derives masses; integer units preserve every digit in peptide sums.
These computed digits do not imply equal experimental accuracy.
"""
from decimal import Decimal

ATOMS = {key: Decimal(value) for key,value in {
    'C':'12', 'H':'1.00782503223', 'N':'14.00307400443',
    'O':'15.99491461957', 'S':'31.9720711744',
}.items()}
SCALE = 10**11

def formula_mass(composition):
    return sum((ATOMS[atom]*count for atom,count in composition.items()), Decimal(0))

WATER_DECIMAL = Decimal('18.010564684')
WATER = float(WATER_DECIMAL)
# Residue formulas: free amino acid minus H2O.
RESIDUE_FORMULAS = {
    'A':{'C':3,'H':5,'N':1,'O':1},
    'R':{'C':6,'H':12,'N':4,'O':1},
    'N':{'C':4,'H':6,'N':2,'O':2},
    'D':{'C':4,'H':5,'N':1,'O':3},
    'C':{'C':3,'H':5,'N':1,'O':1,'S':1},
    'E':{'C':5,'H':7,'N':1,'O':3},
    'Q':{'C':5,'H':8,'N':2,'O':2},
    'G':{'C':2,'H':3,'N':1,'O':1},
    'H':{'C':6,'H':7,'N':3,'O':1},
    'I':{'C':6,'H':11,'N':1,'O':1},
    'L':{'C':6,'H':11,'N':1,'O':1},
    'K':{'C':6,'H':12,'N':2,'O':1},
    'M':{'C':5,'H':9,'N':1,'O':1,'S':1},
    'F':{'C':9,'H':9,'N':1,'O':1},
    'P':{'C':5,'H':7,'N':1,'O':1},
    'S':{'C':3,'H':5,'N':1,'O':2},
    'T':{'C':4,'H':7,'N':1,'O':2},
    'W':{'C':11,'H':10,'N':2,'O':1},
    'Y':{'C':9,'H':9,'N':1,'O':2},
    'V':{'C':5,'H':9,'N':1,'O':1},
}
RESIDUE_DECIMAL = {aa:formula_mass(formula) for aa,formula in RESIDUE_FORMULAS.items()}
FREE_DECIMAL = {aa:mass+WATER_DECIMAL for aa,mass in RESIDUE_DECIMAL.items()}
MONO = {aa:float(mass) for aa,mass in RESIDUE_DECIMAL.items()}
FREE_AA = {aa:float(mass) for aa,mass in FREE_DECIMAL.items()}
RESIDUE_UNITS = {aa:int(mass*SCALE) for aa,mass in RESIDUE_DECIMAL.items()}
