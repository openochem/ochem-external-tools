/* *****************************************************************************
   Copyright (C) 1996-2022 by Eugene V. Radchenko

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
***************************************************************************** */

#include "struct.h"

#include <ctype.h>
#include <iomanip>
#include <math.h>

struct atomData {
    char name[3];
    int val;
};
const atomData atomDataNormal[] =
{
    { "H",1}, {"He",0}, {"Li",1}, {"Be",2}, { "B",3}, { "C",4}, { "N",3}, { "O",2}, { "F",1}, {"Ne",0},
    {"Na",1}, {"Mg",2}, {"Al",3}, {"Si",4}, { "P",3}, { "S",2}, {"Cl",1}, {"Ar",0}, { "K",1}, {"Ca",2},
    {"Sc",3}, {"Ti",4}, { "V",0}, {"Cr",0}, {"Mn",0}, {"Fe",0}, {"Co",0}, {"Ni",0}, {"Cu",0}, {"Zn",0},
    {"Ga",3}, {"Ge",4}, {"As",3}, {"Se",2}, {"Br",1}, {"Kr",0}, {"Rb",1}, {"Sr",2}, { "Y",0}, {"Zr",4},
    {"Nb",0}, {"Mo",0}, {"Tc",0}, {"Ru",0}, {"Rh",0}, {"Pd",0}, {"Ag",1}, {"Cd",0}, {"In",0}, {"Sn",4},
    {"Sb",0}, {"Te",2}, { "I",1}, {"Xe",0}, {"Cs",1}, {"Ba",2}, {"La",0}, {"Ce",0}, {"Pr",0}, {"Nd",0},
    {"Pm",0}, {"Sm",0}, {"Eu",0}, {"Gd",0}, {"Tb",0}, {"Dy",0}, {"Ho",0}, {"Er",0}, {"Tm",0}, {"Yb",0},
    {"Lu",0}, {"Hf",4}, {"Ta",0}, { "W",0}, {"Re",0}, {"Os",0}, {"Ir",0}, {"Pt",0}, {"Au",0}, {"Hg",0},
    {"Tl",0}, {"Pb",0}, {"Bi",0}, {"Po",0}, {"At",1}, {"Rn",0}, {"Fr",0}, {"Ra",2}, {"Ac",0}, {"Th",0},
    {"Pa",0}, { "U",0}, {"Np",0}, {"Pu",0}, {"Am",0}, {"Cm",0}, {"Bk",0}, {"Cf",0}, {"Es",0}, {"Fm",0},
    {"Md",0}, {"No",0}, {"Lr",0}
};

TAtomAttr::TAtomAttr()
    : atno(0), eS(0), eD(0), eT(0), charge(0), mult(0), nH(0), penDif(-1),
    label(0), x(0), y(0), z(0)
{
}


TBondAttr::TBondAttr()
    : order(0), cyc(0), stereo(0), dir(0)
{
}


void TSBond::Invert()
{
    Attr.dir = ~Attr.dir;
}


TStructure::TStructure()
{
}

TStructure::TStructure(const TStructure& str)
    : TAGraph(str), name(str.name)
{
}

int TStructure::GetAtNo(const string& aname)
{
    for (int i = 0; i < sizeof(atomDataNormal) / sizeof(*atomDataNormal); i++)
        if (aname == atomDataNormal[i].name)
            return i + 1;
    return 0;
}

string TStructure::GetAtName(int atno)
{
    if (atno <= 0 || atno > sizeof(atomDataNormal) / sizeof(*atomDataNormal))
        return "Xx";
    return atomDataNormal[atno - 1].name;
}

int TStructure::GetNormVal(int atno)
{
    if (atno <= 0 || atno > sizeof(atomDataNormal) / sizeof(*atomDataNormal))
        return 0;
    return atomDataNormal[atno - 1].val;
}


