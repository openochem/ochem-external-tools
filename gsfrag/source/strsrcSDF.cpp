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

#include "strsrcSDF.h"

#include <string>
#include <sstream>
#include <fstream>
#include "global.h"

#include <iostream>


string StrSourceSDF::getline()
{
    std::getline(file, buf);
    buf = fullstrip(buf, '\r');
    curline++;
    //cout << curline << ": " << buf.length() << "[" << buf << "]\n";
    //  buf = fullstrip(buf);
    return buf;
}


struct SDFAtom {
    float x, y, z;
    string type;
    int mdiff;
    int charge;
    int stereo;
    int atno;
    int nH, eS, eD, eT, nArom;
    int mapNo;

    SDFAtom()
        : x(0), y(0), z(0), mdiff(0), charge(0), stereo(0),
        atno(0), nH(0), eS(0), eD(0), eT(0), nArom(0), mapNo(-1)
    {}
};

struct SDFBond {
    int aorg, atrg;
    int type;
    int stereo;

    SDFBond()
        : aorg(-1), atrg(-1), type(1), stereo(0)
    {}
};

void StrSourceSDF::LoadStructures()
{
    int state = 0;
    TStructureP str;
    string name;
    int natoms, nbonds, nstrat;
    vector<SDFAtom> atoms;

    while (!file.fail() && file.peek() != EOF)
    {
        try
        {
            if (state == 0) // structure header
            {
                name = fullstrip(getline());
                getline();
                getline();
                natoms = 0;
                nbonds = 0;
                state = 1;
                continue;
            }
            else if (state == 1) // count line
            {
                getline();
                natoms = stov<int>(safesstr(buf, 0, 3));
                nbonds = stov<int>(safesstr(buf, 3, 3));
                if (natoms <= 0)
                    throw bad_conv("Invalid atom count");
                if (nbonds < 0)
                    throw bad_conv("Invalid bond count");
                state = 2;
                continue;
            }
            else if (state == 2) // atom block
            {
                str = new TStructure();
                str->name = name;
                nstrat = 0;
                atoms.clear();
                atoms.resize(natoms, SDFAtom());
                for (int i = 0; i < natoms; i++)
                {
                    SDFAtom& atom = atoms[i];
                    getline();
                    atom.x = stov<float>(safesstr(buf, 0, 10));
                    atom.y = -stov<float>(safesstr(buf, 10, 10));
                    atom.z = stov<float>(safesstr(buf, 20, 10));
                    atom.type = fullstrip(safesstr(buf, 30, 4));
                    atom.mdiff = stov<int>(safesstr(buf, 34, 2));
                    atom.charge = stov<int>(safesstr(buf, 36, 3));
                    atom.stereo = stov<int>(safesstr(buf, 39, 3));
                    atom.atno = TStructure::GetAtNo(atom.type);
                    if (atom.atno <= 0)
                    {
                        throw bad_conv("Unknown atom '" + atom.type + '\'');
                    }
                    else if (atom.atno > atnoH)
                    {
                        TSAtom* satom = new TSAtom;
                        str->AddVertex(satom);
                        satom->Attr.atno = atom.atno;
                        atom.mapNo = nstrat++;
                    }
                }
                state = 3;
                continue;
            }
            else if (state == 3) // bond block
            {
                for (int i = 0; i < nbonds; i++)
                {
                    SDFBond bond;
                    getline();
                    bond.aorg = stov<int>(safesstr(buf, 0, 3));
                    bond.atrg = stov<int>(safesstr(buf, 3, 3));
                    bond.type = stov<int>(safesstr(buf, 6, 3));
                    bond.stereo = stov<int>(safesstr(buf, 9, 3));
                    bond.aorg--;
                    bond.atrg--;
                    if (bond.aorg < 0 || bond.atrg < 0 || atoms[bond.aorg].atno == 0 || atoms[bond.atrg].atno == 0)
                        throw bad_conv("Invalid bond");
                    if (atoms[bond.aorg].atno == atnoH)
                    {
                        atoms[bond.atrg].nH++;
                        atoms[bond.atrg].eS++;
                        bond.aorg = -1;
                        bond.atrg = -1;
                    }
                    else if (atoms[bond.atrg].atno == atnoH)
                    {
                        atoms[bond.aorg].nH++;
                        atoms[bond.aorg].eS++;
                        bond.aorg = -1;
                        bond.atrg = -1;
                    }
                    if (bond.aorg >= 0 && bond.atrg >= 0)
                    {
                        TSBond* sbond = new TSBond;
                        str->MakeEdge(atoms[bond.aorg].mapNo, atoms[bond.atrg].mapNo, sbond);
                        sbond->Attr.order = bond.type;
                        switch (bond.type)
                        {
                        case 1:
                            atoms[bond.aorg].eS++;
                            atoms[bond.atrg].eS++;
                            break;
                        case 2:
                            atoms[bond.aorg].eD++;
                            atoms[bond.atrg].eD++;
                            break;
                        case 3:
                            atoms[bond.aorg].eT++;
                            atoms[bond.atrg].eT++;
                            break;
                        case 4:
                            atoms[bond.aorg].nArom++;
                            atoms[bond.atrg].nArom++;
                            break;
                        }
                    }
                }

                for (int i = 0; i < natoms; i++)
                {
                    SDFAtom& matom = atoms[i];
                    if (matom.mapNo < 0)
                        continue;
                    TSAtom& satom = str->Atom(matom.mapNo);
                    satom.Attr.nH = matom.nH;
                    satom.Attr.eS = matom.eS;
                    satom.Attr.eD = matom.eD;
                    if (matom.nArom)
                    {
                        satom.Attr.eD++;
                        satom.Attr.eS += matom.nArom - 1;
                    }
                    satom.Attr.eT = matom.eT;
                    int val = satom.Attr.eS + 2 * satom.Attr.eD + 3 * satom.Attr.eT;
                    int normval = TStructure::GetNormVal(matom.atno);
                    if (normval > val)
                    {
                        val = normval - val;
                        satom.Attr.eS += val;
                        satom.Attr.nH += val;
                    }
                    satom.Attr.x = matom.x;
                    satom.Attr.y = matom.y;
                    satom.Attr.z = matom.z;
                }

                if (!str->name.empty())
                    str->name += ' ';
                str->name += '(' + setname + '#' + vtos(strs.size() + 1) + ')';
                strs.push_back(str);
                cout << "MESSAGE: Molecule " << strs.size() << " is loaded\n";
                cout.flush();
                state = 4;
                continue;
            }
            else if (state == 4) // data header
            {
                getline();
                if (buf == "$$$$")
                {
                    state = 0;
                }
                continue;
            }
        }
        catch (const bad_conv& bc)
        {
            str = new TStructure();
            str->name += "SDF error [" + bc.msg + "] line #" + curline + " str #" + vtos(strs.size() + 1) + " near '" + buf + "'";
            strs.push_back(str);
            state = 4;
        }
    }
}


TStructureP StrSourceSDF::GetStructure(unsigned i)
{
    return strs[i];
}


unsigned StrSourceSDF::GetStructCount()
{
    return strs.size();
}

StrSourceSDF::StrSourceSDF(string sdfile)
    : file(sdfile.c_str()), curline(0)
{
    try {
        LoadStructures();
    }
    catch (const bad_conv& bc)
    {
        strs.clear();
        throw bc;
    }
}
