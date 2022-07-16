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

#ifndef __STRUCT_H
#define __STRUCT_H

#include <string>
using namespace std;

#include "graph.h"
typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef char int8;
typedef short int16;
typedef int int32;

enum atnums
{
    atnoH = 1, atnoHe,
    atnoLi, atnoBe, atnoB, atnoC, atnoN, atnoO, atnoF, atnoNe,
    atnoNa, atnoMg, atnoAl, atnoSi, atnoP, atnoS, atnoCl, atnoAr,
    atnoK, atnoCa, atnoSc, atnoTi, atnoV, atnoCr, atnoMn, atnoFe, atnoCo,
    atnoNi, atnoCu, atnoZn, atnoGa, atnoGe, atnoAs, atnoSe, atnoBr, atnoKr,
    atnoRb, atnoSr, atnoY, atnoZr, atnoNb, atnoMo, atnoTc, atnoRu, atnoRh,
    atnoPd, atnoAg, atnoCd, atnoIn, atnoSn, atnoSb, atnoTe, atnoI, atnoXe,
    atnoCs, atnoBa, atnoLa, atnoCe, atnoPr, atnoNd, atnoPm, atnoSm, atnoEu,
    atnoGd, atnoTb, atnoDy, atnoHo, atnoEr, atnoTm, atnoYb, atnoLu, atnoHf,
    atnoTa, atnoW, atnoRe, atnoOs, atnoIr, atnoPt, atnoAu, atnoHg, atnoTl,
    atnoPb, atnoBi, atnoPo, atnoAt, atnoRn,
    atnoFr, atnoRa, atnoAc, atnoTh, atnoPa, atnoU, atnoNp, atnoPu, atnoAm,
    atnoCm, atnoBk, atnoCf, atnoEs, atnoFm, atnoMd, atnoNo, atnoLr
};


struct TAtomAttr
{
    uint8 atno;
    unsigned eS : 3; //no. of single bonds
    unsigned eD : 3; //no. of double bonds
    unsigned eT : 2; //no. of triple bonds
    int8 charge;
    int8 mult;
    uint8 nH;
    uint8 label;
    int penDif;
    float x, y, z;

    TAtomAttr();

    enum {
        f_atno = 1,
        f_envir = 2,
        f_charge = 4,
        f_mult = 8,
        f_nH = 16,
        f_label = 32,
        f_pen = 64,
        f_coor = 128,
        f_all = f_atno | f_envir | f_charge | f_mult | f_nH | f_label | f_pen | f_coor
    };
};

class TSAtom : public TVertex {
protected:
public:
    TAtomAttr Attr;

    TSAtom()
    {}
    TSAtom(const TSAtom& atom)
        : Attr(atom.Attr) {}
    TSAtom(const TAtomAttr& attr)
        : Attr(attr) {}
    virtual TVertex* Clone() const
    {
        return new TSAtom(*this);
    }

    static TSAtom* FromVertex(TVertex* vert)
    {
        return dynamic_cast<TSAtom*>(vert);
    }
    static const TSAtom* FromVertex(const TVertex* vert)
    {
        return dynamic_cast<const TSAtom*>(vert);
    }

};

struct TBondAttr
{
    int8 order;
    unsigned int cyc : 1;
    unsigned int stereo : 2;   // 0-none, 1-wedge, 2-dash
    unsigned int dir : 1;      // 0-forward, 1-backward

    TBondAttr();

    enum {
        f_order = 1,
        f_cyc = 2,
        f_stereo = 4,
        f_all = f_order | f_cyc | f_stereo
    };
};

class TSBond : public TEdge
{
public:
    TBondAttr Attr;

    TSBond()
    {}
    TSBond(const TSBond& bond)
        : Attr(bond.Attr) {}
    TSBond(const TBondAttr& attr)
        : Attr(attr) {}
    virtual TEdge* Clone() const
    {
        return new TSBond(*this);
    }
    virtual void Invert();

    static TSBond* FromEdge(TEdge* edge)
    {
        return dynamic_cast<TSBond*>(edge);
    }
    static const TSBond* FromEdge(const TEdge* edge)
    {
        return dynamic_cast<const TSBond*>(edge);
    }

};

class TStructure : public TAGraph
{
public:
    TStructure();
    TStructure(const TStructure& str);

    string name;

    virtual void Reset()
    {
        DoReset();
    }
    void DoReset()
    {
        TAGraph::DoReset();
    }

    TSAtom& Atom(int i) const
    {
        return dynamic_cast<TSAtom&>(Vertex(i));
    }

    static int GetAtNo(const string& aname);
    static string GetAtName(int atno);
    static int GetNormVal(int atno);


    virtual TGraph* Clone() const
    {
        return new TStructure(*this);
    }
};
typedef rc_ptr<TStructure> TStructureP;
typedef vector<TStructureP> TStructures;

#endif
