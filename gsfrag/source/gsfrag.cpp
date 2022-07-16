/* ************************************************************************************
   Copyright (C) 2002-2022 by K.S.Fedyaev, E.V.Radchenko, M.I.Skvortsova, V.A.Palyulin

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
************************************************************************************ */

#include "strsrcSDF.h"


#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string.h>

#include "global.h"

#define BLOCKNAME "GSFRAG"


//#include <stdlib.h>
//#include <math.h>
//#include <io.h>
//#include <dos.h>

#define MaxStructSize 100  //Maximal number of atoms in the structures
#define MaxFragSize 16  //Maximal number of atoms in fragments
#define MaxElem 250   //Maximal number of fragments in the database

typedef int FragMatr[MaxFragSize][MaxFragSize];
typedef int FragArr[MaxFragSize];
typedef int StructMatr[MaxStructSize][MaxStructSize];
typedef int StructArr[MaxStructSize];
typedef int BasArr[MaxElem];
typedef float Vector[MaxElem];

struct fragm
{
    int n;
    FragMatr a;
    char nm[16];
    int autom;
};

struct Comp
{
    int in;
    char nm[15];
    char mf[10];
    float y;
    short got;
    int n;
    StructMatr a;
};

typedef fragm FragBase[MaxElem];
typedef char DescName[MaxElem][15];

#include "gsfrag-def.h"

//#include "befunc.cpp"

void substrnum(FragMatr* pa, int na, StructMatr* pb, int nb,
    StructArr* vdegB, FragArr* vdegS,
    StructArr* pused, StructArr* pc, int stp, int* pres)
{
    int i, j;

    for (i = 0; i < nb; i++)
    {
        if ((!pused[0][i]) &&
            (pa[0][stp][stp] == pb[0][i][i]) &&
            (vdegS[0][stp] <= vdegB[0][i]))
        {
            pc[0][stp] = i;
            if (stp)
                for (j = 0; j < stp; j++)
                    if ((pa[0][stp][j] > 0) &&
                        (pb[0][i][pc[0][j]] == 0))
                        goto nextlab;
            if (stp == na - 1)
            {
                pres[0]++;
                goto nextlab;
            }
            pused[0][i] = 1;
            substrnum(pa, na, pb, nb, vdegB, vdegS, pused, pc, stp + 1, pres);
            pused[0][i] = 0;
        }
    nextlab: continue;
    }
}

int SubstCountNum(FragMatr* pa, int na, StructMatr* pb, int nb, int autom)
// This subroutine returnes the number of entrances of fragment "a"
// in matrix "b". "na" and "nb" - dimensions of these matrices.
{
    int i, j;

    StructArr used, vdegr2, c;
    FragArr vdegr1;
    StructArr* vdegB = &vdegr2;
    FragArr* vdegS = &vdegr1;
    StructArr* pused = &used;
    StructArr* pc = &c;

    int res1, * pres;

    for (i = 0; i < na; i++)
    {
        vdegr1[i] = 0;
        for (j = 0; j < na; j++)
            if (i != j)
                vdegr1[i] += pa[0][i][j];
    }
    for (i = 0; i < nb; i++) {
        vdegr2[i] = 0;
        for (j = 0; j < nb; j++)
            if (i != j)
                vdegr2[i] += pb[0][i][j];
    }
    for (i = 0; i < MaxStructSize; i++) used[i] = 0;
    res1 = 0; pres = &res1;
    substrnum(pa, na, pb, nb, vdegB, vdegS, pused, pc, 0, pres);
    return res1 / autom;
}

int GetAdjMatr(TStructure& fs, StructMatr* a)
// Getting Adjacency Matrix of a current structure from *.str-file.
// Returned values: size of matrix or
//                  -number of structure if it is too large or
{

    for (int i = 0; i < MaxStructSize; i++)
        for (int j = 0; j < MaxStructSize; j++)
            a[0][i][j] = 0;

    int n = fs.MaxIndex();  // size of a current structure
    if (n > MaxStructSize)
        return -n;

    for (int ivf = 0; ivf < n; ivf++)
    {
        for (TEdgeIter eiter = fs.Vertex(ivf).first(); eiter != fs.Vertex(ivf).last(); eiter++)
        {
            int ivt = (*eiter)->Term().Index();
            int order = TSBond::FromEdge(fs.Atom(ivf).FindEdge(ivt)->InfoEdge())->Attr.order;
            a[0][ivf][ivt] = order;
        }
    }
    return n;
}

int CalcDescComp(Comp* cmp, vector<float>& res)
// This subprogram calculates fragment descriptors for one compound
//  and printes their values to the output file
// cmp - pointer to the structure of a current compound
// descnum - number of the current (scanning) fragment.

{
    int descnum;
    StructMatr* pb;
    FragMatr* pa;
    int num;

    res.assign(numfrag, 0);

    pb = &cmp[0].a;
    descnum = 0;
    while (descnum < numfrag)
    {
        pa = &frag[descnum].a;
        num = SubstCountNum(pa, frag[descnum].n, pb, cmp[0].n, frag[descnum].autom);
        res[descnum] = num;
        descnum++;
    }
    return descnum;
}

int main(int argc, char* argv[])
{

    if (argc < 2)
        return 0;

    StrSourceSDF strsrc(argv[1]);
    int strnum = strsrc.GetStructCount();
    vector< vector<float> > descr(strnum);


    //	int i,nb;
    Comp cmp, * pcmp = &cmp;
    //	int currx,curry;


    StructMatr* pB;

    for (int s = 0; s < strnum; s++)
    {
        TStructureP str = strsrc.GetStructure(s);
        if (str->MaxIndex() > 0)
        {
            pB = &cmp.a;
            int nb = GetAdjMatr(*str, pB);
            if (nb > 0)
            {
                cmp.n = nb;
                CalcDescComp(&cmp, descr[s]);
            }
            else
            {
                if (nb < 0)
                    str->name = "Structure too large";
            }
        }
        cout << "MESSAGE: Molecule " << (s + 1) << " out of " << strnum << " is calculated\n";
        cout.flush();
    }

    ofstream out("descriptors.txt");
    out << "MOLECULE";
    for (int i = 0; i < numfrag; i++)
        out << '\t' << fullstrip(frag[i].nm);
    out << '\n';

    for (int s = 0; s < strnum; s++)
    {
        out << "MOL" << (s + 1);
        int ndescr = descr[s].size();
        TStructureP str = strsrc.GetStructure(s);
        if (ndescr > 0 && ndescr != numfrag)
        {
            ndescr = 0;
            str->name = "Descriptor number mismatch";
        }
        if (ndescr > 0)
        {
            for (int i = 0; i < ndescr; i++)
            {
                out << '\t' << descr[s][i];
            }
        }
        else
        {
            out << " ERROR: " << BLOCKNAME << ' ' << str->name;
        }
        out << '\n';
    }

    /*	for (i=0;i<numfrag;i++)
        {
            fprintf(fcon,"%3d TOP %s",i+1,frag[i].nm);
            fprintf(fhlp," %-d\t%-d\t%-s\tnumber of fragments ",
                i+1,i+1,frag[i].nm);
            for (int j=0;j<strlen(frag[i].nm);j++)
            {
                if (frag[i].nm[j]=='p')
                {
                    if (j>1) fprintf(fhlp,"_");
                    fprintf(fhlp,"Path");
                }
                if (frag[i].nm[j]=='c')
                {
                    if (j>1) fprintf(fhlp,"_");
                    fprintf(fhlp,"Cyc");
                }
                if ((frag[i].nm[j]>='0')&&(frag[i].nm[j]<='9'))
                    fprintf(fhlp,"%c",frag[i].nm[j]);
                if ((frag[i].nm[j]>='A')&&(frag[i].nm[j]<='Z'))
                {
                    if (!((frag[i].nm[j-1]>='A')&&(frag[i].nm[j-1]<='Z')))
                        fprintf(fhlp,"[");
                    fprintf(fhlp,"%c",frag[i].nm[j]);
                    if (!((frag[i].nm[j+1]>='A')&&(frag[i].nm[j+1]<='Z')))
                        fprintf(fhlp,"]");
                }
            }//for
            fprintf(fhlp,"\ttopological descriptor \n");

            if (i<numfrag-1) fprintf(fcon,"\n");
        }
    */


    return 0;
}




