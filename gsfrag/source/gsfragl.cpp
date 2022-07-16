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

#define BLOCKNAME "GSFRAGL"


//#include <stdio.h>
//#include <string.h>
//#include <stdlib.h>
//#include <math.h>
//#include <io.h>
//#include <dos.h>

#define MaxStructSize 100  //Maximal number of atoms in the structures
#define MaxFragSize 8  //Maximal number of atoms in fragments
#define MaxElem 900   //Maximal number of fragments in the database

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

    int numpos;
    FragArr pos;
    int autom;
    char l[2 * MaxFragSize];
};

struct fragml
{
    int n;
    FragMatr a;
    char nm[16];

    int numpos;
    FragArr pos;
    FragArr autom;
    char l[2 * MaxFragSize];
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

    char l[2 * MaxStructSize];
};

typedef fragm FragBase[MaxElem];
typedef fragml FragBaseL[MaxElem];
typedef char DescName[MaxElem][15];

#include "gsfragl-def.h"

FragBase frag;
int numfrag;

void substrnum(FragMatr* pa, int na, StructMatr* pb, int nb,
    StructArr* vdegB, FragArr* vdegS, char* pla, char* plb,
    StructArr* pused, StructArr* pc, int stp, int* pres)
{
    int i, j;

    for (i = 0; i < nb; i++)
    {
        if ((!pused[0][i]) &&
            (pa[0][stp][stp] == pb[0][i][i]) &&
            (vdegS[0][stp] <= vdegB[0][i]) &&
            (
                (!strncmp(pla + 2 * stp, plb + 2 * i, 2)) ||
                (!strncmp(pla + 2 * stp, "* ", 2))
                ))
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
            substrnum(pa, na, pb, nb, vdegB, vdegS, pla, plb, pused, pc, stp + 1, pres);
            pused[0][i] = 0;
        }
    nextlab: continue;
    }
}

int SubstCountNum(FragMatr* pa, int na, char* pla,
    StructMatr* pb, int nb, char* plb, int autom)
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
    substrnum(pa, na, pb, nb, vdegB, vdegS, pla, plb, pused, pc, 0, pres);
    return res1 / autom;
}

int degree(FragMatr* pa, int na, int k)
{
    int i, deg;

    deg = 0;
    for (i = 0; i < na; i++)
    {
        if (i != k) deg += pa[0][k][i];
    }
    return deg;
}

int GetAdjMatr(TStructure& fs, StructMatr* a, char* l)
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
        string aname = TStructure::GetAtName(fs.Atom(ivf).Attr.atno);
        l[2 * ivf] = aname[0];
        l[2 * ivf + 1] = aname.length() > 1 ? aname[1] : ' ';
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
        num = SubstCountNum(pa, frag[descnum].n, frag[descnum].l, pb,
            cmp[0].n, cmp[0].l, frag[descnum].autom);
        res[descnum] = num;
        descnum++;
    }
    return descnum;
}

int AddLabeledFrags(int numfrag, const char* mark, int val)
{
    FragMatr* pa;

    // numfrag - a number of labeled fragments in the database
    char ssn[10];
    for (int nt = 0; nt < numunlabfrag; nt++)
    {
        //one labeled vertex
        for (int k = 0; k < unlabfrag[nt].numpos; k++)
        {
            pa = &unlabfrag[nt].a;
            if (!(degree(pa, unlabfrag[nt].n, unlabfrag[nt].pos[k]) > val))
            {
                for (int i = 0; i < unlabfrag[nt].n; i++)
                {
                    for (int j = 0; j < unlabfrag[nt].n; j++)
                        frag[numfrag].a[i][j] = unlabfrag[nt].a[i][j];
                }
                strcpy(frag[numfrag].nm, unlabfrag[nt].nm);
                strcat(frag[numfrag].nm, "-");
                sprintf(ssn, "%d", unlabfrag[nt].pos[k] + 1);
                strcat(frag[numfrag].nm, ssn);
                strncat(frag[numfrag].nm, mark, 2);
                //				if (mark[1]=='L') strcat(frag[numfrag].nm,"l"); else
                //					if (mark[1]=='R') strcat(frag[numfrag].nm,"r"); else
                //						strcat(frag[numfrag].nm,mark+1);
                frag[numfrag].n = unlabfrag[nt].n;
                strcpy(frag[numfrag].l, unlabfrag[nt].l);
                frag[numfrag].l[2 * unlabfrag[nt].pos[k]] = mark[0];
                frag[numfrag].l[2 * unlabfrag[nt].pos[k] + 1] = mark[1];
                frag[numfrag].autom = unlabfrag[nt].autom[k];
                /*
                                fprintf(fhlp," %-d\t%-d\t%-s\tnumber of fragments ",
                                    numfrag+1,numfrag+1,frag[numfrag].nm);
                                for (i=0;i<strlen(unlabfrag[nt].nm);i++)
                                {
                                    if (unlabfrag[nt].nm[i]=='p')
                                    {
                                        if (i>1) fprintf(fhlp,"_");
                                        fprintf(fhlp,"Path");
                                    }
                                    if (unlabfrag[nt].nm[i]=='c')
                                    {
                                        if (i>1) fprintf(fhlp,"_");
                                        fprintf(fhlp,"Cyc");
                                    }
                                    if ((unlabfrag[nt].nm[i]>='0')&&(unlabfrag[nt].nm[i]<='9'))
                                        fprintf(fhlp,"%c",unlabfrag[nt].nm[i]);
                                    if ((unlabfrag[nt].nm[i]>='A')&&(unlabfrag[nt].nm[i]<='Z'))
                                    {
                                        if (!((unlabfrag[nt].nm[i-1]>='A')&&(unlabfrag[nt].nm[i-1]<='Z')))
                                            fprintf(fhlp,"[");
                                        fprintf(fhlp,"%c",unlabfrag[nt].nm[i]);
                                        if (!((unlabfrag[nt].nm[i+1]>='A')&&(unlabfrag[nt].nm[i+1]<='Z')))
                                            fprintf(fhlp,"]");
                                    }
                                }//for
                                fprintf(fhlp," with label ");
                                fprintf(fhlp,"%s",frag[numfrag].nm+strlen(frag[numfrag].nm)-2);
                                if (!(frag[numfrag].nm[strlen(frag[numfrag].nm)-1]==' '))
                                    fprintf(fhlp," ");
                                fprintf(fhlp,"on atom %d\t",unlabfrag[nt].pos[k]+1);
                                fprintf(fhlp,"topological descriptor \n");
                */
                numfrag++;
            }
        }
        /*  //two labeled vertices
        if (pFr[0][nt].numpos>1)
        for (k=0;k<pFr[0][nt].numpos;k++)
        for (int l=k+1;l<pFr[0][nt].numpos;l++)
        {
        for (i=0;i<pFr[0][nt].n;i++)
        for (j=0;j<pFr[0][nt].n;j++)
        pFr[0][numfrag].a[i][j]=pFr[0][nt].a[i][j];
        pFr[0][numfrag].n=pFr[0][nt].n;
        strcpy(pFr[0][numfrag].l,pFr[0][nt].l);
        pFr[0][numfrag].l[2*pFr[0][nt].pos[k]]=mark;
        pFr[0][numfrag].l[2*pFr[0][nt].pos[l]]=mark;
        itoa(nt+1,basicfragm[numfrag+1],10);
        strcat(basicfragm[numfrag+1],"-");
        itoa(k+1,st,10);
        strcat(basicfragm[numfrag+1],st);
        strcat(basicfragm[numfrag+1],",");
        itoa(l+1,st,10);
        strcat(basicfragm[numfrag+1],st);
        numfrag++;
        }
        // 3 labeled vertices
        if (pFr[0][nt].numpos>2)
        for (k=0;k<pFr[0][nt].numpos;k++)
        for (int l=k+1;l<pFr[0][nt].numpos;l++)
        for (int m=l+1;m<pFr[0][nt].numpos;m++)
        {
        for (i=0;i<pFr[0][nt].n;i++)
        for (j=0;j<pFr[0][nt].n;j++)
        pFr[0][numfrag].a[i][j]=pFr[0][nt].a[i][j];
        pFr[0][numfrag].n=pFr[0][nt].n;
        strcpy(pFr[0][numfrag].l,pFr[0][nt].l);
        pFr[0][numfrag].l[2*pFr[0][nt].pos[k]]=mark;
        pFr[0][numfrag].l[2*pFr[0][nt].pos[l]]=mark;
        pFr[0][numfrag].l[2*pFr[0][nt].pos[m]]=mark;
        itoa(nt+1,basicfragm[numfrag+1],10);
        strcat(basicfragm[numfrag+1],"-");
        itoa(k+1,st,10);
        strcat(basicfragm[numfrag+1],st);
        strcat(basicfragm[numfrag+1],",");
        itoa(l+1,st,10);
        strcat(basicfragm[numfrag+1],st);
        strcat(basicfragm[numfrag+1],",");
        itoa(m+1,st,10);
        strcat(basicfragm[numfrag+1],st);
        numfrag++;
        }  */
    }
    return numfrag;
}



int main(int argc, char* argv[])
{
    if (argc < 2)
        return 0;

    StrSourceSDF strsrc(argv[1]);
    int strnum = strsrc.GetStructCount();
    vector< vector<float> > descr(strnum);

    for (int i = 0; i < numunlabfrag; i++)
    {
        unlabfrag[i].l[0] = '\0';
        for (int nb = 0; nb < unlabfrag[i].n; nb++)
            strcat(unlabfrag[i].l, "* ");
    }
    numfrag = 0;
    numfrag = AddLabeledFrags(numfrag, "C ", 3);
    numfrag = AddLabeledFrags(numfrag, "S ", 3);
    numfrag = AddLabeledFrags(numfrag, "N ", 3);
    numfrag = AddLabeledFrags(numfrag, "O ", 2);
    numfrag = AddLabeledFrags(numfrag, "Cl", 3);
    numfrag = AddLabeledFrags(numfrag, "Br", 3);
    numfrag = AddLabeledFrags(numfrag, "I ", 3);
    numfrag = AddLabeledFrags(numfrag, "F ", 1);


    Comp cmp;


    StructMatr* pB;
    char* pl;

    for (int s = 0; s < strnum; s++)
    {
        TStructureP str = strsrc.GetStructure(s);
        if (str->MaxIndex() > 0)
        {
            pB = &cmp.a;
            pl = cmp.l;
            int nb = GetAdjMatr(*str, pB, pl);
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




