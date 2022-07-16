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

#include <algorithm>
#include <functional>

#include "graph.h"

int uniq_pack(vector<int>& vect)
{
    int curnum = 0;
    int size = vect.size();
    vector<int> mapval(size, -1);
    for (int i = 0; i < size; i++)
    {
        int j = vect[i];
        if (j >= 0)
        {
            if (mapval[j] < 0)
                mapval[j] = curnum++;
            vect[i] = mapval[j];
        }
    }
    return curnum;
}


TEdge::TEdge()
    : vterm(0)
{
}

TEdge::TEdge(TVertex* termVert)
    : vterm(termVert)
{
}

TEdgeMirror::TEdgeMirror(TVertex* termVert, TEdge* infoEdge)
    : TEdge(termVert), inf(infoEdge)
{
}

TVertex::TVertex()
    : index(-1)
{
    edges.reserve(2);
}

TVertex::TVertex(const TVertex& vert)
    : index(-1)
{
    edges.reserve(vert.edges.size());
}

TVertex::~TVertex()
{
}

class __cmp_term_vert {
public:
    TVertex* vert;
    __cmp_term_vert(TVertex* v)
        : vert(v) {}
    bool operator()(const TEdgeP& edge)
    {
        return &edge->Term() == vert;
    }
};

TEdgeP TVertex::FindEdge(TVertex* vert) const
{
    TEdgeCIter pos = find_if(first(), last(), __cmp_term_vert(vert));
    if (pos != last())
        return *pos;
    else
        return 0;
}

class __cmp_term_vert_index {
public:
    int vert;
    __cmp_term_vert_index(int v)
        : vert(v) {}
    bool operator()(const TEdgeP& edge)
    {
        return edge->Term().Index() == vert;
    }
};

TEdgeP TVertex::FindEdge(int vert) const
{
    TEdgeCIter pos = find_if(first(), last(),
        __cmp_term_vert_index(vert));
    if (pos != last())
        return *pos;
    else
        return 0;
}

static bool __is_info(const TEdgeP& edge)
{
    return edge->IsInfoEdge();
}

int TVertex::CountInfoEdges() const
{
    return count_if(first(), last(), __is_info);
}

void TVertex::Attach(const TEdgeP& edge, TVertex* vert)
{
    TEdgeP edgeThis, edgeVert;
    if (index > vert->index) {
        edge->SetTerm(vert);
        edgeThis = edge;
        edgeVert = new TEdgeMirror(this, edge);
    }
    else {
        edge->Invert();
        edge->SetTerm(this);
        edgeVert = edge;
        edgeThis = new TEdgeMirror(vert, edge);
    }
    edges.push_back(edgeThis);
    vert->edges.push_back(edgeVert);
}

void TVertex::Detach(TVertex* vert)
{
    edges.erase(first() +
        (find_if(first(), last(), __cmp_term_vert(vert)) - first()));
    vert->edges.erase(vert->first() + (find_if(vert->first(),
        vert->last(), __cmp_term_vert(this)) - vert->first()));
}


TGraph::TGraph()
{
}

vector<int> TGraph::Components() const
{
    vector<int> curcomp(VertMask());
    int size = curcomp.size();
    for (int i = 0; i < size; i++)
        if (curcomp[i] >= 0)
            curcomp[i] = i;

    for (TVertCIter viter = first(); viter != last(); viter++)
    {
        TVertex& curvert = **viter;
        int tv = curvert.Index();
        if (!ValidVert(tv))
            continue;
        for (TEdgeIter eiter = curvert.first(); eiter != curvert.last(); eiter++)
        {
            int ta = (*eiter)->Term().Index();
            int iv = curcomp[tv];
            int ia = curcomp[ta];
            if (ia == iv || !ValidEdge(tv, ta))
                continue;
            if (ia < iv)
                swap(ia, iv);
            replace(curcomp.begin() + ia, curcomp.end(), ia, iv);
        }
    }

    return curcomp;
}

int TGraph::NumCompon() const
{
    vector<int> curcomp(Components());
    return uniq_pack(curcomp);
}

int TGraph::IsConnected() const
{
    return NumCompon() == 1;
}

vector<int> TGraph::Components(vector<int>& comp, vector<int>& touch) const
{
    int size = comp.size();
    vector<int> curcomp(comp), mask(size, -1);
    for (int i = 0; i < size; i++)
    {
        int oldcomp = comp[i];
        if (oldcomp >= 0 && touch[i] >= 0 && mask[i] < 0)
        {
            for (int j = 0; j < size; j++)
            {
                if (oldcomp == comp[j])
                {
                    mask[j] = 0;
                    curcomp[j] = j;
                }
            }
        }
    }

    for (TVertCIter viter = first(); viter != last(); viter++)
    {
        TVertex& curvert = **viter;
        int tv = curvert.Index();
        if (mask[tv] < 0 || !ValidVert(tv))
            continue;
        for (TEdgeIter eiter = curvert.first(); eiter != curvert.last(); eiter++)
        {
            int ta = (*eiter)->Term().Index();
            int iv = curcomp[tv];
            int ia = curcomp[ta];
            if (mask[ta] < 0 || ia == iv || !ValidEdge(tv, ta))
                continue;
            if (ia < iv)
                swap(ia, iv);
            replace(curcomp.begin() + ia, curcomp.end(), ia, iv);
        }
    }

    return curcomp;
}


vector<int> TGraph::VertMask() const
{
    vector<int> tmp(MaxIndex(), -1);
    for (TVertCIter iter = first(); iter != last(); iter++)
    {
        int index = (*iter)->Index();
        if (ValidVert(index))
            tmp[index] = index;
    }
    return tmp;
}

class __test_term_vert_valid {
public:
    const TGraph& graph;
    int vf;
    __test_term_vert_valid(const TGraph& g, int vfrom)
        : graph(g), vf(vfrom) {}
    bool operator()(const TEdgeP& edge)
    {
        return graph.ValidEdge(vf, edge->Term().Index());
    }
};

int TGraph::Valence(const TVertex& vert) const
{
    return count_if(vert.first(), vert.last(),
        __test_term_vert_valid(*this, vert.Index()));
}

int TGraph::NumEdges() const
{
    int nEdges = 0;
    for (TVertCIter viter = first(); viter != last(); viter++)
    {
        TVertex& curvert = **viter;
        int tv = curvert.Index();
        if (!ValidVert(tv))
            continue;
        for (TEdgeIter eiter = curvert.first(); eiter != curvert.last(); eiter++)
        {
            int ta = (*eiter)->Term().Index();
            if (!ValidEdge(tv, ta))
                continue;
            if (tv > ta)
                nEdges++;
        }
    }
    return nEdges;
}

TAGraph::TAGraph()
{
    verts.reserve(4);
}

TAGraph::TAGraph(const TAGraph& grp)
{
    verts.reserve(grp.NumVerts());
    for (int i = 0; i < grp.NumVerts(); i++)
    {
        TVertex* curvert = &grp.Vertex(i);
        AddVertex(curvert->Clone());
        for (TEdgeIter iter = curvert->first(); iter != curvert->last(); iter++) {
            TEdge* curedge = *iter;
            if (curedge->IsInfoEdge())
                MakeEdge(i, curedge->Term().Index(), curedge->Clone());
        }
    }
}


TAGraph::~TAGraph()
{
}

vector<int> TAGraph::VertMask() const
{
    return vector<int>(MaxIndex(), 0);
}

bool TAGraph::ValidVert(int vert) const
{
    return vert >= 0 && vert < (int)verts.size();
}

bool TAGraph::ValidEdge(int vf, int vt) const
{
    return ValidVert(vf) && ValidVert(vt) && Vertex(vf).FindEdge(vt);
}


void TAGraph::AddVertex(const TVertexP& vert)
{
    int index = verts.size();
    verts.push_back(vert);
    vert->index = index;
}

void TAGraph::AddVertex(TVertex& org, const TVertexP& vert,
    const TEdgeP& edge)
{
    AddVertex(vert);
    org.Attach(edge, vert);
}

void TAGraph::MakeEdge(TVertex& org, TVertex& term, const TEdgeP& edge)
{
    org.Attach(edge, &term);
}

void TAGraph::DeleteVertex(TVertex* vert)
{
    int pos = -1;
    int i;
    for (i = 0; i < NumVerts(); i++) {
        TVertex* cvert = &Vertex(i);
        if (cvert != vert)
            cvert->Detach(vert);
        else
            pos = i;
    }
    if (pos >= 0) {
        verts.erase(first() + pos);
        for (i = 0; i < NumVerts(); i++)
            Vertex(i).index = i;
    }
}

