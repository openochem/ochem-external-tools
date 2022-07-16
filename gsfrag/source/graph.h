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

#ifndef __GRAPH_H
#define __GRAPH_H

#include <vector>
using namespace std;

#include "rcptr.h"

class TVertex;
class TEdge;
class TEdgeMirror;
class TGraph;
class TAGraph;
class TSubGraph;

class TEdge
{
public:
    TEdge();
    ~TEdge()
    {}
    virtual TEdge* InfoEdge()
    {
        return this;
    }
    virtual const TEdge* InfoEdge() const
    {
        return this;
    }
    virtual bool IsInfoEdge() const
    {
        return true;
    }

    void SetTerm(TVertex* termVert)
    {
        vterm = termVert;
    }
    TVertex& Term() const
    {
        return *vterm;
    }
    virtual void Invert()
    {}
    virtual TEdge* Clone() const
    {
        return new TEdge;
    }
protected:
    TVertex* vterm;
    TEdge(TVertex* termVert);

    friend class TVertex;
};

typedef rc_ptr<TEdge> TEdgeP;
typedef vector<TEdgeP> TEdges;
typedef TEdges::iterator TEdgeIter;
typedef TEdges::const_iterator TEdgeCIter;

class TEdgeMirror : public TEdge {
protected:
    TEdge* inf;
public:
    TEdgeMirror(TVertex* termVert, TEdge* infoEdge);
    ~TEdgeMirror()
    {}
    virtual TEdge* InfoEdge()
    {
        return inf;
    }
    virtual const TEdge* InfoEdge() const
    {
        return inf;
    }
    virtual bool IsInfoEdge()  const
    {
        return false;
    }
    virtual TEdge* Clone() const
    {
        return 0;
    }
};

class TVertex
{
protected:
    TEdges edges;
    int index;
public:
    TVertex();
    TVertex(const TVertex& vert);
    ~TVertex();
    TEdgeIter first()
    {
        return edges.begin();
    }
    TEdgeCIter first() const
    {
        return edges.begin();
    }
    TEdgeIter last()
    {
        return edges.end();
    }
    TEdgeCIter last() const
    {
        return edges.end();
    }
    int Index() const
    {
        return index;
    }
    TEdgeP FindEdge(TVertex* vert) const;
    TEdgeP FindEdge(int vert) const;
    int CountInfoEdges() const;
    int GetNumEdges() const { return edges.size(); }
    const TEdge* GetEdge(int i) const { return edges[i]; }
    void Attach(const TEdgeP& edge, TVertex* vert);
    void Detach(TVertex* vert);

    virtual TVertex* Clone() const
    {
        return new TVertex;
    }
    friend class TAGraph;
};

typedef rc_ptr<TVertex> TVertexP;
typedef vector<TVertexP> TVertexes;
typedef TVertexes::iterator TVertIter;
typedef TVertexes::const_iterator TVertCIter;

class TGraph
{
public:
    TGraph();
    virtual TVertCIter first() const = 0;
    virtual TVertCIter last() const = 0;

    virtual int NumVerts() const = 0;
    vector<int> Components() const;
    vector<int> Components(vector<int>& comp, vector<int>& touch) const;
    int NumCompon() const;
    int IsConnected() const;
    int MaxIndex() const
    {
        return last() - first();
    }
    virtual vector<int> VertMask() const;
    TVertex& operator [](int i)
    {
        return Vertex(i);
    }

    virtual	bool ValidVert(int vert) const = 0;
    bool ValidVert(const TVertex& vert) const
    {
        return ValidVert(vert.Index());
    }
    virtual	bool ValidEdge(int vf, int vt) const = 0;
    int Valence(int vert) const
    {
        return Valence(Vertex(vert));
    }
    virtual int Valence(const TVertex& vert) const;
    int NumEdges() const;
    virtual TGraph* Clone() const = 0;
    virtual TVertex& Vertex(int i) const = 0;
    virtual const TAGraph* BaseGraph() const = 0;
protected:
};

class TAGraph : public TGraph
{
public:
    TAGraph(const TAGraph& grp);
    TAGraph(TSubGraph& grp);
    TAGraph();
    ~TAGraph();
    virtual void Reset()
    {
        DoReset();
    }
    void DoReset()
    {
        verts.erase(first(), last());
    }
    virtual int NumVerts() const
    {
        return verts.size();
    }
    TVertIter first()
    {
        return verts.begin();
    }
    virtual TVertCIter first() const
    {
        return verts.begin();
    }
    TVertIter last()
    {
        return verts.end();
    }
    virtual TVertCIter last() const
    {
        return verts.end();
    }
    virtual vector<int> VertMask() const;
    virtual bool ValidVert(int vert) const;
    virtual	bool ValidEdge(int vf, int vt) const;

    void AddVertex(const TVertexP& vert);
    void AddVertex(TVertex& org, const TVertexP& vert, const TEdgeP& edge);
    void AddVertex(int org, const TVertexP& vert, const TEdgeP& edge)
    {
        AddVertex(Vertex(org), vert, edge);
    }
    void MakeEdge(TVertex& org, TVertex& term, const TEdgeP& edge);
    void MakeEdge(int org, int term, const TEdgeP& edge)
    {
        MakeEdge(Vertex(org), Vertex(term), edge);
    }
    void DeleteVertex(TVertex* vert);

    virtual TVertex& Vertex(int i) const
    {
        return *verts[i];
    }
    int Valence(int vert) const //****using
    {
        return TGraph::Valence(vert);
    }
    virtual int Valence(const TVertex& vert) const
    {
        return vert.edges.size();
    }

    virtual TGraph* Clone() const
    {
        return new TAGraph(*this);
    }

    virtual const TAGraph* BaseGraph() const
    {
        return this;
    }

    friend class TSubGraph;
protected:
    TVertexes verts;
};

typedef rc_ptr<TAGraph> TAGraphP;

typedef vector<pair<int, int> > PVect;


#endif
