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

#ifndef __GLOBAL_H
#define __GLOBAL_H

#pragma warning(disable : 4786)

#include <string>
#include <sstream>
#include <iomanip>
#include <utility>
#include <vector>
#include <cstdlib>
//#include <io.h>
#define _USE_MATH_DEFINES 
#include <math.h>

using namespace std;


inline float square(float x)
{
    return x * x;
}

double roundf(double value, double factor = 1, int dir = 0);
double roundb(double value, double bound = 1, int dir = 0);
double roundp(double value, int dig);
int sign(double x);

const float pi = (float)M_PI;

int getBits(int& packed, int bits);

string safesstr(const string& str, size_t start = 0, size_t len = string::npos);
string fullstrip(const string& str, char c = ' ');
string fullstrip(const string& str, const string& s);

int explode(vector<string>& store, const string& split, const string& str, int merge = 0);
int explode_sep(vector<string>& store, const string& sep, const string& str, int merge = 0);

string toupper(const string& str);
string tolower(const string& str);
string str_replace(const string& pat, const string& rep, const string& str);

class bad_conv {
public:
    string msg;
    bad_conv();
    bad_conv(const string& amsg);
    bad_conv(const bad_conv& bc);
    virtual ~bad_conv();
    virtual void Throw() const { throw* this; }
    friend ostream& operator<< (ostream& str, const bad_conv& bc);
    virtual void put(ostream& str) const;
};

template<class T> T stov(const string& str)
{
    T val;
    istringstream buf(str);
    buf >> val;
    if (buf.fail())
        throw bad_conv(str);
    return val;
}

template<> string stov<string>(const string& str);


string operator+ (const string& str, int val);

template<class T> string vtos(const T& val)
{
    ostringstream buf;
    buf << val;
    return buf.str();
}

string vtosp(double val, int prec);
string vtospx(double val, int prec);

template<class T> string vtosm(const T& val)
{
    ostringstream buff;
    buff.setf(ios::fixed, ios::floatfield);
    buff << val;
    string strf = buff.str();
    ostringstream bufs;
    bufs.setf(ios::scientific, ios::floatfield);
    bufs << val;
    string strs = bufs.str();
    return strf.length() > strs.length() ? strs : strf;
}

void splitpath(const string& path, string& dir, string& name, string& ext);
string getfilepath(const string& pathname);
string getfilename(const string& pathname);
extern string AppPath;

float randunit();
float randuniform(float low, float high);
int randomScaledLinear(int num);
int randomScaledUniform(int num);
int randomScaledCustom(const vector<float>& vect);
int randomScaledCustomW(const vector<float>& vect);

struct fpoint
{
    float x, y;
    fpoint()
        : x(0), y(0)
    {
    }
    fpoint(float vx, float vy)
        : x(vx), y(vy)
    {
    }
    fpoint& operator+= (const fpoint& pt)
    {
        x += pt.x;
        y += pt.y;
        return *this;
    }

    fpoint& operator-= (const fpoint& pt)
    {
        x -= pt.x;
        y -= pt.y;
        return *this;
    }

    static float dist(const fpoint& p1, const fpoint& p2);
};




inline int random(int __num)
{
    return (int)(((unsigned long)rand() * (unsigned long)__num) / (1ul + RAND_MAX));
}

//inline void randomize(void) { srand(clock()); }

template<class value, class limit>
void truncVal(value& val, limit lim)
{
    if (val > lim)
        val = (value)lim;
}



class delim {
    char value;
public:
    delim(char val)
        : value(val)
    {
    }

    inline friend istream& operator>>(istream& stream, const delim& d)
    {
        stream.ignore(256, d.value);
        return stream;
    }
};


class strparam {
protected:
    ios* str;
    strparam* prev;
    int block;

    static int cblock;
    static strparam* cur;

    void reset();
    void restore();
    virtual void set() = 0;
    virtual void clear() = 0;
    ~strparam();
public:
    strparam();
    friend istream& operator>>(istream& stream, strparam& object);
    friend ostream& operator<<(ostream& stream, strparam& object);
    friend ios& strpset(ios&);
    friend ios& strpreset(ios&);
};

class strw : public strparam {
    int oldw;
    int neww;
public:
    void set();
    void clear();
    strw(int width = 0);
    ~strw();
};

class strf : public strparam {
    ios::fmtflags oldf;
    ios::fmtflags newset, newclear;
public:
    void set();
    void clear();
    strf();
    strf(ios::fmtflags fset);
    strf(ios::fmtflags fset, ios::fmtflags fclear);
    ~strf();
};

class strp : public strparam {
    char oldpad;
    char newpad;
public:
    void set();
    void clear();
    strp(char pad = ' ');
    ~strp();
};

class strd : public strparam {
    int oldd;
    int newd;
public:
    void set();
    void clear();
    strd(int dig = 0);
    ~strd();
};


#endif
