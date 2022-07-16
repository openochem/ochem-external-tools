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

#include <iostream>
#include <sstream>
#include <locale>
//#include <io.h>
#include <stdio.h>

#include "global.h"


string safesstr(const string& str, size_t start, size_t len)
{
    size_t length = str.length();
    if (start >= length)
        length = 0;
    else
        length -= start;
    if (length > len)
        length = len;
    if (length)
        return str.substr(start, len);
    else
        return string();
}

string fullstrip(const string& str, char c)
{
    size_t start = str.find_first_not_of(c);
    size_t end = str.find_last_not_of(c);
    if (start == string::npos || end == string::npos)
        return string();
    else
        return str.substr(start, end - start + 1);
}

string fullstrip(const string& str, const string& s)
{
    size_t start = str.find_first_not_of(s);
    size_t end = str.find_last_not_of(s);
    if (start == string::npos || end == string::npos)
        return string();
    else
        return str.substr(start, end - start + 1);
}


int explode(vector<string>& store, const string& split, const string& str, int merge)
{
    store.clear();
    store.reserve(4);
    string buf = str;
    size_t start = 0;
    while (start != string::npos)
    {
        size_t stop = str.find_first_of(split, start);
        size_t len = string::npos;
        if (stop != string::npos)
            len = stop - start;
        string cur = safesstr(buf, start, len);
        if (!cur.empty() || !merge)
            store.push_back(cur);
        if (stop != string::npos)
            stop++;
        start = stop;
    }
    return store.size();
}

int explode_sep(vector<string>& store, const string& sep, const string& str, int merge)
{
    store.clear();
    store.reserve(4);
    string buf = str;
    size_t start = 0;
    size_t seplen = sep.length();
    while (start != string::npos)
    {
        size_t stop = str.find(sep, start);
        size_t len = string::npos;
        if (stop != string::npos)
            len = stop - start;
        string cur = safesstr(buf, start, len);
        if (!cur.empty() || !merge)
            store.push_back(cur);
        if (stop != string::npos)
            stop += seplen;
        start = stop;
    }
    return store.size();
}

string toupper(const string& str)
{
    string tmp;
    tmp.reserve(str.size());
    for (string::const_iterator i = str.begin(); i != str.end(); i++)
    {
        tmp += toupper(*i, locale::classic());
    }
    return tmp;
}

string tolower(const string& str)
{
    string tmp;
    tmp.reserve(str.size());
    for (string::const_iterator i = str.begin(); i != str.end(); i++)
    {
        tmp += tolower(*i, locale::classic());
    }
    return tmp;
}

string str_replace(const string& pat, const string& rep, const string& str)
{
    string tmp = str;
    string::size_type pos = string::npos;
    int len = pat.length();
    while ((pos = tmp.find(pat)) != string::npos)
        tmp.replace(pos, len, rep);
    return tmp;
}

bad_conv::bad_conv()
{
}

bad_conv::bad_conv(const string& amsg)
    : msg(amsg)
{
}

bad_conv::bad_conv(const bad_conv& bc)
{
}

bad_conv::~bad_conv()
{
}

ostream& operator<< (ostream& str, const bad_conv& bc)
{
    bc.put(str);
    str << endl;
    return str;
}

void bad_conv::put(ostream& str) const
{
    str << msg;
}

string operator+ (const string& str, int val)
{
    ostringstream buf;
    buf << val;
    return str + buf.str();
}

void splitpath(const string& path, string& dir, string& name, string& ext)
{
    string tmp = path;
    int start = tmp.find_last_of("\\/");
    if (start != string::npos)
    {
        dir = safesstr(tmp, 0, start + 1);
        tmp = safesstr(tmp, start + 1);
    }
    else
    {
        if (tmp[1] == ':')
        {
            dir = safesstr(tmp, 0, 2);
            tmp = safesstr(tmp, 2);
        }
        else
        {
            dir = string();
        }
    }
    start = tmp.find_last_of('.');
    if (start != string::npos && start > 0)
    {
        ext = safesstr(tmp, start);
        tmp = safesstr(tmp, 0, start);
    }
    else
    {
        ext = string();
    }
    name = tmp;
}

string getfilepath(const string& pathname)
{
    string dir, name, ext;
    splitpath(pathname, dir, name, ext);
    return dir;
}

string getfilename(const string& pathname)
{
    string dir, name, ext;
    splitpath(pathname, dir, name, ext);
    return name + ext;
}

string AppPath;

int getBits(int& packed, int bits)
{
    int mask = ~((~0) << bits);
    int retval = packed & mask;
    packed >>= bits;
    return retval;
}

double roundf(double value, double factor, int dir)
{
    if (factor < 0)
        factor = -factor;
    if (factor < 1e-6)
        factor = -1e-6;
    double temp = value / factor;
    double tempint;
    if (dir > 0)
        tempint = ceil(temp);
    else if (dir < 0)
        tempint = floor(temp);
    else
    {
        temp = modf(temp + 0.5 * sign(value), &tempint);
        // Handle rounding of 0.5 in a special manner
        if (!temp)
        {
            if (modf(tempint / 2, &temp)) // odd number
                tempint -= sign(value);
        }
    }
    return tempint * factor;
}

double roundb(double value, double bound, int dir)
{
    if (bound < 0)
        bound = -bound;
    if (bound < 1e-6)
        bound = 1e-6;
    double bound1 = log10(bound);
    double bound2 = floor(bound1);
    bound2 = pow(10., bound2);
    bound1 = bound / bound2;
    if (bound1 >= 5)
        bound2 *= 5;
    else if (bound1 >= 2)
        bound2 *= 2;
    return roundf(value, bound2, dir);
}


double roundp(double value, int dig)
{
    return roundf(value, pow(10., -dig), 0);
}

int sign(double x)
{
    if (x < 0)
        return -1;
    else if (x > 0)
        return 1;
    return 0;
}


template<> string stov<string>(const string& str)
//template<> string stov<string>(const string &str)
{
    string val;
    val = fullstrip(str);
    return val;
}



string vtosp(double val, int prec)
{
    val = roundp(val, prec);
    char buf[128];
    if (prec <= 0)
    {
        sprintf(buf, "%.10G", val);
    }
    else
    {
        if (val)
        {
            prec += 1 + (int)floor(log10(val * sign(val)));
        }
        char buffmt[16];
        sprintf(buffmt, "%%#.%dG", prec);
        sprintf(buf, buffmt, val);
    }
    return buf;
}

string vtospx(double val, int prec)
{
    if (roundf(val) == val)
        prec = 0;
    return vtosp(val, prec);
}


float randunit()
{
    return rand() / float(RAND_MAX);
}

float randuniform(float low, float high)
{
    return low + (high - low) * randunit();
}


int randomScaledLinear(int num)
{
    float random = randunit();
    float pos = (sqrt(1 + 4 * num * (num + 1) * random) - 1) / 2;
    return min(1 + int(pos), num);
}

int randomScaledUniform(int num)
{
    float pos = num * randunit();
    return min(1 + int(pos), num);
}

int randomScaledCustom(const vector<float>& vect)
{
    float sum = 0;
    int i;
    int num = vect.size();
    for (i = 0; i < num; i++)
        sum += vect[i];
    float pos = sum * randunit();
    for (i = 0; i < num; i++)
    {
        pos -= vect[i];
        if (pos <= 0)
            break;
    }
    return min(1 + i, num);
}

int randomScaledCustomW(const vector<float>& vect)
{
    float sum = 0, minval = 0;
    int i;
    int num = vect.size();
    for (i = 0; i < num; i++)
    {
        float tmp = vect[i];
        sum += tmp;
        if (i == 0 || tmp < minval)
            minval = tmp;
    }
    sum -= num * minval;
    float fact = 1 / (1 + sum);
    float factn = fact / num;
    float pos = randunit();
    for (i = 0; i < num; i++)
    {
        pos -= fact * (vect[i] - minval) + factn;
        if (pos <= 0)
            break;
    }
    return min(1 + i, num);
}


float fpoint::dist(const fpoint& p1, const fpoint& p2)
{
    return sqrt(square(p1.x - p2.x) + square(p1.y - p2.y));
}



strparam* strparam::cur = 0;
int strparam::cblock = 0;

strparam::strparam()
{
    str = 0;
    prev = 0;
    block = 0;
}

strparam::~strparam()
{
    restore();
}

istream& operator>>(istream& stream, strparam& object)
{
    object.str = &stream;
    object.reset();
    return stream;
}

ostream& operator<<(ostream& stream, strparam& object)
{
    object.str = &stream;
    object.reset();
    return stream;
}

ios& strpset(ios& o)
{
    ++strparam::cblock;
    return o;
}

ios& strpreset(ios& o)
{
    if (strparam::cblock)
    {
        if (strparam::cur)
        {
            int curlevel = strparam::cblock;
            while (strparam::cur && strparam::cur->block == curlevel)
                strparam::cur->restore();
        }
        --strparam::cblock;
    }
    return o;
}

void strparam::reset()
{
    prev = cur;
    cur = this;
    block = cblock;
    set();
}

void strparam::restore()
{
    strparam* ptr;
    for (ptr = cur; ptr && ptr != this; ptr = ptr->prev);
    if (ptr)
    {
        do
        {
            ptr = cur;
            if (ptr->str)
            {
                ptr->clear();
                ptr->str = 0;
            }
            cur = ptr->prev;
        } while (cur && ptr != this);
    }
}

strw::strw(int width)
    : neww(width)
{
}

strw::~strw()
{
    restore();
}

void strw::set()
{
    oldw = str->width(neww);
}

void strw::clear()
{
    if (str)
        str->width(oldw);
}

strf::strf()
    : newset(ios::skipws), newclear((ios::fmtflags)~0L)
{
}

strf::strf(ios::fmtflags fset)
    : newset(fset), newclear((ios::fmtflags)0)
{
}

strf::strf(ios::fmtflags fset, ios::fmtflags fclear)
    : newset(fset), newclear(fclear)
{
}

strf::~strf()
{
    restore();
}

void strf::set()
{
    oldf = str->flags();
    str->unsetf(newclear);
    str->setf(newset);
}

void strf::clear()
{
    if (str)
        str->flags(oldf);
}

strp::strp(char pad)
    : newpad(pad)
{
}

strp::~strp()
{
    restore();
}

void strp::set()
{
    oldpad = str->fill(newpad);
}

void strp::clear()
{
    if (str)
        str->fill(oldpad);
}

strd::strd(int dig)
    : newd(dig)
{
}

strd::~strd()
{
    restore();
}

void strd::set()
{
    oldd = str->precision(newd);
}

void strd::clear()
{
    if (str)
        str->precision(oldd);
}

