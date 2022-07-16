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

#ifndef __RCPTR_H
#define __RCPTR_H

class __rc_ptr_base {
protected:
    __rc_ptr_base()
        : prev(0), next(0)
    {
    }

    __rc_ptr_base(const __rc_ptr_base& p)
        : prev(0), next(0)
    {
        slip_in(p);
    }

    void slip_in(const __rc_ptr_base& p)
    {
        prev = &p;
        next = p.next;
        p.next = this;
        if (next)
            next->prev = this;
    }

    int slip_out()
    {
        if (prev || next)
        {
            if (prev)
                prev->next = next;
            if (next)
                next->prev = prev;
            prev = next = 0;
            return 0;
        }
        else
        {
            return 1;
        }
    }

    mutable const __rc_ptr_base* prev, * next;
};


template <class X> class rc_ptr : protected __rc_ptr_base {
public:
    rc_ptr(X* p = 0)
        : data(p), __rc_ptr_base()
    {
    }

    rc_ptr(const rc_ptr& p)
        : data(p.data), __rc_ptr_base(p)
    {
    }

    rc_ptr& operator=(const rc_ptr& p)
    {
        destroy();
        data = p.data;
        slip_in(p);
        return *this;
    }

    rc_ptr& operator=(X* p)
    {
        destroy();
        data = p;
        return *this;
    }

    ~rc_ptr()
    {
        destroy();
    }

    X& operator*() const
    {
        return *data;
    }

    X* operator->() const
    {
        return data;
    }

    operator X* () const
    {
        return data;
    }
protected:
    X* data;
    void destroy()
    {
        if (slip_out())
            delete data;
    }
};

#endif
