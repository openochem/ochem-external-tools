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

FragBaseL unlabfrag = {
    //1
    {2,
    {{0, 1},
     {1, 0}},
    "p1",
    1,{0},
    {1}
    },
    //2
    {3,
    {{0, 1, 0},
     {1, 0, 1},
     {0, 1, 0}},
    "p2",
    2,{0,1},
    {1, 2}
    },
    //3
    {3,
    {{0, 1, 1},
     {1, 0, 1},
     {1, 1, 0}},
    "c3",
    1,{0},
    {2}
    },
    //4
    {4,
    {{0, 1, 0, 0},
     {1, 0, 0, 0},
     {0, 0, 0, 1},
     {0, 0, 1, 0}},
    "p1p1",
    1,{0},
    {2}
    },
    {4,
    {{0, 1, 0, 0},
     {1, 0, 1, 0},
     {0, 1, 0, 1},
     {0, 0, 1, 0}},
    "p3",
    2,{0,1},
    {1,1}
    },
    {4,
    {{0, 1, 1, 1},
     {1, 0, 1, 0},
     {1, 1, 0, 0},
     {1, 0, 0, 0}},
    "c3A",
    3,{0,1,3},
    {2,1,2}
    },
    {4,
    {{0, 1, 1, 0},
     {1, 0, 0, 1},
     {1, 0, 0, 1},
     {0, 1, 1, 0}},
    "c4",
    1,{0},
    {2}
    },
    //8
    {5,
    {{0, 1, 0, 0, 0},
     {1, 0, 0, 0, 0},
     {0, 0, 0, 1, 0},
     {0, 0, 1, 0, 1},
     {0, 0, 0, 1, 0}},
    "p1p2",
    3,{0,2,3},
    {2, 2, 4}
    },
    {5,
    {{0, 1, 0, 0, 0},
     {1, 0, 0, 0, 0},
     {0, 0, 0, 1, 1},
     {0, 0, 1, 0, 1},
     {0, 0, 1, 1, 0}},
    "p1c3",
    2,{0,2},
    {6, 4}
    },
    {5,
    {{0, 1, 1, 1, 0},
     {1, 0, 1, 0, 1},
     {1, 1, 0, 0, 0},
     {1, 0, 0, 0, 0},
     {0, 1, 0, 0, 0}},
    "c3AB",
    3,{0,2,3},
    {1,2,1}
    },
    {5,
    {{0, 1, 0, 1, 1},
     {1, 0, 1, 0, 0},
     {0, 1, 0, 1, 0},
     {1, 0, 1, 0, 0},
     {1, 0, 0, 0, 0}},
    "c4A",
    4,{0,1,2,4},
    {2,1,2,1}
    },
    {5,
    {{0, 1, 0, 0, 1},
     {1, 0, 1, 0, 0},
     {0, 1, 0, 1, 0},
     {0, 0, 1, 0, 1},
     {1, 0, 0, 1, 0}},
    "c5",
    1,{0},
    {2}
    },
    {5,
    {{0, 1, 0, 0, 0},
     {1, 0, 1, 0, 0},
     {0, 1, 0, 1, 0},
     {0, 0, 1, 0, 1},
     {0, 0, 0, 1, 0}},
    "p4",
    3,{0,1,2},
    {1,1,2}
    },
    //14
    {6,
    {{0, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0},
     {0, 0, 0, 1, 0, 0},
     {0, 0, 1, 0, 0, 0},
     {0, 0, 0, 0, 0, 1},
     {0, 0, 0, 0, 1, 0}},
    "p1p1p1",
    1,{0},
    {8}
    },
    {6,
    {{0, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0},
     {0, 0, 0, 1, 0, 0},
     {0, 0, 1, 0, 1, 0},
     {0, 0, 0, 1, 0, 1},
     {0, 0, 0, 0, 1, 0}},
    "p1p3",
    3,{0,2,3},
    {2,2,2}
    },
    {6,
    {{0, 1, 0, 0, 0, 0},
     {1, 0, 1, 0, 0, 0},
     {0, 1, 0, 0, 0, 0},
     {0, 0, 0, 0, 1, 0},
     {0, 0, 0, 1, 0, 1},
     {0, 0, 0, 0, 1, 0}},
    "p2p2",
    2,{0,1},
    {2,4}
    },
    {6,
    {{0, 1, 0, 0, 0, 0},
     {1, 0, 1, 0, 0, 0},
     {0, 1, 0, 0, 0, 0},
     {0, 0, 0, 0, 1, 1},
     {0, 0, 0, 1, 0, 1},
     {0, 0, 0, 1, 1, 0}},
    "p2c3",
    3,{0,1,3},
    {6,12,4}
    },
    {6,
    {{0, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0},
     {0, 0, 0, 1, 1, 1},
     {0, 0, 1, 0, 1, 0},
     {0, 0, 1, 1, 0, 0},
     {0, 0, 1, 0, 0, 0}},
    "p1c3A",
    4,{0,2,3,5},
    {2,4,2,4}
    },
    {6,
    {{0, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0},
     {0, 0, 0, 1, 0, 1},
     {0, 0, 1, 0, 1, 0},
     {0, 0, 0, 1, 0, 1},
     {0, 0, 1, 0, 1, 0}},
    "p1c4",
    2,{0,2},
    {8,4}
    },
    {6,
    {{0, 1, 1, 0, 0, 0},
     {1, 0, 1, 0, 0, 0},
     {1, 1, 0, 0, 0, 0},
     {0, 0, 0, 0, 1, 1},
     {0, 0, 0, 1, 0, 1},
     {0, 0, 0, 1, 1, 0}},
    "c3c3",
    1,{0},
    {12}
    },
    {6,
    {{0, 1, 0, 1, 1, 0},
     {1, 0, 1, 0, 0, 1},
     {0, 1, 0, 1, 0, 0},
     {1, 0, 1, 0, 0, 0},
     {1, 0, 0, 0, 0, 0},
     {0, 1, 0, 0, 0, 0}},
    "c4AB",
    3,{0,2,4},
    {1,1,1}
    },
    {6,
    {{0, 1, 0, 1, 1, 0},
     {1, 0, 1, 0, 0, 0},
     {0, 1, 0, 1, 0, 1},
     {1, 0, 1, 0, 0, 0},
     {1, 0, 0, 0, 0, 0},
     {0, 0, 1, 0, 0, 0}},
    "c4AC",
    3,{0,1,4},
    {2,2,2}
    },
    {6,
    {{0, 1, 1, 1, 0, 0},
     {1, 0, 1, 0, 1, 0},
     {1, 1, 0, 0, 0, 1},
     {1, 0, 0, 0, 0, 0},
     {0, 1, 0, 0, 0, 0},
     {0, 0, 1, 0, 0, 0}},
    "c3ABC",
    2,{0,3},
    {2,2}
    },
    //24
    {6,
    {{0, 1, 0, 0, 0, 1},
     {1, 0, 1, 0, 0, 0},
     {0, 1, 0, 1, 0, 0},
     {0, 0, 1, 0, 1, 0},
     {0, 0, 0, 1, 0, 1},
     {1, 0, 0, 0, 1, 0}},
    "c6",
    1,{0},
    {2}
    },
    {6,
    {{0, 1, 0, 0, 1, 1},
     {1, 0, 1, 0, 0, 0},
     {0, 1, 0, 1, 0, 0},
     {0, 0, 1, 0, 1, 0},
     {1, 0, 0, 1, 0, 0},
     {1, 0, 0, 0, 0, 0}},
    "c5A",
    4,{0,1,2,5},
    {2,1,1,2}
    },
    {6,
    {{0, 1, 0, 0, 0, 0},
     {1, 0, 1, 0, 0, 0},
     {0, 1, 0, 1, 0, 0},
     {0, 0, 1, 0, 1, 0},
     {0, 0, 0, 1, 0, 1},
     {0, 0, 0, 0, 1, 0}},
    "p5",
    3,{0,1,2},
    {1,1,1}
    },
    //27
    {7,
    {{0, 1, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0},
     {0, 0, 0, 1, 0, 0, 0},
     {0, 0, 1, 0, 0, 0, 0},
     {0, 0, 0, 0, 0, 1, 0},
     {0, 0, 0, 0, 1, 0, 1},
     {0, 0, 0, 0, 0, 1, 0}},
    "p1p1p2",
    3,{0,4,5},
    {4,8,16}
    },
    {7,
    {{0, 1, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0},
     {0, 0, 0, 1, 0, 0, 0},
     {0, 0, 1, 0, 1, 0, 0},
     {0, 0, 0, 1, 0, 1, 0},
     {0, 0, 0, 0, 1, 0, 1},
     {0, 0, 0, 0, 0, 1, 0}},
    "p1p4",
    4,{0,2,3,4},
    {2,2,2,4}
    },
    {7,
    {{0, 1, 0, 0, 0, 0, 0},
     {1, 0, 1, 0, 0, 0, 0},
     {0, 1, 0, 0, 0, 0, 0},
     {0, 0, 0, 0, 1, 0, 0},
     {0, 0, 0, 1, 0, 1, 0},
     {0, 0, 0, 0, 1, 0, 1},
     {0, 0, 0, 0, 0, 1, 0}},
    "p2p3",
    4,{0,1,3,4},
    {2,4,2,2}
    },
    //30
    {7,
    {{0, 1, 0, 0, 0, 0, 0},
     {1, 0, 1, 0, 0, 0, 0},
     {0, 1, 0, 1, 0, 0, 0},
     {0, 0, 1, 0, 1, 0, 0},
     {0, 0, 0, 1, 0, 1, 0},
     {0, 0, 0, 0, 1, 0, 1},
     {0, 0, 0, 0, 0, 1, 0}},
    "p6",
    4,{0,1,2,3},
    {1,1,1,2}
    },
    {7,
    {{0, 1, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0},
     {0, 0, 0, 1, 0, 0, 0},
     {0, 0, 1, 0, 0, 0, 0},
     {0, 0, 0, 0, 0, 1, 1},
     {0, 0, 0, 0, 1, 0, 1},
     {0, 0, 0, 0, 1, 1, 0}},
    "p1p1c3",
    2,{0,4},
    {12,16}
    },
    {7,
    {{0, 1, 0, 0, 0, 0, 0},
     {1, 0, 1, 0, 0, 0, 0},
     {0, 1, 0, 1, 0, 0, 0},
     {0, 0, 1, 0, 0, 0, 0},
     {0, 0, 0, 0, 0, 1, 1},
     {0, 0, 0, 0, 1, 0, 1},
     {0, 0, 0, 0, 1, 1, 0}},
    "p3c3",
    3,{0,1,4},
    {6,6,4}
    },
    {7,
    {{0, 1, 1, 0, 0, 0, 0},
     {1, 0, 1, 0, 0, 0, 0},
     {1, 1, 0, 0, 0, 0, 0},
     {0, 0, 0, 0, 1, 1, 1},
     {0, 0, 0, 1, 0, 1, 0},
     {0, 0, 0, 1, 1, 0, 0},
     {0, 0, 0, 1, 0, 0, 0}},
    "c3c3A",
    4,{0,3,4,6},
    {4,12,6,12}
    },
    {7,
    {{0, 1, 1, 0, 0, 0, 0},
     {1, 0, 1, 0, 0, 0, 0},
     {1, 1, 0, 0, 0, 0, 0},
     {0, 0, 0, 0, 1, 0, 1},
     {0, 0, 0, 1, 0, 1, 0},
     {0, 0, 0, 0, 1, 0, 1},
     {0, 0, 0, 1, 0, 1, 0}},
    "c3c4",
    2,{0,3},
    {16,12}
    },
    {7,
    {{0, 1, 0, 0, 0, 0, 0},
     {1, 0, 1, 0, 0, 0, 0},
     {0, 1, 0, 0, 0, 0, 0},
     {0, 0, 0, 0, 1, 0, 1},
     {0, 0, 0, 1, 0, 1, 0},
     {0, 0, 0, 0, 1, 0, 1},
     {0, 0, 0, 1, 0, 1, 0}},
    "p2c4",
    3,{0,1,3},
    {8,16,4}
    },
    {7,
    {{0, 1, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0},
     {0, 0, 0, 1, 0, 0, 1},
     {0, 0, 1, 0, 1, 0, 0},
     {0, 0, 0, 1, 0, 1, 0},
     {0, 0, 0, 0, 1, 0, 1},
     {0, 0, 1, 0, 0, 1, 0}},
    "p1c5",
    2,{0,2},
    {10,4}
    },
    {7,
    {{0, 1, 0, 0, 0, 0, 1},
     {1, 0, 1, 0, 0, 0, 0},
     {0, 1, 0, 1, 0, 0, 0},
     {0, 0, 1, 0, 1, 0, 0},
     {0, 0, 0, 1, 0, 1, 0},
     {0, 0, 0, 0, 1, 0, 1},
     {1, 0, 0, 0, 0, 1, 0}},
    "c7",
    1,{0},
    {2}
    },
    {7,
    {{0, 1, 0, 1, 1, 0, 0},
     {1, 0, 1, 0, 0, 1, 0},
     {0, 1, 0, 1, 0, 0, 1},
     {1, 0, 1, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0},
     {0, 1, 0, 0, 0, 0, 0},
     {0, 0, 1, 0, 0, 0, 0}},
    "c4ABC",
    5,{0,1,3,4,5},
    {1,2,2,1,2}
    },
    {7,
    {{0, 1, 0, 0, 0, 0, 0},
     {1, 0, 1, 0, 0, 0, 0},
     {0, 1, 0, 0, 0, 0, 0},
     {0, 0, 0, 0, 1, 1, 1},
     {0, 0, 0, 1, 0, 1, 0},
     {0, 0, 0, 1, 1, 0, 0},
     {0, 0, 0, 1, 0, 0, 0}},
    "p2c3A",
    5,{0,1,3,4,6},
    {2,4,4,2,4}
    },
    {7,
    {{0, 1, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0},
     {0, 0, 0, 1, 1, 1, 0},
     {0, 0, 1, 0, 1, 0, 1},
     {0, 0, 1, 1, 0, 0, 0},
     {0, 0, 1, 0, 0, 0, 0},
     {0, 0, 0, 1, 0, 0, 0}},
    "p1c3AB",
    4,{0,2,4,5},
    {2,2,4,2}
    },
    {7,
    {{0, 1, 0, 0, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0},
     {0, 0, 0, 1, 0, 1, 1},
     {0, 0, 1, 0, 1, 0, 0},
     {0, 0, 0, 1, 0, 1, 0},
     {0, 0, 1, 0, 1, 0, 0},
     {0, 0, 1, 0, 0, 0, 0}},
    "p1c4A",
    5,{0,2,3,4,6},
    {2,4,2,4,4}
    },
    {7,
    {{0, 1, 0, 0, 1, 1, 0},
     {1, 0, 1, 0, 0, 0, 1},
     {0, 1, 0, 1, 0, 0, 0},
     {0, 0, 1, 0, 1, 0, 0},
     {1, 0, 0, 1, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0},
     {0, 1, 0, 0, 0, 0, 0}},
    "c5AB",
    4,{0,2,3,5},
    {1,1,2,1}
    },
    {7,
    {{0, 1, 0, 0, 1, 1, 0},
     {1, 0, 1, 0, 0, 0, 0},
     {0, 1, 0, 1, 0, 0, 1},
     {0, 0, 1, 0, 1, 0, 0},
     {1, 0, 0, 1, 0, 0, 0},
     {1, 0, 0, 0, 0, 0, 0},
     {0, 0, 1, 0, 0, 0, 0}},
    "c5AC",
    4,{0,1,3,5},
    {1,2,1,1}
    },
    //44
    {7,
    {{0, 1, 0, 0, 0, 1, 1},
     {1, 0, 1, 0, 0, 0, 0},
     {0, 1, 0, 1, 0, 0, 0},
     {0, 0, 1, 0, 1, 0, 0},
     {0, 0, 0, 1, 0, 1, 0},
     {1, 0, 0, 0, 1, 0, 0},
     {1, 0, 0, 0, 0, 0, 0}},
    "c6A",
    5,{0,1,2,3,6},
    {2,1,1,2,2}
    }
};//unlabfrag

int numunlabfrag = sizeof(unlabfrag) / sizeof(unlabfrag[0]); //number of fragment
//int numunlabfrag=238;