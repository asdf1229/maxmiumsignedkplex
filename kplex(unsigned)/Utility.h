#ifndef _UTILITY_H_
#define _UTILITY_H_

#define miv(a,b) ((a)>(b)?(b):(a))
#define mav(a,b) ((a)<(b)?(b):(a))

#include <iostream>
#include <fstream>
#include <vector>
#include <string.h>
#include <cstdlib>
#include <map>
#include <unordered_map>
#include <set>
#include <unordered_set>
#include <queue>
#include <list>
#include <cmath>
#include <ctime>
#include <algorithm>
// #define NDEBUG
#include <assert.h>

#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include <set>

using namespace std;

typedef unsigned int ui;
const int INF = 1000000000;

#include <cassert>

using ept = unsigned long;

#define pb push_back
#define mp make_pair

class Utility {
public:
    static FILE *open_file(const char *file_name, const char *mode) {
        FILE *f = fopen(file_name, mode);
        if(f == nullptr) {
            printf("Can not open file: %s\n", file_name);
            exit(1);
        }

        return f;
    }

    static std::string integer_to_string(long long number) {
        std::vector<ui> sequence;
        if(number == 0) sequence.pb(0);
        while(number > 0) {
            sequence.pb(number%1000);
            number /= 1000;
        }

        char buf[5];
        std::string res;
        for(ui i = sequence.size();i > 0;i --) {
            if(i == sequence.size()) sprintf(buf, "%u", sequence[i-1]);
            else sprintf(buf, ",%03u", sequence[i-1]);
            res += std::string(buf);
        }
        return res;
    }
};

struct Edge{
    ui first, second;
    int sign;
};

static string integer_to_string(long long number) {
    std::vector<ui> sequence;
    if(number == 0) sequence.push_back(0);
    while(number > 0) {
        sequence.push_back(number%1000);
        number /= 1000;
    }
    
    char buf[5];
    std::string res;
    for(unsigned int i = sequence.size();i > 0;i --) {
        if(i == sequence.size()) sprintf(buf, "%u", sequence[i-1]);
        else sprintf(buf, ",%03u", sequence[i-1]);
        res += std::string(buf);
    }
    return res;
}

template <typename T>
int find(T *st, T *ed, T x)
{
    return lower_bound(st, ed, x) - st;
}

#endif /* _UTILITY_H_ */

