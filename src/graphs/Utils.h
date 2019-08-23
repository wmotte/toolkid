#ifndef __UTILS_H
#define __UTILS_H
#include <iostream>

using namespace std;

namespace gsis{

class Utils{
    public:
        static void moveToNextCharacter(istream & in);

        static void skipNewLineCharacter(istream & in);
};
}
#endif
