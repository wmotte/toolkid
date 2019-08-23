#include "Utils.h"

using namespace gsis;

void Utils::moveToNextCharacter(istream & in){
    //skip all white or tab character
    char ch;
    while (true){
        ch = in.get();
        if (ch != ' ' && ch != '\t'){
            in.putback(ch);
            break;
        }
    }
}

void Utils::skipNewLineCharacter(istream & in ){
    char ch;
    while (true){
        ch = in.get();
        if (ch != '\n'){
            in.putback(ch);
            break;
        }
    }
}



