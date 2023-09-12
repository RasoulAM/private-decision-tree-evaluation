#ifndef __PARSER_H
#define __PARSER_H

#include <stdio.h>
#include <string>
#include "utils.h"

class Parser {
    FILE *p = nullptr;
    Node *parseNode();
    uint64_t parseInt();

    public: 
    Parser();
    Node *parseTree(const std::string &fname);
};

#endif 
