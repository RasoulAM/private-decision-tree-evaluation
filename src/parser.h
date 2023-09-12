#ifndef __PARSER_H
#define __PARSER_H

#include <stdio.h>
#include <string>
#include "tree_utils.h"

class Parser {
    FILE *p = nullptr;
    DecisionTreeNode *parseNode();
    uint64_t parseInt();

    public: 
    Parser();
    DecisionTreeNode *parseTree(const std::string &fname);
};

#endif 
