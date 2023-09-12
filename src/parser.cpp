#include "parser.h"

Parser::Parser() { }

inline void consume(FILE *p, int cnt) {
    while(cnt--) {
        getc(p);
    }
}

uint64_t Parser::parseInt() {
    uint64_t ret = 0;
    char ch = 0;
    while(true) {
        ch = getc(p);
        if(ch < '0' || ch > '9') {
            ungetc(ch, p);
            return ret;
        }
        ret = (ret * 10) + (ch - '0');
    }
}

DecisionTreeNode *Parser::parseNode() {
    DecisionTreeNode *ret = NULL;

    consume(p, 10);
    char ch = getc(p);
    
    if(ch == 's') {
        consume(p, 20);
        uint64_t threshold = parseInt();
        consume(p, 13);
        uint64_t feature = parseInt();
        consume(p, 10);
        auto left = parseNode();
        consume(p, 11);
        auto right = parseNode();
        ret = new DecisionTreeNode(threshold, feature);
        ret->left_child = left;
        ret->right_child = right;
        getc(p);
    } else {
        consume(p, 15);
        uint64_t value = parseInt();
        ret = new DecisionTreeNode(value, -1);
        getc(p);
    }

    return ret;
}

DecisionTreeNode *Parser::parseTree(const std::string &fname) {
    p = fopen(fname.c_str(), "r");
    assert(p);

    auto ret = parseNode();

    fclose(p);
    p = nullptr;
    
    return ret;
}
