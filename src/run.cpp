#include "simulator.h"
int main(int argc, char *argv[])
{
    Simulator simulator("../config.cfg", argc, argv);
    simulator.run();
    return 0;
}