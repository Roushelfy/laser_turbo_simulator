#include "simulator.h"
int main(int argc, char *argv[])
{
    Simulator simulator("../config.cfg", argc, argv);
    if (argc > 1 && strcmp(argv[1], "-test_dis") == 0)
    {
        simulator.test_object_distance();
    }
    else
        simulator.run();
    return 0;
}