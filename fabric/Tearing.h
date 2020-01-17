#ifndef TEARING_H
#define TEARING_H

#include <FronTier.h>

#include <vector>


class FabricEdge
{
    public:

        bool sorted;
        double length0;

        double length() const
        {
            return separation(beg,end,3);
        }

    private:

        POINT* beg;
        POINT* end;
};


class FabricTearer
{
    public:

        std::vector<FabricEdge> edges;

        void collectFabricEdges(const INTERFACE* intfc);


    private:
}


#endif
