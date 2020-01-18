#ifndef TEARING_H
#define TEARING_H

#include <FronTier.h>

#include <vector>


class FabricEdge
{
    public:

        double k;
        double tear_threshold;

        double length;
        double tension
        double disp[3];

            //bool visited {false};

        //FabricEdge() = default;

        FabricEdge(POINT* p1, POINT* p2)
            : beg{p1}, end{p2}
        {
            gindex_beg = Gindex(beg);
            gindex_end = Gindex(end);
            length0 = separation(beg,end,3);
        }

        bool checkForTear();

    private:

        double length0;

        //long int gindex;
        long int gindex_beg;
        long int gindex_end;

        POINT* beg;
        POINT* end;

        bool beg_weak {false};
        bool end_weak {false};

        FabricEdge* prev;
        FabricEdge* next;
        
        bool underTension();
        void computeLength();
        void computeDisplacement();
        double computeTension();
};


class FabricTearer
{
    public:

        virtual ~FabricTearer();
        FabricTearer(const Front* front);

        void tearFabric();

        //Testing Functions
        void setEdgeTension(int index, double T);
        void pringEdgeTensions();

    private:
        
        std::vector<int> tear_idx;
        std::vector<FabricEdge*> edges;

        void clearEdges();
        void collectFabricEdges(const Front* front);
        
        void checkForTearingEvents();
        void processTearingEvents();
};


#endif
