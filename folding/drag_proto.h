#ifndef DRAGPROTO_H
#define DRAGPROTO_H

#include "drag.h"
#include "origami.h"
#include "singleton.h"

class DragProtoInit : public Singleton<DragProtoInit> {
    PointDrag pd;
    LineDrag ld;
    GravityDrag gd;
    MultiplePointDrag md;
    FoldDrag sd;
    ZFoldDrag zsd;
    CloseUmbrellaDrag cud;
    RelaxDrag rd;
    GravityBoxDrag gbd;
    CompressDrag cpd;
    SeparateDrag spd;
    RollDrag rold;
    OrigamiFold ogmd;
    AlignDrag ag;
protected:
    friend class Singleton<DragProtoInit>;
    DragProtoInit();
};

#endif
