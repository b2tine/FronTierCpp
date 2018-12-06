#include "drag_proto.h"

std::vector<Drag*> Drag::prototypes;

DragProtoInit::DragProtoInit() {
    Drag::prototypes.push_back(&pd);
    Drag::prototypes.push_back(&ld); 
    Drag::prototypes.push_back(&gd);
    Drag::prototypes.push_back(&md);
    Drag::prototypes.push_back(&sd);
    Drag::prototypes.push_back(&zsd);
    Drag::prototypes.push_back(&cud);
    Drag::prototypes.push_back(&rd);
    Drag::prototypes.push_back(&gbd);
    Drag::prototypes.push_back(&cpd);
    Drag::prototypes.push_back(&spd);
    Drag::prototypes.push_back(&rold);
    Drag::prototypes.push_back(&ogmd);
    Drag::prototypes.push_back(&ag);
}
