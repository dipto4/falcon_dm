#pragma once
#include "falcon.h"

namespace falcon::forces {
    
    class GenericForce {
        protected:
            const int type {} ;
            GenericForce(int init_type): type {init_type}
            {}
        public:
            GenericForce() = default;
            void calculateForce() ;
    };
    
    class NewtonianForce : public GenericForce {
        
        public:
            NewtonianForce() : GenericForce(0)
            {}
    };

    class PN1Force : public GenericForce {
        public:
            PN1Force() : GenericForce(2)
            {}
    };

    class PN2Force : public GenericForce {
        public: 
            PN2Force() : GenericForce(4)
            {}

    };

    class PN25Force : public GenericForce {
        public:
            PN25Force() : GenericForce(5)
            {}
    };
}
