#ifndef SINGLE_H
#define SINGLE_H

// used to perform curiously singleton
template<class T>
class Singleton {
    Singleton(const Singleton&);
    Singleton& operator=(const Singleton&);
protected:
    Singleton() {}
    virtual ~Singleton() {}
public:
    static T& instance() {
        static T theInstance;
        return theInstance;
    }
};

#endif
