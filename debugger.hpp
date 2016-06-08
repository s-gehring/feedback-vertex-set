

#ifndef _H_DEBUGGER
#define _H_DEBUGGER



#include <stdio.h>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <stdint.h>




using namespace std;

#define DEFAULT_LEVEL ALL
typedef unsigned int level;
class Debugger {
  private: 
    static map<string, Debugger*> instances;
    string fname;
    ofstream output;
    level min;
    
    Debugger() {};
    Debugger(string filename, level l);
  public:
  
    static const level ERROR = 1;
    static const level WARNING = 2;
    static const level NOTE = 4;
    static const level DEBUG = 8;
    static const level ALL = 15;
  
  
    static Debugger* get_instance(string filename, level l);
  
    static Debugger* get_instance(string filename);
  
    void log(string s, level l);
  
};













#endif