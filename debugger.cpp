#ifdef __DEBUG
#include "debugger.hpp"

#include <string>

using namespace std;

Debugger::Debugger(string filename, level l) {
  min = l;
  fname = filename;
  output.open(filename);
  if(!output.is_open()) {
    cout << "Can't open file: '"<<filename<<"'."<<endl; 
  }
}

map<string, Debugger*> Debugger::instances;

Debugger* Debugger::get_instance(string filename, level l) {
  if(instances.find(filename) == instances.end()) {
    Debugger* d = new Debugger(filename, l);
    
    instances.insert(make_pair(filename, d));
  }
  return instances[filename];  
}

Debugger* Debugger::get_instance(string filename) {
  return get_instance(filename, DEFAULT_LEVEL); 
}
      

void Debugger::log(string s, level l){
  string pre;
  if(l == ERROR) pre = "ERROR";
  if(l == WARNING) pre="WARNING";
  if(l == NOTE) pre  = "NOTE";
  if(l==DEBUG)pre=   "DEBUG";
  if(this->min & l) {
    time_t raw;
    struct tm * ti;
    time(&raw);
    ti = localtime(&raw);
    char tmp[80];
    strftime(tmp, 80, "%d.%m.%Y %T", ti);
    output << "[" << tmp << "][" << pre <<"]: "<<s<<endl;
  }
}

/*int main() {
  Debugger* x = Debugger::get_instance("./file.log", Debugger::ALL);
  cout << Debugger::ALL << endl;
  x->log("adsjfklasd", Debugger::NOTE);
}*/
#endif