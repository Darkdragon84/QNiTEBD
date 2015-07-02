#ifndef PARSER_H
#define PARSER_H

#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>

using namespace std;

class parser
{
    public:
//        parser(){cout<<"parser created"<<endl;};
        parser() = default;
        parser(int argc, char** argv){parse(argc,argv);};
//        virtual ~parser() {cout<<"parser destroyed"<<endl;};
        virtual ~parser() = default;

        void parse(int argc, char** argv);
        void PrintValues() const;
        template<typename T> bool GetValue(T& val, const string& name, bool abrt=false) const;
        template<typename T> bool GetValue(vector<T>& val, const string& name, bool abrt=false) const;

    protected:
        inline bool is_cfg_file(string arg);
        void parse_from_file(const string& filename);
        void parse_argument(const string& arg);
//        map<string,int> values_;
        map<string,string> values_;

        template<typename T> T sconv(const string& str) const;
    private:
};

template<> int parser::sconv<int>(const string& str) const;
template<> unsigned int parser::sconv<unsigned int>(const string& str) const;
template<> double parser::sconv<double>(const string& str) const;
template<> size_t parser::sconv<size_t>(const string& str) const;
template<> string parser::sconv<string>(const string& str) const;

template<typename T> bool parser::GetValue(T& val, const string& name, bool abrt) const
{
    auto it=values_.find(name); // look for parameter "name"
    bool found=(it!=values_.end()); // see if found
    if (found) val=sconv<T>(it->second); // if found, write value to external variable val
    else
    {
        if (abrt) throw invalid_argument(string("argument \"")+name+string("\" not found"));
        else cout<<"argument \""<<name<<"\" not found, using "<<name<<"="<<val<<endl;
    } // if not found, warn and do nothing to val
    return found;
}

template<typename T> bool parser::GetValue(vector<T>& vec, const string& name, bool abrt) const
{
    string tmp;
    auto it=values_.find(name); // look for parameter "name"
    bool found=(it!=values_.end()); // see if found
    if (found)
    {
        vec.clear();
        tmp = it->second;
        size_t pos1=0,pos2=0;

        while(pos2!=string::npos)
        {
            pos2 = tmp.find(",",pos1);
            vec.push_back(sconv<T>(tmp.substr(pos1,pos2-pos1)));
            pos1 = pos2 + 1;
        }

    }
    else
    {
        if (abrt) throw std::invalid_argument(string("argument \"")+name+string("\" not found"));
        else cout<<"argument \""<<name<<"\" not found"<<endl;
    } // if not found, warn and do nothing to val
    return found;
}

#endif // PARSER_H
