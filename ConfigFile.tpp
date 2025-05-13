#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include "ConfigFile.h"

ConfigFile::ConfigFile(const std::string& fileName)
  : verbosity(1)
{
    std::ifstream file;
    file.open(fileName.c_str());
    if (!file) {
        std::cerr << "[ConfigFile] Impossible d'ouvrir le fichier " << fileName
                  << std::endl;
    } else {
        std::string lineread;
        while (std::getline(file, lineread)) {
            process(lineread);
        }
        file.close();
    }
}

ConfigFile::~ConfigFile() {}

void
ConfigFile::printOut(std::string path) const
{
    std::ofstream outputFile(path.c_str());
    if (outputFile.is_open()) {
        outputFile << toString() << std::endl;
    }
    outputFile.close();
}

std::string
ConfigFile::toString() const
{
    std::string strToReturn;

    for (std::map<std::string, std::string>::const_iterator iter =
           configMap.begin();
         iter != configMap.end();
         ++iter) {
        strToReturn.append(iter->first);
        strToReturn.append("=");
        strToReturn.append(iter->second);
        strToReturn.append("\n");
    }
    return strToReturn;
}

void
ConfigFile::process(const std::string& lineread)
{
    size_t commentPosition = trim(lineread).find('%', 0);
    if (commentPosition != 0 && trim(lineread).length() > 0) {
        size_t equalPosition = lineread.find('=', 1);
        if (equalPosition == std::string::npos) {
            std::cerr << "Ligne sans '=' : " << lineread << std::endl;
        } else {
            std::string key = trim(lineread.substr(0, equalPosition));
            std::string value =
              trim(lineread.substr(equalPosition + 1, lineread.length()));
            std::map<std::string, std::string>::const_iterator val =
              configMap.find(key);
            if (val != configMap.end()) {
                configMap.erase(key);
            }
            configMap.insert(std::pair<std::string, std::string>(key, value));
        }
    }
}

template<typename T>
T
ConfigFile::get(const std::string& key) const
{
    std::map<std::string, std::string>::const_iterator val =
      configMap.find(key);
    T out = T();
    if (val != configMap.end()) {
        std::istringstream iss(val->second);
        iss >> out;
        if (this->verbosity > 0)
            std::cout << "\t" << key << "=" << out << std::endl;
    } else {
        std::cerr << "[ConfigFile] Le parametre suivant est manquant : " << key
                  << std::endl;
    }
    return out;
}

template<>
bool
ConfigFile::get<bool>(const std::string& key) const
{
    std::istringstream iss(configMap.find(key)->second);
    bool result(false);
    iss >> result;
    if (iss.fail()) {
        iss.clear();
        iss >> std::boolalpha >> result;
    }
    if (this->verbosity > 0)
        std::cout << "\t" << key << "=" << result << std::endl;
    return result;
}

void
ConfigFile::setVerbosity(int level)
{
    this->verbosity = level;
}

std::string
ConfigFile::trim(std::string str)
{ // Remove tabs and spaces at the beginning and end of a string
    size_t first = str.find_first_not_of(" \t");
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(" \t");
    return str.substr(first, (last - first + 1));
}
