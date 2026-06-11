
#ifndef NEEDLEMAN_WUNSCH_H
#define NEEDLEMAN_WUNSCH_H
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iterator>
#include <chrono>
#define missMatch 1
#define match 0
#define gap 1
using namespace std;
using Mat = vector<vector<float>>;

void readFastaFile(string archivo, map<string,string>& almacen){
    fstream fastaFile(archivo); 
    bool id_get = false;
    if(!fastaFile.is_open()){
        cout << "NO SE PUEDO ABRIR EL ARCHIVO " << archivo << "\n";
        return;
    }
    string line,ID;
    while(getline(fastaFile,line).good()){
        id_get = false;
        if(line[0] == '>'){
            ID = "";
            for(int i = 1; i < line.size(); i++){
                if(line[i] == ' ' ){
                    id_get = true;
                    break;
                }
                ID += line[i];
            }
        }
        if(!id_get) almacen[ID] += line;
    }
}

void Needleman_Wunsch(Mat& memorization, string& chain1, string& chain2){
    float gap_ = gap;
    for(int i = 1; i < memorization.size(); i++){
        for(int j = 1; j < memorization[0].size(); j++){
            float maximo = 0;
            maximo = memorization[i - 1][j - 1] + (chain1[i - 1] == chain2[j - 1]? match : missMatch);
            maximo = min(memorization[i - 1][j] + gap_, maximo);
            maximo = min(memorization[i][j - 1] + gap_, maximo);
            memorization[i][j] = maximo;
        }
    } 
}

bool equal(float a, float b){
    return abs(a - b) < 0.0001f;
}

pair<string,string> trackeBackNeedleman_Wunsch(Mat& memorization,string chain1,string chain2){
    pair<string, string> ans;
    int i = memorization.size() - 1;
    int j = memorization[0].size() - 1;
    while(i > 0 && j > 0){
            float missOrMatch = (chain1[i - 1] == chain2[j - 1] ? match : missMatch);
            
            // Diagonal
            if(equal(memorization[i][j], memorization[i - 1][j - 1] + missOrMatch)){
                ans.first = chain1[i - 1] + ans.first;
                ans.second = chain2[j - 1] + ans.second;
                i--;
                j--;
            }
            // Condición Arriba (Gap en chain 2)
            else if(equal(memorization[i][j], memorization[i - 1][j] + gap)){ 
                ans.first = chain1[i - 1] + ans.first;
                ans.second = "-" + ans.second;
                i--;
            }
            // Condición Izquierda (Gap en chain 1)
            else if(equal(memorization[i][j], memorization[i][j - 1] + gap)){ 
                ans.first = "-" + ans.first;
                ans.second = chain2[j - 1] + ans.second; 
                j--;

            }
            else{
                cout << "NO ENTRO EN NINGUN CASO\n";
                return ans;
            }
    }
    while(i > 0){
        ans.first  = chain1[i-1] + ans.first;
        ans.second = "-" + ans.second;
        i--;
    }

    while(j > 0){
        ans.first  = "-" + ans.first;
        ans.second = chain2[j-1] + ans.second;
        j--;
    }
   return ans;
}

void initialization(Mat& memorization){
    for(int i = 0; i < memorization.size(); i++){
        memorization[i][0] = i * gap;
    }
    for(int j = 0; j < memorization[0].size(); j++){
        memorization[0][j] = j * gap;
    }
}
#endif