#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iterator>
#include <chrono>

using namespace std;
using Mat = vector<vector<float>>;
using TreeGuia = vector<string*>;
#define missMatch 1
#define match 0
#define gap 1


struct cluster{
    vector<string*> chains;
};

using clusterHash = unordered_map<cluster*,map<cluster*, float>>;
clusterHash clustersMap;


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

/*Mat fillDistance(vector<string> chains){
    Mat distances(chains.size(), vector<float>(chains.size(), 0));
    for(int i = 0; i < chains.size(); i++){
        string chain1 = chains[i];
        for(int j = i + 1; j < chains.size(); j++){
            string chain2 = chains[j];
            Mat memorization(chain1.size() + 1, vector<float>(chain2.size() + 1,0));
            initialization(memorization);
            Needleman_Wunsch(memorization,chain1,chain2);
            distances[i][j] = memorization[memorization.size() - 1][memorization[0].size() - 1];
            distances[j][i] = distances[i][j];
        }
    }
    return distances;
}

pair<int,int> findMinDistance(Mat distances){
    int min_i = 0, min_j = 0;
    float min_distance = numeric_limits<float>::max();

    for(int i = 0; i < distances.size(); i++){
        for(int j = i + 1; j < distances.size(); j++){
            if(distances[i][j] < min_distance){
                min_distance = distances[i][j];
                min_i = i;
                min_j = j;
            }
        }
    }
    return {min_i, min_j};
}*/
void distanceClusters(vector<string>& chains){
    vector<cluster*> clusters; 
    for(int i = 0; i < chains.size(); i++){
        cluster* newCluster = new cluster();
        newCluster->chains.push_back(&chains[i]);
        clusters.push_back(newCluster);
    }

    for(int i = 0; i < clusters.size(); i++){
        for(int j = i + 1; j < clusters.size(); j++){
            Mat memorization(clusters[i]->chains[0]->size() + 1, vector<float>(clusters[j]->chains[0]->size() + 1,0));
            initialization(memorization);
            Needleman_Wunsch(memorization,*clusters[i]->chains[0],*clusters[j]->chains[0]);
            float distance = memorization[memorization.size() - 1][memorization[0].size() - 1];
            clustersMap[clusters[i]][clusters[j]] = distance;
            clustersMap[clusters[j]][clusters[i]] = distance;
        }
    }
}

pair<cluster*, cluster*> findMinDistanceCluster(){
    cluster* minCluster1 = nullptr;
    cluster* minCluster2 = nullptr;
    float minDistance = numeric_limits<float>::max();

    for(const auto& entry : clustersMap){
        cluster* c1 = entry.first;
        for(const auto& neighbor : entry.second){
            cluster* c2 = neighbor.first;
            float distance = neighbor.second;
            if(distance < minDistance){
                minDistance = distance;
                minCluster1 = c1;
                minCluster2 = c2;
            }
        }
    }
    return {minCluster1, minCluster2};
}

void mergeClusters(cluster*c1, cluster* c2){
    cluster* newCluster = new cluster();
    newCluster->chains.insert(newCluster->chains.end(), c1->chains.begin(), c1->chains.end());
    newCluster->chains.insert(newCluster->chains.end(), c2->chains.begin(), c2->chains.end());
    if(clustersMap.size() == 2){
        clustersMap[newCluster] = {};
        clustersMap.erase(c1);
        clustersMap.erase(c2);
        return;
    }
    vector<cluster*> clustersToUpdate;
    for(const auto& sequences: clustersMap[c1]){
        if(sequences.first != c2){
            clustersToUpdate.push_back(sequences.first);
        }
    }

    for(const auto& sequences: clustersToUpdate){
        float distance = 0.f;
        size_t c1Size = c1->chains.size();
        size_t c2Size = c2->chains.size();
        float dstBeforeC1 = clustersMap[sequences][c1];
        float dstBeforeC2 = clustersMap[sequences][c2];
        distance = (c1Size * dstBeforeC1 + c2Size * dstBeforeC2) / (c1Size + c2Size);
        clustersMap[sequences].erase(c1);
        clustersMap[sequences].erase(c2);
        clustersMap[sequences][newCluster] = distance;
        clustersMap[newCluster][sequences] = distance;
    }
    clustersMap.erase(c1);
    clustersMap.erase(c2); 

}


void TreeGuide(vector<string>& chains){
    distanceClusters(chains);
    while(clustersMap.size() > 1){
        pair<cluster*, cluster*> minClusters = findMinDistanceCluster();
        mergeClusters(minClusters.first, minClusters.second);
        cout << "Clusters fusionados. Clusters restantes: " << clustersMap.size() << "\n";
    }
}

int main(){
    vector<string> chains;
    chains.push_back("ATTGGCACCA");
    chains.push_back("ATTTGGACCA");
    chains.push_back("TGGTTCCA");
    chains.push_back("ATTCCACCAC");
    TreeGuide(chains); 
    cout << "Proceso de clustering completado.\n";
    cout << "Numero de clusters finales: " << clustersMap.size() << "\n";
    for(const auto& entry : clustersMap){
        cout << "Cluster con " << entry.first->chains.size() << " cadenas.\n";
        for(const auto& chainPtr : entry.first->chains){
            cout << " - " << *chainPtr << "\n";
        }
    }
    return 0;
}