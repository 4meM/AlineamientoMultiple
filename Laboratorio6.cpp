#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <fstream>
#include <iterator>
#include <chrono>
#include "Needleman_Wunsch.h"

using namespace std;
using Mat = vector<vector<float>>;
using TreeGuia = vector<string*>;



struct cluster{
    vector<string*> chains;
    vector<pair<cluster*,cluster*>> orden;
};

using clusterHash = unordered_map<cluster*,map<cluster*, float>>;
clusterHash clustersMap;




void distanceClusters(vector<string>& chains){
    vector<cluster*> clusters; 
    for(int i = 0; i < chains.size(); i++){
        cluster* newCluster = new cluster();
        newCluster->orden.push_back({nullptr, nullptr});
        newCluster->chains.push_back(&chains[i]);
        clusters.push_back(newCluster);
    }

    for(int i = 0; i < clusters.size(); i++){
        for(int j = i + 1; j < clusters.size(); j++){
            Mat memorization(clusters[i]->chains[0]->size() + 1, vector<float>(clusters[j]->chains[0]->size() + 1,0));
            initialization(memorization, *clusters[i]->chains[0], *clusters[j]->chains[0]);
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
    for(const auto& order : c1->orden){
        if(order.first != nullptr && order.second != nullptr){
            newCluster->orden.push_back({order.first, order.second});
        }
    }
    for(const auto& order : c2->orden){
        if(order.first != nullptr && order.second != nullptr){
            newCluster->orden.push_back({order.first, order.second});
        }
    }
    newCluster->orden.push_back({c1, c2});

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
    }
}


void AlineamientoMultiple() {
    for (int i = 0; i < clustersMap.begin()->first->orden.size(); i++) {

        float distance_min = numeric_limits<float>::max();
        pair<string, string> alignedChains_aux;
        pair<cluster*, cluster*> order = clustersMap.begin()->first->orden[i];

        for (int j = 0; j < order.first->chains.size(); j++) {

            for (int k = 0; k < order.second->chains.size(); k++) {

                Mat memorization(
                    order.first->chains[j]->size() + 1,
                    vector<float>(order.second->chains[k]->size() + 1, 0)
                );

                initialization(memorization, *order.first->chains[j], *order.second->chains[k]);

                Needleman_Wunsch(
                    memorization,
                    *order.first->chains[j],
                    *order.second->chains[k]
                );

                pair<string, string> alignedChains =
                    trackeBackNeedleman_Wunsch(
                        memorization,
                        *order.first->chains[j],
                        *order.second->chains[k]
                    );

                if (distance_min >
                    memorization[memorization.size() - 1][memorization[0].size() - 1]) {

                    distance_min =
                        memorization[memorization.size() - 1][memorization[0].size() - 1];

                    alignedChains_aux = alignedChains;
                }
            }
        }
        cout << "\n";
        cout << "Cadena alineada 1: " << "\n";
        cout << alignedChains_aux.first << "\n";
        cout << "\n";
        cout << "Cadena alineada 2: " << "\n"; 
        cout << alignedChains_aux.second << "\n";
        cout << "\n";
        cout << "\n";
        for (int j = 0; j < order.first->chains.size(); j++) {

            int pointer = 0;

            for (int k = 0; k < alignedChains_aux.first.size(); k++) {

                if (alignedChains_aux.first[k] == '-') {

                    if (pointer <= order.first->chains[j]->size()) {
                        order.first->chains[j]->insert(
                            order.first->chains[j]->begin() + pointer,
                            'X'
                        );
                    } else {
                        order.first->chains[j]->push_back('X');
                    }
                }

                pointer++;
            }
        }

        // Insertar gaps en las cadenas del segundo cluster
        for (int j = 0; j < order.second->chains.size(); j++) {

            int pointer = 0;

            for (int k = 0; k < alignedChains_aux.second.size(); k++) {

                if (alignedChains_aux.second[k] == '-') {

                    if (pointer <= order.second->chains[j]->size()) {
                        order.second->chains[j]->insert(
                            order.second->chains[j]->begin() + pointer,
                            'X'
                        );
                    } else {
                        order.second->chains[j]->push_back('X');
                    }
                }

                pointer++;
            }
        }
    }
}

int main(){
    vector<string> chains;
    chains.push_back("ATTTACGCCT");
    chains.push_back("TTAAGCCAT");
    chains.push_back("TTAATTAACC");
    chains.push_back("ATTTTCCGGA");
    chains.push_back("AATTTACCGCCT");
    TreeGuide(chains); 
    AlineamientoMultiple();
    for(const auto& entry : clustersMap.begin()->first->chains){
        cout << "Cadena alineada: " << *entry << "\n";
    } 
    return 0;
}