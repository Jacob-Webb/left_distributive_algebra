#ifndef EMB_H
#define EMB_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <unistd.h>
#include <vector>

//these are like the "nodes" in our graph, or the terms in our algebra
class EMB {
   public:

    EMB ();
    ~EMB();

    void SetStr(char* nstr);
    int CheckInRng(int ind);
    int CheckInRngFast(int ind);
    int CheckInRngFast(int ind, int start, int end);
    int CheckNotInRngFast(int ind);
    int CheckNotInRngFast(int ind, int start, int end);
    void SortRngPI();
    void BSortRng (int start, int end);
    void BMergeSortRng(int s1, int e1, int e2);
    void BSortARng (int start, int end);
    void BMergeSortARng(int s1, int e1, int e2);
    void SortRng();
    int GetRngPI(int ind);
    int AddRng (int ind, int piind, int s1, int s2);
    void SetRngPI (int rind, int pi);
    int AddNotRng (int ind);
    int AddNotRngSeq (int* inds, int num);
    void SetSqr (int sqrind);
    void SetInRange(int index, int value);
    void SetInRngPI(int index, int value);
    void SetNotInRange(int index, int value);
    void SetMyInd(int value);
    void SetFirstPos(int position);
    void SetScore(int value);
    void SetYScore(float value);
    void SetCreationScore(int scr1, int scr2); //pass the scores of nodes this is created from
    void SetRngChanged(bool TorF);
    void SetSortedPI(bool TorF);

    std::vector<int> GetRngVec();
    std::vector<int> GetRngPIVec();
    std::vector<int> GetNotRngVec();
    char* GetStr();
    int GetNumInRange();
    int GetNotInRange();
    int GetSqr();
    int GetFirstPos();
    int GetScore();
    float GetYScore();
    int GetCreationScore();
    int GetMyInd();
    bool GetRngChanged();

    void WriteToFile(FILE* fi);
    void LoadFromFile(FILE* fi);

protected:
    void Clear ();

private:
    std::vector<int> rng; //a list of nodes which are in the range of this nodes, i.e. this node is a left-divisor
    std::vector<int> rngpi; //the pullback of the corresponding nodes in the range, so thisnode*rngpi[i] = rng[i]
    std::vector<int> altrng; // range sorted by index --used for quickly cheking if something is in range
    std::vector<int> notrng; //the nodes which are not in the range of this node
    std::vector<int> altnotrng; //the nodes which are not in the range of this node, sorted
    std::vector<int> src1; //not used
    std::vector<int> src2; // not used

    char* str; // string of actual term

    int sqr; // the square of this node
    int myind; // what is my index
    float yscore; // the y value in the ordering, based on when the node was created
    int score; // not used, I think
    int firstpos; //not used
    int creationscore; // has the score generated by how the node was generated

    bool rngchanged; //if the list of things in my range has changed
    bool sortedpi; //whether range has been sorted by the preimages (rngpi)
    bool altrngsrt; //not used

    enum results { unsure = -1, no = 0, yes = 1 };
};

#endif // EMB_H
