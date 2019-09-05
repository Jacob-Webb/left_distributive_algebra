#include "emb.h"
#include <math.h>

EMB::EMB ()
{
    sqr = -1;
    myind = -1;
    firstpos = 0;
    str = NULL;
    rngchanged = true;
    sortedpi=false;
    altrngsrt=false;
    creationscore = 1;
}
EMB::~EMB()
{
    Clear();
}

void EMB::SetStr(char* nstr)
{
    //str = nstr;

    if (str)
    {
        if (strlen(str) <= strlen(nstr)) return;
        delete[] str;
    }
    str = new char[strlen(nstr)+1];
    strcpy(str, nstr);

}

int EMB::CheckInRng(int ind)
{
    int inRange = rng.size();
    int notInRange = notrng.size();
    for (int i=0; i < inRange; i++)
    {
        if (rng.at(i) == ind) return yes; //yes in range
    }
    for (int i=0; i<notInRange; i++)
    {
        if (notrng.at(i) == ind) return no; //definitely not in range
    }
    return unsure; //not sure
}

int EMB::CheckInRngFast(int ind)
{
    int inRange = GetNumInRange();
    if (ind == -1) return unsure;
    //if (altrngsrt == 0) SortRngPI();


    int result = CheckInRngFast(ind, 0, inRange-1);
    return result;

    for (int i=0; i<inRange; i++)
    {
        if (rng.at(i) == ind) return yes;
    }
    return unsure;
}

int EMB::CheckInRngFast(int ind, int start, int end)
{
    if (ind == -1 || start > end) return unsure;

    if (start == end){
        if (ind == altrng.at(start)) return yes;
        return unsure;
    }
    int mid = (end-start)/2 +start;
    if (altrng[mid] == ind) {return yes;}
    else if (altrng[mid] > ind){ return CheckInRngFast(ind, start, mid-1);}
    else {return CheckInRngFast(ind, mid+1, end); }
}

int EMB::CheckNotInRngFast(int ind)
{
    int notinRange = GetNotInRange();
    if (ind == -1) return unsure;
    //if (altrngsrt == 0) SortRngPI();


    int result = CheckNotInRngFast(ind, 0, notinRange-1);
    return result;
}

int EMB::CheckNotInRngFast(int ind, int start, int end)
{
    if (ind == -1 || start > end) return unsure;

    if (start == end){
        if (ind == altnotrng.at(start)) return yes;
        return unsure;
    }
    int mid = (end-start)/2 +start;
    if (altnotrng[mid] == ind) {return yes;}
    else if (altnotrng[mid] > ind){ return CheckNotInRngFast(ind, start, mid-1);}
    else {return CheckNotInRngFast(ind, mid+1, end); }
}

void EMB::SortRngPI()
{
    if (sortedpi == true) return;

    sortedpi = true;
    BSortRng(0, rng.size()-1);
   /* for (int i=0; i < nrng-1; i++)
    {
        if (rngpi[i] > rngpi[i+1])
        {
            int t = rng[i];
            rng[i] = rng[i+1];
            rng[i+1] = t;

            t = rngpi[i];
            rngpi[i] = rngpi[i+1];
            rngpi[i+1] = t;
            if (i > 0) i = i-2;
            else i--;
        }
    }*/
    altrngsrt = true;
    for (int i=0; i < GetNumInRange()-1; i++)
    {
        altrng[i] = rng.at(i);
    }
    BSortARng(0, GetNumInRange()-1);
    /*for (int i=0; i < nrng-1; i++)
    {
        if (altrng[i] > altrng[i+1])
        {
            int t = altrng[i];
            altrng[i] = altrng[i+1];
            altrng[i+1] = t;

            if (i > 0) i = i-2;
            else i--;
        }
    }*/
}

void EMB::BSortRng (int start, int end){
    if (start >= end) return;
    else if (start-end == 1)
    {
        if (rngpi[end] < rngpi[start])
        {
            int temp = rng.at(end);
            rng[end] = rng.at(start);
            rng[start] = temp;

            temp = rngpi[end];
            rngpi[end] = rngpi[start];
            rngpi[start] = temp;
        }
    }
    else {
        BSortRng(start, (start+end)/2);
        BSortRng((start+end)/2+1, end);
        BMergeSortRng(start, (start+end)/2, end);
    }
}

void EMB::BMergeSortRng(int s1, int e1, int e2){
    int i1=s1;
    int i2=e1+1;
    int* scratch1 = new int[e2-s1+1];
    int* scratch2 = new int[e2-s1+1];

    for (int i=0; i<e2-s1+1; i++)
    {
        if (i2 > e2 || ( i1 <= e1 && rngpi[i1] < rngpi[i2])) {
            scratch1[i] = rng.at(i1);
            scratch2[i] = rngpi[i1];
            i1++;
        }
        else {
            scratch1[i] = rng.at(i2);
            scratch2[i] = rngpi[i2];
            i2++;
        }
    }
    for ( int i=0; i<e2-s1+1; i++)
    {
        rng[s1+i] = scratch1[i];
        rngpi[s1+i] = scratch2[i];
    }
    delete[] scratch1;
    delete[] scratch2;
}

void EMB::BSortARng (int start, int end){
    if (start >= end) return;
    else if (start-end == 1)
    {
        if (altrng[end] < altrng[start])
        {
            int temp = altrng[end];
            altrng[end] = altrng[start];
            altrng[start] = temp;
        }
    }
    else {
        BSortARng(start, (start+end)/2);
        BSortARng((start+end)/2+1, end);
        BMergeSortARng(start, (start+end)/2, end);
    }
}

void EMB::BMergeSortARng(int s1, int e1, int e2){
    int i1=s1;
    int i2=e1+1;
    int* scratch = new int[e2-s1+1];

    for (int i=0; i<e2-s1+1; i++)
    {
        if (i2 > e2 || ( i1 <= e1 && altrng[i1] < altrng[i2])) {
            scratch[i] = altrng[i1];
            i1++;
        }
        else {
            scratch[i] = altrng[i2];
            i2++;
        }
    }
    for ( int i=0; i<e2-s1+1; i++)
    {
        altrng[s1+i] = scratch[i];
    }
    delete[] scratch;
}

void EMB::SortRng()
{
    int inRange = GetNumInRange();
    sortedpi = false;
    altrngsrt = false;
    for (int i=0; i < inRange-1; i++)
    {
        if (rng[i] > rng[i+1])
        {
            int temp = rng[i];
            rng[i] = rng[i+1];
            rng[i+1] = temp;

            temp = rngpi[i];
            rngpi[i] = rngpi[i+1];
            rngpi[i+1] = temp;
            if (i > 0) i = i-2;
            else i--;
        }
    }
    int notInRange = notrng.size();
    for (int i=0; i < notInRange-1; i++)
    {
        if (notrng[i] >notrng[i+1])
        {
            int temp = notrng[i];
            notrng[i] = notrng[i+1];
            notrng[i+1] = temp;

            if (i > 0) i = i-2;
            else i--;
        }
    }
}

int EMB::GetRngPI(int ind)
{
    if (ind == -1) return -1;
    int inRange = GetNumInRange();
    for (int i=0; i<inRange; i++)
    {
        if (rng[i] == ind) return rngpi[i];
    }
    return -1;
}

int EMB::AddRng (int ind, int piind, int s1, int s2)
{
    sortedpi = false;
    altrngsrt = false;

    std::vector<int> tempAltrng;
    std::vector<int> tempRng;
    std::vector<int> tempRngPI;

    int result = CheckInRng(ind);
    if (result == unsure)
    {
        rngchanged = yes;

        int rind = -1;
        int rind2 = -1;
        int inRange = rng.size();
        tempAltrng.reserve(inRange+1);
        tempRng.reserve(inRange+1);
        tempRngPI.reserve(inRange+1);
        for (int i = 0; i < inRange; ++i)
        {
            tempAltrng.push_back(altrng.at(i));
            tempRng.push_back(rng.at(i));
            tempRngPI.push_back(rngpi.at(i));
        }
        tempAltrng.push_back(0);
        tempRng.push_back(0);
        tempRngPI.push_back(0);
        for (int i=0; i<inRange+1; i++)
        {
            if (i < inRange && altrng[i] < ind) continue;
            else if (rind == -1)
            {
                rind = i;
                tempAltrng[rind] = ind;
                continue;
            }
            else { tempAltrng[i] = altrng.at(i - 1); }
        }
        for (int i=0; i<inRange+1; i++)
        {
            if (i < inRange && rng[i] < ind) continue;
            else if (rind2 == -1)
            {
                rind2 = i;
                tempRng[rind2] = ind;
                tempRngPI[rind2] = piind;
                continue;
            }
            else { tempRng[i] = rng.at(i - 1); tempRngPI[i] = rngpi.at(i-1); }
        }

        //src1.push_back(s1);
        //src2.push_back(s2);

        rng = tempRng;
        rngpi = tempRngPI;
        altrng = tempAltrng;
        //if (rind == -1) altrng.push_back(ind);

        return rng.size();
    }
    else if (ind != -1 && result == no)
    {
        printf("BAD ERROR NOT IN RANGE AND IN RANGE: %d %d\n", myind, ind);
        return -1;
    }
    else if (piind != -1)
    {
        SetRngPI(ind, piind);
        return rng.size();
    }
    return rng.size();
}

void EMB::SetRngPI (int rind, int pi)
{
    int inRange = rng.size();
    sortedpi = false;
    altrngsrt = false;
    for (int i=0; i < inRange; i++)
    {
        if (rng[i] == rind)
        {
            if (rngpi[i]!= -1 && rngpi[i] != pi){
                printf("tried to change preimage of %d from %d to %d\n", rind, rngpi[i], pi); return;}
            if (rngpi[i] != pi) rngchanged = true;
            rngpi[i] = pi;
        }
    }
}

int EMB::AddNotRng (int ind)
{
    if (myind == 519 && ind == 521)
    {
//        printf("adding 519 notinrange 521\n");
    }
    int notInRange = notrng.size();
    if (CheckInRng(ind) != unsure) return notInRange;
    rngchanged = true;


    std::vector<int> tempAltrng;

    int rind = -1;
    int notinRange = notrng.size();
    tempAltrng.reserve(notinRange+1);
    for (int i = 0; i < notinRange; ++i)
    {
        tempAltrng.push_back(altnotrng.at(i));
    }
    tempAltrng.push_back(0);
    for (int i=0; i<notinRange+1; i++)
    {
        if (i < notinRange && altnotrng[i] < ind) continue;
        else if (rind == -1)
        {
            rind = i;
            tempAltrng[rind] = ind;
            continue;
        }
        else { tempAltrng[i] = altnotrng.at(i - 1); }
    }

    notrng.push_back(ind);
    altnotrng = tempAltrng;

    return notrng.size();
}

int EMB::AddNotRngSeq (int* inds, int num)
{
    int notInRange = notrng.size();
    if (num == 0) return notInRange;
    rngchanged = true;

    for (int i=0; i<num; i++)
    {        
        AddNotRng(inds[i]);
    }   
    return notrng.size();
}

void EMB::SetSqr (int sqrind)
{
    if (sqr != -1) return;
    rngchanged = true;
    sqr = sqrind;
    AddRng (sqrind, myind, -1, -1);
}

void EMB::SetInRange(int index, int value) { rng[index] = value; }
void EMB::SetInRngPI(int index, int value) { rngpi[index] = value; }
void EMB::SetNotInRange(int index, int value) { notrng[index] = value; }
void EMB::SetMyInd(int value) { myind = value; }
void EMB::SetFirstPos(int position) { firstpos = position; }
void EMB::SetScore(int value) { score = value; }
void EMB::SetYScore(float value) { yscore = value; }
void EMB::SetRngChanged(bool TorF) { rngchanged = TorF; }
void EMB::SetSortedPI(bool TorF) { sortedpi = TorF; }
void EMB::SetCreationScore(int scr1, int scr2) { creationscore = scr1+scr2+1; }

std::vector<int> EMB::GetRngVec() { return rng; }
std::vector<int> EMB::GetRngPIVec() { return rngpi; }
std::vector<int> EMB::GetNotRngVec() { return notrng; }
char* EMB::GetStr() { return str; }
int EMB::GetNumInRange() { return rng.size(); }
int EMB::GetNotInRange() { return notrng.size(); }
int EMB::GetSqr() { return sqr; }
int EMB::GetFirstPos() { return firstpos; }
int EMB::GetScore() { return score; }
float EMB::GetYScore() { return yscore; }
int EMB::GetCreationScore() { return creationscore; }
int EMB::GetMyInd() { return myind; }
bool EMB::GetRngChanged() { return rngchanged; }

void EMB::WriteToFile(FILE* fi){

    int inRange = rng.size();
    int notInRange = notrng.size();
    fprintf(fi, "%d ", sqr);

    fprintf(fi, "%d ", inRange);
    for (int i=0; i<inRange; i++)
    {
        fprintf(fi, "%d %d ", rng[i], rngpi[i]);
    }

    fprintf(fi, "%d ", notInRange);
    for (int i=0; i<notInRange; i++)
    {
        fprintf(fi, "%d ", notrng[i]);
    }
    fprintf(fi, "\n");

}

void EMB::LoadFromFile(FILE* fi){
    Clear();

    int inRange = rng.size();
    int notInRange = notrng.size();
    fscanf(fi, "%d ", &sqr);

    fscanf(fi, "%d ", &inRange);

    for (int i=0; i<inRange; i++)
    {
        fscanf(fi, "%d %d ", &(rng[i]), &(rngpi[i]));
    }

    fscanf(fi, "%d ", &notInRange);

    for (int i=0; i<notInRange; i++)
    {
        fscanf(fi, "%d ", &(notrng[i]));
    }
    SortRngPI();
}

void EMB::Clear ()
{
    sqr = -1; myind = -1;
    rngchanged = true;
    if (str) delete[] str;
    str = NULL;
    sortedpi=false;
    altrngsrt=false;
}
