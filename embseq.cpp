//
//  EMBSEQ.cpp
//  Graphs LDA
//
//  Created by Scott Cramer on 8/16/16.
//  Copyright © 2016 Scott Cramer. All rights reserved.
//

#include "embseq.h"
#include <math.h>

int outputProgress(int cur, int tot, int opt)
{
    if (cur == 0) return 0;
    if (((int)(100*((float)cur)/(float)tot))/opt > ((int)(100*((float)cur-1)/(float)tot))/opt)
    { printf("."); fflush(stdout); }
    return 0;
}

EMBSEQ::EMBSEQ (int max)
{
len = 0;
invcomp = 1;
MAXEMBS = max;
embs = new EMB[MAXEMBS];
eord = new int[MAXEMBS];
eordi = new int[MAXEMBS];
forceq1 = new int[MAXEMBS];
forceq2 = new int[MAXEMBS];
fscore = new int[MAXEMBS];
forceq1[0] = -1;
forceq2[0] = -1;
numiterations = 0;
ILnum = 1;
}
EMBSEQ::~EMBSEQ()
{
delete[] embs;
delete[] eord;
delete[] eordi;
delete[] forceq1;
delete[] forceq2;
delete[] fscore;
}

EMB* EMBSEQ::GetEmbs() { return embs; }
int* EMBSEQ::GetEord() { return eord; }
int* EMBSEQ::GetEordi() { return eordi; }
int EMBSEQ::GetLength() {return  len; }

//add ind * tind to the list of nodes we want to add
void EMBSEQ::AddToForceQ(int ind, int tind, int score)
{
if (FindImage(ind, tind) != -1 || WillEverCompute(ind, tind) == -1) return; // if we already can compute ind*tind or if we will never be able to (because the answer is too big), then we won't add it to the list
for (int i=0; i < MAXEMBS; i++) //add it to the end of the list, if it's not already there.
{
    if (forceq1[i] == -1)
    {
        forceq1[i] = ind;
        forceq2[i] = tind;
        fscore[i] = score;
//            printf("added: %d of %d\n", ind, tind);
        if (i < MAXEMBS-1) { forceq1[i+1] = -1; forceq2[i+1] = -1; nforce++; }
        return;
    }
    if (forceq1[i] == ind && forceq2[i] == tind) return;
}
}

int EMBSEQ::InsertEmb(int ind, int num)
{
if (num == 0 || num + len > MAXEMBS) return 0;

float lower = 0;
float upper = 200;

//printf("bounds %f %f \n", lower, upper);

if (ind != 0) lower = embs[eord[ind-1]].GetYScore();
if (ind != len) upper = embs[eord[ind]].GetYScore();

if (ind == 0 && len > 0) lower = embs[eord[0]].GetYScore()-20;

//printf("bounds %f %f \n", lower, upper);

int i=len+num-1;
for (i=len+num-1; i>=ind+num; i--)
{
    eord[i] = eord[i-num];
}
for (; i>= ind; i--)
{
    eord[i] = len + i-ind;
}
for (i=len; i<len+num; i++)
{
    embs[i].SetMyInd(i);

}
len = len+num;
invcomp = 0;

for (int i=0; i < num; i++)
{
    if (num > 1) embs[len-num+i].SetYScore(lower + (upper-lower)*log(((float)i+2)*2/((float)num+1) + 1) );
    else embs[len-num+i].SetYScore(lower + (upper-lower)*.75);
    //embs[len-num+i].SetYScore(lower + (upper-lower)*(exp(((float)i)/((float)num) + .1)-1)/(exp(1)-1) );
}
return num;

printf("inserted: %d-many at %d, bounds(%f,%f) with scores:", num, ind, lower, upper);
for (int i=0; i < num; i++)
{
    printf ("%f ", embs[len-num+i].GetYScore());

}
printf("\n");

return num;
}

void EMBSEQ::ComputeOrderInverse()
{
if (invcomp == 1) return;

for (int i=0; i<len; i++)
{
    eordi[i] = FindOrderIndComp(i);
}
invcomp = 1;
}

int EMBSEQ::CreateILSeq ( int num){
    ILnum = num;
InsertEmb (0, num+1);

embs[num].SetCreationScore(1,0);
for (int i=0; i<num; i++)
{
    embs[i].SetSqr(num);
    embs[i].SetCreationScore(i+2,0);
    AddToRngSeq(i, i+1, num-1);
    AddToNotRngSeq(i, 0, i);
}
return num;
}

int EMBSEQ::CreateILSeqModified ( int num, int sqrnum){ //sqrnum should be even

    CreateILSeq(num);

    InsertEmb(0,sqrnum);

    int base = num+1;
    for (int i=0; i<sqrnum; i++)
    {
        embs[base+i].SetSqr(i);
        embs[base+i].SetCreationScore(i+3,0);
        AddToRngSeq(base+i, 0, num-1);
        if(i<sqrnum/2-1) { AddToRngSeq(base+i,base+i+1, base+sqrnum/2-1); AddToNotRngSeq(base+i, base+sqrnum/2, base+sqrnum-1);}
        if (i>=sqrnum/2 && i < sqrnum-1)
        {
            AddToRngSeq(base+i,base+i+1, base+sqrnum-1);
            AddToNotRngSeq(base+i, base, base+sqrnum/2-1);
        }
        AddToNotRngSeq(base+i, base, base+i);
    }
    return num;
}

void EMBSEQ::AddToRngSeq (int ind, int start, int finish)
{
if (finish < start || ind == -1 || start == -1) return;
for (int i=start ; i<= finish; i++)
{
    embs[ind].AddRng(i, -1, -1, -1);
}
}

void EMBSEQ::AddToNotRngSeq (int ind, int start, int finish)
{
if (finish < start || ind == -1 || start == -1) return;
for (int i=start ; i<= finish; i++)
{
    embs[ind].AddNotRng(i);
}
}

void EMBSEQ::AddToNotRngSeqOrd (int ind, int start, int finish)
{
if (finish < start || ind == -1 || start == -1) return;
for (int i=start ; i<= finish; i++)
{
    embs[ind].AddNotRng(eord[i]);
}
}

int EMBSEQ::FindLeastSentAbove(int ind, int rind){
int roid = FindOrderInd(rind);
for (int i=0;  i< len; i++)
{
    int im = FindImage(ind, eord[i]);
    if (FindOrderInd(im) >= roid) return eord[i];
}
return len-1;
}

int EMBSEQ::InsertPullback(int ind, int rind)
{
if (ind == -1 || rind == -1) return 0;
if (embs[ind].GetRngPI(rind) != -1) return 0;
int lb1 = 0;
int lb2 = 0;
int ub = len-1;
lb2 = FindLeastSentAbove(ind, rind); //lb2 = ind*(rind+1), basically
//SortRng(ind);
lb1 = FindMaxInRng(ind, rind); //lb1 is the maximal term which has both ind and rind in its range
//SortRng(ind);
ub = FindMinNotInRng(ind, rind); //ub is the minimal index with either ind or rind not in range
//SortRng(ind);
//  if (FindOrderInd(lb) != FindOrderInd(ub)-1)
//    printf("lower bound and upper bound not strict.\n");

//  if (lb1 > lb2) lb1 = lb2;
if (FindOrderInd(lb1)+1 != FindOrderInd(lb2))
{ //printf("BAD PULLBACK ERROR\n");
    return -1;}

for (int i=0; i<embs[rind].GetNumInRange(); i++)
{
    if (embs[ind].CheckInRng(embs[rind].GetRngVec().at(i))==1)
    {
        int pb = embs[ind].GetRngPI(embs[rind].GetRngVec().at(i));
        if (pb == -1) {
           // printf("BAD PULLBACK ERROR (2nd kind).\n");
            return -1;
        }
        //int pb2 = embs[ind].GetRngPI(embs[rind].rngpi[i]);
        //if (pb2 == -1) {printf("don't know about pullback (not so serious).\n"); }
        //embs[len-1].AddRng(pb, pb2);
    }
}

InsertEmb(FindOrderInd(lb1)+1, 1); // new term index = len-1
embs[ind].SetRngPI(rind, len-1); // sets ind*len-1 = rind

embs[len-1].SetCreationScore(embs[ind].GetCreationScore(), embs[rind].GetCreationScore());

embs[len-1].SetSqr(embs[ind].GetRngPI(embs[rind].GetSqr())); // set square of len-1 to be pullback of square of rind
AddSqrSeqRng(len-1); // add in all squares to the range of len-1

for (int i=0; i<embs[rind].GetNumInRange(); i++) // check everyone in the range of rind
{
    if (embs[ind].CheckInRng(embs[rind].GetRngVec().at(i))==1)// check if each element of the range of rind is in the range of ind
    {
        int pb = embs[ind].GetRngPI(embs[rind].GetRngVec().at(i));// get the pullback, if it exists
        if (pb == -1) {
            printf("don't know about pullback (serious).\n");
            continue;
        }
        int pb2 = embs[ind].GetRngPI(embs[rind].GetRngPIVec().at(i));
        if (pb2 == -1) {//printf("don't know about pullback (not so serious).\n");
        }
        embs[len-1].AddRng(pb, pb2, embs[rind].GetRngVec().at(i), embs[rind].GetRngPIVec().at(i));
    }
}

for (int i=0; i < len; i++)
{
    if (embs[ind].CheckInRng(i)==1)
    {
        for (int i2=0; i2 < embs[i].GetNumInRange(); i2++)
        {
            if (embs[ind].CheckInRng(embs[i].GetRngVec().at(i2)) && embs[i].GetRngPIVec().at(i2)!=-1 && embs[ind].CheckInRng(embs[i].GetRngPIVec().at(i2))==1)
            {
                int ipi = embs[ind].GetRngPI(i);
                int r1 = embs[ind].GetRngPI(embs[i].GetRngVec().at(i2));
                int r2 = embs[ind].GetRngPI(embs[i].GetRngPIVec().at(i2));
                if (ipi != -1 && r1 != -1 && r2 != -1)
                    embs[ipi].AddRng(r1, r2, embs[i].GetRngVec().at(i2), embs[i].GetRngPIVec().at(i2));
            }
        }
    }
}

for (int i=0; i<=FindOrderInd(len-1); i++)
{
    embs[len-1].AddNotRng(eord[i]);
}


for (int i=FindOrderInd(len-1); i < len; i++)
{
    embs[eord[i]].AddNotRng(len-1);
}
for (int i=0; i<embs[rind].GetNotInRange(); i++)
{
    if (embs[ind].CheckInRng(embs[rind].GetNotRngVec().at(i))==1)
    {
        int pb = embs[ind].GetRngPI(embs[rind].GetNotRngVec().at(i));
        if (pb == -1) {continue;}
        embs[len-1].AddNotRng(pb);
    }
}

for (int i=0; i < len; i++)
{
    if (embs[i].CheckInRng(ind)==1 && embs[i].CheckInRng(rind)==1)
    {
        embs[i].AddRng(len-1, -1, -1, -1);
    }
}

//RecheckEmbs();

return 1;
}

void EMBSEQ::ReplaceEmb(int newemb, int oldemb)
{
if (newemb == -1 || oldemb == -1) return;

for (int i=0; i< len ; i++)
{

}
}

void EMBSEQ::ReplaceAllInst(int newi, int oldi){
for (int i=0; i< len; i++)
{
    for (int i2=0; i2< embs[i].GetNumInRange(); i2++)
    {
        if (embs[i].GetRngVec().at(i2) == oldi) embs[i].SetInRange(i2, newi);
        if (embs[i].GetRngPIVec().at(i2) == oldi) embs[i].SetInRngPI(i2, newi);
        //if (embs[i].src1[i2] == oldi) embs[i].src1[i2] = newi;
        //if (embs[i].src2[i2] == oldi) embs[i].src2[i2] = newi;
    }
    for (int i2=0; i2< embs[i].GetNotInRange(); i2++)
    {
        if (embs[i].GetNotRngVec().at(i2) == oldi) embs[i].SetNotInRange(i2, newi);
    }
    if (embs[i].GetSqr() == oldi) embs[i].SetSqr(newi);
    embs[i].SetRngChanged(true);
}
for (int i=0; i<len; i++)
{
    if (eord[i] == oldi) eord[i] = newi;
}
}

void EMBSEQ::RecheckEmbs()
{
/*for (int i=0; i<len; i++)
{
    if (embs[i].rngchanged == 1)
    {
        for (int i2=0; i2 < len; i2++)
        {
            if (embs[i2].firstpos > i) embs[i2].firstpos = i2;
        }
        break;
    }
}*/

/* char **cc = new char*[len];
for (int i=0; i < len; i++)
{
    cc[i] = new char[len];
    for (int i2=0; i2 < len; i2++) cc[i][i2] = 1;
}*/

    printf("0..."); fflush(stdout);

    for (int e=0; false && e<len; e++)
    {
        outputProgress(e,len,2);
        for (int i=0; i< len; i++)
        {
            if (FindOrderInd(i)>=FindOrderInd(e)) continue;
            int img = FindImage(e,i);
            if (img == -1) continue;
            for (int i2=FindOrderInd(i); i2 < FindOrderInd(e); i2++)
            {
                int i2ind = eord[i2];
                for (int i3=FindOrderInd(e)+1;i3<= FindOrderInd(img); i3++)
                {
                    int i3ind = eord[i3];
                    embs[i2ind].AddNotRng(i3ind);
                }
            }
        }
    }

    printf("1...");  fflush(stdout);

for (int e=0; e<len-1; e++)
{
    outputProgress(e,len,2);
    int nontriv = 0;
    for (int i=0; i<len-1; i++)
    {
        if (nontriv == 0 && i > embs[e].GetFirstPos()) embs[e].SetFirstPos(i);

        if (FindOrderInd(e) < FindOrderInd(i))
        {
            if (i == 519 && e == 521)
            {
                printf("519, 521: %d %d\n", FindOrderInd(e), FindOrderInd(i));
            }
            embs[i].AddNotRng(e);
        }
        int im = FindImage(e, i); //check if e*i will be computed
        if (im==-1) { // if not, we continue on to the next i
            if (WillEverCompute(e,i)==1) nontriv = 1;
            continue;
        }
        for (int i2=0; i2 < embs[i].GetNumInRange(); i2++) // for everyone in the range of i
        {
            int i2m = FindImage(e, embs[i].GetRngVec().at(i2)); // find under e of the element in the range and its pullback
            int pim = FindImage(e, embs[i].GetRngPIVec().at(i2));

            if (i2m == -1 && pim != -1) // if the image of the element in range is not computed, but the pullback is
            {
                int i2mnew = FindImage(im, pim); // find (e*i)*(e*pullback) = e*range element
                if (i2mnew != -1) embs[e].AddRng(i2mnew, embs[i].GetRngVec().at(i2), -1, -1 );
            }
            else if (pim == -1 && i2m != -1) // if the image of the element in range is computed but pullback is not
            {
                int pimnew = embs[im].GetRngPI(i2m); // find the pullback of e*range element by (e*i)
                if (pimnew != -1) embs[e].AddRng(pimnew, embs[i].GetRngPIVec().at(i2), -1,-1);
            }

           /* if (embs[i].rngpi[i2] != -1 && )
            {
                int imgpi = embs[im].GetRngPI(i2m);
                if (imgpi != -1)
            }*/
            if (nontriv == 0 && ((i2m == -1 && WillEverCompute(e, embs[i].GetRngVec().at(i2)) == 1) || (pim == -1 && WillEverCompute(e, embs[i].GetRngPIVec().at(i2))==1)))
            {
                nontriv = 1;
            }

        }
    }
}

printf("2...");  fflush(stdout);
for (int e=0; e < len; e++)
{
    outputProgress(e, len, 2);
    for (int i=0; i<len; i++)
    {
        if (embs[e].CheckInRngFast(i) != 1) continue; // make sure i is in the range of e
        for (int i2=0; i2 < embs[i].GetNumInRange(); i2++)
        {
            if (embs[e].CheckInRngFast(embs[i].GetRngVec().at(i2)) == 1 && embs[i].GetRngPIVec().at(i2) != -1)
            {//if e has an element in the range of i in its range, then it should have the pullback in its range
                embs[e].AddRng(embs[i].GetRngPIVec().at(i2), -1, -1, -1);
            }
            else if (embs[e].CheckInRngFast(embs[i].GetRngPIVec().at(i2)) == 1)
            { // on the other hand, if it has the pullback in its range, it should have the image by i in its range
                embs[e].AddRng(embs[i].GetRngVec().at(i2), -1, -1, -1);
            }
            else if (embs[e].CheckNotInRngFast(embs[i].GetRngVec().at(i2))==1)
            {
                embs[e].AddNotRng(embs[i].GetRngPIVec().at(i2));
            }
            else if (embs[e].CheckNotInRngFast(embs[i].GetRngPIVec().at(i2))==1)
            {
                embs[e].AddNotRng(embs[i].GetRngVec().at(i2));
            }

        }
    }
}

printf("3...");  fflush(stdout);

for (int e=0; e< len; e++)
{
    outputProgress(e, len, 2);
    for (int i=0; i<len; i++)
    {
        int ei = FindImage(e, i);
        if (ei == -1) continue;
        for (int i2= 0; i2<embs[i].GetNumInRange(); i2++)
        {
            int ei2 = FindImage(e, embs[i].GetRngVec().at(i2));
            if (ei2 == -1) continue;
            embs[ei].AddRng(ei2, FindImage(e, embs[i].GetRngPIVec().at(i2)), -1, -1);
        }
    }

    for (int i=0; i<len; i++)
    {
        if (embs[e].CheckInRngFast(i) != 1) continue;
        int ei = embs[e].GetRngPI(i);
        if (ei == -1) continue;
        for (int i2= 0; i2<embs[i].GetNumInRange(); i2++)
        {
            if (embs[e].CheckInRngFast(embs[i].GetRngVec().at(i2)) != -1) continue;
            int ei2 = embs[e].GetRngPI(embs[i].GetRngVec().at(i2));
            if (ei2 == -1) continue;
            if (ei == 519 && ei2 == 521)
            {
                printf("e: %d i: %d ei: %d ei2: %d, rng: %d rngi2: %d rngi2pb: %d orders: %d %d \n",
                       e, i, ei, ei2, embs[i].GetRngVec().at(i2), embs[i].GetRngPIVec().at(i2), embs[e].GetRngPI(embs[i].GetRngPIVec().at(i2)),
                       FindOrderInd(i), FindOrderInd(embs[i].GetRngVec().at(i2)));
            }
            embs[ei].AddRng(ei2, embs[e].GetRngPI(embs[i].GetRngPIVec().at(i2)), -1, -1);
        }
    }
}
}

int EMBSEQ::FindImage( int ind, int tind){
if (ind == -1 || tind == -1 || embs[ind].GetNumInRange() == 0) return -1;
embs[ind].SortRngPI();

int res =  FindImage(ind, tind, 0, embs[ind].GetNumInRange()-1);

return res;
/*
for (int i=0; i<embs[ind].nrng; i++)
{
    if (embs[ind].rngpi[i]==tind) {
        return embs[ind].rng[i];}
}
return -1;*/
}
int EMBSEQ::FindImage( int ind, int tind, int start, int end){
if (ind == -1 || tind == -1 || start > end) return -1;

if (start == end)
{
    if (embs[ind].GetRngPIVec().at(start) == tind) return embs[ind].GetRngVec().at(start);
    return -1;
}
int i = (end-start)/2+start;
if (embs[ind].GetRngPIVec().at(i) < tind) return FindImage(ind, tind, i+1, end);
else if (embs[ind].GetRngPIVec().at(i) > tind) return FindImage(ind, tind, start, i-1);
else return embs[ind].GetRngVec().at(i);
}
int EMBSEQ::GetSqrSeq(int ind, int* seq)
{
seq[0] = ind;
for (int i=1; i < MAXSQRS; i++)
{
    int sqr = embs[seq[i-1]].GetSqr();
    if (sqr == -1) return i;
    seq[i] = sqr;
}
return MAXSQRS;
}

int EMBSEQ::GetLeastCommonSqr(int ind, int tind){
int sqr1[MAXSQRS];
int sqr2[MAXSQRS];
int sqr1len;
int sqr2len;
sqr1len = GetSqrSeq(ind, sqr1);
sqr2len = GetSqrSeq(tind, sqr2);

for (int i=0; i < sqr1len; i++)
{
    for (int i2=0; i2 < sqr2len; i2++)
    {
        if (sqr1[i] == sqr2[i2]) return sqr1[i];
    }
}
return -1;
}

int EMBSEQ::GetEmbRank(int ind){
    int lcs = GetLeastCommonSqr(ind, eord[0]);
    if (lcs == -1) return -1;

    int sqr1[MAXSQRS];
    int sqr1len = GetSqrSeq(eord[0],sqr1);

    for (int i=0; i< sqr1len; i++)
    {
        if (sqr1[i] == lcs) return i;
    }
    return -1;
}

int EMBSEQ::WillEverCompute (int ind, int tind)
{
if (ind == -1 || tind == -1) return -1;

int res = FindImage(ind, tind);
if (res != -1) return 1;

if (embs[tind].GetSqr() == -1) return -1;

/*
if (GetEmbRank(ind) >= rank && GetEmbRank(tind) >= rank)
{
    int lcs = GetLeastCommonSqr(ind,tind);
    if (GetLeastCommonSqr(lcs, eord[0]) == lcs) return -1;
}*/

return WillEverCompute(ind, embs[tind].GetSqr());

}

int EMBSEQ::WillEverCompute (int ind, int tind, int rank)
{
if (ind == -1 || tind == -1) return -1;

int res = FindImage(ind, tind);
if (res != -1) {
    if (GetEmbRank(res) > rank) return -1;
    //printf ("foundimage\n");
    return 1;
}

if (embs[tind].GetSqr() == -1 //|| GetEmbRank(embs[tind].GetSqr()) > rank
        ) return -1;

if (GetEmbRank(ind) >= rank && GetEmbRank(tind) >= rank)
{
    int lcs = GetLeastCommonSqr(ind,tind);
    if (GetLeastCommonSqr(lcs, eord[0]) == lcs) return -1;
}

//printf("going to square %d %d\n", ind, tind);
return WillEverCompute(ind, embs[tind].GetSqr(), rank);

}

int EMBSEQ::ComputePushForward(int ind, int tind, int iters)
{
if (ind == -1 || tind == -1 || iters < 0) return -1;

int res = FindImage(ind, tind);
if (res != -1) return res;

int lcs = GetLeastCommonSqr(ind, tind);

if (lcs != tind && FindImage(ind, embs[tind].GetSqr())== -1) return -1;

if (lcs == tind && embs[tind].GetSqr() == -1) return -1;

// if (FindOrderInd(ind) > FindOrderInd(tind))
if (FindOrderInd(ind) > FindOrderInd(tind))
{
   // printf("to compute %d of %d, need to add image\n", ind, tind);
    return -1;
}

int r1 = -1;
int r2 = -1;
for (int i=0; i < len; i++)
{
    if (embs[eord[i]].CheckInRngFast(tind) == 1)
    {
        r2 = embs[eord[i]].GetRngPI(tind);
        if (r2 != -1) {
            r1 = eord[i];

            int pr1 = ComputePushForward(ind, r1, iters-1);
            int pr2 = ComputePushForward(ind, r2, iters-1);
            int res = ComputePushForward(pr1, pr2, iters-1);
            if (res != -1){

                embs[ind].AddRng(res, tind, pr1, pr2);
                return res;
            }

            r1 = r2 = -1;
        }
    }
}
return -1;
}

// this function tries to add the node ind*tind. the score is a measure of how far removed from the original node we wanted.
int EMBSEQ::ForceCompPushForward(int ind, int tind, int score)
{
if (ind == -1 || tind == -1) return -1;
//if (rand()%1000 > 50) return -1;

int res = FindImage(ind, tind); //try to find ind*tind. if it already exists, we return.
if (res != -1) return res;

int lcs = GetLeastCommonSqr(ind, tind); //find the least common node above ind a tind. this means we look at the square sequence above ind and above tind and find the least common node.

if (lcs != tind && FindImage(ind, embs[tind].GetSqr())== -1) { AddToForceQ(ind, embs[tind].GetSqr(), score); return -1; } // if we don't even know what ind* (tind*tind) is, then we try to add that first.

if (lcs == tind && embs[tind].GetSqr() == -1) return -1; //if tind is the highest node we are considering, then we won't succeed (the answer will be too high)

if (FindOrderInd(ind) > FindOrderInd(tind)) //if ind > tind, then we actually try to add the node
{
    return ForcePushForward(ind, tind, score);

}

int r1 = -1;
int r2 = -1;
for (int i=0; i < len; i++) //we go through our nodes and check if any is a left-divisor of tind. If so we can write tind = i*r2. Then we try computing ind*(i*r2) = (ind*i)*(ind*r2) = pr1*pr2. Otherwise, we add all these various products to our queue.
{
    if (embs[eord[i]].CheckInRngFast(tind) == 1)
    {
        r2 = embs[eord[i]].GetRngPI(tind);
        if (r2 != -1) {
            r1 = eord[i];

            int pr1 = FindImage(ind, r1);
            int pr2 = FindImage(ind, r2);
            int res = FindImage(pr1, pr2);
            if (res != -1){

                embs[ind].AddRng(res, tind, pr1, pr2);
                return res;
            }
            if (pr1 == -1 && WillEverCompute(ind, r1) == 1) AddToForceQ(ind, r1, score+1);
            if (pr2 == -1  && WillEverCompute(ind, r2) == 1) AddToForceQ(ind, r2, score +1);
            if(res == -1 && WillEverCompute(pr1, pr2) == 1) AddToForceQ(pr1, pr2, score+1);

            r1 = r2 = -1;
        }
    }
}
return -1;
}

//this is where we check whether we are ok adding the node in the right place. Basically this boils down to whether we are sure that we know the ordering relative to all the other nodes.

int EMBSEQ::ForcePushForward(int ind, int tind, int score)
{
if (ind == -1 || tind == -1) return -1;
if (FindOrderInd(ind) <= FindOrderInd(tind) && embs[ind].CheckInRng(tind) != 0) {//printf("shouldn't have to push forward.\n");
    return -1;} // we don't add the node if ind <= tind and tind might be in the range of ind.

if (FindImage(ind, tind) != -1) { //printf("already found image\n");
    return -1;} //check if we already have this node

if (embs[tind].GetSqr() == -1 || FindImage(ind, GetLeastCommonSqr(ind, tind)) == -1) { //printf("square too high \n");
    return -1;} //check if the node is going to be too high for us to consider
if (FindImage(ind, embs[tind].GetSqr()) == -1) { return -1;} //check if we have not yet added ind*(tind*tind)

int liu = FindOrderInd(FindLeastImgAbove(ind, tind)); //liu = (order index of) least node that ind divides above tind
if (liu > 0 && liu-1 > FindOrderInd(ind)) //if liu is not the least point in the order and liu is not directly above ind
{

    int inrng = embs[ind].CheckInRng(eord[liu-1]); //is the node directly below liu in range of ind?
    if (inrng != 1 || (inrng == 1 && embs[ind].GetRngPI(eord[liu-1]) == -1)) //if it might not be, or it is and I don't know what the pullback is, then I don't know where to place the new node, so I add to the queue ind*(a node just below tind). Computing this might help with figuring this out.
    { //printf("BAD PUSHFORWARD ERROR\n");
        if (FindOrderInd(tind) > 0) AddToForceQ(ind, eord[FindOrderInd(tind)-1], score+1);

        return -1;}
}

InsertEmb(liu,1); //if everything checks out (especially the last possibility above), then I can add in my node right below liu.

embs[len-1].SetCreationScore(embs[ind].GetCreationScore(),embs[tind].GetCreationScore());
embs[ind].AddRng(len-1, tind, -1, -1); //note the fact that ind*tind = new node
embs[len-1].SetSqr(FindImage(ind, embs[tind].GetSqr())); //square of new node is ind*(square of tind)
AddSqrSeqRng(len-1); //add square sequence above new node to its range

//now that we've added a new node, we go through and add as much information about the new node as possible. First we figure out what must be in the range of various nodes, and then we figure out what must not be in the range of various nodes.
for (int i=0; i<embs[tind].GetNumInRange(); i++) //we go through everything that tind left divides, and add the corresponding nodes to what the new nodes divides. If tind divides x then the new nodes divides ind*x.
{
    int img = FindImage(ind, embs[tind].GetRngVec().at(i));

    if (img != -1)
    {
        int piimg = FindImage(ind, embs[tind].GetRngPIVec().at(i));

        embs[len-1].AddRng(img, piimg, embs[tind].GetRngVec().at(i), embs[tind].GetRngPIVec().at(i));
    }
}

// if ind*i = img and i*x = y, then (ind*i)*(ind*x) = ind*y
for (int i=0; i < len; i++)
{
    int img = FindImage(ind, i);
    if (img != -1)
    {
        for (int i2=0; i2 < embs[i].GetNumInRange(); i2++)
        {
            int r1img = FindImage(ind, embs[i].GetRngVec().at(i2));
            int r2img = FindImage(ind, embs[i].GetRngPIVec().at(i2));
            if (r1img != -1 && r2img != -1)
            {
                embs[img].AddRng(r1img, r2img, embs[i].GetRngVec().at(i2), embs[i].GetRngPIVec().at(i2));
            }
        }
    }
}

//nodes that appear below the new node in the order are not in the range of the node
embs[len-1].AddNotRngSeq(eord, FindOrderInd(len-1));

/*  for (int i=0; i<=FindOrderInd(len-1); i++)
{
    embs[len-1].AddNotRng(eord[i]);
}*/
// if y is not in the range of tind, then ind*y is not in the range of ind*tind = new node
for (int i=0; i<embs[tind].GetNotInRange(); i++)
{
    int img = FindImage(ind, embs[tind].GetNotRngVec().at(i));
    if (img != -1)
    {
        embs[len-1].AddNotRng(img);
    }
}

//if new node*x = i and i*pi = i2 and pi < new node then new node does not have i2 in its range.
//if on the other hand img = i*i2 and img <= new node then i2 is not in range of new node
for (int i=0; i< len; i++)
{
    if (embs[len-1].CheckInRngFast(i) == 1)
    {
        for (int i2=0; i2 < len; i2++)
        {
            if (embs[i].CheckInRngFast(i2) == 1)
            {
                int pi = embs[i].GetRngPI(i2);
                if (pi != -1 && FindOrderInd(pi) <= FindOrderInd(len-1))
                {
                    embs[len-1].AddNotRng(i2);
                }
            }
            int img = FindImage(i, i2);
            if (img != -1 &&
                FindOrderInd(img) <= FindOrderInd(len-1))
                embs[len-1].AddNotRng(i2);
        }
    }
}

//if ind and tind are in the range of i, then ind*tind is in the range of i

for (int i=0; i < len; i++)
{
    if (embs[i].CheckInRngFast(ind)==1 && embs[i].CheckInRngFast(tind)==1)
    {
        embs[i].AddRng(len-1, -1, -1, -1);
    }
}

// RecheckEmbs();


return len-1;
}
#define MAXATT 3
#define MAXADDFORCE 40
#define MAXNFORCE 200000

//this is where we really try to add our node -- a lot of the details of our algorithm appear here
//we try to add the node ind * tind
//hullv and hullnum and just keeping track of indices.
int EMBSEQ::ForceCompute(int ind, int tind, int* hullv, int* hullnum)
{
int numatt = 0; // number of attempts so far to add (where nothing happens -- anything changing causes this to go back to 0)
int tatt = 0; // keeps track of how long it's been since we've had a "significant breakthrough"; the higher this gets, the more of the queue we try to go through at each stage
int res = -1;
int omax = 0; //old number of nodes to be forced
forceq1[0] = -1;
forceq2[0] = -1;

int lind = ind;
int ltind = tind;
int least = -1;
int lscore = -1;

int *lochullv = new int[MAXEMBS];
int lochullnum = 0;

AddToForceQ(lind, ltind, 0); // add ind * tind to the list of nodes we want to add
while (res == -1 && numatt < MAXATT && omax < MAXNFORCE) //keep doing this as long as we haven't succeeded (and we have not gotten to the max # of attempts)
{
    tatt++; //increment the time since "significant breakthrough"

    nforce = 0; //this empties the queue of nodes to add (not sure why this is set to zero--you might not want to do this)
    for (int i=0; i< MAXEMBS; i++) //check if there were new nodes added to the queue
    {
        if (forceq1[i] == -1) { if (omax != i) numatt = 0;
            omax = i; break; }
    }

    printf("."); //a . means 1 round has gone by
    fflush(stdout);
    for ( int i=0; i < omax && nforce < MAXADDFORCE*tatt; i++)  // we go through the queue and try to see if we are allowed to add in the nodes
    {
        int fcres = ForceCompPushForward( forceq1[i], forceq2[i], fscore[i]);
        if (fcres != -1) numatt = 0;
    }

    // next, we check to see if there is some node i which is a left divisor of the two nodes forceq1[i2] and forceq2[i2] that we are trying to add. If so, we consider the pullbacks o1 and o2. So i * o1 = forceq1[i2] and i * o2 = forceq2[i2]. Then we try to consider o1*o2. If we can add that node we try to add i*(o1*o2). This is basically using an instance of left-distributivity.
    for (int i = 0; i<len && nforce < MAXADDFORCE*tatt*2; i++)
    {
        for (int i2=0; i2<omax  && nforce < MAXADDFORCE*tatt*2; i2++)
        {
            if (embs[i].CheckInRngFast(forceq1[i2]) == 1 && embs[i].CheckInRngFast(forceq2[i2]) == 1)
            {
                int o1 = embs[i].GetRngPI(forceq1[i2]);
                int o2 = embs[i].GetRngPI(forceq2[i2]);
                int fcres = ForceCompPushForward(o1, o2, fscore[i2]); // add o1*o2
                if (fcres != -1) ForceCompPushForward(i, fcres, fscore[i2]); // add i*(o1*o2)
            }
        }
    }

    // this is the same as the last step, but we start from the end of the queue (these are nodes which were recently added, so maybe it's useful to use them)
    for (int i = 0; i<len && nforce < MAXADDFORCE*tatt*3; i++)
    {
        for (int i2=omax-1; i2>=0  && nforce < MAXADDFORCE*tatt*3; i2--)
        {
            if (embs[i].CheckInRngFast(forceq1[i2]) == 1 && embs[i].CheckInRngFast(forceq2[i2]) == 1)
            {
                int o1 = embs[i].GetRngPI(forceq1[i2]);
                int o2 = embs[i].GetRngPI(forceq2[i2]);
                int fcres = ForceCompPushForward(o1, o2, fscore[i2]);
                if (fcres != -1) ForceCompPushForward(i, fcres, fscore[i2]);
            }
        }
    }

    // here we try to add the nodes at the end of the queue, similar to the last step

    for ( int i=omax-1; i >=0 && nforce < MAXADDFORCE*tatt*4; i--)
    {
        int fcres = ForceCompPushForward( forceq1[i], forceq2[i], fscore[i]);
        if (fcres != -1) numatt = 0;
    }


    //next we try looking at i*forceq1[i2] and i*forceq2[i2]. If we can compute both of these, then we try to compute (i*forceq1[i2]) * (i*forceq2[i2]). This will probably help compute forceq1[i2]*forceq2[i2], since i * (forceq1[i2] * forceq2[i2]) = (i*forceq1[i2]) * (i*forceq2[i2]).
    for (int i=0; i< len && nforce < MAXADDFORCE*tatt*5; i++)
    {
        for (int i2=0; i2<omax && nforce < MAXADDFORCE*tatt*5; i2++)
        {
            int im1 = FindImage(i, forceq1[i2]);
            int im2 = FindImage(i, forceq2[i2]);
            if (im1 != -1 && im2 != -1)
            {
                ForceCompPushForward(im1, im2, fscore[i2]+3);
            }
            if (im1 == -1)
            {
                ForceCompPushForward(i, forceq1[i2], fscore[i2]+3);
            }
            if (im2 == -1)
            {
                ForceCompPushForward(i, forceq2[i2], fscore[i2]+3);
            }
        }
    }

    //try again with the original node we wanted (maybe we can compute it now)
    res = ForceCompPushForward(lind, ltind, 0);
    //WriteCreationAni("creation.txt"); //make an animation, if you want
    SortForceQ(); //sort the queue so that simpler things appear first
    least = -1;

    //most of the rest of this is bookkeeping for taking the hull
    for (int i=0; i<MAXEMBS; i++)
    {
        if (forceq1[i] == -1) break;
        if (FindImage(forceq1[i], forceq2[i]) != -1) {
            least = i;
            lscore = fscore[i];
            hullv[*hullnum] = forceq1[i];
            (*hullnum)++;
            hullv[*hullnum] = forceq2[i];
            (*hullnum)++;
            hullv[*hullnum] = FindImage(forceq1[i], forceq2[i]);
            (*hullnum)++;
            break; }
    }
    if (least != -1 && res == -1 && (least < omax/40))// || omax-least > 5000))
    {
        tatt = 0;
        int *nhullv = new int[MAXEMBS];
        int nhullnum = 0;
        nhullv[0] = lind;
        nhullv[1] = ltind;
        for (int i=0; i<*hullnum; i++)
        {
            nhullv[2+i] = hullv[i];
        }
        for (int i=0; i<= least; i++)
        {
            nhullv[(*hullnum)+2*i+2] = forceq1[i];
            nhullv[(*hullnum)+2*i+1+2] = forceq2[i];
        }
        lochullv[lochullnum] = forceq1[least]; lochullnum++;
        lochullv[lochullnum] = forceq2[least]; lochullnum++;
        lochullv[lochullnum] = FindImage(forceq1[least], forceq2[least]); lochullnum++;
        if (res != -1) { lochullv[lochullnum] = res; lochullnum++; }
        for (int i=0; i<lochullnum; i++)
        {
            nhullv[(*hullnum)+(least+1)*2+2+i] = lochullv[i];
        }
        nhullnum = 2+(*hullnum) + (least+1)*2 + lochullnum;
        if (res != -1)
        {
            nhullnum++;
            nhullv[nhullnum-1] = res;
        }
        forceq1[least] = forceq2[least] = -1;
        BasicRestrict(0, nhullv, &nhullnum, 7);
        for (int i=0; i<*hullnum; i++)
        {
            hullv[i] = nhullv[i+2];
        }
        lind = nhullv[0];
        ltind = nhullv[1];
        if (res != -1) res = nhullv[nhullnum-2];
        for (int i=0; i<lochullnum; i++)
        {
            lochullv[i] = nhullv[(*hullnum)+(least+1)*2+2+i];
        }
        delete[] nhullv;
    }
    else if (least != -1 && res == -1 && least < omax/10)
    {
        lochullv[lochullnum] = forceq1[least]; lochullnum++;
        lochullv[lochullnum] = forceq2[least]; lochullnum++;
        lochullv[lochullnum] = FindImage(forceq1[least], forceq2[least]); lochullnum++;
    }
    WeedForceComp();
   // RecheckEmbs();
    //if (len < 100) DoPullbacks();
    numatt++;

}
delete[] lochullv;
return res;
}

void EMBSEQ::SortForceQ()
{
int num = -1;
for (int i=0; i<MAXEMBS; i++)
{
    if (forceq1[i] == -1) num = i;
}
if (num <= 1) return;
SortForceQ(0, num-1);
}
void EMBSEQ::SortForceQ (int start, int end){
if (start >= end) return;
else if (start-end == 1)
{
    if (fscore[end] < fscore[start])
    {
        int t = forceq1[end];
        forceq1[end] = forceq1[start];
        forceq1[start] = t;

        t = forceq2[end];
        forceq2[end] = forceq2[start];
        forceq2[start] = t;

        t = fscore[end];
        fscore[end] = fscore[start];
        fscore[start] = t;
    }
}
else {
    SortForceQ(start, (start+end)/2);
    SortForceQ((start+end)/2+1, end);
    MergeSortForceQ(start, (start+end)/2, end);
}
}
void EMBSEQ::MergeSortForceQ (int s1, int e1, int e2){
int i1=s1;
int i2=e1+1;
int* scratch1 = new int[e2-s1+1];
int* scratch2 = new int[e2-s1+1];
int* scratch3 = new int[e2-s1+1];

for (int i=0; i<e2-s1+1; i++)
{
    if (i2 > e2 || ( i1 <= e1 && fscore[i1] < fscore[i2])) {
        scratch1[i] = forceq1[i1];
        scratch2[i] = forceq2[i1];
        scratch3[i] = fscore[i1];
        i1++;
    }
    else {
        scratch1[i] = forceq1[i2];
        scratch2[i] = forceq2[i2];
        scratch3[i] = fscore[i2];
        i2++;
    }
}
for ( int i=0; i<e2-s1+1; i++)
{
    forceq1[s1+i] = scratch1[i];
    forceq2[s1+i] = scratch2[i];
    fscore[s1+i] = scratch3[i];
}
delete[] scratch1;
delete[] scratch2;
delete[] scratch3;
}

int EMBSEQ::WeedForceComp()
{
int i2=0;
int i=0;
int omax = 0;
int least = -1;
for (i=0; i<MAXEMBS; i++)
{
    if (forceq1[i]==-1) break;
    omax = i+1;
}
for (i=0; i<MAXEMBS; i++)
{
    if (forceq1[i] == -1) break;
    int res = FindImage(forceq1[i], forceq2[i]);
    if (res == -1 )//&& (i == 0 || i > omax-1000))
    {
        forceq1[i2] = forceq1[i];
        forceq2[i2] = forceq2[i];
        i2++;
    }
    else if (least == -1) least = i;
}
forceq1[i2] = -1;
forceq2[i2] = -1;
//    printf("computed %d\n", i-i2);
return least;
}


void EMBSEQ::ComputeAllPushForwards()
{
for (int i=0; i<len; i++)
{
    outputProgress(i, len, 2);
    for (int i2=0; i2< len; i2++)
    {
        ComputePushForward(i, i2, 2);

    }
}
}


void EMBSEQ::ComputeAllPushForwardsCheck()
{
for (int i=0; i<len; i++)
{
    outputProgress(i, len, 2);
    for (int i2=0; i2< len; i2++)
    {
        int res = ComputePushForward(i, i2, 2);
        if (res == -1)
        {
            int will = WillEverCompute(i, i2);
            if (will == 1)
            {
                printf("%d of %d is not computed\n", i, i2);
            }
        }
    }
}
}

void EMBSEQ::ComputeAllPushForwardsCheck(int rank)
{
for (int i=0; i<len; i++)
{
    outputProgress(i, len, 2);
    if (GetEmbRank(i) > rank) continue;
    for (int i2=0; i2< len; i2++)
    {
        if (GetEmbRank(i2) > rank) continue;
        int res = ComputePushForward(i, i2, 2);
        if (i == 40 && i2 == 6)
        {
            printf("result %d\n", res);
        }
        if (res == -1)
        {
            int will = WillEverCompute(i, i2, rank);

             //   printf("result2 %d\n", will);

            if (will == 1)
            {
                printf("%d of %d is not computed\n", i, i2);
                printf("%d, %d, %d, %d\n", GetEmbRank(i), GetEmbRank(i2), GetLeastCommonSqr(i, i2), GetLeastCommonSqr(GetLeastCommonSqr(i,i2),eord[0]));
            }
        }
    }
}
}

int EMBSEQ::FindLeastImgAbove(int ind, int tind){  //finds the least term u such that ind*v = u for some v > tind
int least = eord[len-1]; // I changed this (error)?
int otind = FindOrderInd(tind);
for (int i=otind+1; i < len; i++){
    int img = FindImage(ind, eord[i]);
    if (img != -1) return img;
}
return least;
}

void EMBSEQ::PerformPushForwardChecks(int ind, int tind){
embs[ind].AddRng(len-1, tind, -1, -1);

embs[len-1].SetCreationScore(embs[ind].GetCreationScore(), embs[tind].GetCreationScore());
embs[len-1].SetSqr(FindImage(ind, embs[tind].GetSqr()));
AddSqrSeqRng(len-1);

for (int i=0; i<embs[tind].GetNumInRange(); i++)
{
    int img = FindImage(ind, embs[tind].GetRngVec().at(i));

    if (img != -1)
    {
        int piimg = FindImage(ind, embs[tind].GetRngPIVec().at(i));

        embs[len-1].AddRng(img, piimg, embs[tind].GetRngVec().at(i), embs[tind].GetRngPIVec().at(i));
    }
}

for (int i=0; i < len; i++)
{
    int img = FindImage(ind, i);
    if (img != -1)
    {
        for (int i2=0; i2 < embs[i].GetNumInRange(); i2++)
        {
            int r1img = FindImage(ind, embs[i].GetRngVec().at(i2));
            int r2img = FindImage(ind, embs[i].GetRngPIVec().at(i2));
            if (r1img != -1 && r2img != -1)
            {
                embs[img].AddRng(r1img, r2img, embs[i].GetRngVec().at(i2), embs[i].GetRngPIVec().at(i2));
            }
        }
    }
}

for (int i=0; i<=FindOrderInd(len-1); i++)
{
    embs[len-1].AddNotRng(eord[i]);
}


for (int i=FindOrderInd(len-1); i < len; i++)
{
    embs[eord[i]].AddNotRng(len-1);
}

for (int i=0; i<embs[tind].GetNotInRange(); i++)
{
    int img = FindImage(ind, embs[tind].GetNotRngVec().at(i));
    if (img != -1)
    {
        embs[len-1].AddNotRng(img);
    }
}

for (int i=0; i< len; i++)
{
    if (embs[len-1].CheckInRng(i) == 1)
    {
        for (int i2=0; i2 < len; i2++)
        {
            if (embs[i].CheckInRng(i2) == 1)
            {
                int pi = embs[i].GetRngPI(i2);
                if (pi != -1 && FindOrderInd(pi) <= FindOrderInd(len-1))
                {
                    embs[len-1].AddNotRng(i2);
                }
            }
            int img = FindImage(i, i2);
            if (img != -1 &&
                FindOrderInd(img) <= FindOrderInd(len-1))
                embs[len-1].AddNotRng(i2);
        }
    }
}

for (int i=0; i < len; i++)
{
    if (embs[i].CheckInRng(ind)==1 && embs[i].CheckInRng(tind)==1)
    {
        embs[i].AddRng(len-1, -1, -1, -1);
    }
}

// RecheckEmbs();
}
int EMBSEQ::PushForward(int ind, int tind)
{
int liu = CheckForPushForward(ind, tind);
if (liu < 0) return -1;

InsertEmb(liu,1);

PerformPushForwardChecks(ind, tind);


return len-1;
}

void EMBSEQ::InsertAllPullbacks()
{
int oind = len-1;
int ind;
while (oind >= 0)
{
    outputProgress(len-oind, len, 2);
    ind = eord[oind];
   // SortRng(ind);
    int poind = len;
    while (poind > oind+1)
    {
        poind--;
        int pind = eord[poind];
        if (embs[ind].CheckInRngFast(pind) != 1) continue;
        if (embs[ind].GetRngPI(pind) != -1) continue;

        //printf("%d pullback of %d gives %d\n", ind, pind, len);
        InsertPullback(ind, pind);
        poind = FindOrderInd(pind);
    }
    oind = FindOrderInd(ind);
    oind--;
}
}

void EMBSEQ::InsertAllPullbacksAbove(int start, int end)
{
int olen = len;
int oind = start;
int ind;
while (oind >= end) //perform pullbacks until we reach end
{
    outputProgress(len-oind, len, 2);
    ind = eord[oind];
    if (embs[ind].GetRngChanged() == true || 1) //no longer used
    {
        //SortRng(ind); //sort the range according to the ordering, so that we can pullback from the top
        int poind = len;
        while (poind > oind+1)
        {
            poind--;
            int pind = eord[poind];
            if (embs[ind].CheckInRngFast(pind) != 1) continue;
            if (embs[ind].GetRngPI(pind) != -1) continue;
            //printf("%d pullback of %d gives %d\n", ind, pind, len);
            InsertPullback(ind, pind);
            poind = FindOrderInd(pind);
        }

        embs[ind].SetRngChanged(false);
//            if (olen < len) {RecheckEmbs(); return;}
    }
    oind = FindOrderInd(ind);
    oind--;
}
}

void EMBSEQ::SetRngChangedBelow(int ind){
for (int i=ind; i>= 0; i--){
    embs[eord[i]].SetRngChanged(true);
}
}

void EMBSEQ::DoPullbacks()
{
while (len < MAXEMBS)
{
    int olen = len;
    InsertAllPullbacksAbove(len-1, 0);
   // RecheckEmbs();
   // ComputeAllPushForwards();
    if (olen == len) return;
}
}
int EMBSEQ::CheckForPushForward(int ind, int tind)
{
if (ind == -1 || tind == -1) return -1;
if (FindOrderInd(ind) <= FindOrderInd(tind) && embs[ind].CheckInRng(tind) != 0) {//printf("shouldn't have to push forward.\n");
    return -2;}

if (FindImage(ind, tind) != -1) { //printf("already found image\n");
    return -3;}

if (embs[tind].GetSqr() == -1 || FindImage(ind, GetLeastCommonSqr(ind, tind)) == -1) { //printf("square too high \n");
    return -4;}
if (FindImage(ind, embs[tind].GetSqr()) == -1) { return -5;}

int liu = FindOrderInd(FindLeastImgAbove(ind, tind));
if (liu > 0 && liu-1 > FindOrderInd(ind))
{

    int inrng = embs[ind].CheckInRng(eord[liu-1]);
    if (inrng != 1 || (inrng == 1 && embs[ind].GetRngPI(eord[liu-1]) == -1))
    { //printf("BAD PUSHFORWARD ERROR\n");
        char str1[10000];
        char str2[10000];
        char str3[10000];
        GetNodeString(str1, ind);
        GetNodeString(str2, tind);
        GetNodeString(str3, FindOrderInd(liu));
        //printf("x\n");
        //printf("(%s*%s) %s\n", str1,str2,str3);
        fflush(stdout);
        embs[ind].CheckInRng(eord[liu-1]);
        return -6;}
}
return liu;
}
// "((j*((j*(j*j))*(j*j)))*(((((j*j)*j)*((j*j)*(j*j)))*(j*j))*(((j*j)*j)*j)))" "((j*(((j*j)*j)*j))*(((j*((j*(j*j))*(j*j)))*j)*((j*j)*j)))"

void EMBSEQ::RunGapFill(int ind, int tind){

int i = ind;
int i2 = tind;

int liu = FindOrderInd(FindLeastImgAbove(i, i2));
int inrng = embs[i].CheckInRng(eord[liu-1]);

while (inrng != 1 || (inrng == 1 && embs[i].GetRngPI(eord[liu-1]) == -1))
{

    char str1[MAXSTRLEN];
    char str2[MAXSTRLEN];
    char str3[MAXSTRLEN];
    char str4[MAXSTRLEN];
    GetNodeString(str1, i);
    GetNodeString(str2, i2);
    GetNodeString(str4, eord[liu-1]);
    sprintf(str3, "(%s*%s)", str1, str2);
    int res = CompareTwoStrings(str3, str4, 7, 7);

    printf("compared %s and %s with result %d\n", str3, str4, res);
    if (res == -1) { liu = -1; break; }
    else if (res == 1 && liu > 1)
    {
        liu--;
        inrng = embs[i].CheckInRng(eord[liu-1]);
    }
    else if (res == 0)
    {
        embs[i].AddRng(eord[liu-1], i2, -1, -1);
        SmartAlgo(1);
        liu = -1; break;
    }
    else if (res == 2)
    {
        break;
    }
    else if (res == 1 && liu <= 1)
    {
        break;
    }
}
if (liu != -1)
{
    InsertEmb(liu,1);

    PerformPushForwardChecks(i, i2);

    SmartAlgo(1);
}
}
void EMBSEQ::FillGaps(int opt, int numthread)
{
nthreads = 0;
pthread_t *threads = new pthread_t[numthread];
for (int i=0; i<len; i++)
{
    for (int i2=0; i2<len; i2++)
    {
        int liu = CheckForPushForward(i, i2);
        if (liu == -6)
        {
            if (rand()%100 > 5) continue;
            THREADARG *ta = new THREADARG(this, i, i2);

            while (nthreads >= numthread)
            {
                sleep(1);
            }
            nthreads++;
            pthread_create(&threads[0], NULL, FillGap, (void *) ta);

          /*  liu = FindOrderInd(FindLeastImgAbove(i, i2));
            int inrng = embs[i].CheckInRng(eord[liu-1]);

            while (inrng != 1 || (inrng == 1 && embs[i].GetRngPI(eord[liu-1]) == -1))
            {

                char str1[MAXSTRLEN];
                char str2[MAXSTRLEN];
                char str3[MAXSTRLEN];
                char str4[MAXSTRLEN];
                GetNodeString(str1, i);
                GetNodeString(str2, i2);
                GetNodeString(str4, eord[liu-1]);
                sprintf(str3, "(%s*%s)", str1, str2);
                int res = CompareTwoStrings(str3, str4, 6, 6);

                printf("compared %s and %s with result %d\n", str3, str4, res);
                if (res == -1) { liu = -1; break; }
                else if (res == 1 && liu > 1)
                {
                    liu--;
                    inrng = embs[i].CheckInRng(eord[liu-1]);
                }
                else if (res == 0)
                {
                    embs[i].AddRng(eord[liu-1], i2, -1, -1);
                    SmartAlgo(1);
                    liu = -1; break;
                }
                else if (res == 2)
                {
                    break;
                }
                else if (res == 1 && liu <= 1)
                {
                    break;
                }
            }
            if (liu != -1)
            {
                InsertEmb(liu,1);

                PerformPushForwardChecks(i, i2);

                SmartAlgo(1);
            }*/
        }
    }
}
while (nthreads > 0)
{
    sleep(1);
}
delete[] threads;
}

void EMBSEQ::SmartAlgo(int opt)
{
int upper = len-1;
while (len < MAXEMBS)
{
    int olen = len;
    int start = len-1;
    int end = upper;
    for (int i=len-1; i>= upper && i >= 0; i--)
    {
        int i2 = eord[i];
        if(embs[i2].GetRngChanged() == true) {
            upper = i-1;
            end = i;
            break;
        }
    }
    if (opt == 1)
    {
        for (end=2; end >=0; end--)  //perform 3 times
        {
            printf(".");
            fflush(stdout);
            InsertAllPullbacksAbove(len-1, 0); if (olen < len) { break;} //insert all pullbacks, starting from the top in the order going down. if there is a new term, break, so that we start back from the top
            PushForwardBelow(len-1, 0); if (olen < len) { break; } //if there were no new pullbacks, start computing push forwards.
            //   ComputeAllPushForwards();

        }
    }
    else if (opt == 3)
    {
        for (end=2; end >=0; end--)  //perform 3 times
        {
            printf(".");
            fflush(stdout);
            PushForwardBelow(len-1, 0); if (olen < len) { break; } //if there were no new pullbacks, start computing push forwards.
            //   ComputeAllPushForwards();

        }
    }
    else if (opt == 2)
    {
        for (end=len-1; end >=0; end--)
        {
            start = end;
            printf(".");
            fflush(stdout);
            InsertAllPullbacksAbove(start, end); if (olen < len) { break;}
            PushForwardSeq(start, end); if (olen < len) { break; }
            //   ComputeAllPushForwards();

        }
    }
    else if (opt == 4)
    {
        int olen = len;
        InsertAllPullbacksAbove(len-1, 0);
        if (true || olen == len)
        {
            for (int i2=0; false&& i2< 400000*numiterations; i2++)
            {
                int t1 = rand()%len;
                int t2 = rand()%len;
                //printf("%d %d %d\n", t1, t2, len);
                PushForward(t1, t2);
            }
            int olen = len;
            for (int i2=0; false && i2< olen; i2++)
            {
                outputProgress(i2,olen,2);
                for (int i3=0; i3<olen; i3++)
                {
                    int img = FindImage(i2,i3);
                    if (img == -1 && rand()%1 == 0)
                    {
                        PushForward(i2, i3);
                    }
                }
            }

            printf("skipnumber: %d\n", 10+10*(1<<(ILnum-7)) );
            //if (ILnum > 7) PushForwardBelowPerc(len-1, 0, 10+10*(1<<(ILnum-7)));
            //else PushForwardBelowPerc(len-1, 0, 10);

            int pb = embs[eord[0]].GetRngPI(1);
            pb = embs[eord[0]].GetRngPI(pb);

            int *hullv = new int[1000000]; //we will keep track of the indexing using this array.
            int hullnum = 0;
            for (int i=FindOrderInd(pb)-1; i >= 0; i--)
            {
                int ind = eord[i];
                ForceCompute(pb, ind, hullv, &hullnum);
            }
        }
        numiterations++;
    }
    upper = end;
    printf("done adding embeddings. (size:%d)\n", len);

    printf("rechecking computations..."); fflush(stdout);
    RecheckEmbs();
    printf ("done.\n");

    printf("computing all push forwards..."); fflush(stdout);
    ComputeAllPushForwards();
    printf("done.\n");

    SearchForPattern();
   // WriteCreationAni("realtimecreation.txt");

    if (opt == 4)
    {
        for (int e=0; false && e< len; e++)
        {
            for (int i=0; i<len; i++)
            {
                int img = FindImage(e, i);
                if (img == -1) continue;
                printf("%d, ", FindOrderInd(img)-FindOrderInd(e));
            }
        }

        printf("\n");
        int *hullv = new int[1000000]; //we will keep track of the indexing using this array.
        int hullnum = 0;
        for (int i3=0; i3< len && i3 < 10; i3++)
        {
            hullv[i3] = i3;
            hullnum++;
        }
        AddSqrRtsToHull(hullv, &hullnum, eord[len-1]);
        int pb = embs[eord[0]].GetRngPI(eord[len-1]);
        AddSqrRtsToHull(hullv, &hullnum, pb);
        pb = embs[eord[0]].GetRngPI(pb);
        AddSqrRtsToHull(hullv, &hullnum, pb);

        for (int i3=0; i3*100< len; i3++)
        {
            hullv[hullnum] = i3*100;
            hullnum++;
        }

        printf("taking hull (size:%d) (num iterations:%d)...\n", len, numiterations);
        BasicRestrict(0, hullv, &hullnum, 1000*(numiterations%2));
        printf("done. (size: %d)\n", len);
    }

    if (opt == 1) WriteSeq("backupseq.txt");
    if (opt == 2) WriteSeq("backupseq2.txt");
    WriteDiagram("diagram.txt", 1);
    if (end < 0) break;
    break; // only do one iteration (this is not the previous behavior)
}
}

void EMBSEQ::AddSqrRtsToHull(int *hullv, int * hullnum, int ind)
{
    for (int i3=0; i3 < len; i3++)
    {
        int s1 = embs[i3].GetSqr();
        if (s1 != -1)
        {
            int s2 = embs[s1].GetSqr();
            if (s2 == ind) {
                hullv[(*hullnum)] = i3;
                (*hullnum)++;
            }
        }
    }

    for (int i3=0; i3 < len; i3++)
    {
        if (embs[i3].GetSqr()== -1 || embs[embs[i3].GetSqr()].GetSqr() != ind) continue;
        for (int i4 = 0; i4< len; i4++)
        {
            if (FindOrderInd(i4) >= FindOrderInd(i3) || embs[i4].GetSqr()== -1 || embs[embs[i4].GetSqr()].GetSqr() != ind) continue;
            int p = FindImage(i3, i4);
            if (p != -1)
            {
                hullv[*hullnum] = p;
                (*hullnum)++;
            }
        }
    }
}

void EMBSEQ::PushForwardBelow()
{
int oind = len-1;
int ind;
while (oind >= 0)
{
    outputProgress(len-oind, len, 2);
    ind = eord[oind];

    for (int imgoind = len-1; imgoind >= 0; imgoind--)
    {
        int imgind = eord[imgoind];
        int res = PushForward(ind, imgind);

        //if ( res != -1)
         //   printf("%d pushforward of %d gives %d\n", ind, imgind, len-1);
    }
    oind = FindOrderInd(ind);
    oind--;
}
}

int EMBSEQ::GetSuccPreimage(int ind){
if (ind == -1) return -1;
int oind = FindOrderInd(ind);
if (oind == len-1) return -1;
int nind = eord[oind+1];

return embs[ind].GetRngPI(nind);
}

void EMBSEQ::PushForwardBelow(int start, int end)
{
int olen = len;
int oind = start;
int ind;
while (oind >= end)
{
    outputProgress(len-oind, len, 2);
    ind = eord[oind];

    int maxpre = GetSuccPreimage(ind);
    for (int imgoind = oind-1;//FindOrderInd(maxpre)-1;
         imgoind >= end; imgoind--)
    {
     //   if (end >= FindOrderInd(maxpre)-1) break;
        int imgind = eord[imgoind];
       // int imgind = eord[end];
        int res = PushForward(ind, imgind);

        //if ( res != -1)
         //   printf("%d pushforward of %d gives %d\n", ind, imgind, len-1);
      //  if (olen < len) return;
//            break;
    }
    oind = FindOrderInd(ind);
    oind--;
}
}

void EMBSEQ::PushForwardBelowPerc(int start, int end, int perc)
{
int olen = len;
int oind = start;
int ind;
while (oind >= end)
{
    if (outputProgress(len-oind, len, 2)) break;
    ind = eord[oind];

    int maxpre = GetSuccPreimage(ind);
    for (int imgoind = oind-1;//FindOrderInd(maxpre)-1;
         imgoind >= end; imgoind--)
    {
     //   if (end >= FindOrderInd(maxpre)-1) break;
        int imgind = eord[imgoind];
       // int imgind = eord[end];
        if (rand() %perc != 0) continue;
        int res = PushForward(ind, imgind);

        //if ( res != -1)
         //   printf("%d pushforward of %d gives %d\n", ind, imgind, len-1);
      //  if (olen < len) return;
//            break;
    }
    oind = FindOrderInd(ind);
    oind--;
}
}

void EMBSEQ::PushForwardSeq(int above, int push)
{
int olen = len;
int oind = len-1;
int ind;
while (oind >= above)
{
    outputProgress(len-oind, olen, 2);
    ind = eord[oind];

    int maxpre = oind;//GetSuccPreimage(ind);
    int imgoind = FindOrderInd(maxpre)-1;
    if (push <= imgoind)
    {
        int imgind = eord[push];

        int res = PushForward(ind, imgind);

       // if ( res != -1)
        //    printf("%d pushforward of %d gives %d\n", ind, imgind, len-1);
    }

    oind = FindOrderInd(ind);
    oind--;
}
}

void EMBSEQ::SortRng(int ind)
{
ComputeOrderInverse();
embs[ind].SetSortedPI(false);
for (int i=0; i < embs[ind].GetNumInRange()-1; i++)
{
    if (FindOrderInd(embs[ind].GetRngVec().at(i)) < FindOrderInd(embs[ind].GetRngVec().at(i+1)))
    {
        int temp = embs[ind].GetRngVec().at(i);
        embs[ind].SetInRange(i, embs[ind].GetRngVec().at(i + 1));
        embs[ind].SetInRange(i+1, temp);

        temp = embs[ind].GetRngPIVec().at(i);
        embs[ind].SetInRngPI(i, embs[ind].GetRngPIVec().at(i+1));//rngpi[i] = embs[ind].rngpi[i+1];
        embs[ind].SetInRngPI(i+1, temp);
        if (i > 0) i = i-2;
        else i--;
    }
}
for (int i=0; i < embs[ind].GetNotInRange()-1; i++)
{
    if (FindOrderInd(embs[ind].GetNotRngVec().at(i)) < FindOrderInd(embs[ind].GetNotRngVec().at(i + 1)))
    {
        int temp = embs[ind].GetNotRngVec().at(i);
        embs[ind].SetNotInRange(i, embs[ind].GetNotRngVec().at(i + 1));
        embs[ind].SetNotInRange(i + 1, temp);

        if (i > 0) i = i-2;
        else i--;
    }
}
}

int EMBSEQ::FindOrderInd( int ind)
{
if (ind == -1) return -1;

ComputeOrderInverse();
return eordi[ind];

for (int i=0; i < len; i++)
{
    if (eord[i] == ind) return i;
}
return -1;
}
int EMBSEQ::FindOrderIndComp( int ind)
{

for (int i=0; i < len; i++)
{
    if (eord[i] == ind) return i;
}
return -1;
}

void EMBSEQ::ComputeScores()
{
for (int i=0; i<len; i++)
{
    int nsqrts =0;
    for (int i2=0; i2< len; i2++)
    {
        if (embs[i2].GetSqr() == i) nsqrts++;
    }
    embs[i].SetScore((int)(10-nsqrts - embs[i].GetNumInRange()/2));//score = (int)(10-nsqrts - embs[i].GetNumInRange()/2);
}

}

void EMBSEQ::AddSqrSeqRng(int ind)
{
int sind = embs[ind].GetSqr();
int sindo = ind;
while (sind != -1)
{
    embs[ind].AddRng(sind, sindo, -1, -1);
    sindo = sind;
    sind = embs[sind].GetSqr();
}
}
int EMBSEQ::FindMaxInRng(int ind1, int ind2){ //find the max index which has BOTH ind1 and ind2 in range
int maxi = -1;
for (int i=0; i < len; i++)
{
    if ((ind1 == -1 || embs[eord[i]].CheckInRng(ind1)==1) && (ind2==-1 || embs[eord[i]].CheckInRng(ind2)==1))
    {
        maxi = eord[i];
    }
}
return maxi;
}
int EMBSEQ::FindMinNotInRng(int ind1, int ind2){ // find the min index with either ind1 or ind2 not in range
int maxi = eord[len-1];
for (int i=len-1; i>=0; i--)
{
    if ((ind1 != -1 && embs[eord[i]].CheckInRng(ind1)==0) || (ind2!=-1 && embs[eord[i]].CheckInRng(ind2)==0))
    {
        maxi = eord[i];
    }
}
return maxi;
}
void EMBSEQ::GetNodeString(char* str, int ind){
for (int i=1; i< 10; i++)
{
    GetNodeString(str, ind, i, 1);
    if (str[0] != 0) return;
}
str[0] = 0;
}
void EMBSEQ::PrepNodeStrings()
{
char str[MAXSTRLEN];

 SetRightSqrStrs();
 SetGenSqrStrs();
for (int i=1; i< 3; i++)
{
    for (int i2=0; i2< len; i2++)
    {
        outputProgress(i2, len, 2);
        GetNodeString(str, i2, i, 2);
    }
}
}

void EMBSEQ::SetGenSqrStrs(){ //set the powers of the generator to j, j^1, j^2,...
    int ei = eord[0];
    int i=0;

    char str[MAXSTRLEN];
    while (ei != -1)
    {

        if (i != 0) { sprintf(str, "j^%d", i); embs[ei].SetStr(str); }
        i++;
        ei = embs[ei].GetSqr();
    }
}

void EMBSEQ::SetRightSqrStrs(){ // set the right powers of the generator to j, j_1, j_2, ...
    int ei = eord[0];
    int i=0;

    char str[MAXSTRLEN];
    while (ei != -1)
    {

        if (i > 1) { sprintf(str, "j_%d", i); embs[ei].SetStr(str); }
        i++;
        ei = ComputePushForward(ei, eord[0], 6);
    }
}

void EMBSEQ::GetNodeString(char* str, int ind, int depth, int opt){
// opt == 1: work with word if it already exists
// opt == 2: try to find a new word no matter what.

if (ind == -1 || depth < 0) {str[0] = 0; return;}
if (embs[ind].GetStr() && opt == 1)
{
    //str = embs[ind].GetStr();
    strcpy(str, embs[ind].GetStr());
    return;
}
if (ind == eord[0])
{
    //str = "j";
    strcpy(str, "j");
    embs[ind].SetStr(str);
    return;
}

str[0] = 0;

int r1 = -1;
int r2 = -1;
for (int i=0; i < len; i++)
{
    if (embs[eord[i]].CheckInRngFast(ind) == 1)
    {
        r2 = embs[eord[i]].GetRngPI(ind);

        if (eord[i] == ind || r2 == ind || r2 == -1) continue;
        char ps1[MAXSTRLEN];
        char ps2[MAXSTRLEN];
        GetNodeString(ps1, eord[i], depth-1, 1);
        GetNodeString(ps2, r2, depth-1, 1);
        if (ps1[0] == 0 || ps2[0] == 0) continue;
        //QByteArray strArray = str.toUtf8();
        //char* charStr = strArray.data();
        sprintf(str, "(%s*%s)", ps1, ps2);
        embs[ind].SetStr(str);
        continue;

    }
}

}

void EMBSEQ::RemoveOuterParenth (char* str)
{
char nstr[MAXSTRLEN];
int i = 0;
while (str[i] == 32) i++;
int i2 = strlen(str)-1;
while (str[i2] == 32) i--;
while (str[i] == 40 && str[i2] == 41)
{
    int i3=i+1, pr = 0;
    while (i3 < i2)
    {
        if (str[i3] == 40) pr++;
        if (str[i3] == 41) pr--;
        if (pr < 0) break;
        i3++;
    }
    if (pr < 0) break;
    i++; i2--;
}
str[i2+1] = 0;
strcpy(nstr, &(str[i]));
strcpy(str, nstr);
}

void EMBSEQ::ProcessString(char* str)
{
    QString qs(str);

    for (int i=20; i>0; i--)
    {
        char nstr[MAXSTRLEN];
        sprintf(nstr, "j_%d", i);

        QString qsr("(j*j)");

        for (int i2 = 1; i2 < i; i2++)
        {
            qsr.replace(i2, 1, "(j*j)");
        }

        qs.replace(nstr, qsr, Qt::CaseSensitive);
    }
    for (int i=20; i>0; i--)
    {
        char nstr[MAXSTRLEN];
        sprintf(nstr, "j^%d", i);

        QString qsr("(j*j)");

        for (int i2 = 1; i2 < i; i2++)
        {
            qsr.replace(qsr.length()-1-i2, 1, "(j*j)");
        }

        qs.replace(nstr, qsr, Qt::CaseSensitive);
    }
    strcpy(str, qs.toStdString().c_str());
}

int EMBSEQ::CompFromNotation (char* str)
{
char nstr[MAXSTRLEN];
strcpy(nstr, str);

ProcessString(nstr);

RemoveOuterParenth(nstr);

if (strlen(nstr) == 1)
{
    if (str[0] == 106) return eord[0];
}
else {
    char nstr1[MAXSTRLEN];
    char nstr2[MAXSTRLEN];

    ParseString(nstr, nstr1, nstr2);
    int res1 = CompFromNotation(nstr1);
    int res2 = CompFromNotation(nstr2);

    if (res1 == -1 || res2 == -1) return -1;

    int res = ComputePushForward(res1, res2, 10);
    return res;
}
return -1;
}
int EMBSEQ::BasicRestrict(int hull1, int* hullv, int* hullnum, int num){
if (hull1 == -1) return -1;
int *nhullv = new int[num+(*hullnum)+1];
int *imgs = new int[num+(*hullnum)+1];
for (int i=0; i< num; i++)
{
    nhullv[i] = i;
}
for (int i=0; i<*hullnum; i++)
{
    nhullv[num + i] = hullv[i];
}
nhullv[num+*hullnum] = hull1;

MakeHull(nhullv, num+(*hullnum)+1, imgs);
int img = imgs[num+ *hullnum];
for (int i=0; i<*hullnum; i++)
{
    hullv[i] = imgs[num+i];
}
hullv[*hullnum] = img;
(*hullnum)++;
delete[] nhullv;
delete[] imgs;
return img;
}
#define RESTRICTNUM 7


// attempt to add the string as a node

int EMBSEQ::ForceCompFromNotation (char* str, int* hullv, int* hullnum)
{
char nstr[MAXSTRLEN];
strcpy(nstr, str);

RemoveOuterParenth(nstr); //remove the outer parentheses i.e. (j*j) -> j*j

// DoPullbacks();

if (strlen(nstr) == 1)
{
    if (str[0] == 106) return eord[0];
}
else {
    char nstr1[MAXSTRLEN];
    char nstr2[MAXSTRLEN];

    ParseString(nstr, nstr1, nstr2); //parse the binary relation, so j*(j*j) -> j, (j*j)
    if (strlen(nstr1) == 0 || strlen(nstr2) == 0) return -1;
    int res1 = ForceCompFromNotation(nstr1, hullv, hullnum); // try to add nstr1 as a node
    res1 = BasicRestrict(res1, hullv, hullnum, RESTRICTNUM); // take the hull, which means we throw away nodes which are not essential.
    int res1num = (*hullnum)-1;
    if (res1 == -1) return -1;
    int res2 = ForceCompFromNotation(nstr2, hullv, hullnum); //try to add nstr2 as a node
    res2 = BasicRestrict(res2, hullv, hullnum, RESTRICTNUM); // take the hull
    int res2num = (*hullnum)-1;

    if (res1 == -1 || res2 == -1) return -1;

    res1 = hullv[res1num];
    res2 = hullv[res2num];
    int res = ForceCompute(res1, res2, hullv, hullnum); //try to put add nstr1 * nstr2
    res = BasicRestrict(res, hullv, hullnum, RESTRICTNUM); // take hull
    //printf("*\n");
    printf("x"); //print an x, which means we succeeded at some level
    fflush(stdout);
    return res;
}
return -1;
}

void EMBSEQ::ParseString(char* str, char* p1, char* p2){
int i=0, pr=0;

while (str[i] != 42 || pr > 0)
{
    if (str[i] == 40) pr++;
    if (str[i] == 41) pr--;
    i++;
}
str[i] = 0;
strcpy(p1, str);
str[i] = 42;
strcpy(p2, &(str[i+1]));
}

void EMBSEQ::ComputeImages (int ind)
{
printf("\n");
for (int i=0; i< len; i++)
{
    int rind = ComputePushForward(ind, i, 5);
    int ycoord = -1;
    if (rind != -1) ycoord = embs[rind].GetYScore();
    else { continue; ycoord = -1; }
    int xcoord = embs[i].GetYScore();
    printf("(%d, %d),", xcoord, ycoord);
}
printf("\n");
}

int EMBSEQ::FindXCoord(int ind)
{
int lsqr =  GetLeastCommonSqr(ind, eord[0]);
if (lsqr == ind) return GetPowerEqual(eord[0], ind);
return 1;
}

int EMBSEQ::GetPowerEqual(int ind, int pwr)
{
int i=0;
int tind = ind;
while (tind != pwr && embs[tind].GetSqr() != -1)
{
    tind = embs[tind].GetSqr();
    i++;
}
if (tind == pwr) return i;
return -1;
}

int EMBSEQ::CompareFromNotation(char* str1, char* str2)
{
int p1 = CompFromNotation(str1);
int p2 = CompFromNotation(str2);
if (p1 == -1 || p2 == -1) return -1;
int i1 = FindOrderInd(p1);
int i2 = FindOrderInd(p2);
if (i1 < i2) return 1;
if (i1 > i2) return 2;
return 0;
}

void EMBSEQ::MakeHull(int *hullv, int vlen, int *imgs){
int *keep = new int[len];

for (int i=0; i<len; i++){
    keep[i] = 0;
    for (int i2=0; i2<vlen; i2++) {
        if (hullv[i2] == i) { keep[i] = 1; break; }
    }
}
int okeep = -1;
int nkeep = 0;
while (okeep != nkeep)
{
    okeep = nkeep;
    for (int i=0; i<len; i++)
    {
        /*if (keep[i] == 0)
        {
            for (int i2=0; i2<embs[i].nrng; i2++)
            {
                if (keep[embs[i].rng[i2]] == 1) { keep[i] = 1; break;}
            }
        }*/
        if (keep[i] == 1)
        {
            if (embs[i].GetSqr() != -1) keep[embs[i].GetSqr()] = 1;

            for (int i2=0; i2<embs[i].GetNumInRange(); i2++)
            {
                if (keep[embs[i].GetRngVec().at(i2)] == 1 && embs[i].GetRngPIVec().at(i2) != -1) {
                    keep[embs[i].GetRngPIVec().at(i2)] = 1; }
            }
        }
    }
    nkeep = 0;
    for (int i=0; i<len; i++) { if (keep[i]==1) nkeep++; }
}
int * kimgs = new int[len];
RestrictEmbs(keep, kimgs);
for (int i=0; i<vlen; i++)
{
    imgs[i] = kimgs[hullv[i]];
}
delete[] keep;
delete[] kimgs;
}

void EMBSEQ::RestrictEmbs(int *keep, int *imgs)
{
int nkeep = 0;
for (int i=0; i<len; i++)
{
    if (keep[i] ==1) nkeep++;
}
if (nkeep == 0) return;

int* neord = new int[MAXEMBS];
EMB* nembs = new EMB[MAXEMBS];

int i=0;
int i2=0;
for (i = 0; i<len; i++)
{
    if (keep[i] == 0) continue;
    imgs[i] = i2;
    i2++;
}

i = i2 = 0;
for (i = 0; i< len; i++)
{
    if (keep[eord[i]] == 0) continue;
    neord[i2] = imgs[eord[i]];
    i2++;
}

i = i2 = 0;
for (i = 0; i<len; i++)
{
    if (keep[i] == 0) continue;
    if (embs[i].GetSqr() != -1 && keep[embs[i].GetSqr()] == 1) nembs[imgs[i]].SetSqr(imgs[embs[i].GetSqr()]);//sqr = imgs[embs[i].GetSqr()];
    for (i2 = 0; i2<embs[i].GetNumInRange(); i2++)
    {
        if (keep[embs[i].GetRngVec().at(i2)] == 1)
        {
            int pi = -1;
            if (embs[i].GetRngPIVec().at(i2) != -1 && keep[embs[i].GetRngPIVec().at(i2)] == 1)
                pi = imgs[embs[i].GetRngPIVec().at(i2)];
            nembs[imgs[i]].AddRng(imgs[embs[i].GetRngVec().at(i2)], pi, -1, -1);
        }
    }
    for (i2 = 0; i2 < embs[i].GetNotInRange(); i2++)
    {
        if (keep[embs[i].GetNotRngVec().at(i2)] == 1)
        {
            nembs[imgs[i]].AddNotRng(imgs[embs[i].GetNotRngVec().at(i2)]);
        }
    }
}

for (int i=0; i<MAXEMBS; i++)
{
    if (forceq1[i] == -1) break;
    if (keep[forceq1[i]] == 0 || keep[forceq2[i]] == 0) { forceq1[i] = -1; break; }
    forceq1[i] = imgs[forceq1[i]];
    forceq2[i] = imgs[forceq2[i]];
}
len = nkeep;
invcomp = 0;

delete[] embs;
delete[] eord;

embs = nembs;
eord = neord;

}

void EMBSEQ::PrintLeftDivision(int ind1, int ind2)
{
int res[MAXLEFTDIV*2+1];

int num = FindLeftDivision(ind1, ind2, res);
if (num <= 0) { printf("could not do left division of %d and %d\n", ind1, ind2); return;}

char str1[MAXSTRLEN];
char str2[MAXSTRLEN];
char str3[MAXSTRLEN];

GetNodeString(str1, ind1);
GetNodeString(str2, ind2);
printf("left division from %d:%s to %d:%s\n", ind1, str1, ind2, str2);
for (int i=0; i< num; i++)
{
    GetNodeString(str3, res[i]);
    printf("%d: %s\n", res[i], str3);
}
}
int EMBSEQ::FindLeftDivision(int ind1, int ind2, int* res){

for (int i=1; i< 10; i++)
{
    int num = FindLeftDivision(ind1, ind2, res, i);
    if (num > 0) return num;
}
return -1;
}

int EMBSEQ::FindLeftDivision(int ind1, int ind2, int* res, int depth)
{
if (depth < 0 || ind1 == -1 || ind2 == -1) return -1;
int o1 = FindOrderInd(ind1);
int o2 = FindOrderInd(ind2);

if (o1 > o2) {res[0] = -1; return 0;}
if (o1 == o2) {res[0] = -1; return 0;}
int nres[MAXLEFTDIV];

int bestnum = MAXLEFTDIV;
if (embs[ind1].CheckInRng(ind2) == 1) {res[0] = embs[ind1].GetRngPI(ind2); res[1] = -1; return 1;}
for (int i=0; i < embs[ind1].GetNumInRange(); i++)
{
    int onew = FindOrderInd(embs[ind1].GetRngVec().at(i));
    if (onew > o2) continue;
    int resnum = FindLeftDivision(embs[ind1].GetRngVec().at(i), ind2, nres, depth-1);
    if (resnum > 0 && resnum < bestnum && resnum < MAXLEFTDIV)
    {
        for (int i2=0; i2<resnum; i2++)
        {
            res[i2+1] = nres[i2];
        }
        res[0] = embs[ind1].GetRngVec().at(i);
        res[resnum+2] = -1;
        bestnum = resnum;
    }
}
if (bestnum >=MAXLEFTDIV) bestnum = 0;
else bestnum++;
return bestnum;
}
int EMBSEQ::WriteDiagram(const char* file, int opt)
{

FILE* fi = fopen(file,  "w");
if (!fi) return 0;

ComputeScores();

for (int i=0; i<len; i++)
{
    float x = (float) 500 *  log(embs[eord[i]].GetCreationScore()*10);

    float y = (float) -1000 * log(200+embs[eord[i]].GetYScore());

    fprintf(fi, "\\node[] (A%d) at (%f,%f) {$%s$};\n", eord[i], x/500, -y/500, embs[eord[i]].GetStr());
    //fprintf(fi, "\\node[] (A%d) at (%f,%f) {$%s$};\n", eord[i], (float)embs[eord[i]].GetScore(), (float)i, embs[eord[i]].GetStr());
}
for (int i=0; i<len; i++)
{
    if (embs[i].GetSqr() != -1)
    {
        fprintf(fi, "\\draw [->] (A%d) -- (A%d) ;\n", i, embs[i].GetSqr());

    }
}
for (int i=0; i<len && opt == 1; i++)
{
    for (int i2=0; i2 < embs[i].GetNumInRange(); i2++)
    {
        if (embs[i].GetRngVec().at(i2)==embs[i].GetSqr()) continue;
        fprintf(fi, "\\draw [dashed] (A%d) --node[auto] {$%d$}  (A%d) ;\n", i, embs[i].GetRngPIVec().at(i2), embs[i].GetRngVec().at(i2));
    }
}
fclose(fi);
return 1;
}
int EMBSEQ::WriteCreationAni(const char* file)
{

FILE* fi = fopen(file,  "w");
if (!fi) return 0;

ComputeScores();
ComputeOrderInverse();

for (int n=0; n < len; n++)
{


    fprintf(fi, "\\begin{tikzpicture}\n\n" );

for (int i=0; i<=n; i++)
{
    fprintf(fi, "\\node[] (A%d) at (%d,%f) {$%d$};\n", i, embs[i].GetScore(), 0.3*(float)FindOrderInd(i), i);
}
for (int i=0; i<=n; i++)
{
    if (embs[i].GetSqr() != -1 && embs[i].GetSqr() <= n)
    {
        fprintf(fi, "\\draw [->] (A%d) -- (A%d) ;\n", i, embs[i].GetSqr());

    }
}

    fprintf(fi, "\n\n\\end{tikzpicture}\n\n\\newpage\n\n");

}

fprintf(fi, "\\end{document}");
fclose(fi);
return 1;
}

int EMBSEQ::WriteSeq(const char* file)
{
FILE* fi = fopen(file,  "w");
if (!fi) return 0;

fprintf(fi, "%d\n", len);

for (int i=0; i<len; i++){
    fprintf(fi, "%d ", eord[i]);
}

fprintf(fi, "\n");

for (int i=0; i<len; i++)
{
    embs[i].WriteToFile(fi);
}

fclose(fi);
return 1;
}

int EMBSEQ::LoadSeq(const char* file){

FILE* fi = fopen(file,  "r");
if (!fi) return 0;


fscanf(fi, "%d\n", &len);

for (int i=0; i<len; i++){
    fscanf(fi, "%d ", &(eord[i]));
}


for (int i=0; i<len; i++)
{
    embs[i].LoadFromFile(fi);
}
invcomp = 0;

fclose(fi);
return 1;
}

void EMBSEQ::SearchForPattern() // searches for a pattern where a^2 = b and c = b * a, but there is a d with a < d < b and c in the rng of d
{
    printf("searching for pattern type 1...");
    for (int i=0; false && i< len; i++)
    {
        for (int i3=0; i3 < embs[i].GetNumInRange(); i3++)
        {
            int rngi = embs[i].GetRngVec().at(i3);
            if (rngi == -1) continue;

            int img = ComputePushForward(rngi, i, 3);
            if (img == -1) continue;
            for (int i2=0; i2 < len; i2++)
            {
                if (FindOrderInd(i2) > FindOrderInd(i) && FindOrderInd(i2) < FindOrderInd(rngi)) {

                   if(embs[i2].CheckInRngFast(img) == 1)
                   {
                        printf("Found example pattern!\n");
                        PrepNodeStrings();
                        printf("%s, %s, %s\n", embs[i].GetStr(), embs[rngi].GetStr(), embs[i2].GetStr());
                    }
                    for (int i4=0; i4 < embs[i2].GetRngVec().size(); i4++)
                    {
                        if (FindOrderInd(embs[i2].GetRngVec().at(i4)) > FindOrderInd(rngi) && FindOrderInd(embs[i2].GetRngVec().at(i4)) < FindOrderInd(img))
                        {
                            printf("Found (weak) example pattern!\n");
                            PrepNodeStrings();
                            printf("%s, %s, %s, %s\n", embs[i].GetStr(), embs[rngi].GetStr(), embs[i2].GetStr(), embs[embs[i2].GetRngVec().at(i4)].GetStr());
                        }
                    }
                }
            }
        }
    }
    printf("done.\n");

    PrepNodeStrings();

    printf("searching for pattern type 2..."); fflush(stdout);
    for (int i=0; false && i< len; i++)
    {
        int rank1 = GetEmbRank(i);
        for (int i3=0; i3 < embs[i].GetNumInRange(); i3++)
        {
            int rngi = embs[i].GetRngVec().at(i3);
            if (rngi == -1) continue;
            int pbi = embs[i].GetRngPI(rngi);
            if (pbi == -1) continue;

            int rankrng = GetEmbRank(rngi);
            int rankpb = GetEmbRank(pbi);

            if ((rankrng < rank1 && rankrng >= rankpb) ||
                    (rankrng >= rank1 && rankrng < rankpb) ||
                    (rankrng < rank1-1 && rankrng < rankpb-1))
            {
                printf("Found example pattern!\n");
                printf("%s, %s, %s\n", embs[i].GetStr(), embs[pbi].GetStr(), embs[rngi].GetStr());
                printf("ranks: %d, %d, %d\n", rank1, rankpb, rankrng);
            }

        }
    }
    printf("done.\n");


    printf("searching for pattern type 3..."); fflush(stdout);

    int ind1 = embs[eord[0]].GetRngPI(0);
    int ind2 = embs[eord[0]].GetRngPI(1);
    int ind3 = embs[eord[0]].GetRngPI(2);
    int ind4 = embs[eord[0]].GetRngPI(3);
    printf("indices: %d %d %d %d\n", ind1, ind2, ind3, ind4);
    for (int i=0; i< len; i++)
    {
        if (embs[i].GetSqr() != ind1) continue;
        for (int i2=0; i2< len; i2++)
        {
            if (embs[i2].GetSqr() != ind2) continue;
            for (int i3=0; i3<len; i3++){
                if (embs[i3].GetSqr() != ind3) continue;
                for (int i4=0; i4<len; i4++){
                    if (embs[i4].GetSqr() != ind4) continue;

                    //printf("%d, %d, %d, %d\n", i, i2, i3, i4);
                    if (FindOrderInd(i) > FindOrderInd(i2)) continue;
                    //if (embs[i].CheckInRng(i2) != 1) continue;

                    if(FindOrderInd(i3) > FindOrderInd(i4)) continue;
                   // if (embs[i3].CheckInRng(i4) != 1) continue;
                    //printf("rng completely checks out\n");
                    //printf("%d, %d, %d, %d\n", i, i2, i3, i4);

                    printf("*");
                    int img1 = FindImage(i2, i); if(img1 == -1) { img1 = PushForward(i2,i); if (img1 == -1) continue;}
                    printf("&");
                    int img2 = FindImage(i4,i3); if (img2 == -1) { img2 = PushForward(i4,i3); if (img2 == -1) continue;}
                    printf("order: %d, %d, %d, %d\n", FindOrderInd(i2), FindOrderInd(i4), FindOrderInd(img1), FindOrderInd(img2));
                    if (FindOrderInd(i2)< FindOrderInd(i4) && FindOrderInd(i4) < FindOrderInd(img1) && FindOrderInd(img1) < FindOrderInd(img2))
                    {
                        printf("Found example pattern (3)!\n");
                        printf("%s, %s, %s, %s\n", embs[i].GetStr(), embs[i2].GetStr(),embs[i3].GetStr(),embs[i4].GetStr());

                    }
                }
            }
        }
    }
    printf("done.\n");
    printf("searching for pattern type 4..."); fflush(stdout);

    for (int j1=0; false && j1 < len; j1++)
    {
        outputProgress(j1, len, 2);
        int ind1 = j1;
        for (int j2=0; j2< len; j2++)
        {
            if (FindOrderInd(j2) <= FindOrderInd(j1) || embs[j1].GetSqr() != embs[j2].GetSqr()) continue;
           int  ind2 = j2;
            for (int j3=0; j3 < len; j3++)
            {
                if (FindOrderInd(j3) <= FindOrderInd(j2)) continue;
               int  ind3 = j3;
                for (int j4=0; j4<len; j4++)
                {

                    if (FindOrderInd(j4) <= FindOrderInd(j3) || embs[j3].GetSqr() != embs[j4].GetSqr()) continue;
                    int ind4 = j4;

                    //int ind1 = embs[eord[0]].GetRngPI(0);
                    //int ind2 = embs[eord[0]].GetRngPI(1);
                    //int ind3 = embs[eord[0]].GetRngPI(2);
                    //int ind4 = embs[eord[0]].GetRngPI(3);
                   // printf("indices: %d %d %d %d\n", ind1, ind2, ind3, ind4);
                    for (int i=0; i< len; i++)
                    {
                        if (embs[i].GetSqr() != ind1) continue;
                        for (int i2=0; i2< len; i2++)
                        {
                            if (embs[i2].GetSqr() != ind2) continue;
                            for (int i3=0; i3<len; i3++){
                                if (embs[i3].GetSqr() != ind3) continue;
                                for (int i4=0; i4<len; i4++){
                                    if (embs[i4].GetSqr() != ind4) continue;

                                    //printf("%d, %d, %d, %d\n", i, i2, i3, i4);
                                    if (FindOrderInd(i) > FindOrderInd(i2)) continue;
                                    //if (embs[i].CheckInRng(i2) != 1) continue;

                                    //printf("*");
                                    if(FindOrderInd(i3) > FindOrderInd(i4)) continue;
                                   // if (embs[i3].CheckInRng(i4) != 1) continue;
                                    //printf("rng completely checks out\n");
                                    //printf("%d, %d, %d, %d\n", i, i2, i3, i4);
                                    int img1 = FindImage(i2, i); if(img1 == -1) continue;
                                    int img2 = FindImage(i4,i3); if (img2 == -1) continue;
                                    //printf("order: %d, %d, %d, %d\n", FindOrderInd(i2], FindOrderInd(i4], FindOrderInd(img1], FindOrderInd(img2]);

                                    if (FindOrderInd(i2)< FindOrderInd(i4) && FindOrderInd(i4) < FindOrderInd(img1))
                                    {
                                        printf("Found example pattern (4)-almost!\n");
                                        printf("%s, %s, %s, %s\n", embs[i].GetStr(), embs[i2].GetStr(),embs[i3].GetStr(),embs[i4].GetStr());

                                    }
                                    if (FindOrderInd(i2)< FindOrderInd(i4) && FindOrderInd(i4) < FindOrderInd(img1) && FindOrderInd(img1) < FindOrderInd(img2))
                                    {
                                        printf("Found example pattern (4)!\n");
                                        printf("%s, %s, %s, %s\n", embs[i].GetStr(), embs[i2].GetStr(),embs[i3].GetStr(),embs[i4].GetStr());

                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    printf("done.\n");
}
