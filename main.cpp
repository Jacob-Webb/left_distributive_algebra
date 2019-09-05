//
//  main.cpp
//  Graphs LDA
//
//  Created by Scott Cramer on 8/16/16.
//  Copyright © 2016 Scott Cramer. All rights reserved.
//

#include <iostream>
#include "embseq.h"
#include "graphwidget.h"

#include <QApplication>
#include <QTime>
#include <QMainWindow>

char* CreateRandomStr(int maxdepth, int maxlength, int num)
{
    char* str = new char[3000];
    int depthside[maxdepth];
    int sqrtvar = jASCII+rand()%num;
    //    if(sqrtvar == oASCII) sqrtvar = pASCII;
    for(int i2=0; i2 < maxdepth; i2++) depthside[i2] = 0;
    str[0] = LEFTP;
    int i=1;
    int depth = 0;
    int extrap = 1;
    int dec;
    while (depth >= 0)
    {
        dec = rand()%2;
        if(dec == 0 || depth >= maxdepth || i >= maxlength)
        {
            str[i] = sqrtvar; i++;
            if (depthside[depth] == 0)
            {
                str[i] = starASCII; i++;
                depthside[depth] = 1;

            }
            else
            {
                while (depthside[depth]==1 && depth >= 0)
                {
                    str[i++] = RIGHTP;
                    depth--;
                }
                if (depth >= 0)
                {
                    depthside[depth] = 1;
                    str[i++] = starASCII;
                }
            }
        }
        else
        {
            str[i] = LEFTP; i++; depth++;
            depthside[depth] = 0;
            continue;
        }


    }
    str[i] = 0;
    return str;

}

void CompareSeqs(EMBSEQ* seq1, EMBSEQ* seq2)
{
    int * f = new int[seq1->len];

    char str[MAXSTRLEN];
    for (int i=0; i< seq1->len; i++) f[i] = -1;
    for (int i=0; i<seq1->len; i++)
    {
        int ind = seq1->eord[i];

        seq1->GetNodeString(str, ind);
        int img = seq2->CompFromNotation(str);
        f[ind] = img;
        if (img == -1)
        {
            printf("%d: %s no image\n", ind, str);
        }
       /* if (i == 0) f[ind] = seq2->eord[0];
        else
        {
            int mx1 = seq1->FindMaxInRng(ind, -1);
            if (mx1 == -1) continue;
            int pi1 = seq1->embs[mx1].GetRngPI(ind);
            if (pi1 == -1 || f[pi1] == -1 || f[mx1] == -1) continue;
            int img = seq2->FindImage(f[pi1], f[mx1]);
            if (img == -1) continue;
            f[ind] = img;
        }

        int e1 = ind;
        int e2 = f[ind];
        while (seq1->embs[e1].sqr != -1 && seq2->embs[e2].sqr != -1)
        {
            e1 = seq1->embs[e1].sqr;
            e2 = seq2->embs[e2].sqr;

            f[e1] = e2;
        }*/
    }
    for (int i=0; i< seq1->len; i++)
    {
        printf("(%d,%d),", i, f[i]);
    }
    printf("\n");
    delete f;
}

//EMBSEQ globalseq(200000);


//compare the two strings to see which is greater, i ranges from nstart to nend, with the value of i determining what left power of the generator we are working below. The larger i is, the more chance of succeeding, but the slower things are going to be.
int CompareTwoStrings(char* str1, char* str2, int nstart, int nend)
{
   // if (globalseq.len == 0)
    //    globalseq.LoadSeq("globseq.txt");

    //we try each i from nstart to nend
    for (int i=nstart; i<= nend; i++)
    {
          EMBSEQ globalseq(1000000); //this contains our term data. We initialize with the max # terms
        if (globalseq.len == 0)
            globalseq.CreateILSeq(i); // this creates our base terms, from which we will build up the rest of the structure
        globalseq.DoPullbacks(); // insert all possible pullbacks. This means that we divide on the left, is possible
        int *hullv = new int[1000000]; //we will keep track of the indexing using this array.
        int hullnum = 0;

        int node1 = globalseq.ForceCompFromNotation(str1, hullv, &hullnum); //attempt to insert the term in str1, node1 has the index for str1
        if (node1 == -1) continue;
        int node1i = 0; // will contain the "hull" index of node1 -- this is not the same as the index. We have to have this "hull" index for when we take a "hull" below (this means we throw away a lot of the nodes to keep things small). The "hull" index allows us to recover the new index below.
        for (int i2=0; i2< hullnum; i2++){
            if (hullv[i2] == node1) node1i = i2; //this is the "hull" index
        }
        int node2 = globalseq.ForceCompFromNotation(str2, hullv, &hullnum); //attempt to add in str2 and put the index for it in node2
        if (node1 == -1 || node2 == -1) continue;
        node1 = hullv[node1i]; // the index of node1 could change, so we get it from the "hull" index

        printf("done thinking. finding simplest terms.\n");

        globalseq.PrepNodeStrings(); //this speeds up getting the string
        char str1[10000];
        char str2[10000];
        globalseq.GetNodeString(str1, node1); //obtain the string from the node
        globalseq.GetNodeString(str2, node2);
        printf("term 1 equivalent form: %s\n", str1); // print the strings
        printf("term 2 equivalent form: %s\n", str2);
        if (node1 == node2) return 0; // equivalent terms
        if (globalseq.FindOrderInd(node1) < globalseq.FindOrderInd(node2)) return 1; // node1 < node2
        else return 2; // node1 > node2
    }
    return -1;
}

void *FillGap(void *ta)
{
    THREADARG* tha = (THREADARG*)ta;
    //tha->seq->nthreads++;
    printf("NEW THREAD: compute %d of %d\n", tha->i, tha->i2);
    tha->seq->RunGapFill(tha->i, tha->i2);
    tha->seq->nthreads--;
    delete tha;
    pthread_exit(NULL);
}

int main(int argc, char * argv[]) {

    QApplication app(argc, argv);

    GraphWidget* widget = new GraphWidget;
    QMainWindow mainWindow;
    mainWindow.setCentralWidget(widget);

    mainWindow.show();
    return app.exec();
}
