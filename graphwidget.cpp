/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the examples of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:BSD$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** BSD License Usage
** Alternatively, you may use this file under the terms of the BSD license
** as follows:
**
** "Redistribution and use in source and binary forms, with or without
** modification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     the documentation and/or other materials provided with the
**     distribution.
**   * Neither the name of The Qt Company Ltd nor the names of its
**     contributors may be used to endorse or promote products derived
**     from this software without specific prior written permission.
**
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
** OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
** LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
** DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
** THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
** (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
** OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE."
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "graphwidget.h"
#include "node.h"
#include "edge.h"

#include <iostream>
#include <math.h>
#include <string>

#include <QMouseEvent>
#include <QKeyEvent>



/*
 * constructor sets up scene information, runs the logic of the graph from EMBSEQ, and creates nodes
 * sequence1 is I think a second, smaller option to sequence2. The program only runs sequence2 right now
*/
GraphWidget::GraphWidget(QWidget* parent)
    : QGraphicsView(parent), sequence1(100000), sequence2(1000000)
{
    //set up scene
    currScene = new QGraphicsScene(this);
    currScene->setItemIndexMethod(QGraphicsScene::BspTreeIndex);
    setScene(currScene);    //this QGraphicsView is currently showing currScene
    setCacheMode(CacheBackground);
    setRenderHint(QPainter::Antialiasing);
    setTransformationAnchor(AnchorUnderMouse);
    scale(qreal(0.8), qreal(0.8));
    setMinimumSize(400,400);
    setWindowTitle(tr("EMBSEQ"));
    maxCreationScore = 200;

    //run underlying logic of graph
    runEMBSEQ();


    //create nodes, edges, and add to scene
    createNodes(currScene);
    connectAllNodes(currScene);

}

void GraphWidget::runEMBSEQ()
{
    qDebug() << "What order would you like to examine?";
    int userInput = 7;
    std::cin >> userInput;
    //All functions taken from original main file
    sequence2.CreateILSeq(userInput);
    sequence2.SmartAlgo(1);
    sequence2.PrepNodeStrings();
}

void GraphWidget::runAlgoIteration()
{
    sequence2.SmartAlgo(1);
    sequence2.PrepNodeStrings();
}

void GraphWidget::checkComputations(int rank){
    qDebug() << "Checking for computations to compute...";
    sequence2.ComputeAllPushForwardsCheck(rank);
    qDebug() << "done.";
}

void GraphWidget::createNodes(QGraphicsScene* scene)
{
    qreal x, y;     //x and y will hold the x and y coordinates of nodes

    //for each node created in the sequence2
    for(int i = 0; i < sequence2.GetLength(); ++i)
    {
        if (sequence2.GetEmbs()[i].GetCreationScore() > maxCreationScore) continue;

        if (sequence2.GetEmbRank(i) > 5) {
                                                continue;}
        createOneNode(scene, i);

        /*QString word;
        word = sequence2.GetEmbs()[i].GetStr();                 //get the embedding's word before creating the node b/c the emb's str will be changed after
        Node* node = new Node(this, sequence2.GetEmbs()[i]);    //create a node
        node->setWord(word);
        addNode(node);                                          //add the node to the nodes vector
        node->setOrderValue(sequence2.GetEordi()[i]);           //set the node's < ordering value
        scene->addItem(node);
        this->connect(node, SIGNAL(nodeClicked(Node*)),  //set the signal between the node and this graphics view
                        this, SLOT(drawAssocNodes(Node*)));
        x = (float) 50 *  sequence2.GetEmbs()[i].GetCreationScore();
        y = (float) -50 * node->getOrderValue();
        node->setPos(x, y);                                     //sets the node's variable on the scene at x and y
*/
    }
    for (int i=0; i < nodes.size(); i++)
    {
        nodes.at(i)->updateRange();
    }
}

int GraphWidget::checkCreatedNode(int ind){
    int numNodes = nodes.size();
    for (int i=0; i<numNodes; i++)
    {
        if (nodes.at(i)->getIndex() == ind) return i;
    }
    return -1;
}

void GraphWidget::createOneNode(QGraphicsScene* scene, int ind)
{
    //return;
    int check = checkCreatedNode(ind);
    if (check != -1) {
        if (scene->items().indexOf(nodes[check]) == -1) scene->addItem(nodes[check]);
        return;
    }

    int numNodes = nodes.size();
    qreal x,y;

    QString word;
    word = sequence2.GetEmbs()[ind].GetStr();                 //get the embedding's word before creating the node b/c the emb's str will be changed after
    Node* node = new Node(this, &sequence2.GetEmbs()[ind]);    //create a node
    node->setWord(word);
    addNode(node);                                          //add the node to the nodes vector
    node->setOrderValue(sequence2.GetEordi()[ind]);           //set the node's < ordering value
    scene->addItem(node);
    this->connect(node, SIGNAL(nodeClicked(Node*)),  //set the signal between the node and this graphics view
                    this, SLOT(drawAssocNodes(Node*)));
    x = (float) 500 *  log(sequence2.GetEmbs()[ind].GetCreationScore()*10);

    int ct = 0;
    for (int i2=0; i2 < numNodes; i2++) if (sequence2.GetEordi()[ind] > sequence2.GetEordi()[nodes.at(i2)->getIndex()]) ct++;
    y = (float) -50 * ct;

    y = (float) -1000 * log(200+sequence2.GetEmbs()[ind].GetYScore());

    //qDebug() << "coords: " << x << " " << y;
    node->setPos(x, y);                                     //sets the node's variable on the scene at x and y

}

void GraphWidget::updateNodeCoords()
{
    int numNodes = nodes.size();
    for (int i=0; i< numNodes; i++)
    {
        qreal x,y;

        Node* node = nodes.at(i);
        int ind = node->getIndex();
        node->setOrderValue(sequence2.GetEordi()[ind]);
        x = (float) 500 *  log(sequence2.GetEmbs()[ind].GetCreationScore()*10);

        int ct = 0;
        for (int i2=0; i2 < numNodes; i2++) if (sequence2.GetEordi()[ind] > sequence2.GetEordi()[nodes.at(i2)->getIndex()]) ct++;
        y = (float) -50 * ct;

        y = (float) -1000 * log(200+sequence2.GetEmbs()[ind].GetYScore());
        //qDebug() << "coords: " << x << " " << y;
        node->setPos(x, y);
    }
   //currScene->update(0, 0, 400, 400);

    readjustCoords();
    redrawEdges();
}

void GraphWidget::readjustCoords()
{
    int numNodes = nodes.size();
    if (numNodes == 0) return;
    qreal maxx = -100000, maxy = -100000, minx = 100000, miny = 100000;

  //  qDebug() << "size" << currScene->items().size();
    for (int i=0; i < numNodes;i++)
    {
        Node* node = nodes.at(i);
        if (currScene->items().indexOf(node) == -1) continue;
        QPointF point = node->pos();
        if (maxx < point.rx()) maxx = point.rx();
        if (minx > point.rx()) minx = point.rx();
        if (maxy < point.ry()) maxy = point.ry();
        if (miny > point.ry()) miny = point.ry();
    }

  //  qDebug() << "old: " << minx << " " << maxx << " " << miny << " " << maxy;

    if (maxx <= minx || maxy <= miny) return;


    qreal newmaxx = 3000;
    qreal newminx = 100;
    qreal newmaxy = 0;
    qreal newminy = -5000;
    for (int i=0; i< numNodes; i++)
    {
        Node* node = nodes.at(i);
        if (currScene->items().indexOf(node) == -1) continue;
        QPointF point = node->pos();
      //  qDebug() << "old " << node->pos();
        node->setPos((point.rx()-minx)/(maxx-minx)*(newmaxx-newminx)+newminx,
                     (point.ry()-miny)/(maxy-miny)*(newmaxy-newminy)+newminy);
      //  qDebug() << "new " << node->pos();
    }
}

void GraphWidget::redrawEdges()
{
    int numEdges = edges.size();
    for(int i = 0; i < numEdges; ++i)
        edges[i]->adjust();
}

void GraphWidget::addNode(Node* node) { nodes.push_back(node); }
void GraphWidget::addEdge(Edge* edge) { edges.push_back(edge); }

QVector<Node*> GraphWidget::getNodes() { return nodes; }
QVector<Edge*> GraphWidget::getEdges() { return edges; }

Node* GraphWidget::getNodeFromIndex(int ind){
    for (int i=0; i < nodes.size(); i++)
    {
        if (nodes.at(i)->getIndex() == ind) return nodes.at(i);
    }
    return NULL;
}

void GraphWidget::drawAssocNodes(Node* coreNode)
{

    int numEdges = edges.size();
    int numNodes = nodes.size();
    int numSceneItems = currScene->items().size();

    //If we are dealing with a graph that has ALL information we erase everything so that we can draw what information we want
    if(numSceneItems == numEdges + numNodes)
    {
        //remove all edges on screen
        for(int i = 0; i < numEdges; ++i)
            currScene->removeItem(edges[i]);

        //remove all nodes on screen
        for(int i = 0; i < numNodes; ++i)
            currScene->removeItem(nodes.at(i));
    }

    //add coreNode and all nodes in range of coreNode back to the screen
    currScene->addItem(coreNode);
    int numInRange = coreNode->getInRange().size();
    int numNotInRange = coreNode->getEmb()->GetNotInRange();
    //for each node in coreNode's range
    for(int i = 0; i < numInRange; ++i)
    {
        //go through every available node on the graph
        for(int j = 0; j < numNodes; ++j)
        {
            //if the node is in coreNode's range
            if(nodes.at(j)->getIndex() == coreNode->getInRange().at(i))
            {
                //add back to the scene
                currScene->addItem(nodes.at(j));
                break;
            }
        }
    }

    //add all edges between coreNode and in range nodes back to screen
    //for each node in the graph
    for(int i = 0; i < numNodes; ++i)
    {
        //and for every edge in the graph
        for(int j = 0; j < numEdges; ++j)
        {
            //if an edge shares another node as it's destination and coreNode as its source, add it
            if(edges.at(j)->destNode() == nodes.at(i) && edges.at(j)->sourceNode() == coreNode)
            {
                currScene->addItem(edges.at(j));
                break;
            }
        }
    }

    //for each node
    for(int i = 0; i < numNodes; ++i)
    {
        //if coreNode is in range, add it
        if (nodes.at(i)->checkInRange(coreNode) != -1)
        {
            currScene->addItem(nodes.at(i));
        }
    }
    //add all edges between some node and the coreNode
    //and for every edge in the graph
        for(int j = 0; j < numEdges; ++j)
        {
            //if an edge shares another node as it's destination and coreNode as its source, add it
            if(edges.at(j)->destNode() == coreNode)
            {
                currScene->addItem(edges.at(j));
            }
        }

    //clear the console in linux.
    system("clear");

    //print information on node and associated nodes to console
    qDebug() << "This node" << coreNode->getIndex() << "has: ";
    for(int i = 0; i < numInRange; ++i)
    {
        printf( "%d ", coreNode->getInRange().at(i));
    }
    printf("\n");
    qDebug() << "in range.";
    /*qDebug() << "This node" << coreNode->getIndex() << "has: ";
    for(int i = 0; i < numNotInRange; ++i)
    {
        qDebug() << coreNode->getEmb()->GetNotRngVec().at(i);
    }
    qDebug() << " not in range.";*/

    qDebug() << "This node's word" << coreNode->getWord();
    qDebug() << "node rank: " << sequence2.GetEmbRank(coreNode->getIndex());
    qDebug() << "node order number: " << coreNode->getOrderValue();
}

void GraphWidget::drawAssocNodesWOErase(Node* coreNode)
{

    int numEdges = edges.size();
    int numNodes = nodes.size();
    int numSceneItems = currScene->items().size();

    //add coreNode and all nodes in range of coreNode back to the screen
    currScene->addItem(coreNode);
    int numInRange = coreNode->getInRange().size();
    int numNotInRange = coreNode->getEmb()->GetNotInRange();
    //for each node in coreNode's range
    for(int i = 0; i < numInRange; ++i)
    {
        //go through every available node on the graph
        for(int j = 0; j < numNodes; ++j)
        {
            //if the node is in coreNode's range
            if(nodes.at(j)->getIndex() == coreNode->getInRange().at(i))
            {
                //add back to the scene
                currScene->addItem(nodes.at(j));
                break;
            }
        }
    }

    //add all edges between coreNode and in range nodes back to screen
    //for each node in the graph
    for(int i = 0; i < numNodes; ++i)
    {
        //and for every edge in the graph
        for(int j = 0; j < numEdges; ++j)
        {
            //if an edge shares another node as it's destination and coreNode as its source, add it
            if(edges.at(j)->destNode() == nodes.at(i) && edges.at(j)->sourceNode() == coreNode)
            {
                currScene->addItem(edges.at(j));
                break;
            }
        }
    }

    //for each node
    for(int i = 0; i < numNodes; ++i)
    {
        //if coreNode is in range, add it
        if (nodes.at(i)->checkInRange(coreNode) != -1)
        {
            currScene->addItem(nodes.at(i));
        }
    }
    //add all edges between some node and the coreNode
    //and for every edge in the graph
        for(int j = 0; j < numEdges; ++j)
        {
            //if an edge shares another node as it's destination and coreNode as its source, add it
            if(edges.at(j)->destNode() == coreNode)
            {
                currScene->addItem(edges.at(j));
            }
        }

}

void GraphWidget::connectAllNodes(QGraphicsScene* scene)
{
    QColor blue(Qt::blue);          //blue will represent an edge between a node in range of another
    QColor red(Qt::red);            //red will represent an edge between a node and its square
    int numNodes = nodes.size();

    //for each node in created in the graph
    for(int i = 0; i < numNodes; ++i)
    {
        if (scene->items().indexOf(nodes.at(i)) == -1) continue;
       int rangeSize = nodes.at(i)->getInRange().size();
       //and for each value in the range of the node
       for(int j = 0; j < rangeSize; ++j)
       {
           if (scene->items().indexOf(nodes.at(j)) == -1) continue;
           //and for every other node
           for(int k = 0; k < numNodes; ++k)
           {
               //if (nodes.at(i)->getIndex() == )

               //if the inner and outer node are not the same and the second node is in range of the first
                if((nodes.at(i) != nodes.at(k)) && (nodes.at(k)->getIndex() == nodes.at(i)->getInRange().at(j)))
                {
                   //and if the node in range is a square of the node
                   if(nodes.at(k)->getIndex() == nodes.at(i)->getSquare())
                   {
                       //add an edge representing a square association
                       createOneEdge(scene, nodes.at(i), nodes.at(k), red);
                       /*Edge* edge = new Edge(nodes.at(i), nodes.at(k), red);
                       addEdge(edge);
                       scene->addItem(edge);*/
                       break;
                   }
                   else
                   {
                      //otherwise add an edge representing an in range assocation
                       createOneEdge(scene, nodes.at(i), nodes.at(k), blue);
                      /*Edge* edge = new Edge(nodes.at(i), nodes.at(k), blue);
                      addEdge(edge);
                      scene->addItem(edge);*/
                      break;
                   }
               }
           }
       }
    }
}

void GraphWidget::createOneEdge(QGraphicsScene* scene, Node *src, Node *dest, QColor col){
    int numEdges = edges.size();
    for (int i=0; i<numEdges; i++)
    {
        if (edges.at(i)->sourceNode() == src && edges.at(i)->destNode() == dest)
        {
            if (scene->items().indexOf(edges.at(i)) == -1) scene->addItem(edges.at(i));
            return;
        }
    }
    Edge* edge = new Edge(src, dest, col);
    addEdge(edge);
    scene->addItem(edge);
}

#ifndef QT_NO_WHEELEVENT
void GraphWidget::wheelEvent(QWheelEvent *event)
{
    //zoom in or out on the scene based on the angle of the mouse wheel movement
    scaleView(pow((double)2, -event->delta() / 240.0));
}
#endif

void GraphWidget::keyPressEvent(QKeyEvent* event)
{
    switch (event->key()) {
    case Qt::Key_Up:
        zoomIn();
        break;
    case Qt::Key_Down:
        zoomOut();
        break;
    case Qt::Key_Right:
        maxCreationScore+=50;
        createNodes(currScene);
        connectAllNodes(currScene);
        updateNodeCoords();
        break;
    case Qt::Key_I:
        inputCommand();
        break;
    case Qt::Key_N:
        runAlgoIteration();
        updateNodeCoords();
        createNodes(currScene);
        connectAllNodes(currScene);
        updateNodeCoords();
        break;
    case Qt::Key_M:
        while(1)
        {
            runAlgoIteration();
            updateNodeCoords();
            createNodes(currScene);
            connectAllNodes(currScene);
            updateNodeCoords();
        }
        break;
    case Qt::Key_C:
        checkComputations(4);
        break;
    case Qt::Key_S:
        sequence2.WriteDiagram("diagram-new.txt", 1);
        break;
    default:
        QGraphicsView::keyPressEvent(event);
    }
}

/*
void GraphWidget::mousePressEvent(QMouseEvent* event)
{
    QGraphicsView::mousePressEvent(event);
    //this->currScene->clear();
}
*/
void GraphWidget::mouseDoubleClickEvent(QMouseEvent* event)
{
    int numNodes = nodes.size();
    int numEdges = edges.size();
    int numSceneItems = currScene->items().size();

    //if we are not dealing with a graph with ALL nodes and edges
    if(numSceneItems != numNodes + numEdges)
    {
        //add every node back to the screen
        for(int i = 0; i < numNodes; ++i)
            currScene->addItem(nodes.at(i));

        //add every edge back to the screen
        for(int i = 0; i < numEdges; ++i)
            currScene->addItem(edges.at(i));

        //clear console of warnings that some nodes and edges are already on screen
        system("clear");
    }
}

void GraphWidget::inputCommand()
{
    std::cin.clear();
    qDebug() << "What would you like to do? (1 for compute term, 2 for compute pullback, 5 for compute product, 3 for select index, 4 for change graph)";
    int userInput;
    std::cin >> userInput;

    qDebug() <<"huh?";

    if (userInput == 1)
    {
        qDebug() << "Input the word you would like to compute:";
        std::string str;
        std::cin >> str;
        int res = sequence2.CompFromNotation((char*)(str.c_str()));
        if (res == -1) qDebug()<< "Could not compute the term.\n";
        else
        {
            qDebug() << "Term has index " << res << "\n";
            qDebug() << "alternate word: " << sequence2.GetEmbs()[res].GetStr() << "\n";

        }
    }
    else if (userInput == 2)
    {
        qDebug() << "To compute a^(-1) b (the pullback of b by a), enter index of a:";
        int aind;
        std::cin >> aind;
        qDebug() << "Enter index of b:";
        int bind;
        std::cin >> bind;
        if (aind < 0 || aind >= sequence2.GetLength()) { qDebug() << "index error."; return; }
        if (bind < 0 || bind >= sequence2.GetLength()) { qDebug() << "index error."; return; }
        int res = sequence2.GetEmbs()[aind].GetRngPI(bind);
        if (res == -1) qDebug()<< "Could not compute the pullback.\n";
        else
        {
            qDebug() << "Pullback has index " << res << "\n";
            qDebug() << "and word: " << sequence2.GetEmbs()[res].GetStr() << "\n";

        }
    }
    else if (userInput == 5)
    {
        qDebug() << "To compute a * b, the product of a and b, enter index of a:";
        int aind;
        std::cin >> aind;
        qDebug() << "Enter index of b:";
        int bind;
        std::cin >> bind;
        if (aind < 0 || aind >= sequence2.GetLength()) { qDebug() << "a index error."; return;}
        if (bind < 0 || bind >= sequence2.GetLength()) { qDebug() << "b index error."; return; }
        int res = sequence2.ComputePushForward(aind, bind, 10);
        if (res == -1) qDebug()<< "Could not compute the product.\n";
        else
        {
            qDebug() << "product has index " << res << "\n";
            qDebug() << "and word: " << sequence2.GetEmbs()[res].GetStr() << "\n";

        }
    }
    else if (userInput == 3)
    {
        qDebug() << "Input the index you would like to select:";
        int aind;
        std::cin >> aind;
        int res = aind;
        if (aind < 0 || aind >= sequence2.GetLength())  qDebug() << "index error.";
        else
        {
            qDebug() <<"index images:";
            sequence2.ComputeImages(aind);
            qDebug() << "Word for index: " << sequence2.GetEmbs()[res].GetStr() << "\n";
            qDebug() << "order index: " << sequence2.FindOrderInd(res);
            qDebug() << "index rank: " << sequence2.GetEmbRank(res);
            Node* node = getNodeFromIndex(aind);
            if (node) drawAssocNodes(node);
            else {
                createOneNode(currScene, aind);
                connectAllNodes(currScene);
                updateNodeCoords();
            }

        }
    }
    else if (userInput == 4)
    {
        qDebug() << "How would you like to change the graph: (1 for add node range, 2 add nodes with term in range)";
        int userInput2;
        std::cin >> userInput2;
        if (userInput2 == 1)
        {
            qDebug() << "enter lower bound for range:";
            int lb;
            std::cin >> lb;
            qDebug() << "enter upper bound for range:";
            int ub;
            std::cin >> ub;
            if (lb < 0 || ub < 0 || lb > sequence2.GetLength() || ub > sequence2.GetLength() || lb > ub) {qDebug() << "bounds error"; return;}

            //remove all edges on screen
          /*  for(int i = 0; i < edges.size(); ++i)
                if (currScene->items().indexOf(edges[i])!= -1) currScene->removeItem(edges[i]);

            //remove all nodes on screen
            for(int i = 0; i < nodes.size(); ++i)
                if (currScene->items().indexOf(nodes[i])!= -1)  currScene->removeItem(nodes.at(i));
*/
            for (int i=lb; i<= ub; i++)
            {
                int ind = sequence2.GetEord()[i];
                drawAssocNodes(nodes.at(ind));
            }
            //connectAllNodes(currScene);
           // updateNodeCoords();

        }
        else if (userInput2 == 2)
        {
            qDebug() << "enter index that must be in range:";
            int rind;
            std::cin >> rind;
            if (rind < 0 || rind > sequence2.GetLength() ) {qDebug() << "index error"; return;}

            for (int i=0; i< sequence2.GetLength(); i++)
            {
                if (sequence2.GetEmbs()[i].CheckInRngFast(rind) == 1)
                    createOneNode(currScene, i);
            }
            connectAllNodes(currScene);
            updateNodeCoords();

        }
    }
}


void GraphWidget::scaleView(qreal scaleFactor)
{
    qreal factor = transform().scale(scaleFactor, scaleFactor).mapRect(QRectF(0, 0, 1, 1)).width();
    if (factor < 0.07 || factor > 100)
        return;

    scale(scaleFactor, scaleFactor);
}

void GraphWidget::zoomIn()
{
    scaleView(qreal(1.2));
}

void GraphWidget::zoomOut()
{
    scaleView(1 / qreal(1.2));
}
