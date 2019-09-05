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

#ifndef GRAPHWIDGET_H
#define GRAPHWIDGET_H

#include <QGraphicsView>
#include <QVector>
#include <QDebug>

#include "embseq.h"


class Node;
class Edge;

class GraphWidget : public QGraphicsView
{
    Q_OBJECT

public:
    /*
     * The constructor sets up the scene, runs the underlying EMBSEQ program,
     * and adds all of the nodes and edges to the graph.
    */
    GraphWidget(QWidget* parent = 0);

    /*
     * runEMBSEQ runs the main functions from EMBSEQ that set up the logic of the graph
    */
    void runEMBSEQ();
    /*
     * runAlgoIteration runs the main functions from EMBSEQ that set up the logic of the graph
    */
    void runAlgoIteration();
    /*
     * checkComputations checks to see if all computations can be made for a certain rank
    */
    void checkComputations(int rank);

    /*
     * createNodes creates the nodes for the graph taken from underlying logic and places them on graph
    */
    void createNodes(QGraphicsScene* scene);
    /*
     * createOneNode adds a specific node to the graph and puts it in the scene, adding edges as well.
    */
    void createOneNode(QGraphicsScene* scene, int ind);
    /*
     * updateNodeCoords updates the coordinates of nodes that have been created
    */
    void updateNodeCoords();
    /*
     * readjustCoords updates the coordinates of nodes that have been created
    */
    void readjustCoords();
    /*
     * redrawEdges removes all edges and adds them so they are redrawn
    */
    void redrawEdges();
    /*
     * checkCreatedNode checks to see if the specific node has already been created
    */
    int checkCreatedNode(int ind);

    /*
     * addNode adds a node to this QGraphicsView's list of nodes
    */
    void addNode(Node* node);

    /*
     * addEdge adds an edge to this QGraphicsView's list of edges
    */
    void addEdge(Edge* edge);

    /*
     * createOneEdge creates a new edge (if it has not already been created)
    */
    void createOneEdge(QGraphicsScene* scene, Node* src, Node* dest, QColor col);

    /*
     * getNodes() returns the list of nodes that have been added to the graph
    */
    QVector<Node*> getNodes();

    /*
     * getNodeFromIndex() returns the node with the given index
    */
    Node* getNodeFromIndex(int ind);

    /*
     * getEdges() returns the list of edges that have been added to the graph
    */
    QVector<Edge*> getEdges();

public slots:
    /*
     * drawAssocNodes erases everything from currScene then adds the given node
     * and all nodes in the node's range along with associated edges back to currScene
    */
    void drawAssocNodes(Node* node);

    //Same, except no erasing to begin with
    void drawAssocNodesWOErase(Node* node);

protected:
    /*
     * connectAllNodes searches all relationships between all nodes and creates edges between them accordingly
    */
    void connectAllNodes(QGraphicsScene* scene);

    /*
     * wheelEvent zooms in or out on the scene with the movement of the mouse wheel
    */
    #ifndef QT_NO_WHEELEVENT
    void wheelEvent(QWheelEvent *event);
    #endif

    /*
     * keyPressEvent zooms in when the up key is pressed and out when the down key is pressed
    */
    void keyPressEvent(QKeyEvent* event);

    /*
     * mousePressEvent isn't used but could have some functionality when the scene is pressed
    */
    //void mousePressEvent(QMouseEvent* event);

    /*
     * mouseDoubleClickEvent redraws all nodes and edges so that the complete graph is on scene
    */
    void mouseDoubleClickEvent(QMouseEvent* event);

    /*
     * inputCommand inputs a command from the command line
    */
    void inputCommand();

    /*
     * scaleView transforms the scene to scale up or down on the scene
    */
    void scaleView(qreal scaleFactor);

    void zoomIn();
    void zoomOut();

private:
    EMBSEQ sequence1;           //one set of sequences of embeddings from logic
    EMBSEQ sequence2;           //right now the main sequence of embeddings from logic

    int maxCreationScore;       // max creation score of node to display

    QGraphicsScene* currScene;  //pointer to current scene filling the QGraphicsView

    QVector<Node*> nodes;       //list of nodes added to currScene
    QVector<Edge*> edges;       //list of edges added to currScene
};

#endif // GRAPHWIDGET_H

