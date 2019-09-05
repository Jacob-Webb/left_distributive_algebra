/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the examples of the Qt Toolkit.
**
*ification, are permitted provided that the following conditions are
** met:
**   * Redistributions of source code must retain the above copyright
**     notice, this list of conditions and the following disclaimer.
**   * Redistributions in binary form must reproduce the above copyright
**     notice, this list of conditions and the following disclaimer in
**     * $QT_BEGIN_LICENSE:BSD$
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
** modthe documentation and/or other materials provided with the
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

#include "node.h"

#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QPainter>
#include <QStyleOption>
#include <iostream>
#include <QString>
#include <QDebug>

using namespace std;

/*
 * this constructor isn't used. Should set up interaction with the Graphwidget
*/
Node::Node(GraphWidget* graphWidget)
    : graph(graphWidget)
{
    setCacheMode(DeviceCoordinateCache);
    setZValue(1);
}

/*
 * construct with an emb object to initialize variables
*/
Node::Node(GraphWidget* graphWidget, EMB* embedding)
    : graph(graphWidget)
{
    //option for caching drawn objects to speed up redrawing when needed
    setCacheMode(DeviceCoordinateCache);
    setZValue(1);      //zValue is the stacking order of drawn objects.

    thisEmb = embedding;
    index = thisEmb->GetMyInd();
    square = thisEmb->GetSqr();
    //use thisEmb's rng vector to populate node's inRange vector
    for(int i = 0; i < thisEmb->GetNumInRange(); ++i)
    {
        inRange.push_back(thisEmb->GetRngVec().at(i));
    }
}

void Node::updateRange()
{
    for(int i = 0; i < thisEmb->GetNumInRange(); ++i)
    {
        int rngSize = inRange.size();
        int embRng = thisEmb->GetRngVec().at(i);
        int added = 0;
        for (int i2=0; i2 < rngSize; i2++) { if (inRange.at(i2) == embRng) { added = 1; break; }}
        if (added == 0) inRange.push_back(embRng);
    }
}

void Node::setEmb(EMB* embedding) { thisEmb = embedding; }
void Node::setOrderValue(int ordering) { orderValue = ordering; }
void Node::setWord(QString embWord) { word = embWord; }

EMB* Node::getEmb() { return thisEmb; }
QVector<int> Node::getInRange() { return inRange; }
int Node::getSquare() { return square; }
int Node::getIndex() { return index; }
int Node::getOrderValue() { return orderValue; }
QString Node::getWord() { return word; }

int Node::checkInRange(Node* nd)
{
    int sz = inRange.size();
    for (int i=0; i< sz; i++)
    {
        if (inRange.at(i) == nd->index)
            return i;
    }
    return -1;
}

QRectF Node::boundingRect() const
{
    qreal adjust = 2;       //adjust for outline stroke
    return QRectF( -10 - adjust, -10 - adjust, 30 + adjust, 30 + adjust);
}

QPainterPath Node::shape() const
{
    QPainterPath path;
    path.addEllipse(-10, -10, 30, 30);
    return path;
}

void Node::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget*)
{
    painter->setPen(Qt::NoPen);
    painter->setBrush(Qt::darkGray);
    painter->drawEllipse(-7, -7, 20, 20);

    QRadialGradient gradient(-3, -3, 10);
    gradient.setColorAt(0, Qt::yellow);
    gradient.setColorAt(1, Qt::darkYellow);
    painter->setBrush(gradient);

    label = QString::number(index);
    // To add the index as a label
    painter->setPen(QPen(Qt::black, 0));
    painter->drawEllipse(-10, -10, 25, 25);     //draw the boundary of the node
    painter->drawText(-5, 5, label);

}

void Node::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    emit nodeClicked(this);
}
