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
