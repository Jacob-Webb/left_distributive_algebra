#ifndef NODE_H
#define NODE_H

#include "emb.h"

#include <QGraphicsObject>
#include <QString>
#include <QVariant>
#include <QVector>
#include <QDebug>

class Edge;

class GraphWidget;
QT_BEGIN_NAMESPACE
class QGraphicsSceneMouseEvent;
QT_END_NAMESPACE

/*
 * Each node will be a graphical representation of an emb
*/
class Node : public QGraphicsObject{

    Q_OBJECT

public:
    /*
     * this constructor isn't used. Should set up interaction with the Graphwidget
    */
    Node(GraphWidget* graphWidget);

    /*
     * construct with an emb object to initialize variables
    */
    Node(GraphWidget* graphWidget, EMB* embedding);

    /*
     * updateRange updates the Range vector by adding indices from the embedding
    */
    void updateRange();

    /*
     * setEmb sets the "thisEmb" member variable to the "embedding" arg
    */
    void setEmb(EMB* embedding);

    /*
     * setOrderValue sets the "orderValue" data member with the "ordering" arg
    */
    void setOrderValue(int ordering);

    /*
     * setWord takes in a QString to set the word of the node
    */
    void setWord(QString embWord);

    /*
     * getEmb returns "thisEmb" data member
    */
    EMB* getEmb();

    /*
     * getInRange returns the inRange vector
    */
    QVector<int> getInRange();

    /*
     * checkInRange checks if a certain node is in the range
    */
    int checkInRange(Node* nd);

    /*
     * getSquare returns the "square" data member
    */
    int getSquare();

    /*
     * getIndex() returns the "index" data member
    */
    int getIndex();

    /*
     * getOrderValue returns the "orderValue" data member
    */
    int getOrderValue();

    /*
     * getWord returns the "word" data member
    */
    QString getWord();

    /*
     * getPosition returns the "position" data member
    */
    QPointF getPosition();

    /*
     * boundingRect returns a QRectF for the bounding rectangle of the graphics item
    */
    QRectF boundingRect() const Q_DECL_OVERRIDE;

    /*
     * shape returns an ellipse of similar dimensions to the graphics item's bounding rectangle
    */
    QPainterPath shape() const Q_DECL_OVERRIDE;

    /*
     * paint draws the shape, colors in, and writes the label to the node
    */
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget) Q_DECL_OVERRIDE;

signals:
    /*
     * nodeClicked signals to a slot when the mouse is pressed on this. Sends this
    */
    void nodeClicked(Node* node);

protected:
    /*
     * mousePressEvent sends a signal to the graph that the node has been clicked by the mouse
    */
    void mousePressEvent(QGraphicsSceneMouseEvent* event) Q_DECL_OVERRIDE;

private:
    GraphWidget* graph;     //graph represents the graph that the node will be displayed on
    QString label;          //should hold the index of the graph as a string to be used to label each node
    QString word;           //represents the string for the embedding i.e. j or (j*j)

    QVector<int> inRange;   //a list to match the indices of other nodes that are in range of this

    int square;             //int representing the index of the node that is the square of this
    int orderValue;         //the orderValue represents the < ordering of the node in relation to other nodes
    int index;              //the int representation of this node
    EMB* thisEmb;            //the emb object that this node is associated with
};
#endif // NODE_H
