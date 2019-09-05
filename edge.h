#ifndef EDGE_H
#define EDGE_H

#include <QtWidgets/qgraphicsitem.h>

class Node;

class Edge : public QGraphicsItem
{
public:
    Edge(Node* sourceNode, Node* destNode);
    Edge(Node* sourceNode, Node* destNode, QColor color);

    Node* sourceNode() const;
    Node* destNode() const;

    void setColor(QColor color);

    void adjust();

    enum { Type = UserType + 2 };
    int type() const Q_DECL_OVERRIDE { return Type; }

protected:
    QRectF boundingRect() const Q_DECL_OVERRIDE;
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget);

private:
    Node* source;
    Node* dest;

    QColor color;

    QPointF sourcePoint;
    QPointF destPoint;
    qreal arrowSize;
};

#endif // EDGE_H
