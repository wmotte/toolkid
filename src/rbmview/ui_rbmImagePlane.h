/********************************************************************************
** Form generated from reading ui file 'rbmImagePlane.ui'
**
** Created: Tue May 26 12:07:27 2009
**      by: Qt User Interface Compiler version 4.5.0
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_RBMIMAGEPLANE_H
#define UI_RBMIMAGEPLANE_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QWidget>
#include "QVTKWidget.h"

QT_BEGIN_NAMESPACE

class Ui_ImagePlaneClass
{
public:
    QHBoxLayout *horizontalLayout;
    QVTKWidget *m_VTK;

    void setupUi(QWidget *ImagePlaneClass)
    {
        if (ImagePlaneClass->objectName().isEmpty())
            ImagePlaneClass->setObjectName(QString::fromUtf8("ImagePlaneClass"));
        ImagePlaneClass->resize(400, 300);
        horizontalLayout = new QHBoxLayout(ImagePlaneClass);
        horizontalLayout->setSpacing(0);
        horizontalLayout->setMargin(0);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        m_VTK = new QVTKWidget(ImagePlaneClass);
        m_VTK->setObjectName(QString::fromUtf8("m_VTK"));

        horizontalLayout->addWidget(m_VTK);


        retranslateUi(ImagePlaneClass);

        QMetaObject::connectSlotsByName(ImagePlaneClass);
    } // setupUi

    void retranslateUi(QWidget *ImagePlaneClass)
    {
        ImagePlaneClass->setWindowTitle(QApplication::translate("ImagePlaneClass", "ImagePlane", 0, QApplication::UnicodeUTF8));
        Q_UNUSED(ImagePlaneClass);
    } // retranslateUi

};

namespace Ui {
    class ImagePlaneClass: public Ui_ImagePlaneClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_RBMIMAGEPLANE_H
