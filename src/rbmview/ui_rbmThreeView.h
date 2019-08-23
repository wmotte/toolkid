/********************************************************************************
** Form generated from reading ui file 'rbmThreeView.ui'
**
** Created: Tue May 26 12:07:27 2009
**      by: Qt User Interface Compiler version 4.5.0
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_RBMTHREEVIEW_H
#define UI_RBMTHREEVIEW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QFrame>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>
#include <QtGui/QWidget>
#include "rbmImagePlane.h"

QT_BEGIN_NAMESPACE

class Ui_ThreeViewClass
{
public:
    QGridLayout *gridLayout;
    QVBoxLayout *verticalLayout;
    QLabel *m_Label1;
    rbm::ImagePlane *m_Widget1;
    QVBoxLayout *verticalLayout_2;
    QLabel *m_Label2;
    rbm::ImagePlane *m_Widget2;
    QVBoxLayout *verticalLayout_3;
    QLabel *m_Label3;
    rbm::ImagePlane *m_Widget3;
    QFrame *m_Frame;
    QVBoxLayout *verticalLayout_4;
    QHBoxLayout *horizontalLayout_2;
    QGridLayout *gridLayout_2;
    QDoubleSpinBox *m_SpinX;
    QDoubleSpinBox *m_SpinY;
    QDoubleSpinBox *m_SpinZ;
    QSpinBox *m_SpinT;
    QLabel *m_LabelX;
    QLabel *m_LabelY;
    QLabel *m_LabelZ;
    QSpacerItem *horizontalSpacer;
    QSpacerItem *horizontalSpacer_2;
    QHBoxLayout *horizontalLayout_6;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QDoubleSpinBox *m_SpinLevel;
    QLabel *label_2;
    QDoubleSpinBox *m_SpinWindow;
    QSpacerItem *horizontalSpacer_5;
    QHBoxLayout *horizontalLayout_5;
    QLabel *label_3;
    QLabel *m_LabelIntensity;
    QSpacerItem *horizontalSpacer_3;
    QHBoxLayout *horizontalLayout_7;
    QPushButton *pushButton;
    QSpacerItem *horizontalSpacer_4;
    QSpacerItem *verticalSpacer;

    void setupUi(QWidget *ThreeViewClass)
    {
        if (ThreeViewClass->objectName().isEmpty())
            ThreeViewClass->setObjectName(QString::fromUtf8("ThreeViewClass"));
        ThreeViewClass->resize(678, 503);
        ThreeViewClass->setAutoFillBackground(true);
        gridLayout = new QGridLayout(ThreeViewClass);
        gridLayout->setSpacing(6);
        gridLayout->setMargin(11);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        verticalLayout = new QVBoxLayout();
        verticalLayout->setSpacing(6);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        m_Label1 = new QLabel(ThreeViewClass);
        m_Label1->setObjectName(QString::fromUtf8("m_Label1"));
        QSizePolicy sizePolicy(QSizePolicy::Ignored, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(m_Label1->sizePolicy().hasHeightForWidth());
        m_Label1->setSizePolicy(sizePolicy);
        QPalette palette;
        QBrush brush(QColor(255, 255, 255, 255));
        brush.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Base, brush);
        QBrush brush1(QColor(255, 85, 0, 255));
        brush1.setStyle(Qt::SolidPattern);
        palette.setBrush(QPalette::Active, QPalette::Window, brush1);
        palette.setBrush(QPalette::Inactive, QPalette::Base, brush);
        palette.setBrush(QPalette::Inactive, QPalette::Window, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette.setBrush(QPalette::Disabled, QPalette::Window, brush1);
        m_Label1->setPalette(palette);
        m_Label1->setAutoFillBackground(true);

        verticalLayout->addWidget(m_Label1);

        m_Widget1 = new rbm::ImagePlane(ThreeViewClass);
        m_Widget1->setObjectName(QString::fromUtf8("m_Widget1"));
        QSizePolicy sizePolicy1(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(m_Widget1->sizePolicy().hasHeightForWidth());
        m_Widget1->setSizePolicy(sizePolicy1);

        verticalLayout->addWidget(m_Widget1);


        gridLayout->addLayout(verticalLayout, 0, 0, 1, 1);

        verticalLayout_2 = new QVBoxLayout();
        verticalLayout_2->setSpacing(6);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        m_Label2 = new QLabel(ThreeViewClass);
        m_Label2->setObjectName(QString::fromUtf8("m_Label2"));
        sizePolicy.setHeightForWidth(m_Label2->sizePolicy().hasHeightForWidth());
        m_Label2->setSizePolicy(sizePolicy);
        QPalette palette1;
        palette1.setBrush(QPalette::Active, QPalette::Base, brush);
        palette1.setBrush(QPalette::Active, QPalette::Window, brush1);
        palette1.setBrush(QPalette::Inactive, QPalette::Base, brush);
        palette1.setBrush(QPalette::Inactive, QPalette::Window, brush1);
        palette1.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette1.setBrush(QPalette::Disabled, QPalette::Window, brush1);
        m_Label2->setPalette(palette1);
        m_Label2->setAutoFillBackground(true);

        verticalLayout_2->addWidget(m_Label2);

        m_Widget2 = new rbm::ImagePlane(ThreeViewClass);
        m_Widget2->setObjectName(QString::fromUtf8("m_Widget2"));
        sizePolicy1.setHeightForWidth(m_Widget2->sizePolicy().hasHeightForWidth());
        m_Widget2->setSizePolicy(sizePolicy1);
        m_Widget2->setMinimumSize(QSize(267, 183));

        verticalLayout_2->addWidget(m_Widget2);


        gridLayout->addLayout(verticalLayout_2, 0, 1, 1, 1);

        verticalLayout_3 = new QVBoxLayout();
        verticalLayout_3->setSpacing(6);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        m_Label3 = new QLabel(ThreeViewClass);
        m_Label3->setObjectName(QString::fromUtf8("m_Label3"));
        sizePolicy.setHeightForWidth(m_Label3->sizePolicy().hasHeightForWidth());
        m_Label3->setSizePolicy(sizePolicy);
        QPalette palette2;
        palette2.setBrush(QPalette::Active, QPalette::Base, brush);
        palette2.setBrush(QPalette::Active, QPalette::Window, brush1);
        palette2.setBrush(QPalette::Inactive, QPalette::Base, brush);
        palette2.setBrush(QPalette::Inactive, QPalette::Window, brush1);
        palette2.setBrush(QPalette::Disabled, QPalette::Base, brush1);
        palette2.setBrush(QPalette::Disabled, QPalette::Window, brush1);
        m_Label3->setPalette(palette2);
        m_Label3->setAutoFillBackground(true);

        verticalLayout_3->addWidget(m_Label3);

        m_Widget3 = new rbm::ImagePlane(ThreeViewClass);
        m_Widget3->setObjectName(QString::fromUtf8("m_Widget3"));
        sizePolicy1.setHeightForWidth(m_Widget3->sizePolicy().hasHeightForWidth());
        m_Widget3->setSizePolicy(sizePolicy1);

        verticalLayout_3->addWidget(m_Widget3);


        gridLayout->addLayout(verticalLayout_3, 1, 0, 1, 1);

        m_Frame = new QFrame(ThreeViewClass);
        m_Frame->setObjectName(QString::fromUtf8("m_Frame"));
        sizePolicy1.setHeightForWidth(m_Frame->sizePolicy().hasHeightForWidth());
        m_Frame->setSizePolicy(sizePolicy1);
        m_Frame->setFrameShape(QFrame::StyledPanel);
        m_Frame->setFrameShadow(QFrame::Raised);
        verticalLayout_4 = new QVBoxLayout(m_Frame);
        verticalLayout_4->setSpacing(6);
        verticalLayout_4->setMargin(11);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setSpacing(6);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        gridLayout_2 = new QGridLayout();
        gridLayout_2->setSpacing(6);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        m_SpinX = new QDoubleSpinBox(m_Frame);
        m_SpinX->setObjectName(QString::fromUtf8("m_SpinX"));
        m_SpinX->setDecimals(2);

        gridLayout_2->addWidget(m_SpinX, 0, 0, 1, 1);

        m_SpinY = new QDoubleSpinBox(m_Frame);
        m_SpinY->setObjectName(QString::fromUtf8("m_SpinY"));
        m_SpinY->setDecimals(2);

        gridLayout_2->addWidget(m_SpinY, 0, 1, 1, 1);

        m_SpinZ = new QDoubleSpinBox(m_Frame);
        m_SpinZ->setObjectName(QString::fromUtf8("m_SpinZ"));
        m_SpinZ->setDecimals(2);

        gridLayout_2->addWidget(m_SpinZ, 0, 2, 1, 1);

        m_SpinT = new QSpinBox(m_Frame);
        m_SpinT->setObjectName(QString::fromUtf8("m_SpinT"));

        gridLayout_2->addWidget(m_SpinT, 0, 3, 1, 1);

        m_LabelX = new QLabel(m_Frame);
        m_LabelX->setObjectName(QString::fromUtf8("m_LabelX"));

        gridLayout_2->addWidget(m_LabelX, 1, 0, 1, 1);

        m_LabelY = new QLabel(m_Frame);
        m_LabelY->setObjectName(QString::fromUtf8("m_LabelY"));

        gridLayout_2->addWidget(m_LabelY, 1, 1, 1, 1);

        m_LabelZ = new QLabel(m_Frame);
        m_LabelZ->setObjectName(QString::fromUtf8("m_LabelZ"));

        gridLayout_2->addWidget(m_LabelZ, 1, 2, 1, 1);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_2->addItem(horizontalSpacer, 1, 3, 1, 1);


        horizontalLayout_2->addLayout(gridLayout_2);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_2->addItem(horizontalSpacer_2);


        verticalLayout_4->addLayout(horizontalLayout_2);

        horizontalLayout_6 = new QHBoxLayout();
        horizontalLayout_6->setSpacing(6);
        horizontalLayout_6->setObjectName(QString::fromUtf8("horizontalLayout_6"));
        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setSpacing(6);
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        label = new QLabel(m_Frame);
        label->setObjectName(QString::fromUtf8("label"));

        horizontalLayout->addWidget(label);

        m_SpinLevel = new QDoubleSpinBox(m_Frame);
        m_SpinLevel->setObjectName(QString::fromUtf8("m_SpinLevel"));
        m_SpinLevel->setDecimals(4);

        horizontalLayout->addWidget(m_SpinLevel);

        label_2 = new QLabel(m_Frame);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        horizontalLayout->addWidget(label_2);

        m_SpinWindow = new QDoubleSpinBox(m_Frame);
        m_SpinWindow->setObjectName(QString::fromUtf8("m_SpinWindow"));
        m_SpinWindow->setDecimals(4);
        m_SpinWindow->setMaximum(99);

        horizontalLayout->addWidget(m_SpinWindow);


        horizontalLayout_6->addLayout(horizontalLayout);

        horizontalSpacer_5 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_6->addItem(horizontalSpacer_5);


        verticalLayout_4->addLayout(horizontalLayout_6);

        horizontalLayout_5 = new QHBoxLayout();
        horizontalLayout_5->setSpacing(6);
        horizontalLayout_5->setObjectName(QString::fromUtf8("horizontalLayout_5"));
        label_3 = new QLabel(m_Frame);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        horizontalLayout_5->addWidget(label_3);

        m_LabelIntensity = new QLabel(m_Frame);
        m_LabelIntensity->setObjectName(QString::fromUtf8("m_LabelIntensity"));

        horizontalLayout_5->addWidget(m_LabelIntensity);

        horizontalSpacer_3 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_5->addItem(horizontalSpacer_3);


        verticalLayout_4->addLayout(horizontalLayout_5);

        horizontalLayout_7 = new QHBoxLayout();
        horizontalLayout_7->setSpacing(6);
        horizontalLayout_7->setObjectName(QString::fromUtf8("horizontalLayout_7"));
        pushButton = new QPushButton(m_Frame);
        pushButton->setObjectName(QString::fromUtf8("pushButton"));

        horizontalLayout_7->addWidget(pushButton);

        horizontalSpacer_4 = new QSpacerItem(40, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        horizontalLayout_7->addItem(horizontalSpacer_4);


        verticalLayout_4->addLayout(horizontalLayout_7);

        verticalSpacer = new QSpacerItem(20, 40, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout_4->addItem(verticalSpacer);


        gridLayout->addWidget(m_Frame, 1, 1, 1, 1);


        retranslateUi(ThreeViewClass);

        QMetaObject::connectSlotsByName(ThreeViewClass);
    } // setupUi

    void retranslateUi(QWidget *ThreeViewClass)
    {
        ThreeViewClass->setWindowTitle(QApplication::translate("ThreeViewClass", "ThreeView", 0, QApplication::UnicodeUTF8));
        m_Label1->setText(QApplication::translate("ThreeViewClass", "Z-Y", 0, QApplication::UnicodeUTF8));
        m_Label2->setText(QApplication::translate("ThreeViewClass", "Z-X", 0, QApplication::UnicodeUTF8));
        m_Label3->setText(QApplication::translate("ThreeViewClass", "X-Y", 0, QApplication::UnicodeUTF8));
        m_LabelX->setText(QApplication::translate("ThreeViewClass", "0.00", 0, QApplication::UnicodeUTF8));
        m_LabelY->setText(QApplication::translate("ThreeViewClass", "0.00", 0, QApplication::UnicodeUTF8));
        m_LabelZ->setText(QApplication::translate("ThreeViewClass", "0.00", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("ThreeViewClass", "Level:", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("ThreeViewClass", "Window:", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("ThreeViewClass", "Intensity:", 0, QApplication::UnicodeUTF8));
        m_LabelIntensity->setText(QApplication::translate("ThreeViewClass", "0.00", 0, QApplication::UnicodeUTF8));
        pushButton->setText(QApplication::translate("ThreeViewClass", "Toggle graph", 0, QApplication::UnicodeUTF8));
        Q_UNUSED(ThreeViewClass);
    } // retranslateUi

};

namespace Ui {
    class ThreeViewClass: public Ui_ThreeViewClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_RBMTHREEVIEW_H
