/********************************************************************************
** Form generated from reading ui file 'rbmView.ui'
**
** Created: Tue May 26 12:07:27 2009
**      by: Qt User Interface Compiler version 4.5.0
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_RBMVIEW_H
#define UI_RBMVIEW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QMainWindow>
#include <QtGui/QMenu>
#include <QtGui/QMenuBar>
#include <QtGui/QStatusBar>
#include <QtGui/QWidget>
#include "rbmThreeView.h"

QT_BEGIN_NAMESPACE

class Ui_ViewClass
{
public:
    QAction *action_Open;
    QAction *action_Quit;
    QAction *action_Three_planes;
    QAction *action_One_frame;
    QWidget *centralwidget;
    QHBoxLayout *horizontalLayout_2;
    rbm::ThreeView *widget;
    QMenuBar *menubar;
    QMenu *menu_File;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *ViewClass)
    {
        if (ViewClass->objectName().isEmpty())
            ViewClass->setObjectName(QString::fromUtf8("ViewClass"));
        ViewClass->resize(673, 523);
        ViewClass->setDockNestingEnabled(true);
        action_Open = new QAction(ViewClass);
        action_Open->setObjectName(QString::fromUtf8("action_Open"));
        action_Quit = new QAction(ViewClass);
        action_Quit->setObjectName(QString::fromUtf8("action_Quit"));
        action_Three_planes = new QAction(ViewClass);
        action_Three_planes->setObjectName(QString::fromUtf8("action_Three_planes"));
        action_One_frame = new QAction(ViewClass);
        action_One_frame->setObjectName(QString::fromUtf8("action_One_frame"));
        centralwidget = new QWidget(ViewClass);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        centralwidget->setGeometry(QRect(0, 26, 673, 474));
        horizontalLayout_2 = new QHBoxLayout(centralwidget);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        widget = new rbm::ThreeView(centralwidget);
        widget->setObjectName(QString::fromUtf8("widget"));

        horizontalLayout_2->addWidget(widget);

        ViewClass->setCentralWidget(centralwidget);
        menubar = new QMenuBar(ViewClass);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 673, 26));
        menu_File = new QMenu(menubar);
        menu_File->setObjectName(QString::fromUtf8("menu_File"));
        ViewClass->setMenuBar(menubar);
        statusbar = new QStatusBar(ViewClass);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        statusbar->setGeometry(QRect(0, 500, 673, 23));
        ViewClass->setStatusBar(statusbar);

        menubar->addAction(menu_File->menuAction());
        menu_File->addAction(action_Open);
        menu_File->addSeparator();
        menu_File->addAction(action_Quit);

        retranslateUi(ViewClass);
        QObject::connect(action_Quit, SIGNAL(triggered()), ViewClass, SLOT(close()));

        QMetaObject::connectSlotsByName(ViewClass);
    } // setupUi

    void retranslateUi(QMainWindow *ViewClass)
    {
        ViewClass->setWindowTitle(QApplication::translate("ViewClass", "rbmView", 0, QApplication::UnicodeUTF8));
        action_Open->setText(QApplication::translate("ViewClass", "&Open", 0, QApplication::UnicodeUTF8));
        action_Quit->setText(QApplication::translate("ViewClass", "&Quit", 0, QApplication::UnicodeUTF8));
        action_Quit->setShortcut(QApplication::translate("ViewClass", "Alt+Q", 0, QApplication::UnicodeUTF8));
        action_Three_planes->setText(QApplication::translate("ViewClass", "&Three planes", 0, QApplication::UnicodeUTF8));
        action_One_frame->setText(QApplication::translate("ViewClass", "&One frame", 0, QApplication::UnicodeUTF8));
        menu_File->setTitle(QApplication::translate("ViewClass", "&File", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class ViewClass: public Ui_ViewClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_RBMVIEW_H
