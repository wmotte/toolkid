/****************************************************************************
** Meta object code from reading C++ file 'rbmThreeView.h'
**
** Created: Tue May 26 12:07:30 2009
**      by: The Qt Meta Object Compiler version 61 (Qt 4.5.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "rbmThreeView.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'rbmThreeView.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 61
#error "This file was generated using the moc from 4.5.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_rbm__ThreeView[] = {

 // content:
       2,       // revision
       0,       // classname
       0,    0, // classinfo
       7,   12, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors

 // slots: signature, parameters, type, tag, flags
      18,   16,   15,   15, 0x0a,
      41,   39,   15,   15, 0x0a,
      64,   62,   15,   15, 0x0a,
      87,   85,   15,   15, 0x0a,
     112,  105,   15,   15, 0x0a,
     139,  133,   15,   15, 0x0a,
     159,   15,   15,   15, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_rbm__ThreeView[] = {
    "rbm::ThreeView\0\0x\0UpdateValueX(double)\0"
    "y\0UpdateValueY(double)\0z\0UpdateValueZ(double)\0"
    "t\0UpdateValueT(int)\0window\0"
    "UpdateWindow(double)\0level\0"
    "UpdateLevel(double)\0ViewGraph()\0"
};

const QMetaObject rbm::ThreeView::staticMetaObject = {
    { &QWidget::staticMetaObject, qt_meta_stringdata_rbm__ThreeView,
      qt_meta_data_rbm__ThreeView, 0 }
};

const QMetaObject *rbm::ThreeView::metaObject() const
{
    return &staticMetaObject;
}

void *rbm::ThreeView::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_rbm__ThreeView))
        return static_cast<void*>(const_cast< ThreeView*>(this));
    return QWidget::qt_metacast(_clname);
}

int rbm::ThreeView::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: UpdateValueX((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 1: UpdateValueY((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 2: UpdateValueZ((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 3: UpdateValueT((*reinterpret_cast< int(*)>(_a[1]))); break;
        case 4: UpdateWindow((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 5: UpdateLevel((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 6: ViewGraph(); break;
        default: ;
        }
        _id -= 7;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
