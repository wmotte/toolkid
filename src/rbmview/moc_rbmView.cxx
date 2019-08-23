/****************************************************************************
** Meta object code from reading C++ file 'rbmView.h'
**
** Created: Tue May 26 12:07:29 2009
**      by: The Qt Meta Object Compiler version 61 (Qt 4.5.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "rbmView.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'rbmView.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 61
#error "This file was generated using the moc from 4.5.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_rbm__View[] = {

 // content:
       2,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   12, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors

 // slots: signature, parameters, type, tag, flags
      17,   11,   10,   10, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_rbm__View[] = {
    "rbm::View\0\0point\0PlotPointSelected(QwtDoublePoint&)\0"
};

const QMetaObject rbm::View::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_rbm__View,
      qt_meta_data_rbm__View, 0 }
};

const QMetaObject *rbm::View::metaObject() const
{
    return &staticMetaObject;
}

void *rbm::View::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_rbm__View))
        return static_cast<void*>(const_cast< View*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int rbm::View::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: PlotPointSelected((*reinterpret_cast< QwtDoublePoint(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 1;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
