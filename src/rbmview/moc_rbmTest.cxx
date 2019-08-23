/****************************************************************************
** Meta object code from reading C++ file 'rbmTest.h'
**
** Created: Tue May 26 12:07:27 2009
**      by: The Qt Meta Object Compiler version 61 (Qt 4.5.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "rbmTest.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'rbmTest.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 61
#error "This file was generated using the moc from 4.5.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_rbm__Test[] = {

 // content:
       2,       // revision
       0,       // classname
       0,    0, // classinfo
       2,   12, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors

 // signals: signature, parameters, type, tag, flags
      13,   11,   10,   10, 0x05,

 // slots: signature, parameters, type, tag, flags
      34,   11,   10,   10, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_rbm__Test[] = {
    "rbm::Test\0\0d\0ValueUpdated(double)\0"
    "SetValue(double)\0"
};

const QMetaObject rbm::Test::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_rbm__Test,
      qt_meta_data_rbm__Test, 0 }
};

const QMetaObject *rbm::Test::metaObject() const
{
    return &staticMetaObject;
}

void *rbm::Test::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_rbm__Test))
        return static_cast<void*>(const_cast< Test*>(this));
    return QObject::qt_metacast(_clname);
}

int rbm::Test::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: ValueUpdated((*reinterpret_cast< double(*)>(_a[1]))); break;
        case 1: SetValue((*reinterpret_cast< double(*)>(_a[1]))); break;
        default: ;
        }
        _id -= 2;
    }
    return _id;
}

// SIGNAL 0
void rbm::Test::ValueUpdated(double _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
