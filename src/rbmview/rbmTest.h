#ifndef RBMTEST_H_
#define RBMTEST_H_

#include <QObject>
#include <iostream>

namespace rbm
{

class Test : public QObject
{
  Q_OBJECT

public:
  Test()
  {
    m_Value = 0;
  }

  virtual ~Test()
  {

  }

public slots:
  void SetValue( double d )
  {
    d = ( double ) ( ( int ) d );
    bool changed = d != m_Value;
    std::cout << "Test: set value d=" << d << " " << ( changed ? "changed" : "" ) << std::endl;

    if ( changed )
      {
      m_Value = d;
      ValueUpdated( m_Value );
      }
  }

signals:
  void ValueUpdated( double d );

private:
  double m_Value;
};

}

#endif /* RBMTEST_H_ */
