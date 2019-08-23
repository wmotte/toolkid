#include "rbmView.h"

int main( int argc, char ** argv )
{
  QApplication app( argc, argv );
  rbm::View view;
  view.show();

  for( int i = 1; i < argc; ++i )
  	{
  	view.LoadImage( argv[ i ] );
  	}

  return app.exec();
}
