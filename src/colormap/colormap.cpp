#include "colormap.h"
#include "itkImageLinearIteratorWithIndex.h"

namespace colormap
{

void ColorMap::Run( const parameters& list )
{
  // [ 1 ] Check input ...
  Checks( list );

  // [ 2 ] Read images ...
  ReadImages( list.inputFileNames );

  // [ 3 ] Read mask ...
  ReadMask( list.maskFileName );

  // [ 4 ] Read background ...
  ReadBackground( list.backgroundFileName );

  // [ 9 ] Resample such that all images have the same dimension, and enforce isotropic in-plane resolution
  if ( m_Background )
    {
    ResampleImages( m_Background, list.sliceDirection, list.scaleFactor, list.interpolationOrder );
    }
  else if ( m_Images.size() > 0 )
    {
    ResampleImages( m_Images[ 0 ], list.sliceDirection, list.scaleFactor, list.interpolationOrder );
    }

  // [ 6 ] Crop ...
  CropAndMask( !list.noMask );

  // [ 7 ] Threshold images ...
  VectorType lowerThreshold( list.overlaysLowerThreshold );
  VectorType upperThreshold( list.overlaysUpperThreshold );

  if ( list.overlaysLowerThreshold.empty() )
  {
    lowerThreshold = GetMinValues();
  }

  if ( list.overlaysUpperThreshold.empty() )
  {
    upperThreshold = GetMaxValues();
  }

  // [ 8 ] Rescale thresholded images and background (if specified) between 0 - 255...
  LevelSetBackground( list.backgroundMin, list.backgroundMax, list.percentile );
  RescaleImages( lowerThreshold, upperThreshold );

  // [ 10 ] Select LUTS ...
  std::vector< std::string > LUTS = list.LUTS;
  if ( list.LUTS.empty() )
    LUTS = SelectLuts( list.inputFileNames.size() );

  // [ 11 ] Generate RGB images ...
  GenerateRGB( LUTS );

  // [ 12 ] Merge background with black empty image ...
  RGBImageType::Pointer merged = MergeBackground( list.backgroundTransparency, list.defaultR, list.defaultG, list.defaultB );

  // [ 13 ] Merge RGB images ...
  merged = MergeRGBImages( merged, list.overlayTransparencies );

  // [ 14 ] Permute to put slices in Z-direction ...
  merged = PermuteImage( merged, list.sliceDirection );

  std::vector< int > selectedSlices = list.selectedSlices;
  if ( selectedSlices.size() == 0 )
  {
    selectedSlices = CalculateSlices( merged, list.skipSlices, list.startSlice, list.endSlice );
  }

  // [ 15 ] Create overlay image ...
  OverlayImageType::Pointer overlay = Overlay( list.numberOfColumns, selectedSlices, merged );

  // [ 16 ] Flip
  overlay = Flip( overlay, true, true );

  // [17] add color bars
  if ( list.addColorBars )
  {
    overlay = AddColorBars( overlay, LUTS, lowerThreshold, upperThreshold, list.colorBarWidth, list.colorBarHorizontal );
  }

  // [ 17 ] Write output to file ...
  WriteOverlay( overlay, list.outputFileName );
}

ColorMap::RGBImageType::Pointer ColorMap::PermuteImage( RGBImageType::Pointer image, unsigned int sliceDirection )
{
  PermuteFilterType::PermuteOrderArrayType permutation;

  switch( sliceDirection )
  {
    case 0:
    {
      permutation[ 0 ] = 2;
      permutation[ 1 ] = 1;
      permutation[ 2 ] = 0;
      break;
    }
    case 1:
    {
      permutation[ 0 ] = 0;
      permutation[ 1 ] = 2;
      permutation[ 2 ] = 1;
      break;
    }

    default:
    {
      return image;
    }
  }

  PermuteFilterType::Pointer filter = PermuteFilterType::New();
  filter->SetInput( image );
  filter->SetOrder( permutation );
  filter->Update();

  return filter->GetOutput();
}


std::vector< int > ColorMap::CalculateSlices( RGBImageType::Pointer image, int skip, int start, int end )
{
  std::vector< int > selectedSlices;

//  image->Print( std::cout );

  int totalNumberOfSlices = image->GetLargestPossibleRegion().GetSize()[ 2 ];
  start = start == -1 ? 0 : start;
  end = end == -1 ? totalNumberOfSlices - 1 : end;

//  std::cout << "Number of slices: " << totalNumberOfSlices << std::endl;
//  std::cout << "Start: " << start << ", end: " << end << std::endl;

  if ( start < 0 || end < 0 || start >= totalNumberOfSlices || end >= totalNumberOfSlices || start > end )
    {
    start = 0;
    end = totalNumberOfSlices - 1;
    }

  if ( skip < 1 )
    {
    skip = 1;
    }

  for( int i = start; i <= end; i += skip )
    {
    selectedSlices.push_back( i );
//    std::cout << "Select with skip=" << skip << ": " << i << std::endl;
    }

  return selectedSlices;
}

void ColorMap::LevelSetBackground( PixelType min, PixelType max, double percentile )
{
  if ( !m_UseBackground )
    {
    return;
    }

  double actualMin, actualMax;
  this->Percentile( m_Background, m_Mask, percentile, actualMin, actualMax );

  if ( min == itk::NumericTraits< double >::NonpositiveMin() )
    {
    min = actualMin;
    }

  if ( max == itk::NumericTraits< double >::max() )
    {
    max = actualMax;
    }

  IntensityWindowingImageFilterType::Pointer filter = IntensityWindowingImageFilterType::New();
  filter->SetInput( m_Background );
  filter->SetOutputMinimum( 0 );
  filter->SetOutputMaximum( 255 );
  filter->SetWindowMinimum( min );
  filter->SetWindowMaximum( max );
  filter->Update();
  m_Background = filter->GetOutput();
}

ColorMap::RGBImageType::Pointer ColorMap::MergeBackground( PixelType transparency, int defaultR, int defaultG, int defaultB )
{
  if( ( transparency < 0 ) || ( transparency > 1 ) )
  {
    std::cerr << "*** ERROR ***: Background transparency should be between 0 and 1!" << std::endl;
    exit( EXIT_FAILURE );
  }


  if( m_RGBImages.empty() && !m_UseBackground )
  {
    std::cerr << "*** ERROR ***: Could not merge background, with no rgb images in memory!" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Allocate ...
  RGBImageType::Pointer merged = RGBImageType::New();
  if ( m_RGBImages.empty() )
    {
    merged->CopyInformation( m_Background );
    merged->SetRegions( m_Background->GetLargestPossibleRegion() );
    }
  else
    {
    merged->CopyInformation( m_RGBImages.at( 0 ) );
    merged->SetRegions( m_RGBImages.at( 0 )->GetLargestPossibleRegion() );
    }

  merged->Allocate();
  RGBPixelType zero;
  
  zero[0] = defaultR; zero[1] = defaultG; zero[2] = defaultB;
  merged->FillBuffer( zero ); // Black background image...

  std::cout<<zero<<std::endl;
  if( m_UseBackground )
  {
    RGBIteratorType it( merged, merged->GetLargestPossibleRegion() );
    IteratorType bit( m_Background, m_Background->GetLargestPossibleRegion() );

    for( it.GoToBegin(), bit.GoToBegin(); ! it.IsAtEnd(), ! bit.IsAtEnd(); ++it, ++bit )
    {
      PixelType backgroundValue = m_Background->GetPixel( bit.GetIndex() );

      if( backgroundValue != 0 )
      {
        RGBPixelType pixel;
        pixel[0] = pixel[1] = pixel[2] = transparency * backgroundValue;
        merged->SetPixel( it.GetIndex(), pixel );
      }
    }
  }

  return merged;
}

void ColorMap::Percentile( ImageType::Pointer input, ImageType::Pointer mask, double percentile, double& lower, double& upper )
{
  itk::ImageRegionConstIterator< ImageType > itInput( input, input->GetLargestPossibleRegion() );
  itk::ImageRegionConstIterator< ImageType > itMask;

  if ( mask )
  {
    itMask = itk::ImageRegionConstIterator< ImageType > ( mask, mask->GetLargestPossibleRegion() );
  }

  std::vector< double > values;
  for( ; !itInput.IsAtEnd(); ++itInput )
  {
    if ( mask )
    {
      if ( itMask.Value() >= 1 )
      {
        values.push_back( itInput.Value() );
      }

      ++itMask;
      continue;
    }
    else
    {
      values.push_back( itInput.Value() );
    }
  }

  std::sort( values.begin(), values.end() );

  int indexLower = ( 1. - percentile / 100. ) * static_cast< double >( values.size() - 1. );
  int indexUpper = ( percentile / 100. ) * static_cast< double >( values.size() - 1. );

  lower = values[ indexLower ];
  upper = values[ indexUpper ];
}

ColorMap::ImageType::Pointer ColorMap::MaskImage( ImageType::Pointer input, ImageType::Pointer mask )
{
  ImageType::Pointer masked = ImageType::New();
  masked->CopyInformation( input );
  masked->SetRegions( input->GetLargestPossibleRegion() );
  masked->Allocate();

  itk::ImageRegionConstIterator< ImageType > itInput( input, input->GetLargestPossibleRegion() );
  itk::ImageRegionConstIterator< ImageType > itMask( mask, input->GetLargestPossibleRegion() );
  itk::ImageRegionIterator< ImageType > itOutput( masked, masked->GetLargestPossibleRegion() );

  for( ; !itInput.IsAtEnd(); ++itInput, ++itMask, ++itOutput )
    {
    if ( itMask.Value() < 1 )
      {
      itOutput.Set( 0 );
      }
    else
      {
      itOutput.Set( itInput.Get() );
      }
    }

  return masked;
}

void ColorMap::CropAndMask( bool applyMask )
{
  if( ! m_UseMask )
    {
    return;
    }

  ImageType::RegionType extractRegion;
  IteratorType it( m_Mask, m_Mask->GetLargestPossibleRegion() );

  ImageType::IndexType minimumIndex;
  ImageType::IndexType maximumIndex;
  minimumIndex[0] = itk::NumericTraits< ImageType::IndexType::IndexValueType >::max();
  minimumIndex[1] = itk::NumericTraits< ImageType::IndexType::IndexValueType >::max();
  minimumIndex[2] = itk::NumericTraits< ImageType::IndexType::IndexValueType >::max();

  maximumIndex[0] = itk::NumericTraits< ImageType::IndexType::IndexValueType >::min();
  maximumIndex[1] = itk::NumericTraits< ImageType::IndexType::IndexValueType >::min();
  maximumIndex[2] = itk::NumericTraits< ImageType::IndexType::IndexValueType >::min();

  for( it.GoToBegin(); ! it.IsAtEnd(); ++it )
  {
    ImageType::IndexType index = it.GetIndex();

    if( it.Value() < 1 )
    {
      continue;
    }

    if( index[0] < minimumIndex[0] )
      minimumIndex[0] = index[0];
    if( index[1] < minimumIndex[1] )
      minimumIndex[1] = index[1];
    if( index[2] < minimumIndex[2] )
      minimumIndex[2] = index[2];

    if( index[0] > maximumIndex[0] )
      maximumIndex[0] = index[0];
    if( index[1] > maximumIndex[1] )
      maximumIndex[1] = index[1];
    if( index[2] > maximumIndex[2] )
      maximumIndex[2] = index[2];
  }

  ImageType::SizeType size;

  size[0] = maximumIndex[0] - minimumIndex[0] + 1;
  size[1] = maximumIndex[1] - minimumIndex[1] + 1;
  size[2] = maximumIndex[2] - minimumIndex[2] + 1;

  extractRegion.SetSize( size );
  extractRegion.SetIndex( minimumIndex );

//  std::cout << extractRegion << std::endl;

  for( unsigned int i = 0; i < m_Images.size(); i++ )
  {
    ExtractImageFilterType::Pointer filter = ExtractImageFilterType::New();
    ImageType::Pointer image = m_Images.at( i );
    if ( applyMask )
    {
      image = MaskImage( image, m_Mask );
    }

    filter->SetInput( image );
    filter->SetExtractionRegion( extractRegion );
    filter->Update();
    m_Images.at( i ) = filter->GetOutput();
  }

  if( m_UseBackground )
  {
    ExtractImageFilterType::Pointer filter = ExtractImageFilterType::New();
    ImageType::Pointer image = m_Background;
    if ( applyMask )
    {
      image = MaskImage( image, m_Mask );
    }
    filter->SetInput( image );
    filter->SetExtractionRegion( extractRegion );
    filter->Update();
    m_Background = filter->GetOutput();
  }

  if( m_Mask )
  {
    ExtractImageFilterType::Pointer filter = ExtractImageFilterType::New();
    filter->SetInput( m_Mask );
    filter->SetExtractionRegion( extractRegion );
    filter->Update();
    m_Mask = filter->GetOutput();
  }
}

ColorMap::OverlayImageType::Pointer ColorMap::Flip( const OverlayImageType::Pointer& overlay, bool x, bool y )
{
  FlipImageFilterType::Pointer filter = FlipImageFilterType::New();
  filter->SetInput( overlay );

  itk::FixedArray< bool, 2 > axes;
  axes[ 0 ] = x;
  axes[ 1 ] = y;

  filter->SetFlipAxes( axes );
  filter->Update();

  return filter->GetOutput();
}

ColorMap::OverlayImageType::Pointer ColorMap::Overlay(
    unsigned int nCols,
    std::vector< int > selectedSlices,
    const RGBImageType::Pointer& image3D
    )
{
  RGBImageType::SizeType size3D = image3D->GetLargestPossibleRegion().GetSize();
  unsigned int numberOfSlices = selectedSlices.size();

  // checks...
  if( nCols > numberOfSlices )
  {
    nCols = numberOfSlices;
  }

  if( nCols < 1 )
  {
    nCols = 1;
  }

  int nRows = ( numberOfSlices / nCols );
  if ( numberOfSlices % nCols > 0 )
  {
    ++nRows;
  }

  // Allocate ...
  OverlayImageType::Pointer overlay = OverlayImageType::New();
  OverlayImageType::RegionType overlayRegion;
  OverlayImageType::SizeType overlaySize;
  overlaySize[ 0 ] = size3D[ 0 ] * nCols + 20 * m_LUTImages.size();
  overlaySize[ 1 ] = size3D[ 1 ] * nRows;
  overlayRegion.SetSize( overlaySize );
  overlay->SetRegions( overlayRegion );
  overlay->Allocate();
  RGBPixelType zero; zero[ 0 ] = 0; zero[ 1 ] = 0; zero[ 2 ] = 0;
  overlay->FillBuffer( zero );

  // ITERATE...
  RGBIteratorType it( image3D, image3D->GetLargestPossibleRegion( ) );

  OverlayImageType::RegionType outputRegion = overlayRegion;
  OverlayImageType::SizeType outputSize = outputRegion.GetSize();
  OverlayImageType::IndexType outputIndex = outputRegion.GetIndex();

  RGBImageType::RegionType inputRegion = image3D->GetLargestPossibleRegion();
  RGBImageType::IndexType inputIndex = inputRegion.GetIndex();
  RGBImageType::SizeType inputSize = inputRegion.GetSize();

  outputSize[ 0 ] = size3D[ 0 ];
  outputSize[ 1 ] = size3D[ 1 ];

  outputIndex[ 0 ] = 0;
  outputIndex[ 1 ] = 0;

  inputSize[ 2 ] = 1;
  int slice = 0;

  int startIndex = inputIndex[ 2 ];
  for( int row = 0; row < nRows; ++row )
  {
    for( unsigned int col = 0; col < nCols && slice < numberOfSlices; ++col )
    {
      inputIndex[ 2 ] = startIndex + selectedSlices[ slice++ ];

      outputIndex[ 0 ] = ( nCols - col - 1 ) * outputSize[ 0 ];
      outputIndex[ 1 ] = ( nRows - row - 1 ) * outputSize[ 1 ];

      inputRegion.SetSize( inputSize );
      inputRegion.SetIndex( inputIndex );

      outputRegion.SetSize( outputSize );
      outputRegion.SetIndex( outputIndex );

      itk::ImageRegionConstIterator< RGBImageType > itInput( image3D, inputRegion );
      itk::ImageRegionIterator< OverlayImageType > itOutput( overlay, outputRegion );

      for( ; !itOutput.IsAtEnd(); ++itInput, ++itOutput )
      {
        itOutput.Set( itInput.Get() );
      }
    }
  }

  return overlay;
}

void ColorMap::WriteOverlay( const OverlayImageType::Pointer& image, const std::string& outputFileName )
{
  OverlayImageFileWriterType::Pointer writer = OverlayImageFileWriterType::New();
  writer->SetFileName( outputFileName );
  writer->SetInput( image );

  try
  {
    writer->Update();
  } catch ( itk::ExceptionObject& e )
  {
    std::cerr << "*** ERROR ***: Could not write " << outputFileName << "!" << std::endl;
    exit( EXIT_FAILURE );
  }
}

ColorMap::RGBImageType::Pointer ColorMap::MergeRGBImages( const RGBImageType::Pointer& background,
    const VectorType& overlayTransparencies )
{
  if ( m_RGBImages.empty() )
  {
    return background;
  }

  VectorType transparencies( overlayTransparencies );

  if( transparencies.empty() )
    transparencies = VectorType( m_RGBImages.size(), 1.0 );

  if( transparencies.size() != m_RGBImages.size() )
  {
    std::cerr << "*** ERROR ***: Number of overlay transparency values is not equal to number of overlay images!" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Allocate ...
  RGBImageType::Pointer merged = background;

  for( unsigned int i = 0; i < m_RGBImages.size(); i++ )
  {
    PixelType trans = transparencies.at( i );
    if ( ( trans < 0 ) || ( trans > 1 ) )
    {
      std::cerr << "*** ERROR ***: Overlay transparency should be between 0 and 1!" << std::endl;
      exit( EXIT_FAILURE );
    }

    RGBImageType::Pointer image = m_RGBImages.at( i );
    RGBIteratorType it( image, image->GetLargestPossibleRegion() );

    for ( it.GoToBegin(); !it.IsAtEnd(); ++it )
    {
      RGBPixelType value = image->GetPixel( it.GetIndex() );
      RGBPixelType backgroundValue = merged->GetPixel( it.GetIndex() );

      if( ( value[0] != 0 ) || ( value[1] != 0 ) || ( value[2] != 0 ) )
      {
        RGBPixelType newValue;

        newValue[0] = backgroundValue[0] *
                 ( 1. - trans )
                 + value[0] * trans;

        newValue[1] = backgroundValue[1] *
                 ( 1. - trans )
                 + value[1] * trans;

        newValue[2] = backgroundValue[2] *
                 ( 1. - trans )
                 + value[2] * trans;

        merged->SetPixel( it.GetIndex(), newValue );
      }
    }
  }

  return merged;
}

void ColorMap::GenerateRGB( const std::vector< std::string >& LUTS )
{
  if ( LUTS.size() != m_Images.size() )
  {
    std::cerr << "*** ERROR ***: Could not generate colormaps, "
      "number of images does not match number of LUTS!" << std::endl;
    exit( EXIT_FAILURE );
  }

  for ( unsigned int i = 0; i < m_Images.size(); i++ )
    {
    m_RGBImages.push_back( GenerateRGB( LUTS.at( i ), m_Images.at( i ) ) );
    }
}

ColorMap::OverlayImageType::Pointer ColorMap::AddColorBars( OverlayImageType::Pointer input, const std::vector< std::string >& LUTS, const VectorType& lowerThreshold, const VectorType& upperThreshold, int width, bool horizontal )
{
  if ( LUTS.size() != m_Images.size() )
  {
    std::cerr << "*** ERROR ***: Could not generate colormaps, "
      "number of images does not match number of LUTS!" << std::endl;
    exit( EXIT_FAILURE );
  }

  OverlayImageType::Pointer output = OverlayImageType::New();
  output->CopyInformation( input );

  OverlayImageType::RegionType inputRegion = input->GetLargestPossibleRegion();
  OverlayImageType::SizeType inputSize = inputRegion.GetSize();
  OverlayImageType::IndexType inputIndex = inputRegion.GetIndex();
  OverlayImageType::RegionType outputRegion = inputRegion;
  OverlayImageType::SizeType outputSize = inputSize;
  OverlayImageType::IndexType outputIndex = outputRegion.GetIndex();

  int firstIndex = horizontal ? 1 : 0;
  int secondIndex = horizontal ? 0 : 1;

  outputSize[ firstIndex ] += m_Images.size() * ( width + 1 ) + 1;

  outputRegion.SetSize( outputSize );
  output->SetRegions( outputRegion );
  output->Allocate();
  output->FillBuffer( itk::NumericTraits< RGBPixelType >::max() );

  itk::ImageRegionConstIterator< OverlayImageType > itInput( input, inputRegion );
  itk::ImageRegionIterator< OverlayImageType > itOutput( output, inputRegion );
  for( ; !itOutput.IsAtEnd(); ++itInput, ++itOutput )
    {
    itOutput.Set( itInput.Get() );
    }

  outputSize[ firstIndex ] = width;
  outputSize[ secondIndex ] -= 2;
  outputIndex[ firstIndex ] = inputSize[ firstIndex ] + 1;
  outputIndex[ secondIndex ] = 1;

  for ( unsigned int i = 0; i < m_Images.size(); i++ )
    {
    ImageType::Pointer image = ImageType::New();
    ImageType::RegionType region;
    ImageType::SizeType size;
    ImageType::IndexType index;

    size.Fill( 1 );
    index.Fill( 0 );

    size[ 0 ] = outputSize[ 0 ];
    size[ 1 ] = outputSize[ 1 ];

    region.SetSize( size );
    region.SetIndex( index );

    image->SetRegions( region );
    image->Allocate();

    itk::ImageLinearIteratorWithIndex< ImageType > it( image, region );
    it.SetDirection( firstIndex );
    int counter = 0;
    for( ; !it.IsAtEnd(); it.NextLine(), ++counter )
      {
      PixelType position = static_cast< PixelType >( counter ) / static_cast< PixelType >( size[ firstIndex ] - 1 );
      if ( !horizontal )
        {
        position = 1. - position;
        }

      for( ; !it.IsAtEndOfLine(); ++it )
        {
        it.Set( position * ( upperThreshold[ i ] - lowerThreshold[ i ] ) + lowerThreshold[ i ] );
        }
      }

    RGBImageType::Pointer lutImage = GenerateRGB( LUTS.at( i ), image );

    outputRegion.SetIndex( outputIndex );
    outputRegion.SetSize( outputSize );

    itk::ImageRegionConstIterator< RGBImageType > itLUT( lutImage, region );
    itk::ImageRegionIterator< OverlayImageType > itOutput( output, outputRegion );

    for( ; !itLUT.IsAtEnd(); ++itLUT, ++itOutput )
      {
      itOutput.Set( itLUT.Get() );
      }

    outputIndex[ firstIndex ] += ( width + 1 );
    }

  return output;
}

ColorMap::ColormapType::Pointer ColorMap::GetCustomColorMap( const std::string& name )
{
    std::stringstream ss;
    ss << "/home1/wim/workspace/toolkid/src/colormap/" << name << ".txt";

  std::ifstream str( ss.str().c_str() );

  std::string line;

    if( str.fail() )
    {
      std::cerr << "*** ERROR ***: Custom colormap, could not read from: " << ss.str() << "!" << std::endl;
      exit( EXIT_FAILURE );
    }

    ColormapType::Pointer colormap = ColormapType::New();

    // Get red values
    {
      std::getline( str, line );
      std::istringstream iss( line );
      float value;
      ColormapType::ChannelType channel;
      while ( iss >> value )
      {
        channel.push_back( value );
      }
      colormap->SetRedChannel( channel );
    }

    // Get green values
    {
      std::getline( str, line );
      std::istringstream iss( line );
      float value;
      ColormapType::ChannelType channel;
      while ( iss >> value )
      {
        channel.push_back( value );
      }
      colormap->SetGreenChannel( channel );
    }

    // Get blue values
    {
      std::getline( str, line );
      std::istringstream iss( line );
      float value;
      ColormapType::ChannelType channel;
      while ( iss >> value )
      {
        channel.push_back( value );
      }
      colormap->SetBlueChannel( channel );
    }

    return colormap;
}

ColorMap::RGBImageType::Pointer ColorMap::GenerateRGB( const std::string& colormapString, ImageType::Pointer image )
{
  RGBFilterType::Pointer rgbfilter = RGBFilterType::New();

  rgbfilter->SetInput( image );

  if ( colormapString == "red" )
  {
    rgbfilter->SetColormap( RGBFilterType::Red );
  } else if ( colormapString == "green" )
  {
    rgbfilter->SetColormap( RGBFilterType::Green );
  } else if ( colormapString == "blue" )
  {
    rgbfilter->SetColormap( RGBFilterType::Blue );
  } else if ( colormapString == "grey" )
  {
    rgbfilter->SetColormap( RGBFilterType::Grey );
  } else if ( colormapString == "cool" )
  {
    rgbfilter->SetColormap( RGBFilterType::Cool );
  } else if ( colormapString == "hot" )
  {
    rgbfilter->SetColormap( RGBFilterType::Hot );
  } else if ( colormapString == "spring" )
  {
    rgbfilter->SetColormap( RGBFilterType::Spring );
  } else if ( colormapString == "autumn" )
  {
    rgbfilter->SetColormap( RGBFilterType::Autumn );
  } else if ( colormapString == "winter" )
  {
    rgbfilter->SetColormap( RGBFilterType::Winter );
  } else if ( colormapString == "copper" )
  {
    rgbfilter->SetColormap( RGBFilterType::Copper );
  } else if ( colormapString == "summer" )
  {
    rgbfilter->SetColormap( RGBFilterType::Summer );
  } else if ( colormapString == "jet" )
  {
    rgbfilter->SetColormap( RGBFilterType::Jet );
  } else if ( colormapString == "hsv" )
  {
    rgbfilter->SetColormap( RGBFilterType::HSV );
  } else if ( colormapString == "overunder" )
  {
    rgbfilter->SetColormap( RGBFilterType::OverUnder );
  }
  else if ( colormapString == "bone" )
  {
    rgbfilter->SetColormap( GetCustomColorMap( "bone" ) );
  }
  else if ( colormapString == "flag" )
  {
    rgbfilter->SetColormap( GetCustomColorMap( "flag" ) );
  }
  else if ( colormapString == "lines" )
  {
    rgbfilter->SetColormap( GetCustomColorMap( "lines" ) );
  }
  else if ( colormapString == "prism" )
  {
    rgbfilter->SetColormap( GetCustomColorMap( "prism" ) );
  }
  else if ( colormapString == "vga" )
  {
    rgbfilter->SetColormap( GetCustomColorMap( "vga" ) );
  }
  else if ( colormapString == "colorcube" )
  {
    rgbfilter->SetColormap( GetCustomColorMap( "colorcube" ) );
  }
  else
  {
    std::cerr << "*** ERROR ***: No valid colormap found for: " << colormapString << "!" << std::endl;
    std::cerr << "Options are: \"grey\", \"red\", \"green\", \"blue\", \"copper\", \"jet\", \"hsv\", \"spring\", \"summer\", \"autumn\", \"winter\", \"hot\", \"cool\", "
    << "\"bone\", \"flag\", \"lines\", \"prism\", \"vga\", \"colorcube\"" << std::endl;

    exit( EXIT_FAILURE );
  }

  rgbfilter->GetColormap()->SetMinimumRGBComponentValue( 0 );
  rgbfilter->GetColormap()->SetMaximumRGBComponentValue( 255 );

  try
  {
    rgbfilter->Update();
  } catch ( ... )
  {
    std::cerr << "*** ERROR ***: Could not create RGB filter!" << std::endl;
    exit( EXIT_FAILURE );
  }

  // Mask pixel with color-value 0 (some LUTS paint zero-value too).
  MaskImageFilterType::Pointer maskFilter = MaskImageFilterType::New();
  maskFilter->SetInput1( rgbfilter->GetOutput() );
  maskFilter->SetInput2( image );
  maskFilter->Update();

  return maskFilter->GetOutput();
}

std::vector< std::string > ColorMap::SelectLuts( unsigned int total )
{
  std::vector< std::string > luts( total );

  std::vector< std::string > options;

  options.push_back( "hot" );
  options.push_back( "cool" );

  options.push_back( "red" );
  options.push_back( "green" );
  options.push_back( "blue" );


  options.push_back( "copper" );
  options.push_back( "jet" );
  options.push_back( "hsv" );

  options.push_back( "spring" );
  options.push_back( "summer" );
  options.push_back( "autumn" );
  options.push_back( "winter" );

  options.push_back( "bone" );
  options.push_back( "flag" );
  options.push_back( "lines" );
  options.push_back( "prism" );
  options.push_back( "vga" );
  options.push_back( "colorcube" );

  for ( unsigned int i = 0; i < luts.size(); i++ )
  {
    if ( i < options.size() )
      luts.at( i ) = options.at( i );
    else
      luts.at( i ) = options.at( i % options.size() );
  }

  return luts;
}

void ColorMap::RescaleImages( VectorType lower, VectorType upper )
{
  // scale between min-max
  for( unsigned int i = 0; i < m_Images.size(); i++ )
  {
    if ( lower[ i ] > upper[ i ] )
      {
      for( itk::ImageRegionIterator< ImageType > it( m_Images[ i ], m_Images[ i ]->GetLargestPossibleRegion() );
           !it.IsAtEnd(); ++it )
        {
        it.Value() *= -1;
        }

      lower[ i ] *= -1;
      upper[ i ] *= -1;
      }

    IntensityWindowingImageFilterType::Pointer filter = IntensityWindowingImageFilterType::New();
    filter->SetInput( m_Images.at( i ) );
    filter->SetOutputMinimum( 0 );
    filter->SetOutputMaximum( 255 );
    filter->SetWindowMinimum( lower[ i ] );
    filter->SetWindowMaximum( upper[ i ] );
    filter->Update();
    m_Images.at( i ) = filter->GetOutput();
  }
}

ColorMap::ImageType::Pointer ColorMap::ResampleImage( ImageType::Pointer image, ImageType::Pointer reference, int interpolationOrder )
{
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetSplineOrder( interpolationOrder );

  FilterType::Pointer filter = FilterType::New();
  filter->SetInput( image );
  filter->SetReferenceImage( reference );
  filter->SetInterpolator( interpolator );
  filter->SetUseReferenceImage( true );
  filter->Update();

  return filter->GetOutput();
}

void ColorMap::ResampleImages( ImageType::Pointer reference, unsigned int sliceDirection, double scaleFactor, int interpolationOrder )
{
  ImageType::SpacingType spacing = reference->GetSpacing();
  ImageType::RegionType region = reference->GetLargestPossibleRegion();
  ImageType::SizeType size = region.GetSize();
  ImageType::IndexType index = region.GetIndex();

  int direction[ 3 ];
  direction[ 0 ] = ( sliceDirection + 1 ) % 3;
  direction[ 1 ] = ( sliceDirection + 2 ) % 3;
  direction[ 2 ] = sliceDirection;

  if ( spacing[ direction[ 0 ] ] != spacing[ direction[ 1 ] ] || scaleFactor != 1 )
    {
    int smallestIndex = spacing[ direction[ 0 ] ] < spacing[ direction[ 1 ] ] ? direction[ 0 ] : direction[ 1 ];
    itk::FixedArray< double, 3 > scale;
    for( int i = 0; i < 2; ++i )
      {
      scale[ direction[ i ] ] = ( spacing[ direction[ i ] ] / spacing[ smallestIndex ] ) * scaleFactor;
      }
    scale[ sliceDirection ] = 1.;

    for( int i = 0; i < 3; ++i )
      {
      spacing[ i ] /= scale[ i ];
      size[ i ] = vcl_ceil( scale[ i ] * size[ i ] );
      }

    region.SetSize( size );

    ImageType::Pointer image = ImageType::New();
    image->CopyInformation( reference );
    image->SetRegions( region );
    image->SetSpacing( spacing );
    image->Allocate();

    reference = image;
    }

  for( unsigned int i = 0; i < m_Images.size(); ++i )
    {
    m_Images[ i ] = ResampleImage( m_Images[ i ], reference, interpolationOrder );
    }

  if ( m_Background )
    {
    m_Background = ResampleImage( m_Background, reference, interpolationOrder );
    }

  if ( m_Mask )
    {
    m_Mask = ResampleImage( m_Mask, reference, 0 );
    }
}

void ColorMap::ReadImages( const std::vector< std::string >& inputFileNames )
{
  for ( unsigned int i = 0; i < inputFileNames.size(); i++ )
  {
    ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
    reader->SetFileName( inputFileNames.at( i ) );
    reader->Update();
    m_Images.push_back( reader->GetOutput() );
  }
}

void ColorMap::ReadMask( const std::string& maskFileName )
{
  if ( maskFileName.empty() )
  {
    m_UseMask = false;
  } else
  {
    ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
    reader->SetFileName( maskFileName );
    reader->Update();
    m_Mask = reader->GetOutput();
    m_UseMask = true;
  }
}

/**
 * Read background.
 */
void ColorMap::ReadBackground( const std::string& backgroundFileName )
{
  if ( backgroundFileName.empty() )
  {
    m_UseBackground = false;
  } else
  {
    ImageFileReaderType::Pointer reader = ImageFileReaderType::New();
    reader->SetFileName( backgroundFileName );
    reader->Update();
    m_Background = reader->GetOutput();
    m_UseBackground = true;
  }
}


ColorMap::VectorType ColorMap::GetMinValues()
{
  VectorType minValues;

  for ( unsigned int i = 0; i < m_Images.size(); i++ )
  {
    MinimumMaximumImageFilterType::Pointer filter = MinimumMaximumImageFilterType::New();
    filter->SetInput( m_Images.at( i ) );
    filter->Update();
    minValues.push_back( filter->GetMinimum() );
  }
  return minValues;
}

ColorMap::VectorType ColorMap::GetMaxValues()
{
  VectorType maxValues;

  for ( unsigned int i = 0; i < m_Images.size(); i++ )
  {
    MinimumMaximumImageFilterType::Pointer filter = MinimumMaximumImageFilterType::New();
    filter->SetInput( m_Images.at( i ) );
    filter->Update();
    maxValues.push_back( filter->GetMaximum() );
  }
  return maxValues;
}

void ColorMap::Checks( const parameters& list )
{
  std::vector< std::string > checkFiles( list.inputFileNames );

  if( ! list.backgroundFileName.empty() )
    checkFiles.push_back( list.backgroundFileName );
  if( ! list.maskFileName.empty() )
    checkFiles.push_back( list.maskFileName );

  for ( unsigned int i = 0; i < checkFiles.size(); i++ )
  {
    std::string inputFileName = checkFiles.at( i );

    if ( graph::Graph< PixelType >::GetImageDimensions( inputFileName ) != 3 )
    {
      std::cout << "*** ERROR ***: Input \"" << inputFileName << "\" is not a 3D image!" << std::endl;
      exit( EXIT_FAILURE );
    }
  }

  if( ! list.overlaysLowerThreshold.empty()
      && ( list.overlaysLowerThreshold.size() != list.inputFileNames.size() ) )
  {
    std::cerr << "*** ERROR ***: Number of lower threshold values is not equal to number of input images!" << std::endl;
    exit( EXIT_FAILURE );
  }

  if( ! list.overlaysUpperThreshold.empty()
      && ( list.overlaysUpperThreshold.size() != list.inputFileNames.size() ) )
  {
    std::cerr << "*** ERROR ***: Number of upper threshold values is not equal to number of input images!" << std::endl;
    exit( EXIT_FAILURE );
  }

  if( ! list.overlayTransparencies.empty()
      && ( list.overlayTransparencies.size() != list.inputFileNames.size() ) )
  {
    std::cerr << "*** ERROR ***: Number of overlay transparencies is not equal to number of input images!" << std::endl;
    exit( EXIT_FAILURE );
  }
}

} // end namespace colormap


