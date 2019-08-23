#include "tkdCmdParser.h"
#include "vnl/algo/vnl_fft_1d.h"
#include "vnl/vnl_matrix_ref.h"
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageRegionIteratorWithIndex.h"

#define logMacro( x ) \
{ \
  ::std::ostringstream ss; \
  ss << x; \
  this->LogMessage( ss.str() ); \
}

class RSFilter
{
public:
  typedef float PixelType;
  typedef itk::Image< PixelType, 4 > ImageType;
  typedef itk::ImageMaskSpatialObject< 3 > MaskType;
  typedef MaskType::ImageType MaskImageType;
  typedef itk::ImageFileReader< ImageType > ReaderType;
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typedef itk::ImageFileReader< MaskImageType > MaskReaderType;
  typedef vnl_matrix_ref< PixelType > ImageMatrixType;
  typedef vnl_vector< PixelType > VectorType;
  typedef vnl_matrix< PixelType > MatrixType;
  typedef itk::RecursiveGaussianImageFilter< ImageType, ImageType > SmoothingType;

  RSFilter( bool verbose = false, const std::string& logFileName = "", bool allVoxels = false )
    : m_Verbose( verbose ), m_AllVoxels( allVoxels )
  {
    if ( logFileName != "" )
      {
      m_Log.open( logFileName.c_str() );
      }
  }

  void LogMessage( const std::string& text )
  {
    if ( this->m_Log.is_open() )
      {
      this->m_Log << text << std::endl;
      this->m_Log.flush();
      }

    if ( this->m_Verbose )
      {
      std::cout << text << std::endl;
      }
  }

  inline void Swap( PixelType& a, PixelType& b ) const
  {
    PixelType tmp = a;
    a = b;
    b = tmp;
  }

  /** (c) Press et al., Numerical Recipes 3rd ed., p.612 **/
  void FFT( PixelType* data, const int& n, const int& isign )
  {
    int nn, mmax, m, j, istep, i;
    PixelType wtemp, wr, wpr, wpi, wi, theta, tempr, tempi, norm;

    nn = n << 1;
    j = 1;

    for( i = 1; i < nn; i += 2 )
      {
      if ( j > i )
        {
        this->Swap( data[ j - 1 ], data[ i - 1 ] );
        this->Swap( data[ j ], data[ i ] );
        }
      m = n;
      while( m >= 2 && j > m )
        {
        j -= m;
        m >>= 1;
        }
      j += m;
      }

    mmax = 2;

    while( nn > mmax )
      {
      istep = mmax << 1;
      theta = static_cast< PixelType >( isign ) * ( 6.28318530717959 / mmax );
      wtemp = sin( 0.5 * theta );
      wpr = -2.0 * wtemp * wtemp;
      wpi = sin( theta );
      wr = 1.0;
      wi = 0.0;

      for( m = 1; m < mmax; m +=2 )
        {
        for( i = m; i <= nn; i += istep )
          {
          j = i + mmax;
          tempr = wr * data[ j - 1 ] - wi * data[ j ];
          tempi = wr * data[ j ] + wi * data[ j - 1 ];
          data[ j - 1 ] = data[ i - 1 ] - tempr;
          data[ j ] = data[ i ] - tempi;
          data[ i - 1 ] += tempr;
          data[ i ] += tempi;
          }
        wtemp = wr;
        wr = wtemp * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
        }
      mmax = istep;
      }

    if ( isign == -1 )
      {
      norm = 1. / n;
      for( i = 0; i < nn; ++i )
        {
        data[ i ] *= norm;
        }
      }
  }

  /** (c) Press et al., Numerical Recipes 3rd ed., p.784 **/
  void FitLine( const VectorType& data, int ignore, PixelType& a, PixelType& b, PixelType& mean )
  {
    int numberOfPoints = data.size();

    PixelType sx = 0;
    PixelType sy = 0;
    PixelType st2 = 0;

    for( int i = ignore; i < ( numberOfPoints - ignore ); ++i )
      {
      sx += static_cast< PixelType >( i );
      sy += data( i );
      }

    PixelType sxoss = sx / static_cast< PixelType >( numberOfPoints - 2 * ignore);

    for( int i = ignore; i < ( numberOfPoints - ignore ); ++i )
      {
      PixelType t = static_cast< PixelType >( i ) - sxoss;
      st2 += t * t;
      a += t * data( i );
      }
    a /= st2;
    b = ( sy - sx * a ) / static_cast< PixelType >( numberOfPoints - 2 * ignore );
    mean = sy / static_cast< PixelType >( numberOfPoints - 2 * ignore );
  }

  void RemoveTrend( VectorType& signal, int ignore = 0, bool addMean = false )
  {
    PixelType a, b, mean;
    this->FitLine( signal, ignore, a, b, mean );

    const int n = signal.size();
    for( int i = 0; i < n; ++i )
      {
      signal( i ) = signal( i ) - a * static_cast< PixelType >( i ) - b;
      if ( addMean )
        {
        signal( i ) += mean;
        }
      }
  }

  void ReadImage( const std::string& filename )
  {
    logMacro( "Read " << filename );
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( filename.c_str() );
    reader->Update();

    m_Image = reader->GetOutput();
    reader = 0;

    m_Region = m_Image->GetLargestPossibleRegion();
    m_Size = m_Region.GetSize();
    m_Spacing = m_Image->GetSpacing();
    m_Origin = m_Image->GetOrigin();

    m_NumberOfVoxels = m_Size[ 0 ] * m_Size[ 1 ] * m_Size[ 2 ];
    m_NumberOfTimepoints = m_Size[ 3 ];

    logMacro( "Number of voxels: " << m_NumberOfVoxels );
    logMacro( "Number of time points: " << m_NumberOfTimepoints );
  }

  void BuildMask( const std::string& filename = "" )
  {
    logMacro( "Build mask. File=" << filename );

    // build mask
    m_MaskImage = MaskImageType::New();
    for( int i = 0; i < 3; ++i )
      {
      m_MaskSize[ i ] = m_Size[ i ];
      m_MaskOrigin[ i ] = m_Origin[ i ];
      m_MaskSpacing[ i ] = m_Spacing[ i ];
      }
    m_MaskImage->SetSpacing( m_MaskSpacing );
    m_MaskImage->SetOrigin( m_MaskOrigin );
    m_MaskRegion.SetSize( m_MaskSize );
    m_MaskImage->SetRegions( m_MaskRegion );
    m_MaskImage->Allocate();
    m_MaskImage->FillBuffer( 0 );

    // read mask
    if ( filename != "" )
      {
      int insideCount = 0;
      int outsideCount = 0;
      MaskReaderType::Pointer maskReader = MaskReaderType::New();
      maskReader->SetFileName( filename.c_str() );
      maskReader->Update();

      MaskType::Pointer mask = MaskType::New();
      mask->SetImage( maskReader->GetOutput() );

      itk::ImageRegionIteratorWithIndex< MaskImageType > it( m_MaskImage, m_MaskRegion );
      MaskImageType::PointType point;
      for( ; !it.IsAtEnd(); ++it )
        {
        m_MaskImage->TransformIndexToPhysicalPoint( it.GetIndex(), point );
        if ( mask->IsInside( point ) )
          {
          ++insideCount;
          it.Set( 255 );
          }
        else
          {
          ++outsideCount;
          }
        }

      logMacro( "Mask voxels: " << insideCount << "/" << ( insideCount + outsideCount ) );
      }
    else
      {
      m_MaskImage->FillBuffer( 255 );
      }

    m_Mask = MaskType::New();
    m_Mask->SetImage( m_MaskImage );
    m_MaskBuffer = m_MaskImage->GetPixelContainer()->GetBufferPointer();
  }

  PixelType* GetImageBuffer()
  {
    return m_Image->GetPixelContainer()->GetBufferPointer();
  }

  void IntensityNormalization( PixelType target = 1000.0 )
  {
    logMacro( "Intensity normalization. Target=" << target );

    const int& numberOfVoxels = m_NumberOfVoxels;
    const int& numberOfTimepoints = m_NumberOfTimepoints;
    ImageMatrixType matrix( numberOfTimepoints, numberOfVoxels, this->GetImageBuffer() );
    const unsigned char* mask = m_MaskBuffer;

    std::vector< PixelType > list;

    for( int i = 0; i < numberOfVoxels; ++i )
      {
      if ( !mask[ i ] )
        {
        continue;
        }

      for( int j = 0; j < numberOfTimepoints; ++j )
        {
        list.push_back( matrix( j, i ) );
        }
      }

    std::sort( list.begin(), list.end() );

    PixelType medianIntensity = list[ list.size() / 2 ];

    logMacro( "Mode intensity: " << medianIntensity );
    logMacro( "Scaling factor: " << ( target / medianIntensity ) );

    for( int i = 0; i < numberOfVoxels; ++i )
      {
      for( int j = 0; j < numberOfTimepoints; ++j )
        {
        matrix( j, i ) *= ( target / medianIntensity );
        }
      }
  }

  void LinearDetrend( int ignorePoints = 0, bool addMean = false )
  {
    logMacro( "Linear de-trend. Ignore " << ignorePoints << " points" << ( addMean ? ", add mean back." : "." ) );

    const int& numberOfVoxels = m_NumberOfVoxels;
    const int& numberOfTimepoints = m_NumberOfTimepoints;
    ImageMatrixType matrix( numberOfTimepoints, numberOfVoxels, this->GetImageBuffer() );

    for( int i = 0; i < numberOfVoxels; ++i )
      {
      VectorType column = matrix.get_column( i );
      this->RemoveTrend( column, ignorePoints, addMean );
      matrix.set_column( i, column );
      }
  }

  void Smooth( PixelType fwhm = 1.0, bool smoothXY = false )
  {
    PixelType sigma = fwhm / vcl_sqrt( 8. * vcl_log( 2. ) );

    logMacro( "Smooth. FWHM=" << fwhm << ", sigma=" << sigma );
    if ( smoothXY )
      {
      logMacro( "Only smooth in X-Y plane" );
      }

    for( int i = 0; i < 3; ++i )
      {
      if ( i == 2 && smoothXY )
        {
        break;
        }

      SmoothingType::Pointer smoother = SmoothingType::New();
      smoother->SetSigma( sigma );
      smoother->SetDirection( i );
      smoother->SetOrder( SmoothingType::ZeroOrder );
      smoother->SetInput( m_Image );
      smoother->Update();

      ImageType::Pointer image = smoother->GetOutput();
      smoother = 0;
      m_Image = image;
      }
  }

  void Filter( const PixelType& highPass, const PixelType& lowPass, float sampleSpacing = 0 )
  {
    const int& numberOfVoxels = m_NumberOfVoxels;
    const int& numberOfTimepoints = m_NumberOfTimepoints;
    ImageMatrixType matrix( numberOfTimepoints, numberOfVoxels, this->GetImageBuffer() );

    if ( sampleSpacing <= 0 )
      {
      sampleSpacing = m_Spacing[ 3 ];
      }

    logMacro( "Filter [" << highPass << ", " << lowPass << "] Hz, sample spacing=" << sampleSpacing << " s" );

    PixelType sampleFactor = static_cast< PixelType >( numberOfTimepoints ) * sampleSpacing;
    int center = round( static_cast< PixelType >( numberOfTimepoints ) * 0.5 );
    int lowerIndex1 = static_cast< int >( vcl_ceil( highPass * sampleFactor ) );
    int higherIndex1 = static_cast< int >( vcl_floor( lowPass * sampleFactor ) );

    if ( higherIndex1 > center )
      {
      higherIndex1 = center;
      }

    if ( lowerIndex1 < 1 )
      {
      lowerIndex1 = 1;
      }

    int lowerIndex2 = numberOfTimepoints - higherIndex1;
    int higherIndex2 = numberOfTimepoints - lowerIndex1;

//    std::cout << "[" << lowerIndex1 << ", " << higherIndex1 << "], [" << lowerIndex2 << ", " << higherIndex2 << "] = [" << highPass << ", " << lowPass << "] Hz" << std::endl;

    vnl_fft_1d< PixelType > fft( numberOfTimepoints );

    vnl_vector< vcl_complex< PixelType > > buffer( numberOfTimepoints );
    buffer.fill( itk::NumericTraits< vcl_complex< PixelType > >::Zero );

    for( int i = 0; i < numberOfVoxels; ++i )
      {
      for( int j = 0; j < numberOfTimepoints; ++j )
        {
        buffer( j ) = vcl_complex< PixelType >( matrix( j, i ), 0 );
        }

      fft.fwd_transform( buffer );

      for( int j = 1; j < numberOfTimepoints; ++j )
        {
        if ( j >= lowerIndex1 && j <= higherIndex1 )
          {
          continue;
          }

        if ( j >= lowerIndex2 && j <= higherIndex2 )
          {
          continue;
          }

        buffer( j ) = 0;
        }

      fft.bwd_transform( buffer );

      for( int j = 0; j < numberOfTimepoints; ++j )
        {
        const vcl_complex< PixelType >& p = buffer( j );
        matrix( j, i ) = p.real();
        }
      }
  }

  void Filter2( const PixelType& highPass, const PixelType& lowPass, float sampleSpacing = 0 )
  {
    const int& numberOfVoxels = m_NumberOfVoxels;
    const int& numberOfTimepoints = m_NumberOfTimepoints;
    ImageMatrixType matrix( numberOfTimepoints, numberOfVoxels, this->GetImageBuffer() );

    if ( sampleSpacing <= 0 )
      {
      sampleSpacing = m_Spacing[ 3 ];
      }

    int bufferLength = pow( 2, vcl_ceil( vcl_log( numberOfTimepoints ) / vcl_log( 2. ) ) );
    PixelType sampleFactor = static_cast< PixelType >( bufferLength ) * sampleSpacing;
    int lowerIndex1 = static_cast< int >( vcl_ceil( highPass * sampleFactor ) );
    int higherIndex1 = static_cast< int >( vcl_floor( lowPass * sampleFactor ) );
    int lowerIndex2 = bufferLength - higherIndex1 - 1;
    int higherIndex2 = bufferLength - lowerIndex1 - 1;

//    std::cout << "[" << lowerIndex1 << ", " << higherIndex1 << "], [" << lowerIndex2 << ", " << higherIndex2 << "] = [" << highPass << ", " << lowPass << "] Hz" << std::endl;

    vnl_vector< vcl_complex< PixelType > > buffer( bufferLength );
    buffer.fill( itk::NumericTraits< vcl_complex< PixelType > >::Zero );

    for( int i = 0; i < numberOfVoxels; ++i )
      {
      for( int j = 0; j < numberOfTimepoints; ++j )
        {
        buffer( j ) = vcl_complex< PixelType >( matrix( j, i ), 0 );
        }

      for( int j = numberOfTimepoints; j < bufferLength; ++j )
        {
        buffer( j ) = 0;
        }

      this->FFT( reinterpret_cast< PixelType* >( buffer.data_block() ), bufferLength, 1 );

      for( int j = 0; j < bufferLength; ++j )
        {
        if ( j >= lowerIndex1 && j <= higherIndex1 )
          {
          continue;
          }

        if ( j >= lowerIndex2 && j <= higherIndex2 )
          {
          continue;
          }

        buffer( j ) = 0;
        }

      this->FFT( reinterpret_cast< PixelType* >( buffer.data_block() ), bufferLength, 1 );

      for( int j = 0; j < numberOfTimepoints; ++j )
        {
        const vcl_complex< PixelType >& p = buffer( j );
        matrix( j, i ) = p.real();
        }
      }
  }

  void GlobalSignalRegression()
  {
    logMacro( "Global signal regression" );

    const int& numberOfVoxels = m_NumberOfVoxels;
    const int& numberOfTimepoints = m_NumberOfTimepoints;
    ImageMatrixType matrix( numberOfTimepoints, numberOfVoxels, this->GetImageBuffer() );
    const unsigned char* mask = m_MaskBuffer;

    // determine mean and remove linear trend
    VectorType mean( numberOfTimepoints );
    mean.fill( 0 );
    int count = 0;
    for( int i = 0; i < numberOfVoxels; ++i )
      {
      if ( mask[ i ] == 0 )
        {
        continue;
        }

      mean += matrix.get_column( i );
      ++count;
      }

    // calculate mean signal
    mean /= static_cast< PixelType >( count );

    // regression on all voxels
    if ( m_AllVoxels )
      {
      count = numberOfVoxels;
      }

    // build masked input matrix
    MatrixType B( numberOfTimepoints, count );
    count = 0;
    for( int i = 0; i < numberOfVoxels; ++i )
      {
      if ( mask[ i ] == 0 )
        {
        if ( m_AllVoxels )
          {
          ++count;
          }
        continue;
        }

      B.set_column( count++, matrix.get_column( i ) );
      }

    // regression (cf. Fox 2009, Appendix)
    ImageMatrixType g( numberOfTimepoints, 1, mean.data_block() );
    MatrixType gPlus = g.transpose() * ( 1.0 / dot_product( mean, mean ) );
    MatrixType BetaG = gPlus * B;
    MatrixType Bg = B - ( g * BetaG );

    // reconstruct
    count = 0;
    matrix.fill( 0 );
    for( int i = 0; i < numberOfVoxels; ++i )
      {
      if ( !m_AllVoxels && mask[ i ] == 0 )
        {
        continue;
        }

      matrix.set_column( i, Bg.get_column( count++ ) );
      }
  }

  void Write( const std::string& filename )
  {
    logMacro( "Write " << filename );
    WriterType::Pointer writer = WriterType::New();
    writer->SetInput( m_Image );
    writer->SetFileName( filename.c_str() );
    writer->Update();
  }

protected:
  ImageType::Pointer m_Image;
  int m_NumberOfVoxels;
  int m_NumberOfTimepoints;
  ImageType::RegionType m_Region;
  ImageType::SizeType m_Size;
  ImageType::SpacingType m_Spacing;
  ImageType::PointType m_Origin;
  MaskImageType::RegionType m_MaskRegion;
  MaskImageType::SizeType m_MaskSize;
  MaskImageType::SpacingType m_MaskSpacing;
  MaskImageType::PointType m_MaskOrigin;
  MaskType::Pointer m_Mask;
  MaskImageType::Pointer m_MaskImage;
  unsigned char* m_MaskBuffer;
  std::ofstream m_Log;
  bool m_Verbose;
  bool m_AllVoxels;
};

int main( int argc, char ** argv )
{
  typedef RSFilter::PixelType PixelType;

  tkd::CmdParser p( "rsfilter", "Intensity normalization, smoothing, filtering, global signal regression for resting state fMRI" );

  std::string inputFileName, outputFileName, maskFileName;
  PixelType lowPass = 0;
  PixelType highPass = 0;
  PixelType sampleSpacing = 0;
  PixelType fwhm = 0;
  PixelType intensityTarget = 0;
  bool addMean = false;
  bool smoothXY = false;
  bool noDetrend = false;
  bool noGlobalSignalRegression = false;
  int ignorePoints = 0;
  bool verbose = false;
  std::string logFileName;
  bool allVoxels = false;

  p.AddArgument( inputFileName, "input" )
    ->AddAlias( "i" )
    ->SetDescription( "4D input image" )
    ->SetRequired( true );

  p.AddArgument( outputFileName, "output" )
    ->AddAlias( "o" )
    ->SetDescription( "4D filtered output image" )
    ->SetRequired( true );

  p.AddArgument( maskFileName, "mask" )
    ->AddAlias( "m" )
    ->SetDescription( "3D mask image" );

  p.AddArgument( lowPass, "lowpass" )
    ->AddAlias( "lp" )
    ->SetDescription( "Low pass (Hz)" );

  p.AddArgument( highPass, "highpass" )
    ->AddAlias( "hp" )
    ->SetDescription( "High pass (Hz)" );

  p.AddArgument( sampleSpacing, "sample-spacing" )
    ->AddAlias( "ss" )
    ->SetDescription( "Sample spacing (in seconds; default: 4th voxel dimension)" );

  p.AddArgument( addMean, "add-mean" )
    ->AddAlias( "am" )
    ->SetDescription( "Add mean back after linear de-trending" );

  p.AddArgument( fwhm, "fwhm" )
    ->SetDescription( "Gaussian smoothing Full Width at Half Maximum" );

  p.AddArgument( intensityTarget, "ing" )
    ->SetDescription( "Global intensity normalization, target value" );

  p.AddArgument( smoothXY, "smooth-xy" )
    ->AddAlias( "sxy" )
    ->SetDescription( "Perform smoothin in X-Y plane only" );

  p.AddArgument( noDetrend, "no-detrend" )
    ->AddAlias( "nd" )
    ->SetDescription( "No linear de-trending" );

  p.AddArgument( noGlobalSignalRegression, "no-global-regression" )
    ->AddAlias( "ng" )
    ->SetDescription( "No global signal regression" );

  p.AddArgument( ignorePoints, "ignore" )
    ->AddAlias( "ig" )
    ->SetDescription( "Ignore samples at start and end during trend estimation (default: 0)" );

  p.AddArgument( verbose, "verbose" )
    ->AddAlias( "v" )
    ->SetDescription( "Log to screen" );

  p.AddArgument( logFileName, "log" )
    ->AddAlias( "l" )
    ->SetDescription( "Log to file" );

  p.AddArgument( allVoxels, "all-voxels" )
    ->AddAlias( "all" )
    ->SetDescription( "Apply transformations to all voxels" );

  if ( !p.Parse( argc, argv ) )
    {
    p.PrintUsage( std::cout );
    return -1;
    }

  RSFilter filter( verbose, logFileName, allVoxels );

  // read

  filter.ReadImage( inputFileName );
  filter.BuildMask( maskFileName );

  // intensity normalization

  if ( intensityTarget > 0 )
    {
    filter.IntensityNormalization( intensityTarget );
    }

  // linear de-trend

  if ( !noDetrend )
    {
    filter.LinearDetrend( ignorePoints, addMean );
    }

  // smooth

  if ( fwhm > 0 )
    {
    filter.Smooth( fwhm, smoothXY );
    }

  // filter

  if ( lowPass > 0 )
    {
    filter.Filter( highPass, lowPass, sampleSpacing );
    }

  // global signal regression

  if ( !noGlobalSignalRegression )
    {
    filter.GlobalSignalRegression();
    }

  // write

  filter.Write( outputFileName );

  return 0;
}
