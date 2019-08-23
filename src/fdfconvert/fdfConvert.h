#ifndef __fdfConvert_h__
#define __fdfConvert_h__
#include "fdfDictionary.h"
#include <string>
#include <set>
#include "itkImage.h"

namespace fdf
{

/**
 * \class FDFConvert
 * \brief Read directory of multi-slice (multi-echo, arrayed) reconstructed Varian images (FDF format) and export to time-series data set
 *
 * Data acquired with Varian pulse sequences are automatically reconstructed by VNMRJ into FDF images.
 * For 2D acquisitions, reconstructed slices are stored in a directory with extension .img.
 * FDFConvert scans a directory for FDF-files according to the following naming convention:
 *   slice[xyz]image[xyz]echo[xyz].fdf
 * A series of 3D image volumes is obtained from the FDF-files, finally ordered by echo number.
 * Thus, the output of FDFConvert is a 4D series of image volumes image1echo1,image2echo1,...imageNecho1,image1echo2,...,imageNechoM.
 * Currently, FDFConvert does not support conversion of 3D acquisitions (e.g. GE3D or FSE3D).
 */
class FDFConvert
{
public:
  typedef itk::Image< float, 2 > SliceType;
  typedef itk::Image< float, 3 > ImageType;
  typedef itk::Image< float, 4 > SeriesType;
  
  typedef SliceType::Pointer SlicePointer;
  typedef ImageType::Pointer ImagePointer;
  typedef SeriesType::Pointer SeriesPointer;
  
  FDFConvert();

  /**
   * Read FDF files from a reconstruction directory.
   * From the filenames the number of images, slices, and echo's is calculated.
   */  
  void Read( const std::string& path );
  
  /**
   * Write a single 4D volume to file
   */
  void Write( const std::string& output );
  
  /**
   * Increase verbosity level
   */
  void SetVerbose( bool verbose );

protected:
  /**
   * Read a single slice FDF file. Upon first call, it will initialize m_Data and set appropriate spacing and origins.
   */
  void LoadSlice( const std::string& filename, int slice, int volume );
  
  /**
   * Read a header line from a 2D FDF file
   */
  bool ReadLine( const std::string& line, Dictionary& dictionary );

  int m_NumberOfSlices;
  int m_NumberOfEchos;
  int m_NumberOfImages;
  double m_MinZ;
  double m_MaxZ;  
  bool m_Verbose;
  SeriesPointer m_Data;
};

} // end namespace fdf

#endif /*__fdfConvert_h__*/
