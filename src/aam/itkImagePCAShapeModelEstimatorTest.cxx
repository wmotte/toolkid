/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkImagePCAShapeModelEstimatorTest.cxx,v $
  Language:  C++
  Date:      $Date: 2008/02/03 04:05:34 $
  Version:   $Revision: 1.11 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif
// Insight classes
#include "itkImage.h"
#include "itkVector.h"
#include "vnl/vnl_matrix_fixed.h"
#include "vnl/vnl_math.h"
#include "itkImageRegionIterator.h"
#include "itkLightProcessObject.h"
#include "itkTextOutput.h"

#include "itkImagePCAShapeModelEstimator.h"

//Data definitions 
#define   IMGWIDTH            2
#define   IMGHEIGHT           2
#define   NDIMENSION          2
#define   NUMTRAINIMAGES      3
#define   NUMLARGESTPC        2

// class to support progress feeback


class ShowProgressObject
{
public:
  ShowProgressObject(itk::LightProcessObject * o)
    {m_Process = o;}
  void ShowProgress()
    {std::cout << "Progress " << m_Process->GetProgress() << std::endl;}
  itk::LightProcessObject::Pointer m_Process;
};

/**
 * Main.
 */
int main(int, char* [] )
{

  itk::OutputWindow::SetInstance(itk::TextOutput::New().GetPointer());

  //------------------------------------------------------
  //Create 3 simple test images with
  //------------------------------------------------------
  typedef itk::Image<double,NDIMENSION> InputImageType; 
  typedef itk::Image<double,NDIMENSION> OutputImageType; 
  typedef itk::Image<double,NDIMENSION> MeanImageType;


  typedef InputImageType::PixelType ImagePixelType;

  typedef InputImageType::PixelType InputImagePixelType;

  typedef
    itk::ImageRegionIterator< InputImageType > InputImageIterator;

  typedef
    itk::ImageRegionIterator< OutputImageType > OutputImageIterator;
  
  InputImageType::Pointer image1 = InputImageType::New();

  InputImageType::Pointer image2 = InputImageType::New();

  InputImageType::Pointer image3 = InputImageType::New();

  InputImageType::SizeType inputImageSize = {{ IMGWIDTH, IMGHEIGHT }};

  InputImageType::IndexType index;
  index.Fill(0);
  InputImageType::RegionType region;

  region.SetSize( inputImageSize );
  region.SetIndex( index );

  //--------------------------------------------------------------------------
  // Set up Image 1 first
  //--------------------------------------------------------------------------

  image1->SetLargestPossibleRegion( region );
  image1->SetBufferedRegion( region );
  image1->Allocate();

  // setup the iterators
  InputImageIterator image1It( image1, image1->GetBufferedRegion() );

  //--------------------------------------------------------------------------
  // Set up Image 2 first
  //--------------------------------------------------------------------------

  image2->SetLargestPossibleRegion( region );
  image2->SetBufferedRegion( region );
  image2->Allocate();

  // setup the iterators
  InputImageIterator image2It( image2, image2->GetBufferedRegion() );

  //--------------------------------------------------------------------------
  // Set up Image 3 first
  //--------------------------------------------------------------------------

  image3->SetLargestPossibleRegion( region );
  image3->SetBufferedRegion( region );
  image3->Allocate();

  // setup the iterators
  InputImageIterator image3It( image3, image3->GetBufferedRegion() );

  //--------------------------------------------------------------------------
  //Manually create and store each vector
  //--------------------------------------------------------------------------
  //Image no. 1
  for( int i = 0; i< 4; i++ )
    {
    image1It.Set( 1 ); ++image1It;
    }
  //Image no. 2
  image2It.Set( 2 ); ++image2It;
  image2It.Set( 0 ); ++image2It;
  image2It.Set( 0 ); ++image2It;
  image2It.Set( 2 ); ++image2It;

  //Image no. 3
  image3It.Set( 0 ); ++image3It;
  image3It.Set( 3 ); ++image3It;
  image3It.Set( 3 ); ++image3It;
  image3It.Set( 0 ); ++image3It;

  //----------------------------------------------------------------------
  // Test code for the Shape model estimator
  //----------------------------------------------------------------------

  //----------------------------------------------------------------------
  //Set the image model estimator
  //----------------------------------------------------------------------
  typedef itk::ImagePCAShapeModelEstimator<InputImageType, OutputImageType> 
    ImagePCAShapeModelEstimatorType;

  ImagePCAShapeModelEstimatorType::Pointer 
    applyPCAShapeEstimator = ImagePCAShapeModelEstimatorType::New();

  //----------------------------------------------------------------------
  //Set the parameters of the clusterer
  //----------------------------------------------------------------------
  applyPCAShapeEstimator->SetNumberOfTrainingImages( NUMTRAINIMAGES );
  applyPCAShapeEstimator->SetNumberOfPrincipalComponentsRequired( NUMLARGESTPC + 1 );
  applyPCAShapeEstimator->SetNumberOfPrincipalComponentsRequired( NUMLARGESTPC );
  applyPCAShapeEstimator->SetInput(0, image1);
  applyPCAShapeEstimator->SetInput(1, image2);
  applyPCAShapeEstimator->SetInput(2, image3);

  applyPCAShapeEstimator->Update();

  //Test the printself function to increase coverage
  applyPCAShapeEstimator->Print(std::cout);

  //Exercise TypeMacro in superclass
  typedef ImagePCAShapeModelEstimatorType::Superclass GenericEstimatorType;
  std::cout << applyPCAShapeEstimator->GenericEstimatorType::GetNameOfClass() << std::endl;

  //Print out the number of training images and the number of principal 
  //components
  std::cout << "The number of training images are: " <<
    applyPCAShapeEstimator->GetNumberOfTrainingImages() << std::endl;

  std::cout << "The number of principal components desired are: " <<
    applyPCAShapeEstimator->GetNumberOfPrincipalComponentsRequired() << std::endl;

  //Print the eigen vectors
  vnl_vector<double> eigenValues = 
    applyPCAShapeEstimator->GetEigenValues();
  unsigned int numEigVal =  eigenValues.size();
  std::cout << "Number of returned eign-values: " << numEigVal << std::endl;

  std::cout << "The " << 
    applyPCAShapeEstimator->GetNumberOfPrincipalComponentsRequired() << 
    " largest eigen values are:" << std::endl;

  for(unsigned int i= 0; i< vnl_math_min( numEigVal, (unsigned int)NUMLARGESTPC ); i++ )
    {
    std::cout << eigenValues[ i ] << std::endl; 
    }  
  std::cout << "" << std::endl;
  std::cout << "" << std::endl;

  
  //Print the MeanImage
  OutputImageType::Pointer outImage = applyPCAShapeEstimator->GetOutput( 0 );
  OutputImageIterator outImageIt( outImage, outImage->GetBufferedRegion() );
  outImageIt.GoToBegin();

  std::cout << "The mean image is:" << std::endl;
  while(!outImageIt.IsAtEnd() )
    {
    std::cout << (double)(outImageIt.Get()) << ";"  << std::endl;  
    ++outImageIt; 
    } 
  std::cout << "  " << std::endl;

  //Print the largest two eigen vectors
  for (unsigned int j=1; j< NUMLARGESTPC + 1; j++ )
    {
    OutputImageType::Pointer outImage2 = applyPCAShapeEstimator->GetOutput( j );
    OutputImageIterator outImage2It( outImage2, outImage2->GetBufferedRegion() );
    outImage2It.GoToBegin();

    std::cout << "" << std::endl;
    std::cout << "The eigen vector number: " << j << " is:" << std::endl;
    while(!outImage2It.IsAtEnd() )
      {
      std::cout << (double) (outImage2It.Get()) << ";"  << std::endl;  
      ++outImage2It; 
      } 
    std::cout << "  " << std::endl;

    }

  //Test for the eigen values for the test case precomputed using Matlab/Splus
  std::cout << "" << std::endl;
  if( (eigenValues[2] < 6 || eigenValues[2] > 6.1) || (eigenValues[1] >0.1) )
    {
    std::cout<< "Test Passed" << std::endl;
    }
  else
    {
    std::cout<< "Test failed" << std::endl;
    return EXIT_FAILURE;
    }


  return EXIT_SUCCESS;
}
