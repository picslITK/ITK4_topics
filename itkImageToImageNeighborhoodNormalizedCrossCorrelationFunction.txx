/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkImageToImageNeighborhoodNormalizedCrossCorrelationFunction_txx
#define __itkImageToImageNeighborhoodNormalizedCrossCorrelationFunction_txx

#include "itkImageToImageNeighborhoodNormalizedCrossCorrelationFunction.h"
#include "itkNumericTraits.h"

namespace itk
{
/**
 * Constructor
 */
template< class TFixedImage, class TMovingImage >
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::ImageToImageNeighborhoodNormalizedCrossCorrelationFunction()
{
  this->m_Parameters.Fill( NumericTraits< ParametersValueType >::Zero );

  m_FixedImageGradientCalculator = FixedImageGradientCalculatorType::New();
  m_MovingImageGradientCalculator = MovingImageGradientCalculatorType::New();

}

template< class TFixedImage, class TMovingImage >
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::~ImageToImageNeighborhoodNormalizedCrossCorrelationFunction()
{
}

/**
 * PrintSelf
 */
template< class TFixedImage, class TMovingImage >
void
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Parameters: " << m_Parameters << std::endl;
}

template< class TFixedImage, class TMovingImage >
unsigned int
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::GetNumberOfParameters() const { return 5; }

template< class TFixedImage, class TMovingImage >
typename ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::MeasureType
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::GetValue( const ParametersType & parameters ) const
{
  this->m_Parameters = parameters;
  return 1.0;
}


template< class TFixedImage, class TMovingImage >
void
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::GetDerivative(const ParametersType &,
                           DerivativeType & derivative) const { derivative.Fill(0.0); }

template< class TFixedImage, class TMovingImage >
void
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::Initialize(void) throw ( itk::ExceptionObject ) {}



template< class TFixedImage, class TMovingImage >
void
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::IntializeGradientCalculator()
 {
    // GS: test: don't know if this is thread safe?
    // setup gradient calculator


    m_FixedImageGradientCalculator->SetInputImage( m_FixedImage );
    m_MovingImageGradientCalculator->SetInputImage( m_MovingImage  );

 }




template< class TFixedImage, class TMovingImage >
void
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::IntializeScanning(const ImageRegionType &scan_region, ScanningIteratorType &scan_it, ScanMemType &scan_mem, ScanParaType &scan_para)
 {
    scan_para.scan_region = scan_region;
    scan_para.I = this->m_FixedImage;
    scan_para.J = this->m_MovingImage;
    scan_para.r = this->GetRadius();

    int nb_fill_zero = scan_para.I->GetBufferedRegion().GetIndex(0) - (scan_region.GetIndex(0) - scan_para.r[0]);
    if (nb_fill_zero<0) nb_fill_zero=0;
    scan_para.number_of_fill_zero = nb_fill_zero;

    scan_it = ScanningIteratorType(scan_para.r, scan_para.I, scan_region);
    scan_para.window_length = scan_it.Size();
    scan_para.scan_region_begin_index_dim0 = scan_it.GetBeginIndex()[0];



    // GS: test: don't know if this is thread safe?
    // setup gradient calculator
    // m_FixedImageGradientCalculator->SetInputImage( Superclass::m_FixedImage );
    // m_MovingImageGradientCalculator->SetInputImage( Superclass::m_MovingImage  );
 }


template< class TFixedImage, class TMovingImage >
void
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::UpdateQueuesAtBeginingOfLine(const ScanningIteratorType &scan_it, ScanMemType &scan_mem, const ScanParaType &scan_para)
 {

    int nb_fill_zero = scan_para.number_of_fill_zero;
    unsigned int hoodlen = scan_para.window_length;
//    RadiusType r = this->GetRadius();

    scan_mem.Qsuma2 = SumQueueType( nb_fill_zero, 0.0 );
    scan_mem.Qsumb2 = SumQueueType( nb_fill_zero, 0.0 );
    scan_mem.Qsuma = SumQueueType( nb_fill_zero, 0.0 );
    scan_mem.Qsumb = SumQueueType( nb_fill_zero, 0.0 );
    scan_mem.Qsumab = SumQueueType( nb_fill_zero, 0.0 );
    scan_mem.Qcount = SumQueueType( nb_fill_zero, 0.0 );


    typename FixedImageType::IndexType oindex = scan_it.GetIndex();

    // Now add the rest of the values from each hyperplane
    for( unsigned int i = nb_fill_zero; i < ( 2*scan_para.r[0] + 1 ); i++ )
    {
        float suma2 = 0.0;
        float sumb2 = 0.0;
        float suma = 0.0;
        float sumb = 0.0;
        float sumab = 0.0;
        float count = 0.0;

        for( unsigned int indct = i; indct < hoodlen; indct += ( 2*scan_para.r[0] + 1 ) )
        {
            bool isInBounds = true;
            scan_it.GetPixel( indct, isInBounds );

            typename FixedImageType::IndexType index = scan_it.GetIndex( indct );

            // if ( !isInBounds || ( this->m_FixedImageMask && this->m_FixedImageMask->GetPixel( index ) < 0.25 ) )
            if ( !isInBounds)
            {
                // std::cout << "DEBUG: error" << std::endl;
                continue;
            }

            float a = scan_para.I->GetPixel( index );
            float b = scan_para.J->GetPixel( index );

            suma2 += a*a;
            sumb2 += b*b;
            suma += a;
            sumb += b;
            sumab += a*b;
            count += 1.0;
        }

        scan_mem.Qsuma2.push_back( suma2 );
        scan_mem.Qsumb2.push_back( sumb2 );
        scan_mem.Qsuma.push_back( suma );
        scan_mem.Qsumb.push_back( sumb );
        scan_mem.Qsumab.push_back( sumab );
        scan_mem.Qcount.push_back( count );

    }

 }


template< class TFixedImage, class TMovingImage >
void
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::UpdateQueuesToNextScanWindow(const ScanningIteratorType &scan_it, ScanMemType &scan_mem, const ScanParaType &scan_para)
 {
    const unsigned int hoodlen = scan_para.window_length;

    // Increment the iterator and check to see if we're at the end of the
    // line.  If so, go to the next line.  Otherwise, add the
    // the values for the next hyperplane.

    float suma2 = 0.0;
    float sumb2 = 0.0;
    float suma = 0.0;
    float sumb = 0.0;
    float sumab = 0.0;
    float count = 0.0;

    for( unsigned int indct = 2*scan_para.r[0]; indct < hoodlen; indct += ( 2*scan_para.r[0] + 1 ) )
    {
        bool isInBounds = true;
        scan_it.GetPixel( indct, isInBounds );
        typename FixedImageType::IndexType index = scan_it.GetIndex( indct );

        //GS: this is a weird notation: force prob mask < 0.25 to be zero
        // if ( !isInBounds || ( this->m_FixedImageMask && this->m_FixedImageMask->GetPixel( index ) < 0.25 ) )
        if ( !isInBounds )
        {
            continue;
        }

        float a = scan_para.I->GetPixel( index );
        float b = scan_para.J->GetPixel( index );

        suma2 += a*a;
        sumb2 += b*b;
        suma += a;
        sumb += b;
        sumab += a*b;
        count += 1.0;
    }

    scan_mem.Qsuma2.push_back( suma2 );
    scan_mem.Qsumb2.push_back( sumb2 );
    scan_mem.Qsuma.push_back( suma );
    scan_mem.Qsumb.push_back( sumb );
    scan_mem.Qsumab.push_back( sumab );
    scan_mem.Qcount.push_back( count );


    scan_mem.Qsuma2.pop_front();
    scan_mem.Qsumb2.pop_front();
    scan_mem.Qsuma.pop_front();
    scan_mem.Qsumb.pop_front();
    scan_mem.Qsumab.pop_front();
    scan_mem.Qcount.pop_front();

 }


template< class TFixedImage, class TMovingImage >
void
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::UpdateQueues(const ScanningIteratorType &scan_it, ScanMemType &scan_mem, const ScanParaType &scan_para)
 {

    if ( scan_it.GetIndex()[0] == scan_para.scan_region_begin_index_dim0) {
        UpdateQueuesAtBeginingOfLine(scan_it, scan_mem, scan_para);
    } else {
        UpdateQueuesToNextScanWindow(scan_it, scan_mem, scan_para);
    }

 }

template< class TFixedImage, class TMovingImage >
void
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::ComputeInformationFromQueues(const ScanningIteratorType &scan_it, ScanMemType &scan_mem, const ScanParaType &scan_para)
 {
    // Test to see if there are any voxels we need to handle in the current
    // window.

    float suma2 = 0.0;
    float sumb2 = 0.0;
    float suma = 0.0;
    float sumb = 0.0;
    float sumab = 0.0;
    float count = 0.0;

    typename SumQueueType::iterator itcount = scan_mem.Qcount.begin();
    while( itcount != scan_mem.Qcount.end() )
    {
        count += *itcount;
        ++itcount;
    }

    // If there are values, we need to calculate the different quantities
    if( count > 0 )
    {

        typename SumQueueType::iterator ita2 = scan_mem.Qsuma2.begin();
        typename SumQueueType::iterator itb2 = scan_mem.Qsumb2.begin();
        typename SumQueueType::iterator ita = scan_mem.Qsuma.begin();
        typename SumQueueType::iterator itb = scan_mem.Qsumb.begin();
        typename SumQueueType::iterator itab = scan_mem.Qsumab.begin();

        while( ita2 != scan_mem.Qsuma2.end() )
        {
            suma2 += *ita2;
            sumb2 += *itb2;
            suma += *ita;
            sumb += *itb;
            sumab += *itab;

            ++ita2;
            ++itb2;
            ++ita;
            ++itb;
            ++itab;
        }

        float fixedMean = suma / count;
        float movingMean = sumb / count;

        float sff = suma2 - fixedMean*suma - fixedMean*suma + count*fixedMean*fixedMean;
        float smm = sumb2 - movingMean*sumb - movingMean*sumb + count*movingMean*movingMean;
        float sfm = sumab - movingMean*suma - fixedMean*sumb + count*movingMean*fixedMean;

        typename FixedImageType::IndexType oindex = scan_it.GetIndex();

        float val = scan_para.I->GetPixel( oindex ) - fixedMean;
        float val1 = scan_para.J->GetPixel( oindex ) - movingMean;

        scan_mem.Ia = val;
        scan_mem.Ja = val1;
        scan_mem.sfm = sfm;
        scan_mem.sff = sff;
        scan_mem.smm = smm;
    }
 }



template< class TFixedImage, class TMovingImage >
void
ImageToImageNeighborhoodNormalizedCrossCorrelationFunction< TFixedImage, TMovingImage >
::ComputeUpdateBothDirection(const ScanningIteratorType &scan_it,  ScanMemType &scan_mem, const ScanParaType &scan_para,
        VectorType &deriv, VectorType &derivInv, InternalComputationValueType &local_cc)
 {

    //typename TDeformationField::PixelType deriv, derivInv;
    deriv.Fill(0.0);
    derivInv.Fill(0.0);
    local_cc=1.0;

    typedef InternalComputationValueType TMPREALTYPE;

    IndexType index = scan_it.GetIndex();

//    if (this->m_FixedImageMask) if (this->m_FixedImageMask->GetPixel( index ) < 0.25 )
//        return;

    TMPREALTYPE sff=scan_mem.sff;
    TMPREALTYPE smm=scan_mem.smm;
    TMPREALTYPE sfm=scan_mem.sfm;
    TMPREALTYPE Ji=scan_mem.Ja; //finitediffimages[1]->GetPixel(index);
    TMPREALTYPE Ii=scan_mem.Ia; //finitediffimages[0]->GetPixel(index);

    CovariantVectorType gradI,gradJ;

    if ( sff == 0.0 || smm == 0.0) return;

    gradI = m_FixedImageGradientCalculator->EvaluateAtIndex( index );
    gradJ = m_MovingImageGradientCalculator->EvaluateAtIndex( index );

    for (int qq=0; qq<ImageDimension; qq++)
    {
        deriv[qq]   -=2.0*sfm/(sff*smm)*( Ji - sfm/sff*Ii )*gradI[qq];
        derivInv[qq]-=2.0*sfm/(sff*smm)*( Ii - sfm/smm*Ji )*gradJ[qq];
    }

    if (fabs(sff*smm) > 0) local_cc = sfm*sfm / ( sff * smm );


    if (deriv.GetNorm() > 1e20 || derivInv.GetNorm() > 1e20 ) {
      std::cout << " ------------- DEBUG-----------------------------" << std::endl
            << "deriv.GetNorm() = " << deriv.GetNorm() << std::endl
          << " derivInv.GetNorm() = " << derivInv.GetNorm() << std::endl
          << "sff=" << sff << std::endl
          << "smm=" << smm << std::endl
          << "sfm=" << sfm << std::endl;

      exit(-1);
    }



    return;

 }


} // end namespace itk

#endif
