/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    $RCSfile: itkDemonsImageToImageObjectMetric.hxx,v $
 Language:  C++
 Date:      $Date: $
 Version:   $Revision: $

 Copyright (c) Insight Software Consortium. All rights reserved.
 See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notices for more information.

 =========================================================================*/
#ifndef __itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetric_hxx
#define __itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetric_hxx

#include "itkANTSSparseNeighborhoodCorrelationImageToImageObjectMetric.h"
#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "vnl/vnl_math.h"

namespace itk {

/**
 * Constructor
 */
template<class TFixedImage, class TMovingImage, class TVirtualImage>
ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric() {

    //modify the callback function.
    this->m_ValueAndDerivativeThreader->SetThreadedGenerateData(
            Self::SparseSamplingGetValueAndDerivativeMultiThreadedCallback);

    this->m_NumberOfSampling = 1;

}

template<class TFixedImage, class TMovingImage, class TVirtualImage>
ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::~ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric() {
}

/*
 * GetValueAndDerivative
 */
template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::GetValueAndDerivative(MeasureType & value,
        DerivativeType & derivative) {
    // This starts threading, and will iterate over virtual image region and
    // call GetValueAndDerivativeProcessPoint.
    this->GetValueAndDerivativeMultiThreadedInitiate(derivative);

    // Sums up results from each thread, and optionally averages them.
    // Derivative results are written directly to \c derivative.
    this->GetValueAndDerivativeMultiThreadedPostProcess(true /*doAverage*/);

    value = this->GetValueResult();
}

/** This function computes the local voxel-wise contribution of
 *  the metric to the global integral of the metric/derivative.
 */
template<class TFixedImage, class TMovingImage, class TVirtualImage>
bool ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::SparseGetValueAndDerivativeProcessPoint(
        const VirtualImageRegionType &scan_region,
        const VirtualIndexType & virtualIndex,
        const VirtualPointType & virtualPoint,
        const FixedImagePointType & mappedFixedPoint,
        const FixedImagePixelType & fixedImageValue,
        const FixedImageDerivativesType & fixedImageDerivatives,
        const MovingImagePointType & mappedMovingPoint,
        const MovingImagePixelType & movingImageValue,
        const MovingImageDerivativesType & movingImageDerivatives,
        MeasureType & metricValueReturn, DerivativeType & localDerivativeReturn,
        ThreadIdType threadID) {
    bool pointIsValid;

    ScanningIteratorType scan_it;
    ScanParaType scan_para;
    ScanMemType scan_mem;

    /* Create an iterator over the virtual sub region */
    this->InitializeScanning(scan_region, virtualIndex, scan_it, scan_mem,
            scan_para, threadID);

    Superclass::UpdateQueuesAtBeginingOfLine(scan_it, scan_mem, scan_para,
            threadID);

    pointIsValid = Superclass::ComputeInformationFromQueues(scan_it, scan_mem,
            scan_para, threadID);

    if (pointIsValid) {
        Superclass::ComputeMovingTransformDerivative(scan_it, scan_mem,
                scan_para, localDerivativeReturn, metricValueReturn, threadID);
    }

    return pointIsValid;
}

/**
 * Print out internal information about this class
 */
template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::PrintSelf(std::ostream& os,
        Indent indent) const {
    Superclass::PrintSelf(os, indent);
    os << indent << "Number of sampling: " << m_NumberOfSampling << std::endl;
}

/**
 * Initialize
 */
template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::Initialize(void) throw (ExceptionObject) {
    this->Superclass::Initialize();
}

template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::SparseSamplingGetValueAndDerivativeMultiThreadedCallback(
        const ThreaderInputObjectType& virtualImageSubRegion,
        ThreadIdType threadID, MetricBaseclass * self) {

    Self *dataHolder = dynamic_cast<Self *>(self);

    VirtualIndexType virtualIndex;
    VirtualPointType virtualPoint;
    FixedOutputPointType mappedFixedPoint;
    FixedImagePixelType fixedImageValue;
    FixedImageDerivativesType fixedImageDerivatives;
    MovingOutputPointType mappedMovingPoint;
    MovingImagePixelType movingImageValue;
    MovingImageDerivativesType movingImageDerivatives;
    bool pointIsValid;
    MeasureType metricValueResult;
    MeasureType metricValueSum = 0;



    /* Get pre-allocated local results object. This actually provides very
     * little benefit, since this only gets called once for each thread. However
     * if we get up to the hundres of threads, it might have an impact */
    DerivativeType & localDerivativeResult =
            dataHolder->m_LocalDerivativesPerThread[threadID];


    // get random sampling routine here!!
    VirtualImageSampleContainer sampleLocations;

    if (dataHolder->m_NumberOfSampling > 0) {

        unsigned int numberOfSamplingInThisThread =
                dataHolder->m_NumberOfSampling
                        / dataHolder->GetNumberOfThreads();
        sampleLocations.resize(numberOfSamplingInThisThread);
        dataHolder->SampleVirtualImageRegion(virtualImageSubRegion,
                sampleLocations);
    } else { // use the full image region
        dataHolder->FullSampleVirtualImageRegion(virtualImageSubRegion,
                sampleLocations);
    }



    for (unsigned int i = 0; i < sampleLocations.size(); i++) {

        virtualIndex = sampleLocations[i];

        dataHolder->m_VirtualDomainImage->TransformIndexToPhysicalPoint(
                virtualIndex, virtualPoint);

        /* Transform the point into fixed and moving spaces, and evaluate.
         * These methods will check that the point lies within the mask if
         * one has been set, and then verify they lie in the fixed or moving
         * space as appropriate.
         * If both tests pass, the point is evaluated and pointIsValid is
         * returned as \c true.
         * Do this in a try block to catch exceptions and print more useful info
         * then we otherwise get when exceptions are caught in MultiThreader. */
        try {
            dataHolder->TransformAndEvaluateFixedPoint(
                    virtualIndex, virtualPoint,
                    mappedFixedPoint, pointIsValid, fixedImageValue,
                    false /*compute gradient*/, fixedImageDerivatives,
                    threadID);
            if (pointIsValid) {
                dataHolder->TransformAndEvaluateMovingPoint(
                        virtualIndex, virtualPoint,
                        mappedMovingPoint, pointIsValid, movingImageValue,
                        false /*compute gradient*/, movingImageDerivatives,
                        threadID);
            }
        } catch (ExceptionObject & exc) {
            //NOTE: there must be a cleaner way to do this:
            std::string msg("Caught exception: \n");
            msg += exc.what();
            ExceptionObject err(__FILE__, __LINE__, msg);
            throw err;
        }

        /* Call the user method in derived classes to do the specific
         * calculations for value and derivative. */
        try {
            if (pointIsValid) {
                pointIsValid =
                        dataHolder->SparseGetValueAndDerivativeProcessPoint(
                                virtualImageSubRegion, virtualIndex,
                                virtualPoint, mappedFixedPoint, fixedImageValue,
                                fixedImageDerivatives, mappedMovingPoint,
                                movingImageValue, movingImageDerivatives,
                                metricValueResult, localDerivativeResult,
                                threadID);
            }
        } catch (ExceptionObject & exc) {
            //NOTE: there must be a cleaner way to do this:
            std::string msg(
                    "Exception in GetValueAndDerivativeProcessPoint:\n");
            msg += exc.what();
            ExceptionObject err(__FILE__, __LINE__, msg);
            throw err;
        }

        /* Assign the results */
        if (pointIsValid) {
            dataHolder->m_NumberOfValidPointsPerThread[threadID]++;
            metricValueSum += metricValueResult;
            /* Store the result. This depends on what type of
             * transform is being used. */
            dataHolder->StoreDerivativeResult(localDerivativeResult,
                    virtualIndex, threadID);
        }
    }


    /* Store metric value result for this thread. */
    dataHolder->m_MeasurePerThread[threadID] = metricValueSum;
}


template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::InitializeScanning(
        const VirtualImageRegionType &scan_region,
        const VirtualIndexType &start_index, ScanningIteratorType &scan_it,
        ScanMemType &scan_mem, ScanParaType &scan_para,
        const ThreadIdType threadID) {

    scan_para.scan_region = scan_region;
    scan_para.I = this->m_FixedImage;
    scan_para.J = this->m_MovingImage;
    scan_para.r = this->GetRadius();

    int nb_fill_zero = scan_para.I->GetBufferedRegion().GetIndex(0)
            - (start_index[0] - scan_para.r[0]);
    if (nb_fill_zero < 0)
        nb_fill_zero = 0;

    scan_para.number_of_fill_zero = nb_fill_zero;

    scan_it = ScanningIteratorType(scan_para.r, scan_para.I, scan_region);
    scan_it.SetLocation(start_index);

    scan_para.window_length = scan_it.Size();
    scan_para.scan_region_begin_index_dim0 = scan_it.GetBeginIndex()[0];

}


template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::SampleVirtualImageRegion(
                const VirtualImageRegionType &region,
                VirtualImageSampleContainer & samples)
                {

//    if (samples.size() != m_NumberOfSampling) {
//        throw ExceptionObject(__FILE__, __LINE__,
//                "Sample size does not match desired number of samples");
//    }

    // Set up a random iterator within the user specified fixed image region.
    typedef ImageRandomConstIteratorWithIndex<VirtualImageType> RandomIterator;
    RandomIterator randIter(this->m_VirtualDomainImage, region);

    typename VirtualImageSampleContainer::iterator iter;
    typename VirtualImageSampleContainer::const_iterator end = samples.end();

    randIter.SetNumberOfSamples(m_NumberOfSampling);
    randIter.GoToBegin();
    for (iter = samples.begin(); iter != end; ++iter) {
        // Get sampled index
        VirtualIndexType index = randIter.GetIndex();
        // Translate index to point
        *iter = index;

        // Jump to random position
        ++randIter;
    }
}


template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSSparseNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::FullSampleVirtualImageRegion(
                const VirtualImageRegionType &region,
                VirtualImageSampleContainer & samples)
                {

//    if (samples.size() != m_NumberOfSampling) {
//        throw ExceptionObject(__FILE__, __LINE__,
//                "Sample size does not match desired number of samples");
//    }


    samples.clear();
    samples.resize(region.GetNumberOfPixels());

    // Set up a random iterator within the user specified fixed image region.
    //typedef ImageRandomConstIteratorWithIndex<VirtualImageType> RandomIterator;
    typedef ImageRegionConstIteratorWithIndex< FixedImageType > RegionIterator;

    RegionIterator regionIter(this->m_VirtualDomainImage, region);

    typename VirtualImageSampleContainer::iterator iter;
    typename VirtualImageSampleContainer::const_iterator end = samples.end();

    regionIter.GoToBegin();
    for (iter = samples.begin(); iter != end; ++iter) {
        // Get sampled index
        VirtualIndexType index = regionIter.GetIndex();
        // Translate index to point
        *iter = index;

        // Jump to random position
        ++regionIter;
    }
}

} // end namespace itk

#endif
