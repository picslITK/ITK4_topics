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
#ifndef __itkANTSNeighborhoodCorrelationImageToImageObjectMetric_hxx
#define __itkANTSNeighborhoodCorrelationImageToImageObjectMetric_hxx

#include "itkANTSNeighborhoodCorrelationImageToImageObjectMetric.h"
#include "itkNumericTraits.h"

namespace itk {
/**
 * Constructor
 */
template<class TFixedImage, class TMovingImage, class TVirtualImage>
ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage, TMovingImage,
        TVirtualImage>::ANTSNeighborhoodCorrelationImageToImageObjectMetric() {

    //modify the callback function.
    this->m_ValueAndDerivativeThreader->SetThreadedGenerateData(
            Self::NeighborhoodScanningWindowGetValueAndDerivativeMultiThreadedCallback);

}

template<class TFixedImage, class TMovingImage, class TVirtualImage>
ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage, TMovingImage,
        TVirtualImage>::~ANTSNeighborhoodCorrelationImageToImageObjectMetric() {
}

/**
 * Initialize
 */
template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::Initialize(void) throw (ExceptionObject) {
    this->Superclass::Initialize();
}

template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
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

template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::PrintSelf(std::ostream & os,
        Indent indent) const {
    Superclass::PrintSelf(os, indent);
    os << indent << "Correlation window radius: " << m_Radius << std::endl;
}

template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::NeighborhoodScanningWindowGetValueAndDerivativeMultiThreadedCallback(
        const ThreaderInputObjectType& virtualImageSubRegion,
        ThreadIdType threadID, Superclass * dataHolderA) {

    Self *dataHolder = dynamic_cast<Self *>(dataHolderA);

    VirtualPointType virtualPoint;
    FixedOutputPointType mappedFixedPoint;
    // FixedImagePixelType fixedImageValue;
    FixedImageDerivativesType fixedImageDerivatives;
    MovingOutputPointType mappedMovingPoint;
    // MovingImagePixelType movingImageValue;
    MovingImageDerivativesType movingImageDerivatives;
    bool pointIsValid;
    MeasureType metricValueResult;
    MeasureType metricValueSum = 0;

//    DerivativeType localDerivativeResult(
//            dataHolder->GetNumberOfLocalParameters());

    DerivativeType & localDerivativeResult = dataHolder->m_LocalDerivativesPerThread[threadID];

    typedef Self MetricType;

    ScanningIteratorType scan_it;
    ScanParaType scan_para;
    ScanMemType scan_mem;

    /* Create an iterator over the virtual sub region */
    dataHolder->InitializeScanning(virtualImageSubRegion, scan_it, scan_mem,
            scan_para, threadID);

    /* Iterate over the sub region */
    scan_it.GoToBegin();
    while (!scan_it.IsAtEnd()) {
        /* Get the virtual point */
        dataHolder->m_VirtualDomainImage->TransformIndexToPhysicalPoint(
                scan_it.GetIndex(), virtualPoint);

        /* Call the user method in derived classes to do the specific
         * calculations for value and derivative. */
        try {
            dataHolder->UpdateQueues(scan_it, scan_mem, scan_para, threadID);
            pointIsValid = dataHolder->ComputeInformationFromQueues(scan_it,
                    scan_mem, scan_para, threadID);
            dataHolder->ComputeMovingTransformDerivative(scan_it, scan_mem,
                    scan_para, localDerivativeResult, metricValueResult,
                    threadID);

        } catch (ExceptionObject & exc) {
            //NOTE: there must be a cleaner way to do this:
            std::string msg("Caught exception: \n");
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
                    scan_it.GetIndex(), threadID);
        }

        //next index
        ++scan_it;
    }

    /* Store metric value result for this thread. */
    dataHolder->m_MeasurePerThread[threadID] = metricValueSum;

}

template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::InitializeScanning(
        const ImageRegionType &scan_region, ScanningIteratorType &scan_it,
        ScanMemType &scan_mem, ScanParaType &scan_para,
        const ThreadIdType threadID) {

    scan_para.scan_region = scan_region;
    scan_para.I = this->m_FixedImage;
    scan_para.J = this->m_MovingImage;
    scan_para.V = this->m_VirtualDomainImage;
    scan_para.r = this->GetRadius();

    int nb_fill_zero = scan_para.V->GetBufferedRegion().GetIndex(0)
            - (scan_region.GetIndex(0) - scan_para.r[0]);
    if (nb_fill_zero < 0)
        nb_fill_zero = 0;
    scan_para.number_of_fill_zero = nb_fill_zero;

    scan_it = ScanningIteratorType(scan_para.r, scan_para.I, scan_region);
    scan_para.window_length = scan_it.Size();
    scan_para.scan_region_begin_index_dim0 = scan_it.GetBeginIndex()[0];

}

template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::UpdateQueuesAtBeginingOfLine(
        const ScanningIteratorType &scan_it, ScanMemType &scan_mem,
        const ScanParaType &scan_para, const ThreadIdType threadID) {

    int nb_fill_zero = scan_para.number_of_fill_zero;
    unsigned int hoodlen = scan_para.window_length;

    scan_mem.Qsuma2 = SumQueueType(nb_fill_zero, 0.0);
    scan_mem.Qsumb2 = SumQueueType(nb_fill_zero, 0.0);
    scan_mem.Qsuma = SumQueueType(nb_fill_zero, 0.0);
    scan_mem.Qsumb = SumQueueType(nb_fill_zero, 0.0);
    scan_mem.Qsumab = SumQueueType(nb_fill_zero, 0.0);
    scan_mem.Qcount = SumQueueType(nb_fill_zero, 0.0);

    typename FixedImageType::IndexType oindex = scan_it.GetIndex();

    typedef InternalComputationValueType LocalRealType;

    // Now add the rest of the values from each hyperplane
    for (unsigned int i = nb_fill_zero; i < (2 * scan_para.r[0] + 1); i++) {
        LocalRealType suma2 = 0.0;
        LocalRealType sumb2 = 0.0;
        LocalRealType suma = 0.0;
        LocalRealType sumb = 0.0;
        LocalRealType sumab = 0.0;
        LocalRealType count = 0.0;

        for (unsigned int indct = i; indct < hoodlen;
                indct += (2 * scan_para.r[0] + 1)) {

            bool isInBounds = true;
            scan_it.GetPixel(indct, isInBounds);

            typename VirtualImageType::IndexType index = scan_it.GetIndex(indct);

//            if ( !isInBounds || ( this->m_FixedImageMask && this->m_FixedImageMask->GetPixel( index ) < 0.25 ) )
            if (!isInBounds) {
                // std::cout << "DEBUG: error" << std::endl;
                continue;
            }

            VirtualPointType virtualPoint;
            FixedOutputPointType mappedFixedPoint;
            FixedImagePixelType fixedImageValue;
            FixedImageDerivativesType fixedImageDerivatives;
            MovingOutputPointType mappedMovingPoint;
            MovingImagePixelType movingImageValue;
            MovingImageDerivativesType movingImageDerivatives;
            bool pointIsValid;

            this->m_VirtualDomainImage->TransformIndexToPhysicalPoint(index,
                    virtualPoint);

//            if (0) {
             if (this->m_PreWarpImages) {
                try {
                    this->EvaluatePreWarpedImagesAtIndex(index, virtualPoint,
                            false, /* compute gradient */
                            fixedImageValue, movingImageValue,
                            fixedImageDerivatives, movingImageDerivatives,
                            pointIsValid, threadID);
//                    /* Get the point in moving and fixed space for use below */
//                    mappedFixedPoint = self->m_FixedTransform->TransformPoint(
//                            virtualPoint);
//                    mappedMovingPoint = self->m_MovingTransform->TransformPoint(
//                            virtualPoint);
                } catch (ExceptionObject & exc) {
                    //NOTE: there must be a cleaner way to do this. We want to add
                    // the this filename and line number to give user more useful
                    // information about where the exception was generated.
                    std::string msg("Caught exception: \n");
                    msg += exc.what();
                    ExceptionObject err(__FILE__, __LINE__, msg);
                    throw err;
                }
            } else {

                try {
                    this->TransformAndEvaluateFixedPoint( //
                            index,
                            virtualPoint, //
                            mappedFixedPoint, //
                            pointIsValid, //
                            fixedImageValue, //
                            false/*compute gradient*/, // true /*compute gradient*/,
                            fixedImageDerivatives, //
                            threadID);
                    if (pointIsValid) {
                        this->TransformAndEvaluateMovingPoint( //
                                index,
                                virtualPoint, //
                                mappedMovingPoint, //
                                pointIsValid, //
                                movingImageValue, //
                                false/*compute gradient*/, // true /*compute gradient*/,
                                movingImageDerivatives, //
                                threadID);
                    }

                } catch (ExceptionObject & exc) {
                    //NOTE: there must be a cleaner way to do this:
                    std::string msg("Caught exception: \n");
                    msg += exc.what();
                    ExceptionObject err(__FILE__, __LINE__, msg);
                    throw err;
                }
            }

            if (pointIsValid) {

                LocalRealType a = fixedImageValue; //scan_para.I->GetPixel(index);
                LocalRealType b = movingImageValue; // scan_para.J->GetPixel(index);

                suma2 += a * a;
                sumb2 += b * b;
                suma += a;
                sumb += b;
                sumab += a * b;
                count += 1.0;
            }

        }

        scan_mem.Qsuma2.push_back(suma2);
        scan_mem.Qsumb2.push_back(sumb2);
        scan_mem.Qsuma.push_back(suma);
        scan_mem.Qsumb.push_back(sumb);
        scan_mem.Qsumab.push_back(sumab);
        scan_mem.Qcount.push_back(count);

    }

}

template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::UpdateQueuesToNextScanWindow(
        const ScanningIteratorType &scan_it, ScanMemType &scan_mem,
        const ScanParaType &scan_para, const ThreadIdType threadID) {

    const unsigned int hoodlen = scan_para.window_length;

    typedef InternalComputationValueType LocalRealType;

    LocalRealType suma2 = 0.0;
    LocalRealType sumb2 = 0.0;
    LocalRealType suma = 0.0;
    LocalRealType sumb = 0.0;
    LocalRealType sumab = 0.0;
    LocalRealType count = 0.0;

    for (unsigned int indct = 2 * scan_para.r[0]; indct < hoodlen;
            indct += (2 * scan_para.r[0] + 1)) {
        bool isInBounds = true;

        scan_it.GetPixel(indct, isInBounds);
        typename VirtualImageType::IndexType index = scan_it.GetIndex(indct);

        //GS: this is a weird notation: force prob mask < 0.25 to be zero
        // if ( !isInBounds || ( this->m_FixedImageMask && this->m_FixedImageMask->GetPixel( index ) < 0.25 ) )
        if (!isInBounds) {
            continue;
        }

        VirtualPointType virtualPoint;
        FixedOutputPointType mappedFixedPoint;
        FixedImagePixelType fixedImageValue;
        FixedImageDerivativesType fixedImageDerivatives;
        MovingOutputPointType mappedMovingPoint;
        MovingImagePixelType movingImageValue;
        MovingImageDerivativesType movingImageDerivatives;
        bool pointIsValid;

        this->m_VirtualDomainImage->TransformIndexToPhysicalPoint(index,
                virtualPoint);

//        if (0) {
        if (this->m_PreWarpImages) {
            try {
                this->EvaluatePreWarpedImagesAtIndex(index, virtualPoint, false, /* compute gradient */
                fixedImageValue, movingImageValue, fixedImageDerivatives,
                        movingImageDerivatives, pointIsValid, threadID);
                //                    /* Get the point in moving and fixed space for use below */
                //                    mappedFixedPoint = self->m_FixedTransform->TransformPoint(
                //                            virtualPoint);
                //                    mappedMovingPoint = self->m_MovingTransform->TransformPoint(
                //                            virtualPoint);
            } catch (ExceptionObject & exc) {
                //NOTE: there must be a cleaner way to do this. We want to add
                // the this filename and line number to give user more useful
                // information about where the exception was generated.
                std::string msg("Caught exception: \n");
                msg += exc.what();
                ExceptionObject err(__FILE__, __LINE__, msg);
                throw err;
            }
        } else {
            try {
                this->TransformAndEvaluateFixedPoint( //
                        index,
                        virtualPoint, //
                        mappedFixedPoint, //
                        pointIsValid, //
                        fixedImageValue, //
                        false/*compute gradient*/, // true /*compute gradient*/,
                        fixedImageDerivatives, //
                        threadID);
                if (pointIsValid) {
                    this->TransformAndEvaluateMovingPoint( //
                            index,
                            virtualPoint, //
                            mappedMovingPoint, //
                            pointIsValid, //
                            movingImageValue, //
                            false/*compute gradient*/, // true /*compute gradient*/,
                            movingImageDerivatives, //
                            threadID);
                }
            } catch (ExceptionObject & exc) {
                //NOTE: there must be a cleaner way to do this:
                std::string msg("Caught exception: \n");
                msg += exc.what();
                ExceptionObject err(__FILE__, __LINE__, msg);
                throw err;
            }

        }

        if (pointIsValid) {

            LocalRealType a = fixedImageValue; //scan_para.I->GetPixel(index);
            LocalRealType b = movingImageValue; // scan_para.J->GetPixel(index);

            suma2 += a * a;
            sumb2 += b * b;
            suma += a;
            sumb += b;
            sumab += a * b;
            count += 1.0;
        }

    }

    scan_mem.Qsuma2.push_back(suma2);
    scan_mem.Qsumb2.push_back(sumb2);
    scan_mem.Qsuma.push_back(suma);
    scan_mem.Qsumb.push_back(sumb);
    scan_mem.Qsumab.push_back(sumab);
    scan_mem.Qcount.push_back(count);

    scan_mem.Qsuma2.pop_front();
    scan_mem.Qsumb2.pop_front();
    scan_mem.Qsuma.pop_front();
    scan_mem.Qsumb.pop_front();
    scan_mem.Qsumab.pop_front();
    scan_mem.Qcount.pop_front();

}

template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::UpdateQueues(
        const ScanningIteratorType &scan_it, ScanMemType &scan_mem,
        const ScanParaType &scan_para, const ThreadIdType threadID) {

    if (scan_it.GetIndex()[0] == scan_para.scan_region_begin_index_dim0) {
        UpdateQueuesAtBeginingOfLine(scan_it, scan_mem, scan_para, threadID);
    } else {
        UpdateQueuesToNextScanWindow(scan_it, scan_mem, scan_para, threadID);
    }

}

template<class TFixedImage, class TMovingImage, class TVirtualImage>
bool ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::ComputeInformationFromQueues(
        const ScanningIteratorType &scan_it, ScanMemType &scan_mem,
        const ScanParaType &scan_para, const ThreadIdType threadID) {

    // Test to see if there are any voxels we need to handle in the current
    // window.

    typedef InternalComputationValueType LocalRealType;

    LocalRealType suma2 = 0.0;
    LocalRealType sumb2 = 0.0;
    LocalRealType suma = 0.0;
    LocalRealType sumb = 0.0;
    LocalRealType sumab = 0.0;
    LocalRealType count = 0.0;

    typename SumQueueType::iterator itcount = scan_mem.Qcount.begin();
    while (itcount != scan_mem.Qcount.end()) {
        count += *itcount;
        ++itcount;
    }

    if (count < 0) {
        // no points avaiable in the queue, perhapse out of image region
        return false;
    }

    // If there are values, we need to calculate the different quantities
    typename SumQueueType::iterator ita2 = scan_mem.Qsuma2.begin();
    typename SumQueueType::iterator itb2 = scan_mem.Qsumb2.begin();
    typename SumQueueType::iterator ita = scan_mem.Qsuma.begin();
    typename SumQueueType::iterator itb = scan_mem.Qsumb.begin();
    typename SumQueueType::iterator itab = scan_mem.Qsumab.begin();

    while (ita2 != scan_mem.Qsuma2.end()) {
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

    LocalRealType fixedMean = suma / count;
    LocalRealType movingMean = sumb / count;

    LocalRealType sff = suma2 - fixedMean * suma - fixedMean * suma
            + count * fixedMean * fixedMean;
    LocalRealType smm = sumb2 - movingMean * sumb - movingMean * sumb
            + count * movingMean * movingMean;
    LocalRealType sfm = sumab - movingMean * suma - fixedMean * sumb
            + count * movingMean * fixedMean;

    typename VirtualImageType::IndexType oindex = scan_it.GetIndex();

    VirtualPointType virtualPoint;
    FixedOutputPointType mappedFixedPoint;
    FixedImagePixelType fixedImageValue;
    FixedImageDerivativesType fixedImageDerivatives;
    MovingOutputPointType mappedMovingPoint;
    MovingImagePixelType movingImageValue;
    MovingImageDerivativesType movingImageDerivatives;
    bool pointIsValid;

    this->m_VirtualDomainImage->TransformIndexToPhysicalPoint(oindex,
            virtualPoint);

//    if (0) {
    if (this->m_PreWarpImages) {
        try {
            this->EvaluatePreWarpedImagesAtIndex(oindex, virtualPoint, true, /* compute gradient */
            fixedImageValue, movingImageValue, fixedImageDerivatives,
                    movingImageDerivatives, pointIsValid, threadID);
            //                    /* Get the point in moving and fixed space for use below */
//                       mappedFixedPoint = self->m_FixedTransform->TransformPoint(
//                               virtualPoint);
            mappedMovingPoint = this->m_MovingTransform->TransformPoint(
                    virtualPoint);
        } catch (ExceptionObject & exc) {
            //NOTE: there must be a cleaner way to do this. We want to add
            // the this filename and line number to give user more useful
            // information about where the exception was generated.
            std::string msg("Caught exception: \n");
            msg += exc.what();
            ExceptionObject err(__FILE__, __LINE__, msg);
            throw err;
        }
    } else {
        try {
            this->TransformAndEvaluateFixedPoint(
                    oindex,
                    virtualPoint, //
                    mappedFixedPoint, //
                    pointIsValid, //
                    fixedImageValue, //
                    true/*compute gradient*/, // true /*compute gradient*/,
                    fixedImageDerivatives, threadID);
            if (pointIsValid) {
                this->TransformAndEvaluateMovingPoint(
                        oindex,
                        virtualPoint, //
                        mappedMovingPoint, //
                        pointIsValid, //
                        movingImageValue, //
                        true/*compute gradient*/, // true /*compute gradient*/,
                        movingImageDerivatives, threadID);
            }
        } catch (ExceptionObject & exc) {
            //NOTE: there must be a cleaner way to do this:
            std::string msg("Caught exception: \n");
            msg += exc.what();
            ExceptionObject err(__FILE__, __LINE__, msg);
            throw err;
        }

    }

    if (pointIsValid) {

        float val = fixedImageValue - fixedMean; // scan_para.I->GetPixel(oindex) - fixedMean;
        float val1 = movingImageValue - movingMean; // scan_para.J->GetPixel(oindex) - movingMean;
        scan_mem.Ia = val;
        scan_mem.Ja = val1;
        scan_mem.sfm = sfm;
        scan_mem.sff = sff;
        scan_mem.smm = smm;

        scan_mem.gradI = fixedImageDerivatives;
        scan_mem.gradJ = movingImageDerivatives;

        scan_mem.mappedMovingPoint = mappedMovingPoint;

    }

    return true;

}

template<class TFixedImage, class TMovingImage, class TVirtualImage>
void ANTSNeighborhoodCorrelationImageToImageObjectMetric<TFixedImage,
        TMovingImage, TVirtualImage>::ComputeMovingTransformDerivative(
        const ScanningIteratorType &scan_it, ScanMemType &scan_mem,
        const ScanParaType &scan_para, DerivativeType &deriv,
        MeasureType &local_cc, const ThreadIdType threadID) {

    MovingImageDerivativesType deriv_wrt_image;

    //    MovingImageDerivativesType deriv_wrt_image(0.0);
    //    DerivativeType deriv_wrt_image(VirtualImageDimension);
    //    deriv_wrt_image.Fill(0.0);

    //    deriv.Fill(0.0);
    //    derivInv.Fill(0.0);
    local_cc = 1.0;

    typedef InternalComputationValueType LocalRealType;

    //    IndexType index = scan_it.GetIndex();

    LocalRealType sff = scan_mem.sff;
    LocalRealType smm = scan_mem.smm;
    LocalRealType sfm = scan_mem.sfm;
    LocalRealType Ji = scan_mem.Ja;
    LocalRealType Ii = scan_mem.Ia;

    FixedImageDerivativesType gradI = scan_mem.gradI;
    MovingImageDerivativesType gradJ = scan_mem.gradJ;

    if (sff == 0.0 || smm == 0.0) {
        deriv.Fill(0.0);
        return;
    }

    for (unsigned int qq = 0; qq < VirtualImageDimension; qq++) {
//        deriv_wrt_image[qq] -= 2.0 * sfm / (sff * smm) * (Ji - sfm / sff * Ii)
//        *gradI[qq];
        deriv_wrt_image[qq] = 2.0 * sfm / (sff * smm) * (Ii - sfm / smm * Ji)
                * gradJ[qq];
    }

    if (fabs(sff * smm) > 0)
        local_cc = sfm * sfm / (sff * smm);

    /* Use a pre-allocated jacobian object for efficiency */
    MovingTransformJacobianType & jacobian =
            this->m_MovingTransformJacobianPerThread[threadID];

    /** For dense transforms, this returns identity */
    this->m_MovingTransform->GetJacobianWithRespectToParameters(
            scan_mem.mappedMovingPoint, jacobian);

    unsigned int numberOfLocalParameters = this->m_MovingTransform->GetNumberOfLocalParameters();

    for (unsigned int par = 0; par < numberOfLocalParameters; par++) {

        double sum = 0.0;
        for (unsigned int dim = 0; dim < this->MovingImageDimension; dim++) {
            sum += deriv_wrt_image[dim] * jacobian(dim, par);
        }

        deriv[par] = sum;
    }

    return;

}

} // end namespace itk

#endif
