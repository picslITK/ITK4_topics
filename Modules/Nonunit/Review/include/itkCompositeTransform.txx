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
#ifndef __itkCompositeTransform_txx
#define __itkCompositeTransform_txx

#include "itkCompositeTransform.h"
#include <string.h> // for memcpy on some platforms

namespace itk
{

/**
 * Constructor
 */
template
<class TScalar, unsigned int NDimensions>
CompositeTransform<TScalar, NDimensions>::
CompositeTransform() : Superclass( NDimensions, 0 )
{
    this->m_TransformQueue.clear();
    this->m_TransformsToOptimizeFlags.clear();
    this->m_TransformsToOptimizeQueue.clear();
    this->m_PreviousTransformsToOptimizeUpdateTime = 0;
}

/**
 * Destructor
 */
template
<class TScalar, unsigned int NDimensions>
CompositeTransform<TScalar, NDimensions>::
~CompositeTransform()
{
}

template
<class TScalar, unsigned int NDimensions>
bool CompositeTransform<TScalar, NDimensions>
::IsLinear() const
 {
    typename TransformQueueType::const_iterator it;

    for( it = this->m_TransformQueue.begin();
            it != this->m_TransformQueue.end(); ++it )
    {
        if ( !(*it)->IsLinear() )
        {
            return false;
        }
    }
    return true;
 }

/**
 * Transform point
 */
template
<class TScalar, unsigned int NDimensions>
typename CompositeTransform<TScalar, NDimensions>
::OutputPointType
 CompositeTransform<TScalar, NDimensions>
::TransformPoint( const InputPointType& inputPoint ) const
 {
    OutputPointType outputPoint( inputPoint );
    typename TransformQueueType::const_iterator it;
    /* Apply in reverse queue order.  */
    it = this->m_TransformQueue.end();
    do
    {
        it--;
        outputPoint = (*it)->TransformPoint( outputPoint );
    }
    while (it != this->m_TransformQueue.begin() );

    return outputPoint;
 }

/**
 * return an inverse transformation
 */
template
<class TScalar, unsigned int NDimensions>
bool
CompositeTransform<TScalar, NDimensions>
::GetInverse( Self *inverse ) const
 {
    typename TransformQueueType::const_iterator it;

    inverse->ClearTransformQueue();

    for( it = this->m_TransformQueue.begin();
            it != this->m_TransformQueue.end(); ++it )
    {
        TransformTypePointer inverseTransform = dynamic_cast<Superclass *>(
                ( ( *it )->GetInverseTransform() ).GetPointer() );
        if ( !inverseTransform )
        {
            inverse->ClearTransformQueue();
            return false;
        }
        else
        {
            /* This also sets TransformToOptimizeFlags list, but it's reset below. */
            inverse->PushFrontTransform( inverseTransform );
        }
    }

    /* Copy the optimization flags */
    inverse->m_TransformsToOptimizeFlags.clear();
    for( TransformsToOptimizeFlagsType::iterator
            ofit = this->m_TransformsToOptimizeFlags.begin();
            ofit != this->m_TransformsToOptimizeFlags.end(); ofit++ )
    {
        inverse->m_TransformsToOptimizeFlags.push_front( *ofit );
    }

    return true;
 }

/**
 * Return an inverse of this transform
 */
template
<class TScalar, unsigned int NDimensions>
typename CompositeTransform<TScalar, NDimensions>
::InverseTransformBasePointer
 CompositeTransform<TScalar, NDimensions>
::GetInverseTransform() const
 {
    Pointer inverseTransform = New();
    if( this->GetInverse( inverseTransform ) )
    {
        return inverseTransform.GetPointer();
    }
    else
    {
        return NULL;
    }
 }

template
<class TScalar, unsigned int NDimensions>
const typename CompositeTransform< TScalar, NDimensions >::JacobianType &
CompositeTransform<TScalar, NDimensions>
::GetJacobian( const InputPointType & p ) const
 {
    /* Returns a concatenated MxN array, holding the Jacobian of each sub
     * transform that is selected for optimization. The order is the same
     * as that in which they're applied, i.e. reverse order.
     * M rows = dimensionality of the transforms
     * N cols = total number of parameters in the selected sub transforms. */

    this->m_Jacobian.SetSize( NDimensions, this->GetNumberOfParameters() );
    unsigned int offset = 0;
    OutputPointType transformedPoint( p );

    for( signed long tind = (signed long) this->GetNumberOfTransforms()-1;
            tind >= 0; tind-- )
    {
        TransformTypePointer transform = this->GetNthTransform( tind );
        if( this->GetNthTransformToOptimize( tind ) )
        {
            /* Copy from another matrix, element-by-element */
            /* The matrices are row-major, so block copy is less obviously better */
            this->m_Jacobian.update(
                    transform->GetJacobian( transformedPoint ), 0, offset );
            offset += transform->GetParameters().Size();
        }
        /* Transform the point so it's ready for next transform's Jacobian */
        transformedPoint = transform->TransformPoint( transformedPoint );
    }

    return this->m_Jacobian;
 }


template
<class TScalar, unsigned int NDimensions>
void
CompositeTransform<TScalar, NDimensions>
::GetJacobianWithRespectToParameters( const InputPointType & p, JacobianType &j ) const
 {
    /* Returns a concatenated MxN array, holding the Jacobian of each sub
     * transform that is selected for optimization. The order is the same
     * as that in which they're applied, i.e. reverse order.
     * M rows = dimensionality of the transforms
     * N cols = total number of parameters in the selected sub transforms. */
    j.SetSize( NDimensions, this->GetNumberOfLocalParameters() );
    unsigned int offset = 0;
//    unsigned int offset_previous = -1;
    OutputPointType transformedPoint( p );

    for( signed long tind = (signed long) this->GetNumberOfTransforms()-1;
            tind >= 0; tind-- )
    {
        TransformTypePointer transform = this->GetNthTransform( tind );
        if( this->GetNthTransformToOptimize( tind ) )
        {
            /* Copy from another matrix, element-by-element */
            /* The matrices are row-major, so block copy is less
             * obviously better */

            // to do: why parameters are listed from N-1 to 1???
            typename TransformType::JacobianType current_jacobian;
            current_jacobian.SetSize(
              NDimensions, transform->GetNumberOfLocalParameters());
            transform->GetJacobianWithRespectToParameters(
              transformedPoint, current_jacobian );

            // debug: force only the closes transform to update!!
//            if (offset > 0){
//                current_jacobian.Fill(0.0);
//            }

            j.update( current_jacobian, 0, offset );
//            std::cout << "cur_tr=" << transform->GetParameters() << std::endl;
//            std::cout << "new_j=" << std::endl << j << std::endl;

            offset += transform->GetNumberOfLocalParameters();

        }

        /** The composite transform needs to compose previous jacobians
         *  (those closer to the originating point) with the current
         *  transform's jacobian.  We therefore update the previous
         *  jacobian by multiplying the current matrix jumping over the
         *  first transform. The matrix here refers to  dT/dx at the point.
         *  For example, in the affine transform, this is the affine matrix.
         *  TODO1: for general transform, there should be something like
         *  GetPartialDerivativeOfPointCoordinates
         *
         *  Also, noted the multiplication contains all the affine matrix from
         *  all transforms no matter they are going to be optimized or not
         *
         */

        if (offset > 0){
            JacobianType old_j = j.extract(NDimensions,offset,0,0);
            j.update( transform->GetJacobianWithRespectToPosition() * old_j, 0, 0);

        }

        /* Transform the point so it's ready for next transform's Jacobian */
        transformedPoint = transform->TransformPoint( transformedPoint );
    }


    return;
 }


template
<class TScalar, unsigned int NDimensions>
const typename CompositeTransform< TScalar, NDimensions >::ParametersType &
CompositeTransform<TScalar, NDimensions>
::GetParameters( ) const
 {
    TransformQueueType transforms = this->GetTransformsToOptimizeQueue();
    if( transforms.size() == 1 )
    {
      // Return directly to avoid copying. Most often we'll have only a single
      // active transform, so we'll end up here.
      return transforms[0]->GetParameters();
    }
    else
    {
        /* Resize destructively. But if it's already this size, nothing is done so
         * it's efficient. */
        this->m_Parameters.SetSize( this->GetNumberOfParameters() );

        unsigned int      offset = 0;
        typename TransformQueueType::const_iterator it;

        it = transforms.end();
        do
        {
            it--;
            const ParametersType & subParameters = (*it)->GetParameters();
            /* use vnl_vector data_block() to get data ptr */
            memcpy( &(this->m_Parameters.data_block())[offset],
                    subParameters.data_block(),
                    subParameters.Size()
                    * sizeof( ParametersValueType ) );
            offset += subParameters.Size();

        } while (it != transforms.begin() );
    }

    return this->m_Parameters;
 }

template
<class TScalar, unsigned int NDimensions>
void
CompositeTransform<TScalar, NDimensions>
::SetParameters(const ParametersType & inputParameters)
 {
    /* We do not copy inputParameters into m_Parameters,
     * to avoid unnecessary copying. */

    /* Assumes input params are concatenation of the parameters of the
     sub transforms currently selected for optimization, in
     the order of the queue from begin() to end(). */
    TransformQueueType transforms = this->GetTransformsToOptimizeQueue();

    /* Verify proper input size. */
    if( inputParameters.Size() != this->GetNumberOfParameters() )
    {
        itkExceptionMacro(<< "Input parameter list size is not expected size. "
                << inputParameters.Size() << " instead of "
                << this->GetNumberOfParameters() << ".");
    }

    if( transforms.size() == 1 )
      {
      /* Avoid unnecessary copying. See comments below */
      if( &inputParameters == &this->m_Parameters )
        {
        transforms[0]->SetParameters( transforms[0]->GetParameters() );
        }
      else
        {
        transforms[0]->SetParameters(inputParameters);
        }
      }
    else
      {
      unsigned int      offset = 0;
      typename TransformQueueType::const_iterator it;

      it = transforms.end();
      do
        {
        it--;

        /* If inputParams is same object as m_Parameters, we just pass
         * each sub-transforms own m_Parameters in. This is needed to
         * avoid unnecessary copying of parameters in the sub-transforms,
         * while still allowing SetParameters to do any oeprations on the
         * parameters to update member variable states. A hack. */
        ParametersType & subParameters =
          const_cast<ParametersType&>( (*it)->GetParameters() );
        if( &inputParameters == &this->m_Parameters )
          {
          (*it)->SetParameters( subParameters );
          }
        else
          {
          /* New parameter data, so copy it in */
          /* Use vnl_vector data_block() to get data ptr */
          memcpy( subParameters.data_block(),
                  &(inputParameters.data_block())[offset],
                  subParameters.Size()
                  * sizeof( ParametersValueType ) );
          /* Call SetParameters explicitly to include anything extra it does */
          (*it)->SetParameters(subParameters);
          offset += subParameters.Size();
          }
        } while (it != transforms.begin() );
      }
    return;
 }

template
<class TScalar, unsigned int NDimensions>
const typename CompositeTransform< TScalar, NDimensions >::ParametersType &
CompositeTransform<TScalar, NDimensions>
::GetFixedParameters(void) const
 {
    TransformQueueType transforms = this->GetTransformsToOptimizeQueue();
    /* Resize destructively. But if it's already this size, nothing is done so
     * it's efficient. */
    this->m_FixedParameters.SetSize( this->GetNumberOfFixedParameters() );

    unsigned int      offset = 0;
    typename TransformQueueType::const_iterator it;

    it = transforms.end();
    do
    {
        it--;
        const ParametersType & subFixedParameters = (*it)->GetFixedParameters();
        /* use vnl_vector data_block() to get data ptr */
        memcpy( &(this->m_FixedParameters.data_block())[offset],
                subFixedParameters.data_block(),
                subFixedParameters.Size()
                * sizeof( ParametersValueType ) );
        offset += subFixedParameters.Size();
    } while (it != transforms.begin() );

    return this->m_FixedParameters;
 }

template
<class TScalar, unsigned int NDimensions>
void
CompositeTransform<TScalar, NDimensions>
::SetFixedParameters(const ParametersType & inputParameters)
 {
    /* Assumes input params are concatenation of the parameters of the
     sub transforms currently selected for optimization. */
    TransformQueueType transforms = this->GetTransformsToOptimizeQueue();

    unsigned int      offset = 0;
    typename TransformQueueType::const_iterator it;

    /* Verify proper input size. */
    if( inputParameters.Size() != this->GetNumberOfFixedParameters() )
    {
        std::cerr << "CompositeTransform::SetFixedParameters error: object dump: "
                << std::endl << this;
        itkExceptionMacro(<< "Input parameter list size is not expected size. "
                << inputParameters.Size() << " instead of "
                << this->GetNumberOfFixedParameters() << ".");
    }
    this->m_FixedParameters = inputParameters;

    it = transforms.end();
    do
    {
        it--;
        ParametersType & subFixedParameters =
          const_cast<ParametersType&>( (*it)->GetFixedParameters() );
        /* Use vnl_vector data_block() to get data ptr */
        memcpy( subFixedParameters.data_block(),
                &(this->m_FixedParameters.data_block())[offset],
                subFixedParameters.Size()
                * sizeof( ParametersValueType ) );
        /* Call SetParameters explicitly to include anything extra it does */
        (*it)->SetFixedParameters(subFixedParameters);
        offset += subFixedParameters.Size();
    } while (it != transforms.begin() );

    return;
 }


template
<class TScalar, unsigned int NDimensions>
unsigned int
CompositeTransform<TScalar, NDimensions>
::GetNumberOfParameters(void) const
 {
    /* Returns to total number of params in all transforms currently
     * set to be used for optimized.
     * NOTE: We might want to optimize this only to store the result and
     * only re-calc when the composite object has been modified.
     * However, it seems that number of parameter might change for dense
     * field transfroms (deformation, bspline) during processing and
     * we wouldn't know that in this class, so this is safest. */
    unsigned int result = 0;
    TransformTypePointer transform;

    for( signed long tind = (signed long) this->GetNumberOfTransforms()-1;
            tind >= 0; tind-- )
      {
      if( this->GetNthTransformToOptimize( tind ) )
        {
        transform = this->GetNthTransform( tind );
        result += transform->GetNumberOfParameters();
        }
      }
    return result;
 }

template
<class TScalar, unsigned int NDimensions>
unsigned int
CompositeTransform<TScalar, NDimensions>
::GetNumberOfLocalParameters(void) const
 {
    /* Returns to total number of *local* params in all transforms currently
     * set to be used for optimized.
     * NOTE: We might want to optimize this only to store the result and
     * only re-calc when the composite object has been modified.
     * However, it seems that number of parameter might change for dense
     * field transfroms (deformation, bspline) during processing and
     * we wouldn't know that in this class, so this is safest. */
    unsigned int result = 0;
    TransformTypePointer transform;

    for( signed long tind = (signed long) this->GetNumberOfTransforms()-1;
            tind >= 0; tind-- )
      {
      if( this->GetNthTransformToOptimize( tind ) )
        {
        transform = this->GetNthTransform( tind );
        result += transform->GetNumberOfLocalParameters();
        }
      }
    return result;
 }

template
<class TScalar, unsigned int NDimensions>
unsigned int
CompositeTransform<TScalar, NDimensions>
::GetNumberOfFixedParameters(void) const
 {
    /* Returns to total number of params in all transforms currently
     * set to be used for optimized.
     * NOTE: We might want to optimize this only to store the result and
     * only re-calc when the composite object has been modified. */
    unsigned int result = 0;
    TransformTypePointer transform;

    for( signed long tind = (signed long) this->GetNumberOfTransforms()-1;
            tind >= 0; tind-- )
      {
      if( this->GetNthTransformToOptimize( tind ) )
        {
        transform = this->GetNthTransform( tind );
        result += transform->GetFixedParameters().Size();
        }
      }
    return result;
 }

template
<class TScalar, unsigned int NDimensions>
void
CompositeTransform<TScalar, NDimensions>
::UpdateTransformParameters(  DerivativeType & update,
                              ScalarType  factor )
{
  /* Update parameters within the sub-transforms set to be optimized. */
  /* NOTE: We might want to thread this over each sub-transform, if we
   * find we're working with longer lists of sub-transforms that do
   * not implement any threading of their own for UpdateTransformParameters.
   * Since the plan is for a UpdateTransformParameters functor that is
   * user-assignable, we would need a method in the
   * functor to return whether or not it does therading. If all sub-transforms
   * return that they don't thread, we could do each sub-transform in its
   * own thread from here. */
  unsigned int numberOfParameters = this->GetNumberOfParameters();

  if( update.Size() != numberOfParameters )
    {
    itkExceptionMacro("Parameter update size, " << update.Size() << ", must "
                      " be same as transform parameter size, "
                      << numberOfParameters << std::endl);
    }

  unsigned int offset = 0;
  TransformTypePointer subtransform;

  for( signed long tind = (signed long) this->GetNumberOfTransforms()-1;
          tind >= 0; tind-- )
    {
    if( this->GetNthTransformToOptimize( tind ) )
      {
      transform = this->GetNthTransform( tind );
      /* The input values are in a monolithic block, so we have to point
       * to the subregion corresponding to the individual subtransform.
       * This simply creates an Array object with data pointer, no
       * memory is allocated or copied. */
      DerivativeType subUpdate( &((update.data_block())[offset]),
                                subtransform->GetNumberOfParameters(), false );
      /* This call will also call SetParameters, so don't need to call it
       * expliclity here. */
      transform->UpdateTransformParameters( subUpdate, factor );
      offset += subtransform->GetNumberOfParameters();
      }
    }
  this->Modified();
}

/* Get local support flag */
template
<class TScalar, unsigned int NDimensions>
bool
CompositeTransform<TScalar, NDimensions>
::HasLocalSupport() const
{
  /* We only return true if all subtransforms return true for
   * HasLocalSupport. */
  bool result = true;
  for( signed long tind = (signed long) this->GetNumberOfTransforms()-1;
          tind >= 0; tind-- )
    {
    if( this->GetNthTransformToOptimize( tind ) )
      {
      if( ! this->GetNthTransform( tind )->HasLocalSupport() )
        {
        result = false;
        }
      }
    }
  return result;
}


template
<class TScalar, unsigned int NDimensions>
typename CompositeTransform< TScalar, NDimensions >::TransformQueueType &
CompositeTransform<TScalar, NDimensions>
::GetTransformsToOptimizeQueue() const
{
  /* Update the list of transforms to use for optimization only if
   the selection of transforms to optimize may have changed */
  if( this->GetMTime() > this->m_PreviousTransformsToOptimizeUpdateTime )
  {
      this->m_TransformsToOptimizeQueue.clear();
      for( size_t n=0; n < this->m_TransformQueue.size(); n++ )
      {
          /* Return them in the same order as they're found in the main list */
          if( this->GetNthTransformToOptimize( n ) )
          {
              this->m_TransformsToOptimizeQueue.push_back( m_TransformQueue[n] );
          }
      }
      this->m_PreviousTransformsToOptimizeUpdateTime = this->GetMTime();
  }
  return this->m_TransformsToOptimizeQueue;
}

template <class TScalarType, unsigned int NDimensions>
void
CompositeTransform<TScalarType, NDimensions>::
PrintSelf( std::ostream& os, Indent indent ) const
{
    Superclass::PrintSelf( os,indent );

    if( this->m_TransformQueue.empty() )
    {
        os << indent << "Transform queue is empty." << std::endl;
        return;
    }

    os  << indent << "TransformsToOptimizeFlags, begin() to end(): "
            << std::endl << indent << indent;
    for(  TransformsToOptimizeFlagsType::iterator
            it = this->m_TransformsToOptimizeFlags.begin();
            it != this->m_TransformsToOptimizeFlags.end(); it++ )
    {
        os << *it << " ";
    }
    os << std::endl;

    os << indent <<  "Transforms in queue, from begin to end:" << std::endl;
    typename TransformQueueType::const_iterator cit;
    for( cit = this->m_TransformQueue.begin();
            cit != this->m_TransformQueue.end(); ++cit )
    {
        os << indent << ">>>>>>>>>" << std::endl;
        (*cit)->Print( os, indent );
    }
    os << indent <<  "End of Transforms." << std::endl << "<<<<<<<<<<" << std::endl;

    os << indent <<  "TransformsToOptimize in queue, from begin to end:" << std::endl;
    for( cit = this->m_TransformsToOptimizeQueue.begin();
            cit != this->m_TransformsToOptimizeQueue.end(); ++cit )
    {
        os << indent << ">>>>>>>>>" << std::endl;
        (*cit)->Print( os, indent );
    }
    os << indent <<  "End of Transforms." << std::endl << "<<<<<<<<<<" << std::endl;

    os << indent << "PreviousTransformsToOptimizeUpdateTime: "
            <<  m_PreviousTransformsToOptimizeUpdateTime << std::endl;
    os << indent <<  "End of CompositeTransform." << std::endl << "<<<<<<<<<<" << std::endl;
}

} // namespace itk

#endif
