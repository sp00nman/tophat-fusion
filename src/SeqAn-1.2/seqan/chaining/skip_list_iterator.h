 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: skip_list_iterator.h 3420 2009-02-12 12:10:09Z doering@PCPOOL.MI.FU-BERLIN.DE $
 ==========================================================================*/

/*	2006 Hendrik Woehrle
*
*	Skip List Iterator
*
*/

//SEQAN_NO_DDDOC: do not generate documentation for this file


#ifndef SEQAN_HEADER_SKIP_LIST_ITER_H
#define SEQAN_HEADER_SKIP_LIST_ITER_H

namespace seqan
{

//////////////////////////////////////////////////////////////////////////////
// Iterator:
// 
//	bidirectional iterator for dynamic skip list (base layer is a list)
//	random access iterator for static skip list (base layer is an array)
//
//////////////////////////////////////////////////////////////////////////////

		// the standard iterator is a pointer to the entrys in the base layer
	template < typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Iterator< SkipList< TObject, TModus, TSpec, TStructuring >, Standard >
	{
		typedef SkipBaseElement< TObject, TModus, TSpec, TStructuring > * Type;
	};

		// metatypes of the iterator
	template < typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct Reference< Iter< SkipList< TObject, TModus, TSpec, TStructuring >, Standard > >
	{
		typedef TObject & Type;
	};

	template < typename TObject, typename TModus, typename TSpec, typename TStructuring >
	struct GetValue< Iter< SkipList< TObject, TModus, TSpec, TStructuring >, Standard > >
	{
		typedef TObject Type;
	};

	
	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	goNext( SkipBaseElement< TObject, TModus, TSpec, TStructuring > *& it );

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring > inline
	void
	goPrevious( SkipBaseElement< TObject, TModus, TSpec, TStructuring > *& it );

	
//////////////////////////////////////////////////////////////////////////////
//
//	generic iterator functions and overloading of functions
//
//////////////////////////////////////////////////////////////////////////////

	template< typename TType, typename TIterType > inline
	typename Key< typename Value< Iter< TType, TIterType > >::Type >::Type
	key( Iter< TType, TIterType > & it )
	{
		return key( value( it ) );
	}

	template< typename TType, typename TIterType, typename TParam > inline
	typename Key< typename Value< Iter< TType, TIterType > >::Type >::Type
	key( Iter< TType, TIterType > & it,
			TParam & param )
	{
		return key( value( it ), param );
	}

	template < typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TValue>
	inline void
	assignValue( SkipBaseElement< TObject, TModus, TSpec, TStructuring > * me,
					TValue const & _value)
	{
	SEQAN_CHECKPOINT
		// do nothing
	} 

	//const version for iterators as targets
	template < typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TValue>
	inline void
	assignValue( SkipBaseElement< TObject, TModus, TSpec, TStructuring > const * me,
					TValue const & _value)
	{
	SEQAN_CHECKPOINT
		// do nothing
	} 


	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TValue >
	inline void
	moveValue( SkipBaseElement< TObject, TModus, TSpec, TStructuring > *& it,
				TValue const & _value)
	{
	SEQAN_CHECKPOINT
		// do nothing
	}

	template< typename TObject, typename TModus, typename TSpec, typename TStructuring, typename TValue >
	inline void
	moveValue( SkipBaseElement< TObject, TModus, TSpec, TStructuring > const *& it,
				TValue const & _value)
	{
	SEQAN_CHECKPOINT
		// do nothing
	}



	template <typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec>
	inline typename Position<Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const>::Type 
	position(Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & me)
	{
	SEQAN_CHECKPOINT
		typename Position<Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const>::Type pos = 0;
		Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > it = begin( container(me) );
		while( it != me )
		{
			goNext( it );
			++pos;
		}
		return pos;
	}

	//____________________________________________________________________________

	template <typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec, typename TContainer2>
	inline typename Position<Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const>::Type 
	position(Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & me,
			 TContainer2 const &)
	{
	SEQAN_CHECKPOINT
		return position( me );
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator +
	//////////////////////////////////////////////////////////////////////////////

	template < typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec, typename TIntegral>
	inline Iter< SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> >  
	operator + (Iter< SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & left,
				TIntegral right)
	{
	SEQAN_CHECKPOINT
		Iter< SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > buffer(container(left), hostIterator(left) );
		while( right > 0 )
		{
			goNext( buffer );
			--right;
		}
		return buffer;
	}
	template < typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec, typename TIntegral>
	inline Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> >  
	operator + (TIntegral left,
				Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & right)
	{
	SEQAN_CHECKPOINT
		Iter< SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > buffer(container(left), hostIterator(left) );
		while( right > 0 )
		{
			goNext( buffer );
			--right;
		}
		return buffer;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator +=
	//////////////////////////////////////////////////////////////////////////////

	template < typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec, typename TIntegral>
	inline Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > &
	operator += (Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > & left,
				 TIntegral right)
	{
	SEQAN_CHECKPOINT
		while( right > 0 )
		{
			goNext( left );
			--right;
		}
		return left;
	}

	//////////////////////////////////////////////////////////////////////////////
	// operator -
	//////////////////////////////////////////////////////////////////////////////

	template < typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec, typename TIntegral>
	inline Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> >  
	operator - (Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & left,
				TIntegral right)
	{
	SEQAN_CHECKPOINT
		Iter< SkipList< TObject, SkipListDynamic, TSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > buffer(container(left), hostIterator(left) );
		while( right > 0 )
		{
			goPrevious( buffer );
			--right;
		}
		return buffer;
	}

	//____________________________________________________________________________

	// ???
	//template < typename TObject, typename SkipListDynamic, typename TSpec, typename TStructuring, typename TIterator, typename TSpec>
	//inline typename Position<Iter<TContainer, AdaptorIterator<TIterator, TSpec> > >::Type  
	//operator - (Iter<SkipList< TObject, SkipListDynamic, TSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & left,
	//			Iter<SkipList< TObject, SkipListDynamic, TSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > const & right)
	//{
	//SEQAN_CHECKPOINT
	//	return hostIterator(left) - hostIterator(right);
	//}

	//////////////////////////////////////////////////////////////////////////////
	// operator -=
	//////////////////////////////////////////////////////////////////////////////

	template < typename TObject, typename TSLSpec, typename TStructuring, typename TIterator, typename TSpec, typename TIntegral>
	inline Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > &
	operator -= (Iter<SkipList< TObject, SkipListDynamic, TSLSpec, TStructuring >, AdaptorIterator<TIterator, TSpec> > & left,
				 TIntegral right)
	{
	SEQAN_CHECKPOINT
		while( right > 0 )
		{
			goPrevious( left );
			--right;
		}
		return left;;
	}

} // namespace ...

//////////////////////////////////////////////////////////////////////////////


#endif //SEQAN_HEADER_SKIP_LIST_ITER_H
