#include "MEMF.h"

#ifdef _MSC_VER
#pragma warning(disable: 4661)
#endif

// Instantiations: <nodeid, labelid, captype>
// IMPORTANT: 
//    nodeid must be integer type, can be signed or unsigned
//	  labelid must be a signed integer type 
//    captype must be a signed data type

template class MEMF<int,short,int>;
template class MEMF<int,int,int>;
template class MEMF<unsigned int,short,int>;
template class MEMF<unsigned int,unsigned short,int>;
template class MEMF<unsigned int,int,int>;
template class MEMF<unsigned long,int,int>;
template class MEMF<int,short,float>;
template class MEMF<unsigned int,short,float>;
template class MEMF<unsigned int,int,float>;
template class MEMF<unsigned long,int,double>;