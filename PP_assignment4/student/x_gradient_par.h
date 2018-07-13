#ifndef _X_GRADIENT_PAR
#define _X_GRADIENT_PAR

#include <boost/gil/extension/io/jpeg_dynamic_io.hpp>
#include <omp.h>
#include"x_gradient.h"

using namespace boost::gil;

template <typename Out> struct halfdiff_cast_channels; // forward declaration

template <typename SrcView, typename DstView>
void x_gradient(const SrcView& src, const DstView& dst, int num_threads) {

	//TODO implement your solution in here.
      typedef typename channel_type<DstView>::type dst_channel_t;
      #pragma omp parallel
      {
        #pragma omp for schedule(dynamic)
        for (int y=0; y<src.height(); ++y) {
          typename SrcView::x_iterator src_it = src.row_begin(y);
          typename DstView::x_iterator dst_it = dst.row_begin(y);

          for (int x=1; x<src.width()-1; ++x) {
              static_transform(src_it[x-1], src_it[x+1], dst_it[x],
                               halfdiff_cast_channels<dst_channel_t>());
          }
      }
    }

}

#endif // !_X_GRADIENT_PAR_
