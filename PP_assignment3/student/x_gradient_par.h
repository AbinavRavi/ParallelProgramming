#ifndef _X_GRADIENT_PAR
#define _X_GRADIENT_PAR

#include <boost/gil/extension/io/jpeg_dynamic_io.hpp>
#include <omp.h>
#include"x_gradient.h"
#include <thread>
#include <iostream>
#include <vector>

using namespace boost::gil;

template <typename Out> struct halfdiff_cast_channels; // forward declaration
template <typename SrcView, typename DstView>

void parallel_x_gradient(const SrcView& src, const DstView& dst, int i, int start, int end){
  typedef typename channel_type<DstView>::type dst_channel_t;

  for (int y = start*i; y < end; ++y)
  {
      typename SrcView::x_iterator src_it = src.row_begin(y);
      typename DstView::x_iterator dst_it = dst.row_begin(y);

      for (int x = 1; x < src.width() - 1; ++x)
      {
          static_transform(src_it[x - 1], src_it[x + 1], dst_it[x],
                           halfdiff_cast_channels<dst_channel_t>());
      }
  }

}
template <typename SrcView, typename DstView>
void x_gradient(const SrcView& src, const DstView& dst, int num_threads) {

  std::vector<std::thread> t;
	int end = 0;
	int start = src.height()/num_threads;


	// TODO put your solution in here.
//create thread
  for (int i = 0; i< num_threads; i++){
    end = (i < num_threads-1)?(src.height()/num_threads)*(i+1):(src.height()/num_threads)*(i+1)+(src.height()%num_threads);
    t.push_back(std::thread(parallel_x_gradient<SrcView, DstView>, src, dst, i, start, end));
  }
  for (auto &th : t )
	{
		th.join();
	}
}

#endif // !_X_GRADIENT_PAR_
