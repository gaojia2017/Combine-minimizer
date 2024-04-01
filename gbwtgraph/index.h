#ifndef GBWTGRAPH_CONSTRUCTION_H
#define GBWTGRAPH_CONSTRUCTION_H

#include <cstdlib>
#include <functional>

#include <omp.h>

//#include <boost/functional/hash.hpp>
#include <gbwtgraph/gbwtgraph.h>
#include <gbwtgraph/minimizer.h>

/*
  index.h: Minimizer index construction from GBWTGraph.
*/

namespace gbwtgraph
{
// inline size_t 
// combineHash(size_t a, size_t b){
//     size_t seed = 0;
//     boost::hash_combine(seed, a);
//     boost::hash_combine(seed, b);
//     return seed;
// }
//------------------------------------------------------------------------------

/*
  Index the haplotypes in the graph. Insert the minimizers into the provided index.
  Function argument get_payload is used to generate the payload for each position
  stored in the index.
  The number of threads can be set through OMP.
*/
template<class KeyType>
void
index_haplotypes(const GBWTGraph& graph, MinimizerIndex<KeyType>& index,
                 const std::function<payload_type(const pos_t&)>& get_payload)
{
  typedef typename MinimizerIndex<KeyType>::minimizer_type minimizer_type;
  std::cerr <<"index_haplotype****************\n";
  int threads = omp_get_max_threads();

  // Minimizer caching. We only generate the payloads after we have removed duplicate positions.
  std::vector<std::vector<std::pair<minimizer_type, pos_t>>> cache(threads);
  constexpr size_t MINIMIZER_CACHE_SIZE = 1024;
  
  auto flush_cache = [&](int thread_id)
  {
    std::vector<std::pair<minimizer_type, pos_t>>& current_cache = cache[thread_id];
    // int gap = 1;
    // for(size_t i = 0; i < current_cache.size() - gap; i++){
    //     std::get<0>(current_cache[i]).hash = combineHash(std::get<0>(current_cache[i]).hash, std::get<0>(current_cache[i + gap]).hash);
    // }
    gbwt::removeDuplicates(current_cache, false);
    
    std::vector<payload_type> payload;
    payload.reserve(current_cache.size());
    for(size_t i = 0; i < current_cache.size(); i++) { payload.push_back(get_payload(current_cache[i].second)); }
    #pragma omp critical (minimizer_index)
    {
      for(size_t i = 0; i < current_cache.size(); i++)
      {
        // if( current_cache[i].first.hash == 14967316640526472684){
        //   std::cerr<<"debug "<<current_cache[i].first.hash <<"pos"<<current_cache[i].second <<std::endl;
        // }
        index.insert(current_cache[i].first, current_cache[i].second, payload[i]);
      }
    }
    cache[thread_id].clear();
  };

  // Minimizer finding.
  auto  find_minimizers = [&](const std::vector<handle_t>& traversal, const std::string& seq)
  {
    //std::vector<minimizer_type> minimizers = index.minimizers_gap(seq, index.w()); // Calls syncmers() when appropriate.
    std::vector<minimizer_type> minimizers = index.minimizers_combine_gap(seq, index.window_bp());//稀疏+组合一件套
    auto iter = traversal.begin();
    size_t node_start = 0;
    int thread_id = omp_get_thread_num();
    for(minimizer_type& minimizer : minimizers)
    {
      if(minimizer.empty()) { continue; }

      // Find the node covering minimizer starting position.
      size_t node_length = graph.get_length(*iter);
      while(node_start + node_length <= minimizer.offset)
      {
        node_start += node_length;
        ++iter;
        node_length = graph.get_length(*iter);
      }
      pos_t pos { graph.get_id(*iter), graph.get_is_reverse(*iter), minimizer.offset - node_start };
      if(minimizer.is_reverse) { pos = reverse_base_pos(pos, node_length); }
      if(!Position::valid_offset(pos))
      {
        #pragma omp critical (cerr)
        {
          std::cerr << "index_haplotypes(): Node offset " << offset(pos) << " is too large" << std::endl;
        }
        std::exit(EXIT_FAILURE);
      }
      //因为key增加了一个hash字段，所以需要添加赋值语句
      //minimizer.key.hash = minimizer.hash;

      cache[thread_id].emplace_back(minimizer, pos);
    }
    if(cache[thread_id].size() >= MINIMIZER_CACHE_SIZE) { flush_cache(thread_id); }
  };
//组合索引
   auto find_minimizers1 = [&](const std::vector<handle_t>& traversal, const std::string& seq)
  {
    //std::cerr<<"seq.length"<<seq.size()<<std::endl;
    //std::vector<minimizer_type> minimizers = index.minimizers_gap(seq, 8); // Calls syncmers() when appropriate.
    //std::vector<minimizer_type> minimizers = index.minimizers(seq);
    std::vector<minimizer_type> minimizers = index.minimizers_e2(seq, index.window_bp());
    //assert(minimizers.size() >= 2);
    auto iter = traversal.begin();
    size_t id = graph.get_id(*iter);
    //std::cerr<<"**********minimizers.size()"<<minimizers.size()<<"and id "<<graph.get_id(*iter) <<std::endl;
    // if(id == 37797){
    //   std::cerr<<"seq = "<<seq<<std::endl;
    // }
    // if(minimizers.size() == 1){
        
    //     std::cerr<<"size = 1 and id "<<graph.get_id(*iter)<<"seq length = "<<seq.size()<<"node length = "<<graph.get_length(*iter)<<"\n";
    // }
    size_t node_start = 0;
    int thread_id = omp_get_thread_num();
    int gap = 1;
    //gap = 1
    for(int i = 0; i < minimizers.size() - gap; i++)
    {
      minimizer_type& minimizer = minimizers[i];
      if(minimizer.empty()) { continue; }
      // //即最后一个minimzer和自己组合
      // int index = (i < minimizers.size() - gap) ? i + gap : i;
      // minimizer.hash = combineHash(minimizer.hash, minimizers[index].hash);

      //单个minimzer的情况不再组合
      // if(i < minimizers.size() - gap){
      //     minimizer.hash = combineHash(minimizer.hash, minimizers[i + gap].hash);
      //     //minimizer.key.hash = minimizer.hash;
      // }
      // Find the node covering minimizer starting position.
      size_t node_length = graph.get_length(*iter);
      while(node_start + node_length <= minimizer.offset)
      {
        node_start += node_length;
        ++iter;
        node_length = graph.get_length(*iter);
      }
      //std::cerr<<"id = "<<graph.get_id(*iter)<<std::endl;
      pos_t pos { graph.get_id(*iter), graph.get_is_reverse(*iter), minimizer.offset - node_start };
      //std::cerr<<"id = "<<graph.get_id(*iter)<<"offset = "<<minimizer.offset - node_start<<"hash"<<minimizer.hash<<std::endl;
      if(minimizer.is_reverse) { pos = reverse_base_pos(pos, node_length); }
      if(!Position::valid_offset(pos))
      {
        #pragma omp critical (cerr)
        {
          std::cerr << "index_haplotypes(): Node offset " << offset(pos) << " is too large" << std::endl;
        }
        std::exit(EXIT_FAILURE);
      }
      cache[thread_id].emplace_back(minimizer, pos);
    }
    if(cache[thread_id].size() >= MINIMIZER_CACHE_SIZE) { 
      
      flush_cache(thread_id); 
    }
  };

  // 组合稀疏索引
  auto  find_minimizers2 = [&](const std::vector<handle_t>& traversal, const std::string& seq)
  {
    
    //std::vector<minimizer_type> minimizers = index.minimizers_e2(seq, index.window_bp());//稀疏
    //std::vector<minimizer_type> minimizers = index.minimizers(seq);//
    std::vector<minimizer_type> minimizers = index.minimizers_combine_gap(seq, index.w());//稀疏加组合
    auto iter = traversal.begin();
    size_t node_start = 0;
    int thread_id = omp_get_thread_num();
    for(minimizer_type& minimizer : minimizers)
    {
      if(minimizer.empty()) { continue; }
      // Find the node covering minimizer starting position.
      size_t node_length = graph.get_length(*iter);
      while(node_start + node_length <= minimizer.offset)
      {
        node_start += node_length;
        ++iter;
        node_length = graph.get_length(*iter);
      }
      pos_t pos { graph.get_id(*iter), graph.get_is_reverse(*iter), minimizer.offset - node_start };
      if(minimizer.is_reverse) { pos = reverse_base_pos(pos, node_length); }
      if(!Position::valid_offset(pos))
      {
        #pragma omp critical (cerr)
        {
          std::cerr << "index_haplotypes(): Node offset " << offset(pos) << " is too large" << std::endl;
        }
        std::exit(EXIT_FAILURE);
      }
      //因为key增加了一个hash字段，所以需要添加赋值语句
      //minimizer.key.hash = minimizer.hash;

      cache[thread_id].emplace_back(minimizer, pos);
    }
    if(cache[thread_id].size() >= MINIMIZER_CACHE_SIZE) { flush_cache(thread_id); }
  };

  
  /*
    Index the minimizers.
    We do a lot of redundant work by traversing both orientations and finding almost the same minimizers
    in each orientation. If we consider only the windows starting in forward (reverse) orientation,
    we may skip windows that cross from a reverse node to a forward node (from a forward node to a
    reverse node).
  */
  for_each_haplotype_window(graph, index.window_bp(),index.w(), find_minimizers2, (threads > 1));
  //for_each_haplotype_window(graph, index.window_bp(),index.w(), find_minimizers1, false);
  
  for(int thread_id = 0; thread_id < threads; thread_id++) { flush_cache(thread_id); }
}
  
//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_CONSTRUCTION_H
