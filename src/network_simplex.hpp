#include <type_traits>
#include <vector>
#include "full_bipartitegraph.h"
#include "network_simplex_simple.h"

auto network_simplex(const auto &p, const auto &f, auto dist) {
	size_t n = f.size();
	assert(n == p.size());
	typedef typename std::remove_reference_t<decltype(p)>::value_type PointType;
	typedef typename std::remove_reference_t<decltype(f)>::value_type FlowType;
	typedef std::result_of_t<decltype(dist)(const PointType&, const PointType&)> CostType;
	std::vector <FlowType> supply_list, demand_list;
	supply_list.reserve(n); demand_list.reserve(n);
	for (size_t i = 0; i < n; ++i)
		if (f[i] > 0)
			supply_list.push_back(f[i]);
		else if (f[i] < 0)
			demand_list.push_back(f[i]);
	lemon::FullBipartiteDigraph di(supply_list.size(), demand_list.size());
	lemon::NetworkSimplexSimple<lemon::FullBipartiteDigraph, FlowType, CostType, ssize_t> net(di, true, supply_list.size() + demand_list.size(), supply_list.size() * demand_list.size());
	for (size_t i = 0, ID = 0; i < n; ++i)
		if (f[i] > 0)
			for (size_t j = 0; j < n; ++j)
				if (f[j] < 0)
					net.setCost(di.arcFromId(ID++), dist(p[i], p[j]));
	net.supplyMap(supply_list.data(), supply_list.size(), demand_list.data(), demand_list.size());
	net.run();
	return net.totalCost();
}
