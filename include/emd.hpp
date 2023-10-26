#include <concepts>
#include <random>
#include <vector>
#include <cassert>
#include <algorithm>
#include <unordered_map>
#include <memory>
#include <numeric>
#include <limits>
#include <functional>
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/priority_queue.hpp>
#include "lct.hpp"
#include "ett.hpp"
#include "network_simplex/full_bipartitegraph.h"
#include "network_simplex/network_simplex_simple.h"
#ifdef TIME_LOG
#include <fstream>
#include <chrono>
#endif
#if defined(DEBUG) or defined(PYBIND11_DETAILED_ERROR_MESSAGES)
#include <iostream>
#endif

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral  _CostType, std::integral IndexType>
class emd {
	typedef std::common_type_t<_FlowType, _CostType> _RetType;
	ett<_CostType, IndexType> cost_tree;
	std::vector <std::unique_ptr <lct::node<_FlowType> > > lctvertex;
	std::unordered_map <lct::node<_FlowType>*, std::tuple <IndexType, IndexType, typename ett<_CostType, IndexType>::pnode*, typename ett<_CostType, IndexType>::pnode*> > edges;
	std::vector <_PointType> points;
	std::vector <_FlowType> flows;
	std::function <_CostType(const _PointType&, const _PointType&)> dist_func;
	_RetType total_cost;
	std::vector <IndexType> erased_nodes;
	[[nodiscard("Check returns the minimum of cost_tree and has no side effects")]]
	std::tuple <_CostType, IndexType, IndexType> check();
public:
#ifndef NDEBUG
	static std::mt19937_64 rng;
#endif
	emd(emd && val);
	void modify(IndexType idx, const _PointType &new_position);
	void send_flow(IndexType from, IndexType to, _FlowType Delta);
	void erase(IndexType idx);
	IndexType insert_supply(const _PointType &p);
	IndexType insert_demand(const _PointType &p);
	void reserve(size_t capacity);
	template <bool flag = true>
	IndexType insert(const _PointType &p);
	emd(const auto &p, const auto &f, const auto &dist);
	emd operator = (const emd &val) = delete;
	emd&& operator = (emd && val);
	_RetType cost() { return total_cost; }
	size_t size() const;
	~emd();
};


#ifndef NDEBUG
template <typename _PointType, std::signed_integral _FlowType, std::signed_integral _CostType, std::integral IndexType>
std::mt19937_64 emd<_PointType, _FlowType, _CostType, IndexType>::rng((std::random_device())());
#endif

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral _CostType, std::integral IndexType>
void emd<_PointType, _FlowType, _CostType, IndexType>::reserve(size_t capacity) { cost_tree.reserve(capacity); }

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral _CostType, std::integral IndexType>
void emd<_PointType, _FlowType, _CostType, IndexType>::send_flow(IndexType from, IndexType to, _FlowType Delta) {
	if (Delta < 0) {
		Delta = -Delta;
		std::swap(from, to);
	}
	else if (Delta == 0)
		return;
	_CostType dist = cost_tree(from, to);
	while (true) {
		auto [e, f] = lct::sendflow(lctvertex[from].get(), lctvertex[to].get(), Delta);
		total_cost += dist * f;
		Delta -= f;
		if (Delta == 0)
			return;
		auto it = edges.find(e);
		assert(it != edges.end());
		auto [u, v, puu, pvv] = it->second;
		lct::cut(lctvertex[u].get(), lctvertex[v].get());
		edges.erase(it);
		auto [minv, U, V] = cost_tree.template cut<1>(puu, pvv);
		dist += minv;
		edges.emplace(lct::link(lctvertex[U].get(), lctvertex[V].get(), _FlowType()), std::make_tuple(U, V, pvv, puu));
	}
}

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral _CostType, std::integral IndexType>
emd<_PointType, _FlowType, _CostType, IndexType>&& emd<_PointType, _FlowType, _CostType, IndexType>::operator = (emd <_PointType, _FlowType, _CostType, IndexType> && val) {
	cost_tree = std::move(val.cost_tree);
	lctvertex = std::move(val.lctvertex);
	edges = std::move(val.edges);
	points = std::move(val.points);
	flows = std::move(val.flows);
	dist_func = std::move(val.dist_func);
	erased_nodes = std::move(val.erased_nodes);
	total_cost = val.total_cost;
	return std::move(*this);
}

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral _CostType, std::integral IndexType>
emd<_PointType, _FlowType, _CostType, IndexType>::~emd() {
	for (auto [key, val] : edges)
		delete key;
}

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral _CostType, std::integral IndexType>
emd<_PointType, _FlowType, _CostType, IndexType>::emd(const auto &p, const auto &f, const auto &dist): cost_tree(), lctvertex(f.size()), edges(), points(p), flows(f), dist_func(dist), total_cost() {
#ifdef NDEBUG
	static std::mt19937_64 rng((std::random_device())());
#endif
	assert(points.size() == flows.size());
	assert(std::accumulate(flows.begin(), flows.end(), _FlowType()) == _FlowType());
	const size_t n = f.size(), m = 2 * n - 2;
	assert(n > 0);
	assert(n < std::numeric_limits<IndexType>::max());
	cost_tree.potential = (_CostType*)std::malloc(n * sizeof(_CostType));
	cost_tree.distance = (_CostType*)std::malloc(n * (n - 1) / 2 * sizeof(_CostType));
	cost_tree.alloc->_Cap = cost_tree.alloc->_Size = n;
	for (size_t i = 1; i < n; ++i) {
		for (size_t j = 0; j < i; ++j) {
			cost_tree.distance[i * (i - 1) / 2 + j] = dist(points[i], points[j]);
		}
	}
#ifdef SINKHORN
	std::vector <IndexType> a, b, newidx(n);
	a.reserve(n); b.reserve(n);
	for (size_t i = 0; i < flows.size(); ++i)
		if (flows[i] > 0) {
			newidx[i] = a.size();
			a.push_back(i);
		}
		else {
			newidx[i] = b.size();
			b.push_back(i);
		}
	std::vector <double> sk(a.size() * b.size()), seq(a.size() * b.size());
	for (size_t i = 0; i < a.size(); ++i)
		for (size_t j = 0; j < b.size(); ++j)
			sk[i * b.size() + j] = cost_tree[a[i]][b[j]];
	std::copy(sk.begin(), sk.end(), seq.begin());
	std::nth_element(seq.begin(), seq.begin() + a.size() + b.size() - 2, seq.end(), std::greater <double>());
	double lambda = 60 / seq[a.size() + b.size() - 2];
	std::vector <double> ().swap(seq);
	std::for_each(sk.begin(), sk.end(), [lambda](auto &x) { x = std::exp(-x * lambda); });
	for (bool flag = 1; flag; ) {
		flag = 0;
#pragma omp parallel for
		for (size_t i = 0; i < a.size(); ++i) {
			double sum = 0;
			for (size_t j = 0; j < b.size(); ++j)
				sum += sk[i * b.size() + j];
			sum = f[a[i]] / sum;
			if (std::abs(sum - 1) > 1e-6) {
				flag = 1;
				for (size_t j = 0; j < b.size(); ++j)
					sk[i * b.size() + j] *= sum;
			}
		}
#pragma omp parallel for
		for (size_t j = 0; j < b.size(); ++j) {
			double sum = 0;
			for (size_t i = 0; i < a.size(); ++i)
				sum += sk[i * b.size() + j];
			sum = -f[b[j]] / sum;
			if (std::abs(sum - 1) > 1e-6) {
				flag = 1;
				for (size_t i = 0; i < a.size(); ++i)
					sk[i * b.size() + j] *= sum;
			}
		}
	}
	std::vector <typename __gnu_pbds::priority_queue <std::pair <double, IndexType> >::point_iterator> it(n);
	std::vector <uint8_t> vis(n);
	std::vector <IndexType> deg(n), prec(n);
	std::fill(prec.begin(), prec.end(), ~IndexType());
	__gnu_pbds::priority_queue <std::pair <double, IndexType> > pq;
	it[0] = pq.push(std::make_pair(0, 0)); prec[0] = 0;
	for (size_t i = 0; i < n; ++i) {
		assert(!pq.empty());
		IndexType u = pq.top().second; pq.pop(); vis[u] = true;
		for (IndexType v = 0; v < n; ++v)
			if (flows[u] <= 0 && flows[v] > 0) {
				double w = sk[newidx[v] * b.size() + newidx[u]];
				if (prec[v] == ~IndexType()) {
					it[v] = pq.push(std::make_pair(w, v));
					prec[v] = u;
				}
				else if (vis[v] == false && it[v]->first < w) {
					pq.modify(it[v], std::make_pair(w, v));
					prec[v] = u;
				}
			}
			else if (flows[v] <= 0 && flows[u] > 0) {
				double w = sk[newidx[u] * b.size() + newidx[v]];
				if (prec[v] == ~IndexType()) {
					it[v] = pq.push(std::make_pair(w, v));
					prec[v] = u;
				}
				else if (vis[v] == false && it[v]->first < w) {
					pq.modify(it[v], std::make_pair(w, v));
					prec[v] = u;
				}
			}
	}
	assert(pq.empty());
	std::vector <IndexType> stk;
	std::vector <uint8_t>().swap(vis);
#else
	std::vector <IndexType> deg(n), prec(n, -1), newidx(n);
	std::vector <_FlowType> supply_list, demand_list;
	supply_list.reserve(n); demand_list.reserve(n);
	for (size_t i = 0; i < n; ++i)
		if (flows[i] > 0) {
			newidx[i] = supply_list.size();
			supply_list.push_back(flows[i]);
		}
		else {
			newidx[i] = demand_list.size();
			demand_list.push_back(flows[i]);
		}
	if (supply_list.empty()) {
		supply_list.push_back(demand_list.back());
		demand_list.pop_back();
	}
	const auto original_num_threads = omp_get_max_threads();
	omp_set_num_threads(1);
	lemon::FullBipartiteDigraph di(supply_list.size(), demand_list.size());
	lemon::NetworkSimplexSimple <lemon::FullBipartiteDigraph, _FlowType, _CostType, ssize_t> net(di, true, flows.size(), supply_list.size() * demand_list.size());
	for (size_t i = 0, ID = 0; i < n; ++i)
		if (flows[i] > 0) {
			for (size_t j = 0; j < n; ++j)
				if (flows[j] <= 0) {
					net.setCost(di.arcFromId(ID++), cost_tree[i][j]);
				}
		}
	net.supplyMap(&supply_list[0], supply_list.size(), &demand_list[0], demand_list.size());
	net.run();
	omp_set_num_threads(original_num_threads);
	std::vector <IndexType> stk; stk.reserve(n);
	std::vector <IndexType> label(n, -1);
	std::vector <std::tuple <_CostType, IndexType, IndexType> > ga(n * n, std::tuple <_CostType, IndexType, IndexType> (std::numeric_limits<_CostType>::max(), -1, -1));
	IndexType num_clusters = 0;
	for (IndexType i = 0; i < n; ++i)
		if (label[i] == IndexType(-1)) {
			ga[num_clusters * n + num_clusters] = std::make_tuple(0, num_clusters, num_clusters);
			stk.push_back(i); label[i] = num_clusters;
			cost_tree.potential[i] = 0;
			while (!stk.empty()) {
				IndexType u = stk.back(); stk.pop_back();
				for (IndexType v = 0; v < n; ++v)
					if (label[v] == IndexType(-1)) {
						if (flows[u] > 0 && flows[v] <= 0 && net.flow(di.arcFromId(newidx[u] * demand_list.size() + newidx[v])) > 0) {
							assert(flows[v] < 0);
							label[v] = num_clusters;
							stk.emplace_back(v);
							cost_tree.potential[v] = cost_tree.potential[u] + cost_tree[u][v];
						}
						else if (flows[v] > 0 && flows[u] <= 0 && net.flow(di.arcFromId(newidx[v] * demand_list.size() + newidx[u])) > 0) {
							assert(flows[u] < 0);
							label[v] = num_clusters;
							stk.emplace_back(v);
							cost_tree.potential[v] = cost_tree.potential[u] - cost_tree[u][v];
						}
					}
					else if (label[v] != num_clusters) {
						if (flows[u] >= 0 && flows[v] <= 0)
							ga[num_clusters * n + label[v]] = std::min(ga[num_clusters * n + label[v]], std::make_tuple(cost_tree.get_value({ .i = u, .j = v }), u, v));
						if (flows[v] >= 0 && flows[u] <= 0)
							ga[label[v] * n + num_clusters] = std::min(ga[label[v] * n + num_clusters], std::make_tuple(cost_tree.get_value({ .i = v, .j = u }), v, u));
					}
			}
			++num_clusters;
		}
	std::vector <_CostType> newp(num_clusters, std::numeric_limits<_CostType>::max());
	std::vector <size_t> prec_clusters(num_clusters, -1);
	newp[0] = 0;
	std::list <IndexType> que;
	std::vector <typename std::list <IndexType>::iterator> it(num_clusters, que.end());
	que.emplace_back(0); it[0] = que.begin(); prec_clusters[0] = 0;
	while (!que.empty()) {
		IndexType u = que.front(); que.pop_front(); it[u] = que.end();
		for (IndexType v = 0; v < num_clusters; ++v) {
			const size_t e = u * n + v;
			assert(std::get<0>(ga[e]) != std::numeric_limits<_CostType>::max());
			if (newp[v] > newp[u] + std::get<0>(ga[e])) {
				newp[v] = newp[u] + std::get<0>(ga[e]);
				prec_clusters[v] = e;
				if (it[v] != que.end())
					que.erase(it[v]);
				if (!que.empty() && newp[v] <= newp[que.front()]) {
					que.emplace_front(v);
					it[v] = que.begin();
				}
				else {
					que.emplace_back(v);
					it[v] = --que.end();
				}
			}
		}
	}
	for (IndexType i = 0; i < num_clusters; ++i) {
		IndexType s = std::get<2>(ga[prec_clusters[i]]);
		prec[s] = std::get<1>(ga[prec_clusters[i]]);
		assert(stk.empty());
		stk.push_back(s);
		while (!stk.empty()) {
			IndexType u = stk.back(); stk.pop_back();
			for (IndexType v = 0; v < n; ++v)
				if (prec[v] == IndexType(-1) && ((flows[u] > 0 && flows[v] <= 0 && net.flow(di.arcFromId(newidx[u] * demand_list.size() + newidx[v])) > 0) || (flows[v] > 0 && flows[u] <= 0 && net.flow(di.arcFromId(newidx[v] * demand_list.size() + newidx[u])) > 0))) {
					prec[v] = u;
					stk.push_back(v);
				}
		}
	}
#endif
	std::vector <std::tuple <IndexType, IndexType, _FlowType, size_t, size_t> > ge(n);
	std::vector <IndexType> tour;
	std::vector <std::vector <IndexType> > gl(n);
	std::vector <_FlowType> tmpflow(n);
	tour.reserve(m);
	stk.reserve(n);
	std::copy(flows.begin(), flows.end(), tmpflow.begin());
	for (size_t i = 1; i < n; ++i) {
		++deg[prec[i]];
		gl[prec[i]].push_back(i);
	}
	for (size_t i = 0; i < n; ++i)
		if (deg[i] == 0)
			stk.push_back(i);
	for (size_t _ = 1; _ < n; ++_) {
		assert(!stk.empty());
		IndexType u = stk.back(); stk.pop_back();
		IndexType v = prec[u];
		std::get<0>(ge[u]) = v;
		std::get<1>(ge[u]) = u;
		std::get<2>(ge[u]) = -tmpflow[u];
		tmpflow[v] += tmpflow[u];
		if (!--deg[v])
			stk.push_back(v);
	}
	assert(stk.size() == 1);
	std::vector <IndexType>().swap(stk);
	std::vector <std::pair <IndexType, IndexType> > dfsstk;
	dfsstk.reserve(n);
	dfsstk.emplace_back(0, 0);
	cost_tree.potential[0] = 0;
	for (size_t _ = 2; _ < 2 * n; ++_) {
		auto &[u, i] = dfsstk.back();
		if (i == gl[u].size()) {
			std::get<4>(ge[u]) = tour.size();
			dfsstk.pop_back();
			tour.emplace_back(dfsstk.back().first);
		}
		else {
			auto v = gl[u][i++];
			assert(std::get<0>(ge[v]) == u);
			assert(std::get<1>(ge[v]) == v);
			auto w = std::get<2>(ge[v]);
			total_cost += std::abs(w) * cost_tree[u][v];
			cost_tree.potential[v] = cost_tree.potential[u] + (w > 0 || (w == 0 && flows[u] > flows[v]) ? cost_tree[u][v] : -cost_tree[u][v]);
			dfsstk.emplace_back(v, 0);
			std::get<3>(ge[v]) = tour.size();
			tour.emplace_back(v);
		}
	}
	assert(dfsstk.size() == 1);
	assert(dfsstk.back().first == 0 && dfsstk.back().second == gl[0].size());
	ge.erase(ge.begin());
	std::for_each(lctvertex.begin(), lctvertex.end(), [](auto &x) { x.reset(lct::add_node<_FlowType>()); });
	edges.reserve(n);
	typename ett<_CostType, IndexType>::pnode *PNodeChunk = (typename ett<_CostType, IndexType>::pnode*)malloc(m * sizeof(typename ett<_CostType, IndexType>::pnode));
	cost_tree.alloc->top = PNodeChunk;
	cost_tree.alloc->_AllocatedPNodePool.emplace_back(PNodeChunk);
	for (auto [u, v, w, p1, p2] : ge) {
		if (w > 0 || (w == 0 && flows[u] >= flows[v]))
			edges.emplace(lct::link(lctvertex[u].get(), lctvertex[v].get(), w), std::make_tuple(u, v, &PNodeChunk[p1], &PNodeChunk[p2]));
		else
			edges.emplace(lct::link(lctvertex[v].get(), lctvertex[u].get(), _FlowType(-w)), std::make_tuple(v, u, &PNodeChunk[p2], &PNodeChunk[p1]));
	}
	std::vector <typename ett<_CostType, IndexType>::HeightType> ht(m);
	std::for_each(ht.begin(), ht.end(), [](auto &x) { while ((x = std::countl_zero(rng()) + 1) > 64); });
	size_t maxit = std::max_element(ht.begin(), ht.end()) - ht.begin();
	ht.insert(ht.end(), ht.begin(), ht.begin() + maxit);
	ht.erase(ht.begin(), ht.begin() + maxit);
	assert(ht.size() == m);
	assert(*std::max_element(ht.begin(), ht.end()) == ht.front());
	typename ett<_CostType, IndexType>::HeightType maxht = ht.front();
	typename ett<_CostType, IndexType>::node *left[64], *above[64], *const basis = cost_tree.alloc->node_malloc(maxht);
	std::tuple <_CostType, IndexType, IndexType> minv[65];
	PNodeChunk[0].ptr = basis;
	PNodeChunk[m - 1].right = &PNodeChunk[0];
	std::iota(left, left + maxht, basis);
	PNodeChunk[0].ht = maxht;
	std::fill(minv, minv + maxht + 1, std::tuple <_CostType, IndexType, IndexType> (0, tour.front(), tour.front()));
	for (size_t i = 1; i < m; ++i) {
		typename ett<_CostType, IndexType>::node *const ptr = cost_tree.alloc->node_malloc(ht[i]);
		std::tuple <_CostType, IndexType, IndexType> val;
		std::get<0>(val) = cost_tree.get_value({ .i = std::get<1>(val) = tour.front(), .j = std::get<2>(val) = tour[i] });
		for (typename ett<_CostType, IndexType>::HeightType j = 0; j < ht[i]; ++j) {
			if (minv[j + 1] > minv[j])
				minv[j + 1] = minv[j];
			left[j]->data = { .i = std::get<1>(minv[j]), .j = std::get<2>(minv[j]) };
			left[j] = left[j]->right = &ptr[j];
			minv[j] = val;
		}
	}
	for (typename ett<_CostType, IndexType>::HeightType i = 0; i < maxht; ++i) {
		if (minv[i + 1] > minv[i])
			minv[i + 1] = minv[i];
		left[i]->right = above[i] = &basis[i];
		left[i]->data = { .i = std::get<1>(minv[i]), .j = std::get<2>(minv[i]) };
		std::get<1>(minv[i]) = basis[i].data.i;
		std::get<2>(minv[i]) = basis[i].data.j;
		std::get<0>(minv[i]) = cost_tree.get_value(basis[i].data);
	}
	for (size_t i = 1; i < m; ++i) {
		PNodeChunk[i - 1].right = &PNodeChunk[i];
		PNodeChunk[i].ptr = cost_tree.alloc->node_malloc(ht[i]);
		PNodeChunk[i].ht = ht[i];
		const IndexType label = tour[i];
		std::iota(left, left + ht[i], PNodeChunk[i].ptr);
		for (typename ett<_CostType, IndexType>::HeightType j = 0; j < ht[i]; ++j) {
			if (minv[j + 1] > minv[j])
				minv[j + 1] = minv[j];
			above[j]->data = { .i = std::get<1>(minv[j]), .j = std::get<2>(minv[j]) };
			above[j] = above[j]->right;
			above[j]->down = &PNodeChunk[i].ptr[j];
			minv[j] = std::tuple <_CostType, IndexType, IndexType> (0, label, label);
		}
		for (size_t j = i + 1; (j = j == m ? 0 : j) != i; ++j) {
			const typename ett<_CostType, IndexType>::HeightType h = std::min(ht[i], ht[j]);
			typename ett<_CostType, IndexType>::node *const ptr = cost_tree.alloc->node_malloc(h);
			std::tuple <_CostType, IndexType, IndexType> val;
			std::get<1>(val) = label, std::get<2>(val) = tour[j],
			std::get<0>(val) = cost_tree.get_value({ .i = label, .j = tour[j] });
			for (typename ett<_CostType, IndexType>::HeightType k = 0; k < h; ++k) {
				if (minv[k + 1] > minv[k])
					minv[k + 1] = minv[k];
				above[k] = above[k]->right;
				left[k]->data = { .i = std::get<1>(minv[k]), .j = std::get<2>(minv[k]) };
				left[k] = left[k]->right = above[k]->down = &ptr[k];
				minv[k] = val;
			}
			for (typename ett<_CostType, IndexType>::HeightType k = h; k < ht[j]; ++k) {
				if (minv[k + 1] > minv[k])
					minv[k + 1] = minv[k];
				above[k]->data = { .i = std::get<1>(minv[k]), .j = std::get<2>(minv[k]) };
				above[k] = above[k]->right;
				std::get<1>(minv[k]) = above[k]->data.i,
				std::get<2>(minv[k]) = above[k]->data.j,
				std::get<0>(minv[k]) = cost_tree.get_value(above[k]->data);
			}
		}
		for (typename ett<_CostType, IndexType>::HeightType j = 0; j < ht[i]; ++j) {
			if (minv[j + 1] > minv[j])
				minv[j + 1] = minv[j];
			left[j]->right = &PNodeChunk[i].ptr[j];
			left[j]->data = { .i = std::get<1>(minv[j]), .j = std::get<2>(minv[j]) };
			above[j] = &PNodeChunk[i].ptr[j];
			std::get<1>(minv[j]) = PNodeChunk[i].ptr[j].data.i,
			std::get<2>(minv[j]) = PNodeChunk[i].ptr[j].data.j,
			std::get<0>(minv[j]) = cost_tree.get_value(PNodeChunk[i].ptr[j].data);
		}
	}
	for (typename ett<_CostType, IndexType>::HeightType j = 0; j < maxht; ++j) {
		if (minv[j + 1] > minv[j])
			minv[j + 1] = minv[j];
		above[j]->data = { .i = std::get<1>(minv[j]), .j = std::get<2>(minv[j]) };
	}
	for (auto [i, p] = std::make_pair(0ul, basis); i < m; ++i, p = p->right)
		for (typename ett<_CostType, IndexType>::HeightType j = 0; j < ht[i]; ++j)
			above[j] = above[j]->right, above[j]->down = &p[j];
#ifndef NDEBUG
	{
		const typename ett<_CostType, IndexType>::node *p = basis;
		for (size_t i = 0; i < m; ++i, p = p->down) {
			const typename ett<_CostType, IndexType>::node *q = p;
			for (size_t j = 0; j < m; ++j, q = q->right)
				for (typename ett<_CostType, IndexType>::HeightType k = 0; k < std::min(ht[i], ht[j]); ++k) {
					std::tuple <_CostType, IndexType, IndexType> minv(std::numeric_limits<_CostType>::max(), -1, -1);
					for (size_t I = i; I < m && (I == i || ht[I] < k); ++I)
						for (size_t J = j; J < m && (J == j || ht[J] < k); ++J) 
							minv = std::min(minv, std::make_tuple(cost_tree.get_value({ .i = tour[i], .j = tour[j] }), tour[i], tour[j]));
					assert(std::get<1>(minv) == q->data.i);
					assert(std::get<2>(minv) == q->data.j);
				}
			assert(q == p);
		}
		assert(p == basis);
	}
	for (auto [e, t] : edges) {
		auto [u, v, puu, pvv] = t;
		assert(cost_tree.get_value({ .i = u, .j = v }) == 0);
	}
#endif
#if defined(SINKHORN) or not defined(NDEBUG)
	for (auto [minval, from, to] = minv[maxht]; minval < 0; ) {
		auto [e, f] = lct::sendflow(lctvertex[to].get(), lctvertex[from].get());
		total_cost += minval * f;
		auto it = edges.find(e);
		assert(it != edges.end());
		auto [u, v, puu, pvv] = it->second;
#ifndef SINKHORN
		assert(0);
#endif
		lct::cut(lctvertex[u].get(), lctvertex[v].get());
		edges.erase(it);
		edges.emplace(lct::link(lctvertex[from].get(), lctvertex[to].get(), f), std::make_tuple(from, to, puu, pvv));
		std::tie(minval, from, to) = cost_tree.template cut<0>(puu, pvv);
	}
#endif
#ifndef NDEBUG
	for (IndexType i = 0; i < n; ++i)
		for (IndexType j = 0; j < n; ++j)
			assert(cost_tree.get_value({ .i = i, .j = j }) >= 0);
	for (auto [e, t] : edges) {
		auto [u, v, puu, pvv] = t;
		assert(cost_tree.get_value({ .i = u, .j = v }) == 0);
	}
#endif
#ifdef DEBUG
	cost_tree.check();
#endif
}

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral _CostType, std::integral IndexType>
size_t emd<_PointType, _FlowType, _CostType, IndexType>::size() const { return cost_tree.size(); }

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral _CostType, std::integral IndexType>
void emd<_PointType, _FlowType, _CostType, IndexType>::modify(IndexType i, const _PointType &p) {
	const size_t n = size();
	points[i] = p;
	for (size_t j = 0; j < n; ++j)
		if (i != j)
			cost_tree[i][j] = dist_func(points[i], points[j]);
	cost_tree.update_node(i);
	for (auto [le, d] : edges) {
		auto [u, v, puu, pvv] = d;
		if (_CostType old_val = cost_tree.get_value({ .i = u, .j = v }); old_val) {
			total_cost += lct::get_flow(lctvertex[u].get(), lctvertex[v].get()) * old_val;
			assert((u != i) != (v != i));
			cost_tree.update_edge(puu, pvv, old_val);
			assert(cost_tree.get_value({ .i = u, .j = v }) == 0);
		}
	}
#ifdef TIME_LOG
	static std::ofstream fout("log.txt");
	size_t cnt = 0;
	uint64_t t = 0;
#endif
#ifdef DEBUG
	cost_tree.check();
#endif
	for (std::tuple <_CostType, IndexType, IndexType> minv = cost_tree.min(); assert(minv == check()), std::get<0>(minv) < 0; ) {
		auto [minval, from, to] = minv;
		auto [e, f] = lct::sendflow(lctvertex[to].get(), lctvertex[from].get());
		total_cost += minval * f;
		auto it = edges.find(e);
		assert(it != edges.end());
		auto [u, v, puu, pvv] = it->second;
		lct::cut(lctvertex[u].get(), lctvertex[v].get());
		edges.erase(it);
		edges.emplace(lct::link(lctvertex[from].get(), lctvertex[to].get(), f), std::make_tuple(from, to, puu, pvv));
#ifdef TIME_LOG
		t -= std::chrono::system_clock::now().time_since_epoch().count();
		++cnt;
#endif
		minv = cost_tree.template cut<0>(puu, pvv);
#ifdef TIME_LOG
		t += std::chrono::system_clock::now().time_since_epoch().count();
#endif
	}
#ifdef TIME_LOG
	fout << cnt << ' ' << t << '\n';
#endif
	assert(check() == (std::tuple <_CostType, IndexType, IndexType>(0, 0, 0)));
}

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral _CostType, std::integral IndexType>
std::tuple <_CostType, IndexType, IndexType> emd<_PointType, _FlowType, _CostType, IndexType>::check() {
	const size_t n = size();
	for (auto [e, t] : edges) {
		auto [u, v, puu, pvv] = t;
		assert(cost_tree.get_value({ .i = u, .j = v }) == 0);
	}
	std::tuple <_CostType, IndexType, IndexType> minv(std::numeric_limits<_CostType>::max(), -1, -1);
	for (IndexType i = 0; i < n; ++i)
		for (IndexType j = 0; j < n; ++j)
			minv = std::min(minv, std::make_tuple(cost_tree.get_value({ .i = i, .j = j }), i, j));
	return minv;
}

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral  _CostType, std::integral IndexType>
emd<_PointType, _FlowType, _CostType, IndexType>::emd(emd <_PointType, _FlowType, _CostType, IndexType> && val): cost_tree(), lctvertex(), edges(), points(), flows(), dist_func(), total_cost(), erased_nodes() { *this = std::move(val); }

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral  _CostType, std::integral IndexType>
void emd<_PointType, _FlowType, _CostType, IndexType>::erase(IndexType idx) {
	if (std::find(erased_nodes.begin(), erased_nodes.end(), idx) != erased_nodes.end())
		throw std::runtime_error("Point already erased");
	erased_nodes.push_back(idx);
}

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral  _CostType, std::integral IndexType>
IndexType emd<_PointType, _FlowType, _CostType, IndexType>::insert_demand(const _PointType &p) { return insert<1>(p); }

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral  _CostType, std::integral IndexType>
IndexType emd<_PointType, _FlowType, _CostType, IndexType>::insert_supply(const _PointType &p) { return insert<0>(p); }

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral  _CostType, std::integral IndexType>
template <bool flag>
IndexType emd<_PointType, _FlowType, _CostType, IndexType>::insert(const _PointType &p) {
	if (erased_nodes.empty()) {
		size_t n = points.size();
		assert(lctvertex.size() == n);
		assert(flows.size() == n);
		lctvertex.emplace_back(lct::add_node<_FlowType>());
		std::vector <_CostType> distarr(n + 1);
		flows.push_back(0);
		points.push_back(p);
		for (size_t i = 0; i < n; ++i)
			distarr[i] = dist_func(points[i], p);
		auto [from, to, puu, pvv] = cost_tree.template insert<flag>(distarr);
		edges.emplace(lct::link(lctvertex[from].get(), lctvertex[to].get(), _FlowType()), std::make_tuple(from, to, puu, pvv));
		if (flag)
			assert(to == n);
		else
			assert(from == n);
		return n;
	}
	else {
		IndexType ret = erased_nodes.back();
		erased_nodes.pop_back();
		modify(ret, p);
		return ret;
	}
}

template <std::integral IndexType = uint32_t>
auto make_emd(const auto &p, const auto &f, const auto &dist) {
	typedef typename std::remove_reference_t<decltype(p)>::value_type PointType;
	typedef typename std::remove_reference_t<decltype(f)>::value_type FlowType;
	typedef std::result_of_t<decltype(dist)(const PointType&, const PointType&)> CostType;
	return emd<PointType, FlowType, CostType, IndexType>(p, f, dist);
}
