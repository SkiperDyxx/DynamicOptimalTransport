#pragma once
#include <initializer_list>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include <ctime>
#include <vector>
#include <deque>
#include <cassert>
#include <cstdint>
#include <random>
#include <cstring>
#include <limits>
#include <unordered_map>
#include <cstdlib>
#include <algorithm>
#include <bit>
#include <list>
#include <iostream>
#include <concepts>

#ifdef PUSH_UP_LOG
#include <chrono>
size_t push_up_funcs = 0;
size_t push_up_cnt = 0;
size_t push_up_n = 0;
uint64_t push_up_time = 0;
#endif

template <typename _PointType, std::signed_integral _FlowType, std::signed_integral _CostType, std::integral _IndexType>
class emd;

template <typename T = int64_t, typename IndexType = uint32_t>
class ett {
	typedef uint8_t HeightType;
	class node {
	public:
#ifndef NDEBUG
		bool _isTop;
		static bool isTop(const node* ptr);
#endif
		node *right, *down;
		struct DataType {
			IndexType i, j;
			constexpr bool operator == (const DataType &a) const { return i == a.i && j == a.j; }
		} data;
	};
	class pnode {
	public:
		pnode *right;
		node *ptr;
		HeightType ht;
		IndexType data() const;
	};
	T *potential, *distance;
	struct allocator {
		size_t _Cap, _Size, unallocatedNodeChunkSize, unallocatedPNodeChunkSize;
		static constexpr size_t CHUNK = std::max(0x80000000ul, 0x4000000ul);
		pnode *top;
		node* unallocatedNodeChunk;
		pnode* unallocatedPNodeChunk;
		std::vector <node*> _AllocatedNodePool;
		std::vector <pnode*> _AllocatedPNodePool;
		void reserve(size_t capacity);
		node* node_malloc(size_t size);
		pnode* pnode_malloc();
		allocator();
		~allocator();
	} *alloc;
	struct distance_temporary_class {
		T *dptr;
		IndexType i;
		T& operator [](IndexType j);
	};
	T get_value(typename node::DataType val);
	inline T get_value(IndexType i, IndexType j);
	// There is no need for lazy tags & push downs,
	// as the data of each node is its real value plus the potential difference
	// and the potential have already been stored
	void _reserve(size_t capacity);
	inline void push_up_node(node *ptr);
	inline void push_up_node(node *cr, node *cd);
#ifndef DEBUG
	void check();
#endif
	pnode* get_pnode(IndexType u, std::vector <HeightType> &qht, node* head[]);
public:
#ifndef NDEBUG
	static std::mt19937_64 rng;
#endif
#ifdef DEBUG
	void check();
#endif
	template <bool k>
	std::tuple <T, IndexType, IndexType> cut(pnode *u, pnode *v);
	template <bool flag = true>
	std::tuple <IndexType, IndexType, pnode*, pnode*> insert(const auto &arr);
	T operator ()(IndexType u, IndexType v) const;
	distance_temporary_class operator [](IndexType i);
	void reserve(size_t capacity);
	size_t size() const;
//	std::tuple <T, IndexType, IndexType> update_edges(const auto &u, IndexType v);
	void update_edge(const pnode *const puu, const pnode *const pvv, T out_dated_value);
	void update_node(IndexType u);
	[[nodiscard("min has no side effects")]]
	std::tuple <T, IndexType, IndexType> min();
	ett operator = (const ett &val) = delete;
	ett&& operator = (ett && val);
	ett(ett && val);
	template <typename _PointType, std::signed_integral _FlowType, std::signed_integral _CostType, std::integral _IndexType>
	friend class emd;
	ett();
	~ett();
};

template <typename T, typename IndexType>
ett<T, IndexType>&& ett<T, IndexType>::operator = (ett<T, IndexType> && val) {
	std::swap(potential, val.potential);
	std::swap(distance, val.distance);
	std::swap(alloc, val.alloc);
	return std::move(*this);
}

template <typename T, typename IndexType>
void ett<T, IndexType>::update_node(IndexType u) {
#ifdef PUSH_UP_LOG
	++push_up_funcs;
#endif
	static node *uu[64];
	static std::vector <HeightType> hu; hu.reserve(2 * size() - 2);
	pnode *ptr = get_pnode(u, hu, uu);
	assert(hu.size() == size() * 2 - 2);
	const HeightType maxhu = alloc->top->ht;
	while (true) {
		assert(hu.back() == ptr->ht);
		static std::pair <node*, node*> iter[64];
		for (HeightType i = 1; i < maxhu; ++i)
			iter[i] = std::make_pair(uu[i], uu[i]);
		for (HeightType ht : hu)
			for (HeightType i = 1; i < ht; ++i) {
				if (iter[i].first != iter[i].second)
					push_up_node(iter[i].first, iter[i].second);
				iter[i].first = iter[i].first->right;
				iter[i].second = iter[i].second->down;
			}
		for (HeightType i = 1; i < maxhu; ++i)
			push_up_node(uu[i]);
		int R = 0;
		do {
			ptr = ptr->right; ++R;
			std::iota(uu, uu + ptr->ht, &ptr->ptr[0]);
			if (ptr == alloc->top) {
				check();
				return;
			}
		} while (ptr->data() != u);
		hu.insert(hu.end(), hu.begin(), hu.begin() + R);
		hu.erase(hu.begin(), hu.begin() + R);
	}
}

template <typename T, typename IndexType>
void ett<T, IndexType>::update_edge(const typename ett<T, IndexType>::pnode *const puu, const typename ett<T, IndexType>::pnode* const pvv, T out_dated_value) {
#ifdef PUSH_UP_LOG
	++push_up_funcs;
#endif
	static node *uu[64], *vv[64];
	static std::pair <node*, node*> iter[64];
	std::fill(uu, uu + 64, nullptr);
	std::fill(vv, vv + 64, nullptr);
	static std::vector <HeightType> hu, hv;
	hu.clear(); hv.clear();
	const size_t N = size(), M = 2 * N - 2;
	hu.reserve(M); hv.reserve(M);
	std::vector <uint8_t> vis(N);
	for (const pnode *p = pvv; p != puu; p = p->right) {
		std::iota(uu, uu + p->ht, p->ptr);
		hu.push_back(p->ht);
		if (auto u = p->data(); !vis[u]) {
			potential[u] -= out_dated_value;
			vis[u] = true;
		}
	}
	for (const pnode *p = puu; p != pvv; p = p->right) {
		std::iota(vv, vv + p->ht, p->ptr);
		hv.push_back(p->ht);
	}
	assert(hu.size() + hv.size() == M);
	T minval = *std::min_element(potential, potential + N);
	std::for_each(potential, potential + N, [minval](T &x) { x -= minval; });
	const HeightType maxhu = std::find(uu, uu + 64, nullptr) - uu, maxhv = std::find(vv, vv + 64, nullptr) - vv;
	hu.insert(hu.end(), hv.begin(), hv.end());
	hv.insert(hv.end(), hu.begin(), hu.begin() + hu.size() - hv.size());
	for (HeightType i = 1; i < maxhu; ++i)
		iter[i] = std::make_pair(uu[i], uu[i]);
	for (const HeightType ht : hv)
		for (HeightType h = 1; h < std::min(maxhu, ht); ++h) {
			if (iter[h].first != iter[h].second)
				push_up_node(iter[h].first, iter[h].second);
			iter[h].first = iter[h].first->right;
			iter[h].second = iter[h].second->down;
		}
#ifndef NDEBUG
	for (size_t i = 1; i < maxhu; ++i) {
		assert(iter[i].first == uu[i]);
		assert(iter[i].second == uu[i]);
	}
#endif
	for (HeightType i = 1; i < maxhv; ++i)
		iter[i] = std::make_pair(vv[i], vv[i]);
	for (HeightType i = maxhv; i < maxhu; ++i)
		iter[i] = std::make_pair(uu[i], uu[i]);
	for (const HeightType ht : hu)
		for (HeightType h = 1; h < ht; ++h) {
			if (iter[h].first != iter[h].second)
				push_up_node(iter[h].first, iter[h].second);
			iter[h].first = iter[h].first->right;
			iter[h].second = iter[h].second->down;
		}
#ifndef NDEBUG
	for (size_t i = 1; i < maxhv; ++i) {
		assert(iter[i].first == vv[i]);
		assert(iter[i].second == vv[i]);
	}
	for (HeightType i = maxhv; i < maxhu; ++i) {
		assert(iter[i].first == uu[i]);
		assert(iter[i].second == uu[i]);
	}
#endif
	for (HeightType i = 1; i < std::max(maxhu, maxhv); ++i) {
		if (i < maxhu)
			push_up_node(uu[i]);
		if (i < maxhv)
			push_up_node(vv[i]);
	}
	check();
}

template <typename T, typename IndexType>
std::tuple <T, IndexType, IndexType> ett<T, IndexType>::min() {
	const node *p = &alloc->top->ptr[alloc->top->ht - 1];
	size_t n = 1;
	std::tuple <T, IndexType, IndexType> ret;
	std::get<1>(ret) = p->data.i;
	std::get<2>(ret) = p->data.j;
	std::get<0>(ret) = get_value(p->data);
	for (const node *q = p->right; q != p; q = q->right, ++n)
		if (auto val = std::make_tuple(get_value(q->data), q->data.i, q->data.j); val < ret)
			ret = val;
#ifndef NDEBUG
	const node *const ptr = p;
#endif
	for (size_t i = 1; i < n; ++i) {
		p = p->down;
		const node *q = p;
		for (size_t j = 0; j < n; ++j, q = q->right)
			if (auto val = std::make_tuple(get_value(q->data), q->data.i, q->data.j); val < ret)
				ret = val;
		assert(q == p);
	}
	assert(p->down == ptr);
	return ret;
}

template <typename T, typename IndexType> template <bool flag>
std::tuple <T, IndexType, IndexType> ett<T, IndexType>::cut(pnode *puu, pnode *pvv) {
#ifdef PUSH_UP_LOG
	++push_up_funcs;
#endif
	const size_t n = size(), m = 2 * n - 2;
	if (n == 2) {
		assert(puu->right == pvv);
		assert(pvv->right == puu);
		if (flag) {
			std::swap(potential[0], potential[1]);
			if (puu->data() == 1) {
				check();
				assert(potential[0] >= potential[1]);
				return std::tuple <T, IndexType, IndexType>(2 * (potential[0] - potential[1]), 1, 0);
			}
			else {
				check();
				assert(potential[0] <= potential[1]);
				return std::tuple <T, IndexType, IndexType>(2 * (potential[1] - potential[0]), 0, 1);
			}
		}
		else {
			check();
			return std::tuple <T, IndexType, IndexType>(0, 0, 0);
		}
	}
	pnode *pij = pvv, *pji = puu;
	auto iter_right_func = [](node* &v) { v = v->right; }, iter_down_func = [](node* &v) { v = v->down; };
	HeightType huu = puu->ht, hvv = pvv->ht;
	static node *uu[64], *vv[64], *uv[64], *vu[64];
	std::fill(uu, uu + 64, nullptr);
	std::fill(vv, vv + 64, nullptr);
	static std::vector <HeightType> hu, hv;
	hu.clear(); hv.clear(); hu.reserve(2 * m); hv.reserve(2 * m);
	while (pij->right != puu) {
		pij = pij->right;
		hu.push_back(pij->ht);
		std::iota(uu, uu + pij->ht, pij->ptr);
	}
	const HeightType maxhu = std::find(uu, uu + 64, nullptr) - uu;
	std::copy(uu, uu + maxhu, uv);
	std::copy(uu, uu + maxhu, vu);
	while (pji->right != pvv) {
		pji = pji->right;
		hv.push_back(pji->ht);
		const HeightType h = std::min(pji->ht, maxhu);
		std::iota(vv, vv + pji->ht, pji->ptr);
		std::for_each(uv, uv + h, iter_right_func);
		std::for_each(vu, vu + h, iter_down_func);
	}
	assert(hu.size() + hv.size() == m - 2);
	const HeightType maxhv = std::find(vv, vv + 64, nullptr) - vv, minh = std::min(maxhu, maxhv);
	static node *iter_above[64], *iter_below[64], *iter_left[64], *iter_right[64];
	static node *quud[64], *quur[64], *quvd[64], *quvr[64], *qvud[64], *qvur[64], *qvvr[64], *qvvd[64];
	if (maxhu == 0 || maxhv == 0) {
		if (maxhu == 0) {
			std::swap(puu, pvv);
			std::swap(huu, hvv);
			std::swap(pij, pji);
			std::swap(hu, hv);
			std::copy(vv, vv + maxhv, uu);
		}
		node *const quu = puu->ptr, *const qvv = pvv->ptr, *const quv = quu->right, *const qvu = quu->down;
		assert(quv->down == qvv && qvu->right == qvv);
		const HeightType maxh = std::max(maxhu, maxhv);
		std::copy(uu, uu + maxh, quur);
		std::copy(uu, uu + maxh, quud);
		std::for_each(quur, quur + std::min(maxh, huu), iter_right_func);
		std::for_each(quud, quud + std::min(maxh, huu), iter_down_func);
		std::copy(quur, quur + maxh, qvvr);
		std::copy(quud, quud + maxh, qvvd);
		pij->right = pvv->right;
		std::for_each(qvvr, qvvr + std::min(maxh, hvv), iter_right_func);
		std::for_each(qvvd, qvvd + std::min(maxh, hvv), iter_down_func);
		std::copy(qvvd, qvvd + maxh, iter_left);
		std::copy(qvvr, qvvr + maxh, iter_above);
		for (HeightType i = 0; i < maxh; ++i) {
			if (i < huu) {
				quur[i]->down = quur[i]->down->down;
				quud[i]->right = quud[i]->right->right;
			}
			if (i < hvv) {
				qvvr[i]->down = qvvr[i]->down->down;
				qvvd[i]->right = qvvd[i]->right->right;
			}
			if (i < huu && i < hvv) {
				qvvr[i]->down = qvvr[i]->down->down;
				qvvd[i]->right = qvvd[i]->right->right;
				quur[i]->down = quur[i]->down->down;
				quud[i]->right = quud[i]->right->right;
			}
		}
		for (const HeightType h : hu) {
			iter_above[0] = iter_above[0]->right;
			iter_left[0] = iter_left[0]->down;
			iter_above[0]->down = iter_above[0]->down->down->down;
			iter_left[0]->right = iter_left[0]->right->right->right;
			for (HeightType i = 1; i < h; ++i) {
				assert((iter_above[i] == qvvr[i]) == (iter_left[i] == qvvd[i]));
				if (iter_above[i] != qvvr[i])
					push_up_node(iter_above[i], iter_left[i]);
				iter_above[i] = iter_above[i]->right;
				iter_left[i] = iter_left[i]->down;
				if (i < huu) {
					iter_above[i]->down = iter_above[i]->down->down;
					iter_left[i]->right = iter_left[i]->right->right;
				}
				if (i < hvv) {
					iter_above[i]->down = iter_above[i]->down->down;
					iter_left[i]->right = iter_left[i]->right->right;
				}
			}
		}
		for (HeightType i = 1; i < maxh; ++i) {
			assert(iter_above[i] == uu[i]);
			assert(iter_left[i] == uu[i]);
			push_up_node(uu[i]);
		}
		const IndexType v = puu->data();
		std::tuple <T, IndexType, IndexType> ret(std::numeric_limits<T>::max(), -1, -1);
		IndexType u;
		const bool f = flag != (maxhv == 0);
		if (f) {
			for (IndexType u = 0; u < v; ++u)
				ret = std::min(ret, std::tuple <T, IndexType, IndexType> (get_value(u, v), u, v));
			for (IndexType u = v + 1; u < n; ++u)
				ret = std::min(ret, std::tuple <T, IndexType, IndexType> (get_value(u, v), u, v));
			u = std::get<1>(ret);
			potential[v] += std::get<0>(ret);
		}
		else {
			for (IndexType u = 0; u < v; ++u)
				ret = std::min(ret, std::tuple <T, IndexType, IndexType> (get_value(v, u), v, u));
			for (IndexType u = v + 1; u < n; ++u)
				ret = std::min(ret, std::tuple <T, IndexType, IndexType> (get_value(v, u), v, u));
			u = std::get<2>(ret);
			potential[v] -= std::get<0>(ret);
		}
		size_t R = 0;
		while (pij->data() != u) {
			pij = pij->right;
			std::iota(uu, uu + pij->ht, pij->ptr);
			std::iota(quur, quur + std::min(huu, pij->ht), quur[0]->down);
			std::iota(quud, quud + std::min(huu, pij->ht), quud[0]->right);
			std::iota(qvvr, qvvr + std::min(hvv, pij->ht), qvvr[0]->down);
			std::iota(qvvd, qvvd + std::min(hvv, pij->ht), qvvd[0]->right);
			++R;
		}
		hu.insert(hu.end(), hu.begin(), hu.begin() + R);
		hu.erase(hu.begin(), hu.begin() + R);
		pvv->right = pij->right;
		pij->right = puu;
		quu->data.i = quu->data.j = quv->data.i = qvu->data.j = v;
		qvv->data.i = qvv->data.j = qvu->data.i = quv->data.j = u;
		for (HeightType i = 0; i < std::min({ maxh, huu, hvv }); ++i) {
			quv[i].right = quud[i]->right;
			qvu[i].down  = quur[i]->down;
			quur[i]->down = quud[i]->right = &quu[i];
			qvv[i].right = qvvd[i]->right;
			qvv[i].down  = qvvr[i]->down;
			qvvd[i]->right = &qvu[i];
			qvvr[i]->down = &quv[i];
			quud[i] = &quv[i];
			quur[i] = &qvu[i];
			qvvr[i] = qvvd[i] = &qvv[i];
		}
		for (HeightType i = hvv; i < std::min(huu, maxh); ++i) {
			quu[i].right = quud[i]->right;
			quu[i].down = quur[i]->down;
			quur[i] = quud[i] = quur[i]->down = quud[i]->right = &quu[i];
		}
		for (HeightType i = huu; i < std::min(hvv, maxh); ++i) {
			qvv[i].right = qvvd[i]->right;
			qvv[i].down = qvvr[i]->down;
			qvvd[i] = qvvr[i] = qvvd[i]->right = qvvr[i]->down = &qvv[i];
		}
		node *euud = quud[0], *euur = quur[0], *evvd = qvvd[0], *evvr = qvvr[0];
		std::copy(uu, uu + maxh, iter_above);
		std::copy(uu, uu + maxh, iter_left);
		for (const HeightType h : hu) {
			assert(euud->down == evvd);
			assert(euur->right == evvr);
			iter_above[0] = iter_above[0]->right;
			iter_left[0] = iter_left[0]->down;
			euud = euud->right;
			evvd = evvd->right;
			euur = euur->down;
			evvr = evvr->down;
			euud[0].data.i = euur[0].data.j = v;
			evvd[0].data.i = evvr[0].data.j = u;
			assert(iter_above[0]->data.j == iter_left[0]->data.i);
			euud[0].data.j = euur[0].data.i = evvd[0].data.j = evvr[0].data.i = iter_above[0]->data.j;
			evvd[0].down = iter_above[0]->down;
			evvr[0].right = iter_left[0]->right;
			iter_above[0]->down = euud;
			iter_left[0]->right = euur;
			for (HeightType i = 1; i < h; ++i) {
				assert(i >= huu || i >= hvv || quud[i]->down == qvvd[i]);
				assert(i >= huu || i >= hvv || quur[i]->right == qvvr[i]);
				node *r = iter_above[i], *d = iter_left[i];
				if (r != d) {
					push_up_node(r, d);
					if (i < huu)
						push_up_node(r = r->down, d = d->right);
					if (i < hvv)
						push_up_node(r = r->down, d = d->right);
				}
				iter_above[i] = iter_above[i]->right;
				iter_left[i] = iter_left[i]->down;
				if (i < hvv) {
					evvd[i].down = iter_above[i]->down;
					evvr[i].right = iter_left[i]->right;
				}
				else if (i < huu) {
					euud[i].down = iter_above[i]->down;
					euur[i].right = iter_left[i]->right;
				}
				if (i < huu) {
					iter_above[i]->down = &euud[i];
					iter_left[i]->right = &euur[i];
				}
				else if (i < hvv) {
					iter_above[i]->down = &evvd[i];
					iter_left[i]->right = &evvr[i];
				}
			}
		}
		for (HeightType i = 1; i < std::max({ maxh, huu, hvv }); ++i) {
			if (i < maxh) {
				push_up_node(uu[i]);
				node *r = uu[i], *d = uu[i];
				if (i < huu)
					push_up_node(r = r->right, d = d->down);
				if (i < hvv)
					push_up_node(r = r->right, d = d->down);
			}
			if (i < huu) {
				push_up_node(&quu[i]);
				if (i < hvv)
					push_up_node(&quv[i], &qvu[i]);
			}
			if (i < hvv)
				push_up_node(&qvv[i]);
		}
		check();
		return flag ? ret : min();
	}
	std::for_each(uv, uv + std::min(maxhu, huu), iter_right_func);
	std::for_each(vu, vu + std::min(maxhu, huu), iter_down_func);
	node *const qvv = pvv->ptr, *const quu = puu->ptr, *const quv = uv[0]->down->right, *const qvu = vu[0]->down->right;
	for (HeightType i = 0; i < huu; ++i) {
		if (i < maxhu) {
			quud[i] = uu[i]->down;
			quur[i] = uu[i]->right;
			if (i < maxhv) {
				quvd[i] = uv[i]->down;
				qvur[i] = vu[i]->right;
				if (i < hvv) {
					quud[i]->right = quv[i].right;
					quur[i]->down  = qvu[i].down;
				}
				else {
					quud[i]->right = quvd[i]->right;
					quur[i]->down  = qvur[i]->down;
				}
				quvd[i]->right = quu[i].right;
				qvur[i]->down  = quu[i].down;
			}
			else if (i < hvv) {
				quud[i]->right = quv[i].right;
				quur[i]->down  = qvu[i].down;
			}
			else {
				quud[i]->right = quu[i].right;
				quur[i]->down  = quu[i].down;
			}
		}
		else if (i < maxhv) {
			quvd[i] = vv[i]->down;
			qvur[i] = vv[i]->right;
			if (i < hvv) {
				quvd[i] = quvd[i]->down;
				qvur[i] = qvur[i]->right;
			}
			quvd[i]->right = quu[i].right;
			qvur[i]->down  = quu[i].down;
		}
	}
	for (HeightType i = 0; i < hvv; ++i) {
		if (i < maxhv) {
			qvvd[i] = vv[i]->down;
			qvvr[i] = vv[i]->right;
			if (i < maxhu) {
				quvr[i] = uv[i]->right;
				qvud[i] = vu[i]->down;
				if (i < huu) {
					qvvd[i]->right = qvu[i].right;
					qvvr[i]->down  = quv[i].down;
				}
				else {
					qvvd[i]->right = qvud[i]->right;
					qvvr[i]->down  = quvr[i]->down;
				}
				quvr[i]->down  = qvv[i].down;
				qvud[i]->right = qvv[i].right;
			}
			else if (i < huu) {
				qvvr[i]->down  = quv[i].down;
				qvvd[i]->right = qvu[i].right;
			}
			else {
				qvvr[i]->down  = qvv[i].down;
				qvvd[i]->right = qvv[i].right;
			}
		}
		else if (i < maxhu) {
			qvud[i] = uu[i]->down;
			quvr[i] = uu[i]->right;
			if (i < huu) {
				qvud[i] = qvud[i]->down;
				quvr[i] = quvr[i]->right;
			}
			qvud[i]->right = qvv[i].right;
			quvr[i]->down  = qvv[i].down;
		}
	}
	std::copy(uu, uu + maxhu, iter_left);
	std::copy(quud, quud + huu, iter_left);
	std::copy(uu, uu + maxhu, iter_above);
	std::copy(quur, quur + huu, iter_above);
	std::copy(vv, vv + maxhv, iter_right);
	std::copy(vv, vv + maxhv, iter_below);
	std::copy(qvvr, qvvr + std::min(hvv, maxhv), iter_below);
	std::copy(qvvd, qvvd + std::min(hvv, maxhv), iter_right);
	std::copy(vu, vu + std::min(maxhv, maxhu), iter_below);
	std::copy(uv, uv + std::min(maxhv, maxhu), iter_right);
	std::copy(qvur, qvur + std::min(huu, maxhv), iter_below);
	std::copy(quvd, quvd + std::min(huu, maxhv), iter_right);
	uint64_t enable = 0;
	for (auto h : hv) {
		const HeightType mh = std::min(maxhu, h);
		for (HeightType i = 0; i < h; ++i) {
			if (i != 0 && (enable >> i & 1)) {
				push_up_node(iter_below[i], iter_right[i]);
				if (i < mh)
					push_up_node(iter_above[i], iter_left[i]);
			}
			iter_below[i] = iter_below[i]->right;
			iter_right[i] = iter_right[i]->down;
			if (i < maxhu) {
				iter_above[i] = iter_above[i]->right;
				iter_left [i] = iter_left [i]->down;
			}
		}
		enable |= (1ull << h) - 1;
		for (HeightType i = std::max(huu, hvv); i < mh; ++i) {
			std::swap(iter_above[i]->down, iter_below[i]->down);
			std::swap(iter_left[i]->right, iter_right[i]->right);
		}
		for (HeightType i = hvv; i < std::min(mh, huu); ++i) {
			std::tie(iter_above[i]->down, iter_below[i]->down) = std::make_pair(iter_below[i]->down, iter_above[i]->down->down);
			std::tie(iter_left[i]->right, iter_right[i]->right) = std::make_pair(iter_right[i]->right, iter_left[i]->right->right);
		}
		for (HeightType i = huu; i < std::min(mh, hvv); ++i) {
			std::tie(iter_above[i]->down, iter_below[i]->down) = std::make_pair(iter_below[i]->down->down, iter_above[i]->down);
			std::tie(iter_left[i]->right, iter_right[i]->right) = std::make_pair(iter_right[i]->right->right, iter_left[i]->right);
		}
		for (HeightType i = 0; i < std::min({ mh, huu, hvv }); ++i) {
			std::tie(iter_above[i]->down, iter_below[i]->down) = std::make_pair(iter_below[i]->down->down, iter_above[i]->down->down);
			std::tie(iter_left[i]->right, iter_right[i]->right) = std::make_pair(iter_right[i]->right->right, iter_left[i]->right->right);
		}
		for (HeightType i = mh; i < h; ++i) {
			if (i < hvv) {
				iter_down_func(iter_below[i]->down);
				iter_right_func(iter_right[i]->right);
			}
			if (i < huu) {
				iter_down_func(iter_below[i]->down);
				iter_right_func(iter_right[i]->right);
			}
		}
	}
#ifndef NDEBUG
	for (HeightType i = 0; i < minh; ++i) {
		assert(iter_above[i] == uv[i]);
		assert(iter_below[i] == vv[i]);
		assert(iter_left [i] == vu[i]);
		assert(iter_right[i] == vv[i]);
	}
#endif
	for (HeightType i = 1; i < maxhv; ++i)
		push_up_node(vv[i]);
	std::copy(iter_left, iter_left + maxhu, iter_right);
	std::for_each(iter_right, iter_right + std::min(maxhu, hvv), iter_down_func);
	std::copy(iter_above, iter_above + maxhu, iter_below);
	std::for_each(iter_below, iter_below + std::min(maxhu, hvv), iter_right_func);
	std::swap_ranges(iter_left, iter_left + maxhu, iter_above);
	enable = 0;
	for (HeightType h : hu) {
		const HeightType mh = std::min(h, maxhv);
		for (HeightType i = 0; i < h; ++i) {
			assert((iter_above[i] == vu[i]) == (iter_left[i] == uv[i]));
			if (i != 0 && (enable >> i & 1)) {
				push_up_node(iter_below[i], iter_right[i]);
				if (i < mh)
					push_up_node(iter_above[i], iter_left[i]);
			}
			iter_below[i] = iter_below[i]->right;
			iter_right[i] = iter_right[i]->down;
			if (i < mh) {
				iter_above[i] = iter_above[i]->right;
				iter_left [i] = iter_left [i]->down;
			}
		}
		enable |= (1ull << h) - 1;
		for (HeightType i = std::max(huu, hvv); i < mh; ++i) {
			std::swap(iter_above[i]->down, iter_below[i]->down);
			std::swap(iter_left[i]->right, iter_right[i]->right);
		}
		for (HeightType i = hvv; i < std::min(mh, huu); ++i) {
			std::tie(iter_above[i]->down, iter_below[i]->down) = std::make_pair(iter_below[i]->down->down, iter_above[i]->down);
			std::tie(iter_left[i]->right, iter_right[i]->right) = std::make_pair(iter_right[i]->right->right, iter_left[i]->right);
		}
		for (HeightType i = huu; i < std::min(mh, hvv); ++i) {
			std::tie(iter_above[i]->down, iter_below[i]->down) = std::make_pair(iter_below[i]->down, iter_above[i]->down->down);
			std::tie(iter_left[i]->right, iter_right[i]->right) = std::make_pair(iter_right[i]->right, iter_left[i]->right->right);
		}
		for (HeightType i = 0; i < std::min({ mh, huu, hvv }); ++i) {
			std::tie(iter_above[i]->down, iter_below[i]->down) = std::make_pair(iter_below[i]->down->down, iter_above[i]->down->down);
			std::tie(iter_left[i]->right, iter_right[i]->right) = std::make_pair(iter_right[i]->right->right, iter_left[i]->right->right);
		}
		for (HeightType i = mh; i < h; ++i) {
			if (i < huu) {
				iter_down_func(iter_below[i]->down);
				iter_right_func(iter_right[i]->right);
			}
			if (i < hvv) {
				iter_down_func(iter_below[i]->down);
				iter_right_func(iter_right[i]->right);
			}
		}
	}
#ifndef NDEBUG
	for (HeightType i = 0; i < minh; ++i) {
		assert(vu[i] == iter_above[i]);
		assert(uu[i] == iter_below[i]);
		assert(uv[i] == iter_left [i]);
		assert(uu[i] == iter_right [i]);
	}
#endif
	for (HeightType i = 1; i < maxhu; ++i)
		push_up_node(uu[i]);
	for (HeightType i = 1; i < std::min(maxhu, maxhv); ++i)
		push_up_node(vu[i], uv[i]);
	T minval, retval;
	IndexType newu, newv;
	typename node::DataType ret; {
		const node *const basis = flag ? iter_above[minh - 1] : iter_left[minh - 1];
		ret = basis->data;
		minval = get_value(ret);
		size_t n = 1;
		for (const node *p = basis->right; p != basis; ++n, p = p->right) {
			assert(node::isTop(p));
			if (auto nval = get_value(p->data); std::make_tuple(nval, p->data.i, p->data.j) < std::make_tuple(minval, ret.i, ret.j))
				minval = nval, ret = p->data;
		}
		for (const node *q = basis->down; q != basis; q = q->down) {
			const node *p = q;
			for (size_t j = 0; j < n; ++j, p = p->right) {
				assert(node::isTop(p));
				if (auto nval = get_value(p->data); std::make_tuple(nval, p->data.i, p->data.j) < std::make_tuple(minval, ret.i, ret.j))
					minval = nval, ret = p->data;
			}
			assert(p == q);
		}
	}
	if (flag)
		newu = ret.j, newv = ret.i, retval = minval;
	else {
		newu = ret.i, newv = ret.j;
		minval = -minval;
	}
	std::vector <uint8_t> vis(n);
	for (pnode *q = pvv->right; q != puu; q = q->right) {
		assert(q->data() >= 0 && q->data() < n);
		if (!vis[q->data()]) {
			potential[q->data()] += minval;
			vis[q->data()] = 1;
		}
	}
	minval = *std::min_element(potential, potential + n);
	std::for_each(potential, potential + n, [minval](T &x) { x -= minval; });
	pij->right = pvv->right;
	pji->right = puu->right;
	size_t R = 0;
	while (pij->data() != newu) {
		pij = pij->right;
		std::iota(uu, uu + pij->ht, pij->ptr);
		assert(pij->ht == hu[R]);
		++R;
		const HeightType h = pij->ht;
		std::for_each(vu, vu + std::min(h, minh), iter_right_func);
		std::for_each(uv, uv + std::min(h, minh), iter_down_func);
		std::iota(quud, quud + std::min(h, huu), quud[0]->right);
		std::iota(quur, quur + std::min(h, huu), quur[0]->down);
		std::iota(qvud, qvud + std::min(h, hvv), qvud[0]->right);
		std::iota(quvr, quvr + std::min(h, hvv), quvr[0]->down);
	};
	hu.insert(hu.end(), hu.begin(), hu.begin() + R);
	hu.erase(hu.begin(), hu.begin() + R);
	R = 0;
	while (pji->data() != newv) {
		pji = pji->right;
		std::iota(vv, vv + pji->ht, pji->ptr);
		assert(pji->ht == hv[R]);
		++R;
		const HeightType h = pji->ht;
		std::for_each(vu, vu + std::min(h, minh), iter_down_func);
		std::for_each(uv, uv + std::min(h, minh), iter_right_func);
		std::iota(qvvr, qvvr + std::min(h, hvv), qvvr[0]->down);
		std::iota(qvur, qvur + std::min(h, huu), qvur[0]->down);
		std::iota(qvvd, qvvd + std::min(h, hvv), qvvd[0]->right);
		std::iota(quvd, quvd + std::min(h, huu), quvd[0]->right);
	}
	hv.insert(hv.end(), hv.begin(), hv.begin() + R);
	hv.erase(hv.begin(), hv.begin() + R);
	puu->right = pji->right;
	pvv->right = pij->right;
	pij->right = puu;
	pji->right = pvv;
	for (HeightType i = 0; i < huu; ++i) {
		if (i < maxhu) {
			if (i < maxhv) {
				quu[i].right = quvd[i]->right;
				quu[i].down  = qvur[i]->down;
				if (i < hvv) {
					quv[i].right = quud[i]->right;
					qvu[i].down  = quur[i]->down;
					quvd[i]->right = &quv[i];
					qvur[i]->down  = &qvu[i];
				}
				else {
					quvd[i]->right = quud[i]->right;
					qvur[i]->down  = quur[i]->down;
				}
			}
			else if (i < hvv) {
				quv[i].right = quud[i]->right;
				qvu[i].down  = quur[i]->down;
			}
			else {
				quu[i].right = quud[i]->right;
				quu[i].down  = quur[i]->down;
			}
			quud[i]->right = &quu[i];
			quur[i]->down  = &quu[i];
		}
		else if (i < maxhv) {
			quu[i].right = quvd[i]->right;
			quu[i].down  = qvur[i]->down;
			if (i < hvv) {
				quvd[i]->right = &quv[i];
				qvur[i]->down  = &qvu[i];
			}
			else {
				quvd[i]->right = &quu[i];
				qvur[i]->down = &quu[i];
			}
		}
	}
	for (HeightType i = 0; i < hvv; ++i) {
		if (i < maxhv) {
			if (i < maxhu) {
				qvv[i].right = qvud[i]->right;
				qvv[i].down  = quvr[i]->down;
				if (i < huu) {
					qvu[i].right = qvvd[i]->right;
					quv[i].down  = qvvr[i]->down;
					qvud[i]->right = &qvu[i];
					quvr[i]->down  = &quv[i];
				}
				else {
					qvud[i]->right = qvvd[i]->right;
					quvr[i]->down  = qvvr[i]->down;
				}
			}
			else if (i < huu) {
				qvu[i].right = qvvd[i]->right;
				quv[i].down  = qvvr[i]->down;
			}
			else {
				qvv[i].right = qvvd[i]->right;
				qvv[i].down  = qvvr[i]->down;
			}
			qvvd[i]->right = &qvv[i];
			qvvr[i]->down  = &qvv[i];
		}
		else if (i < maxhu) {
			qvv[i].right = qvud[i]->right;
			qvv[i].down  = quvr[i]->down;
			if (i < huu) {
				qvud[i]->right = &qvu[i];
				quvr[i]->down  = &quv[i];
			}
			else {
				qvud[i]->right = &qvv[i];
				quvr[i]->down  = &qvv[i];
			}
		}
	}
	memcpy(iter_above, vu, sizeof(iter_above));
	memcpy(iter_below, uu, sizeof(iter_below));
	memcpy(iter_left,  uv, sizeof(iter_left));
	memcpy(iter_right, uu, sizeof(iter_right));
	node *edge_above = qvv, *edge_below = quv, *edge_left = qvv, *edge_right = qvu;
	quu->data.i = quu->data.j = quv->data.i = qvu->data.j = newv;
	qvv->data.i = qvv->data.j = qvu->data.i = quv->data.j = newu;
	enable = 0;
	for (HeightType h : hu) {
		for (HeightType i = 0; i < h; ++i) {
			if (i != 0 && (enable >> i & 1)) {
				node *r = iter_below[i], *d = iter_right[i];
				push_up_node(r, d); r = r->down; d = d->right;
				if (i < huu) {
					push_up_node(r, d);
					r = r->down; d = d->right;
				}
				if (i < maxhv) {
					r = iter_above[i]; d = iter_left[i];
					push_up_node(r, d);
					r = r->down; d = d->right;
				}
				if (i < hvv)
					push_up_node(r, d);
			}
			iter_below[i] = iter_below[i]->right;
			iter_right[i] = iter_right[i]->down;
			if (i < maxhv) {
				iter_above[i] = iter_above[i]->right;
				iter_left[i] = iter_left[i]->down;
			}
		}
		enable |= (1 << h) - 1;
		edge_above = edge_above->right;
		edge_below = edge_below->right;
		edge_left  = edge_left->down;
		edge_right = edge_right->down;
		assert(iter_above[0]->data.j == iter_below[0]->data.j && iter_left[0]->data.i == iter_right[0]->data.i && iter_above[0]->data.j == iter_left[0]->data.i);
		edge_above->data.i = edge_left->data.j = newu;
		edge_below->data.i = edge_right->data.j = newv;
		edge_above->data.j = edge_left->data.i = edge_below->data.j = edge_right->data.i = iter_above[0]->data.j;
		for (HeightType i = 0; i < std::min(h, hvv); ++i) {
			edge_above[i].down = iter_below[i]->down;
			edge_left[i].right = iter_right[i]->right;
		}
		for (HeightType i = std::max(maxhv, hvv); i < std::min(h, huu); ++i) {
			edge_below[i].down  = iter_below[i]->down;
			edge_right[i].right = iter_right[i]->right;
		}
		for (HeightType i = 0; i < std::min({ h, huu, maxhv }); ++i) {
			edge_below[i].down = iter_above[i]->down;
			edge_right[i].right = iter_left[i]->right;
		}
		for (HeightType i = std::max(huu, hvv); i < std::min(h, maxhv); ++i) {
			std::swap(iter_above[i]->down, iter_below[i]->down);
			std::swap(iter_left[i]->right, iter_right[i]->right);
		}
		for (HeightType i = huu; i < std::min({ h, hvv, maxhv }); ++i) {
			iter_below[i]->down = iter_above[i]->down;
			iter_right[i]->right = iter_left[i]->right;
		}
		for (HeightType i = hvv; i < std::min({ h, huu, maxhv }); ++i) {
			iter_above[i]->down = iter_below[i]->down;
			iter_left[i]->right = iter_right[i]->right;
		}
		for (HeightType i = 0; i < std::min(huu, h); ++i) {
			iter_below[i]->down = &edge_below[i];
			iter_right[i]->right = &edge_right[i];
		}
		for (HeightType i = std::max(huu, maxhv); i < std::min(h, hvv); ++i) {
			iter_below[i]->down = &edge_above[i];
			iter_right[i]->right = &edge_left[i];
		}
		for (HeightType i = 0; i < std::min({ h, hvv, maxhv }); ++i) {
			iter_above[i]->down = &edge_above[i];
			iter_left[i]->right = &edge_left[i];
		}
	}
	if (huu < maxhu) {
		for (HeightType i = 0; i < huu; ++i) {
			iter_above[i] = iter_below[i]->right;
			iter_left[i] = iter_right[i]->down;
		}
		std::copy(iter_below + huu, iter_below + maxhu, iter_above + huu);
		std::copy(iter_right + huu, iter_right + maxhu, iter_left + huu);
	}
	else
		for (HeightType i = 0; i < maxhu; ++i) {
			iter_above[i] = iter_below[i]->right;
			iter_left[i] = iter_right[i]->down;
		}
	std::copy(vv, vv + maxhv, iter_below);
	std::copy(vv, vv + maxhv, iter_right);
	edge_above = quu; edge_left  = quu;
	edge_right = quv; edge_below = qvu;
	enable = 0;
	for (HeightType h: hv) {
		for (HeightType i = 0; i < h; ++i) {
			if (i != 0 && (enable >> i & 1)) {
				node *r = iter_below[i], *d = iter_right[i];
				push_up_node(r, d); r = r->down; d = d->right;
				if (i < hvv) {
					push_up_node(r, d);
					r = r->down; d = d->right;
				}
				if (i < maxhu) {
					r = iter_above[i]; d = iter_left[i];
					push_up_node(r, d);
					r = r->down; d = d->right;
				}
				if (i < huu)
					push_up_node(r, d);
			}
			iter_below[i] = iter_below[i]->right;
			iter_right[i] = iter_right[i]->down;
			if (i < maxhu) {
				iter_above[i] = iter_above[i]->right;
				iter_left[i] = iter_left[i]->down;
			}
		}
		enable |= (1 << h) - 1;
		edge_above = edge_above->right;
		edge_below = edge_below->right;
		edge_left  = edge_left->down;
		edge_right = edge_right->down;
		assert(iter_above[0]->data.j == iter_below[0]->data.j && iter_left[0]->data.i == iter_right[0]->data.i && iter_above[0]->data.j == iter_left[0]->data.i);
		edge_above->data.i = edge_left->data.j = newv;
		edge_below->data.i = edge_right->data.j = newu;
		edge_above->data.j = edge_left->data.i = edge_below->data.j = edge_right->data.i = iter_above[0]->data.j;
		for (HeightType i = 0; i < std::min(h, huu); ++i) {
			edge_above[i].down = iter_below[i]->down;
			edge_left[i].right = iter_right[i]->right;
		}
		for (HeightType i = std::max(huu, maxhu); i < std::min(h, hvv); ++i) {
			edge_below[i].down  = iter_below[i]->down;
			edge_right[i].right = iter_right[i]->right;
		}
		for (HeightType i = 0; i < std::min({ maxhu, h, hvv }); ++i) {
			edge_below[i].down = iter_above[i]->down;
			edge_right[i].right = iter_left[i]->right;
		}
		for (HeightType i = std::max(huu, hvv); i < std::min(h, maxhu); ++i) {
			std::swap(iter_above[i]->down, iter_below[i]->down);
			std::swap(iter_left[i]->right, iter_right[i]->right);
		}
		for (HeightType i = hvv; i < std::min({ h, huu, maxhu }); ++i) {
			iter_below[i]->down = iter_above[i]->down;
			iter_right[i]->right = iter_left[i]->right;
		}
		for (HeightType i = huu; i < std::min({ h, hvv, maxhu }); ++i) {
			iter_above[i]->down = iter_below[i]->down;
			iter_left[i]->right = iter_right[i]->right;
		}
		for (HeightType i = 0; i < std::min(hvv, h); ++i) {
			iter_below[i]->down = &edge_below[i];
			iter_right[i]->right = &edge_right[i];
		}
		for (HeightType i = std::max(hvv, maxhu); i < std::min(h, huu); ++i) {
			iter_below[i]->down = &edge_above[i];
			iter_right[i]->right = &edge_left[i];
		}
		for (HeightType i = 0; i < std::min({ h, huu, maxhu }); ++i) {
			iter_above[i]->down = &edge_above[i];
			iter_left[i]->right = &edge_left[i];
		}
	}
	for (HeightType i = 1; i < std::max({ huu, hvv, maxhu, maxhv }); ++i) {
		if (i < maxhu) {
			push_up_node(uu[i]);
			node *r = uu[i], *d = uu[i];
			if (i < huu)
				push_up_node(r = r->down, d = d->right);
			if (i < maxhv)
				push_up_node(r = vu[i], d = uv[i]);
			if (i < hvv)
				push_up_node(r = r->down, d = d->right);
		}
		if (i < maxhv) {
			push_up_node(vv[i]);
			node *r = vv[i], *d = vv[i];
			if (i < hvv)
				push_up_node(r = r->down, d = d->right);
			if (i < maxhu)
				r = uv[i], d = vu[i];
			if (i < huu)
				push_up_node(r = r->down, d = d->right);
		}
		if (i < huu)
			push_up_node(&quu[i]);
		if (i < hvv)
			push_up_node(&qvv[i]);
		if (i < std::min(huu, hvv))
			push_up_node(&qvu[i], &quv[i]);
	}
	check();
	return flag ? std::make_tuple(retval, ret.i, ret.j) : min();
}

#ifndef NDEBUG
template <typename T, typename IndexType>
std::mt19937_64 ett<T, IndexType>::rng((std::random_device())());
#endif

template <typename T, typename IndexType> template <bool flag>
std::tuple <IndexType, IndexType, typename ett<T, IndexType>::pnode*, typename ett<T, IndexType>::pnode*> ett<T, IndexType>::insert(const auto &_arr) {
#ifdef PUSH_UP_LOG
	++push_up_funcs;
#endif
	const size_t n = alloc->_Size++;
	if (n == std::numeric_limits<T>::max())
		throw std::length_error("Maximum node number: " + std::to_string(std::numeric_limits<T>::max()));
	const size_t m = 2 * n - 2;
	assert(_arr.size() == n || (_arr.size() == n + 1 && _arr.back() == T()));
	reserve(n + 1);
	pnode *const q[2] = { alloc->pnode_malloc(), alloc->pnode_malloc() };
	std::copy(_arr.begin(), _arr.begin() + n, distance + n * (n - 1) / 2);
	const HeightType HT[2] = { q[0]->ht, q[1]->ht };
	if (n == 1) {
		if (alloc->top->ht < q[0]->ht)
			alloc->top = q[0];
		if (alloc->top->ht < q[1]->ht)
			alloc->top = q[1];
		q[0]->right = q[1];
		q[1]->right = q[0];
		potential[!flag] = 0;
		potential[flag] = _arr[0];
		node *const p[2][2] = { alloc->node_malloc(HT[0]), alloc->node_malloc(std::min(HT[0], HT[1])), alloc->node_malloc(std::min(HT[0], HT[1])), alloc->node_malloc(HT[1]) };
		for (HeightType h = 0; h < std::min(HT[0], HT[1]); ++h)
			for (int i = 0; i < 2; ++i)
				for (int j = 0; j < 2; ++j) {
					p[i][j ^ 1][h].right = p[i ^ 1][j][h].down = &p[i][j][h];
					p[i][j][h].data.i = i;
					p[i][j][h].data.j = j;
				}
		for (HeightType h = HT[0]; h < HT[1]; ++h) {
			p[0][0][h].data.i = p[0][0][h].data.j = 0;
			p[0][0][h].right = p[0][0][h].down = &p[0][0][h];
		}
		for (HeightType h = HT[1]; h < HT[0]; ++h) {
			p[1][1][h].data.i = p[1][1][h].data.j = 0;
			p[1][1][h].right = p[1][1][h].down = &p[1][1][h];
		}
		check();
		return std::tuple<IndexType, IndexType, pnode*, pnode*>(!flag, flag, q[flag], q[!flag]);
	}
	size_t u = 0;
	if (flag) {
		T pn = potential[0] + _arr[0];
		for (size_t i = 1; i < n; ++i)
			if (T p = potential[i] + _arr[i]; pn > p)
				u = i, pn = p;
		potential[n] = pn;
	}
	else {
		T pn = potential[0] - _arr[0];
		for (size_t i = 1; i < n; ++i)
			if (T p = potential[i] - _arr[i]; pn < p)
				u = i, pn = p;
		potential[n] = pn;
	}
	static node* head[64];
	static std::vector <HeightType> qht; qht.reserve(m);
	pnode* const quu = get_pnode(u, qht, head);
	if (alloc->top->ht < q[0]->ht)
		alloc->top = q[0];
	if (alloc->top->ht < q[1]->ht)
		alloc->top = q[1];
	const HeightType maxht = std::find(head, head + 64, nullptr) - head;
	static node *rnode[64], *dnode[64], *rtail[2][64], *dtail[2][64];
	std::fill(rnode, rnode + 64, nullptr);
	std::fill(dnode, dnode + 64, nullptr);
	for (HeightType i = 0; i < maxht; ++i) {
		rnode[i] = head[i];
		dnode[i] = head[i];
	}
	node *p[2][2];
	for (int i = 1; i >= 0; --i)
		for (int j = 1; j >= 0; --j) {
			HeightType ht = std::min(HT[i], HT[j]);
			p[i][j] = alloc->node_malloc(ht);
            for (HeightType k = 0; k < ht; ++k) {
				if (rtail[i][k] != nullptr)
					rtail[i][k]->right = &p[i][j][k];
				if (dtail[j][k] != nullptr)
					dtail[j][k]->down = &p[i][j][k];
				rtail[i][k] = &p[i][j][k];
				dtail[j][k] = &p[i][j][k];
			}
		}
	p[1][1][0].data.i = p[1][1][0].data.j = p[1][0][0].data.i = p[0][1][0].data.j = n;
	p[0][0][0].data.i = p[0][0][0].data.j = p[0][1][0].data.i = p[1][0][0].data.j = u;
	q[0]->ptr = p[0][0];
	q[1]->ptr = p[1][1];
	q[0]->right = quu->right;
	q[1]->right = q[0];
	quu->right = q[1];
	const std::tuple <IndexType, IndexType, pnode*, pnode*> ret = flag ? std::make_tuple(u, n, q[1], q[0]) : std::make_tuple(n, u, q[0], q[1]);
	assert(qht.size() == m);
	for (HeightType h : qht) {
		const HeightType ht[2] = { std::min(HT[0], h), std::min(HT[1], h) };
		node *const PR[2] = { alloc->node_malloc(ht[0]), alloc->node_malloc(ht[1]) }, *const PD[2] = { alloc->node_malloc(ht[0]), alloc->node_malloc(ht[1]) };
		rnode[0] = rnode[0]->right;
		dnode[0] = dnode[0]->down;
		assert(rnode[0]->data.j == dnode[0]->data.i);
		const IndexType j = rnode[0]->data.j;
		PR[1]->data.j = PD[1]->data.i = PR[0]->data.j = PD[0]->data.i = j;
		PR[1]->data.i = PD[1]->data.j = n; PR[0]->data.i = PD[0]->data.j = u;
		PR[0]->down = rnode[0]->down;
		PR[1]->down = PR[0];
		rnode[0]->down = PR[1];
		PD[0]->right = dnode[0]->right;
		PD[1]->right = PD[0];
		dnode[0]->right = PD[1];
		for (int _ = 0; _ < 2; ++_) {
			rtail[_][0] = rtail[_][0]->right = PR[_];
			dtail[_][0] = dtail[_][0]->down = PD[_];
		}
		for (HeightType j = 1; j < h; ++j) {
			if (rnode[j] != dnode[j])
				push_up_node(rnode[j], dnode[j]);
			rnode[j] = rnode[j]->right;
			dnode[j] = dnode[j]->down;
			for (int _ = 0; _ < 2; ++_)
				if (j < HT[_]) {
					rtail[_][j]->right = &PR[_][j];
					PR[_][j].down = rnode[j]->down;
					rnode[j]->down = &PR[_][j];
					dtail[_][j]->down = &PD[_][j];
					PD[_][j].right = dnode[j]->right;
					dnode[j]->right = &PD[_][j];
					if (rtail[_][j] != dtail[_][j])
						push_up_node(rtail[_][j], dtail[_][j]);
					else
						push_up_node(rtail[_][j]);
					rtail[_][j] = &PR[_][j];
					dtail[_][j] = &PD[_][j];
				}
		}
	}
	for (int i = 0; i < 2; ++i) {
		rtail[i][0]->right = &p[i][1][0];
		dtail[i][0]->down = &p[1][i][0];
	}
	for (HeightType k = 1; k < std::max({ maxht, HT[0], HT[1] }); ++k) {
		for (int i = 0; i < 2; ++i)
			if (k < HT[i]) {
				rtail[i][k]->right = dtail[i][k]->down = &p[i][i][k];
				push_up_node(&p[i][i][k]);
			}
		if (k < std::min(HT[0], HT[1])) {
			rtail[0][k]->right = &p[0][1][k];
			dtail[0][k]->down = &p[1][0][k];
			push_up_node(&p[1][0][k], &p[0][1][k]);
		}
		if (k < maxht) {
			push_up_node(head[k]);
			for (int i = 1; i >= 0; --i)
				if (k < HT[i])
					push_up_node(rtail[i][k], dtail[i][k]);
		}
	}
	check();
	return ret;
}

template <typename T, typename IndexType>
void ett<T, IndexType>::push_up_node(node *ptr) {
	assert(!node::isTop(&ptr[-1]));
	typename node::DataType mindata = ptr[-1].data;
	T minv = get_value(mindata);
	size_t n = 1;
	for (const node *p = ptr[-1].right; p != &ptr->right[-1]; p = p->right, ++n) {
		assert(node::isTop(p));
		if (auto nval = get_value(p->data); std::make_tuple(nval, p->data.i, p->data.j) < std::make_tuple(minv, mindata.i, mindata.j))
			minv = nval, mindata = p->data;
	}
#ifdef PUSH_UP_LOG
	++push_up_cnt;
	push_up_n += n * n;
	push_up_time -= std::chrono::system_clock::now().time_since_epoch().count();
#endif
	const node *q = ptr[-1].down;
	for (size_t i = 1; i < n; ++i, q = q->down) {
		const node *p = q;
		for (size_t j = 0; j < n; ++j, p = p->right) {
			assert(node::isTop(p));
			if (auto nval = get_value(p->data); std::make_tuple(nval, p->data.i, p->data.j) < std::make_tuple(minv, mindata.i, mindata.j))
				minv = nval, mindata = p->data;
		}
	}
#ifdef PUSH_UP_LOG
	push_up_time += std::chrono::system_clock::now().time_since_epoch().count();
#endif
	assert(!node::isTop(q));
	ptr->data = mindata;
}

template <typename T, typename IndexType>
void ett<T, IndexType>::push_up_node(node *cr, node *cd) {
#ifdef PUSH_UP_LOG
	++push_up_cnt;
	push_up_time -= std::chrono::system_clock::now().time_since_epoch().count();
#endif
	std::tuple <T, IndexType, IndexType> minr, mind;
	minr = mind = std::tuple <T, IndexType, IndexType> (std::numeric_limits<T>::max(), -1, -1);
	std::pair p(&cr[-1], &cd[-1]);
	size_t n = 0;
	auto push_up = [&minr, &mind, this](std::pair <node*, node*> p) {
		if (std::tuple nval(get_value(p.first->data), p.first->data.i, p.first->data.j); nval < minr)
			minr = nval;
		if (std::tuple nval(get_value(p.second->data), p.second->data.i, p.second->data.j); nval < mind)
			mind = nval;
	};
	do {
		push_up(p);
		p = std::make_pair(p.first->right, p.second->down);
		assert((p.first == &cr->right[-1]) == (p.second == &cd->down[-1]));
		++n;
	} while (p.first != &cr->right[-1]);
	for (std::pair q(cr[-1].down, cd[-1].right);
			assert((q.first == &cr->down[-1]) == (q.second == &cd->right[-1])),
			q.first != &cr->down[-1]; q = std::make_pair(q.first->down, q.second->right)) {
		p = q;
		for (size_t i = 0; i < n; ++i, p = std::make_pair(p.first->right, p.second->down))
			push_up(p);
	}
	cr->data = { .i = std::get<1>(minr), .j = std::get<2>(minr) };
	cd->data = { .i = std::get<1>(mind), .j = std::get<2>(mind) };
#ifdef PUSH_UP_LOG
	push_up_time += std::chrono::system_clock::now().time_since_epoch().count();
#endif
}

template <typename T, typename IndexType>
typename ett<T, IndexType>::node* ett<T, IndexType>::allocator::node_malloc(size_t size) {
	if (unallocatedNodeChunkSize < size)
		_AllocatedNodePool.push_back(unallocatedNodeChunk = (node*)std::malloc((unallocatedNodeChunkSize = CHUNK / sizeof(node)) * sizeof(node)));
	node *ptr = &unallocatedNodeChunk[unallocatedNodeChunkSize -= size];
	--size;
#ifndef NDEBUG
	for (size_t k = 0; k < size; ++k)
		ptr[k]._isTop = false;
	ptr[size]._isTop = true;
	ptr->data.i = ptr->data.j = -1;
#endif
	return ptr;
}

template <typename T, typename IndexType>
typename ett<T, IndexType>::pnode* ett<T, IndexType>::allocator::pnode_malloc() {
	if (unallocatedPNodeChunkSize == 0)
		_AllocatedPNodePool.push_back(unallocatedPNodeChunk = (pnode*)std::malloc((unallocatedPNodeChunkSize = CHUNK / sizeof(pnode)) * sizeof(pnode)));
#ifdef NDEBUG
	static std::mt19937_64 rng((std::random_device())());
#endif
	--unallocatedPNodeChunkSize;
	while ((unallocatedPNodeChunk[unallocatedPNodeChunkSize].ht = std::countl_zero(rng()) + 1) > 64);
	return &unallocatedPNodeChunk[unallocatedPNodeChunkSize];
}

template <typename T, typename IndexType>
T ett<T, IndexType>::operator ()(IndexType u, IndexType v) const { return potential[v] - potential[u]; }

template <typename T, typename IndexType>
typename ett<T, IndexType>::distance_temporary_class ett<T, IndexType>::operator [](IndexType i) {
	distance_temporary_class ret = { .dptr = distance, .i = i };
	return ret;
}

template <typename T, typename IndexType>
T& ett<T, IndexType>::distance_temporary_class::operator [](IndexType j) {
	size_t di, dj;
	if (i > j)
		di = i, dj = j;
	else if (i < j)
		di = j, dj = i;
	else
		throw std::out_of_range("Two subscripts cannot be equal");
	return dptr[di * (di - 1) / 2 + dj];
}

template <typename T, typename IndexType>
void ett<T, IndexType>::_reserve(size_t capacity) {
	reserve(capacity);
}

template <typename T, typename IndexType>
void ett<T, IndexType>::reserve(size_t capacity) {
	assert(capacity <= std::numeric_limits<IndexType>::max());
	const size_t _Cap = alloc->_Cap, newN = capacity;
	if (capacity <= _Cap) return;
	distance = (T*)realloc(distance, newN * (newN + 1) / 2 * sizeof(T));
	potential = (T*)realloc(potential, newN * sizeof(T));
	alloc->reserve(capacity);
}

template <typename T, typename IndexType>
void ett<T, IndexType>::allocator::reserve(size_t capacity) { _Cap = capacity; }

template <typename T, typename IndexType>
T ett<T, IndexType>::get_value(IndexType i, IndexType j) { return get_value({ .i = i, .j = j }); }

template <typename T, typename IndexType>
T ett<T, IndexType>::get_value(typename ett<T, IndexType>::node::DataType data) {
	assert(data.i < alloc->_Size);
	assert(data.j < alloc->_Size);
	return data.i == data.j ? 0 : (*this)[data.i][data.j] + potential[data.i] - potential[data.j];
}

template <typename T, typename IndexType>
ett <T, IndexType>::ett(): potential(nullptr), distance(nullptr), alloc(new allocator()) {}

template <typename T, typename IndexType>
ett <T, IndexType>::allocator::allocator(): _Cap(0), _Size(1), unallocatedNodeChunkSize(0), unallocatedPNodeChunkSize(0), top(nullptr), _AllocatedNodePool(), _AllocatedPNodePool() {}

#ifndef NDEBUG
template <typename T, typename IndexType>
bool ett <T, IndexType>::node::isTop(const ett <T, IndexType>::node *ptr) {
	return ptr->_isTop;
}
#endif

template <typename T, typename IndexType>
IndexType ett<T, IndexType>::pnode::data() const {
	assert(ptr->data.i == ptr->data.j);
	return ptr->data.i;
}
template <typename T, typename IndexType>
size_t ett<T, IndexType>::size() const { return alloc->_Size; }

template <typename T, typename IndexType>
ett <T, IndexType>::~ett() {
	free(potential);
	free(distance);
	delete alloc;
}

template <typename T, typename IndexType>
ett<T, IndexType>::allocator::~allocator() {
	for (void* ptr: _AllocatedNodePool)
		free(ptr);
	for (void* ptr: _AllocatedPNodePool)
		free(ptr);
}

template <typename T, typename IndexType>
ett<T, IndexType>::pnode* ett<T, IndexType>::get_pnode(IndexType u, std::vector <typename ett<T, IndexType>::HeightType> &qht, ett<T, IndexType>::node* head[64]) {
	std::fill(head, head + 64, nullptr);
	static std::vector <HeightType> qend; qend.clear(); qht.clear();
	pnode* p = alloc->top;
	while (true) {
		qend.push_back(p->ht);
		std::iota(&head[0], &head[p->ht], &p->ptr[0]);
		if (p->data() == u)
			break;
		p = p->right;
	}
	pnode* ret = p;
	while ((p = p->right) != alloc->top)
		qht.push_back(p->ht);
	qht.insert(qht.end(), qend.begin(), qend.end());
	return ret;
}

template <typename T, typename IndexType>
void ett<T, IndexType>::check() {
#ifdef DEBUG
	std::deque <IndexType> seq;
	pnode *pn = alloc->top;
	assert(size() > 1 || (size() == 1 && pn == nullptr));
	const size_t m = 2 * size() - 2;
	size_t cnt[64]; cnt[0] = m; std::fill(cnt + 1, cnt + 64, std::numeric_limits<size_t>::max());
	for (size_t i = 0; i < m; ++i, pn = pn->right) {
		assert(i == 0 || pn != alloc->top);
		for (HeightType k = 0; k < pn->ht; ++k)
			assert((k == pn->ht - 1) == node::isTop(&pn->ptr[k]));
		for (HeightType k = 0; k < pn->ht; ++k) {
			node *r = pn->ptr[k].right, *d = pn->ptr[k].down;
			if (k > 0) {
				auto olddata = pn->ptr[k].data;
				push_up_node(&pn->ptr[k]);
				assert(pn->ptr[k].data == olddata);
			}
			size_t j = 1;
			while (r != d) {
				assert(j < cnt[k]);
				if (k > 0) {
					auto datar = r->data, datad = d->data;
					push_up_node(r, d);
					assert(r->data == datar);
					assert(d->data == datad);
				}
				else {
					assert(r->data.i == pn->data());
					assert(d->data.j == pn->data());
				}
				r = r->right, d = d->down, ++j;
			}
			assert(r == &pn->ptr[k]);
			if (cnt[k] == std::numeric_limits<size_t>::max())
				cnt[k] = j;
			assert(cnt[k] == j);
		}
	}
#endif
}
